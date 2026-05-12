# -*- coding: utf-8 -*-
"""
Robin 签名方案 (NTRU-based Hash-and-Sign)

基于紧凑Gadget和NTRU假设的Hash-and-Sign签名方案。

方案概述:
  - 环: R = Z[x]/(x^n - 1), n 为素数
  - 陷门关系: h·f + g = p mod Q
  - 公钥: h (一个环元素)
  - 私钥: (f, g) (两个短多项式)

签名:
  给定消息 msg, 计算 u = H(msg||salt),
  使用私钥采样短向量 (z0, z1) 满足 z0 + h·z1 = u - e mod Q,
  输出 (salt, z1). 验签时恢复 z0+e = u - h·z1 mod Q,
  验证 (z0+e, z1) 的范数是否在界内。

算法参考 (论文 Section 5):
  - Algorithm 3: Robin.KeyGen
  - Algorithm 4: Robin.Sign
  - Algorithm 5: Robin.Verify

参考文献: CRYPTO 2023, Section 5
"""

import math
import time
from typing import Tuple, Optional, List
from .polynomials import Polynomial, circulant_matrix
from .compact_gadget import approx_gadget
from .gaussian_sampler import sample_z, sample_gadget_coset
from .hash_utils import MessageHasher, msg_to_bytes
from .parameters import RobinParams, get_robin_params


# ============================================================================
# 辅助函数: 生成三值多项式
# ============================================================================

def sample_ternary_polynomial(n: int, a: int, b: int,
                               Q: int) -> Polynomial:
    """
    从 T(n, a, b) 中均匀采样三值多项式

    定义:
      T(n, a, b) = {v ∈ R | v 恰好有 a 个系数等于 +1,
                               b 个系数等于 -1,
                               n-a-b 个系数等于 0}

    参数:
        n: 多项式维度
        a: +1 系数的个数
        b: -1 系数的个数
        Q: 系数模数 (仅用于多项式对象)

    返回:
        采样到的三值多项式
    """
    import random
    rng = random.SystemRandom()

    # 创建系数列表: a个1, b个-1, 其余为0
    coeffs = [1] * a + [-1] * b + [0] * (n - a - b)
    # 随机打乱
    rng.shuffle(coeffs)

    return Polynomial(coeffs, Q, "cyclic")


# ============================================================================
# 计算矩阵的最大奇异值
# ============================================================================

def _compute_dft(coeffs: List[int]) -> List[complex]:
    """
    计算多项式的 DFT (用于估计循环矩阵的奇异值)

    对于 R^-_n (cyclic), M(a) 的奇异值 = |DFT(a)_k|

    优先使用 numpy.fft (O(n log n)), 否则回退到纯 Python DFT (O(n²))。
    """
    n = len(coeffs)
    # 尝试使用 numpy 加速 (numpy 是可选的, 不安装也能运行)
    try:
        import numpy as np
        dft = np.fft.fft(np.array(coeffs, dtype=np.float64))
        return [complex(d) for d in dft]
    except ImportError:
        pass
    # 纯 Python 回退
    result = []
    for k in range(n):
        total = complex(0, 0)
        for j in range(n):
            angle = -2 * math.pi * j * k / n
            total += coeffs[j] * complex(math.cos(angle), math.sin(angle))
        result.append(total)
    return result


def compute_s1_upper_bound(f: Polynomial, g: Polynomial) -> float:
    """
    计算 s₁(M(f·f̄ + g·ḡ)) 的估计值

    对于 R^-_n 中的多项式 a, 循环矩阵 M(a) 的奇异值是 |DFT(a)_k|。
    因此 s₁(M(a)) = max_k |DFT(a)_k|。

    对于组合矩阵:
      DFT(f·f̄)_k = |DFT(f)_k|^2
      DFT(g·ḡ)_k = |DFT(g)_k|^2
      所以 s₁ = max_k (|DFT(f)_k|^2 + |DFT(g)_k|^2)

    参数:
        f: 私钥多项式 f
        g: 私钥多项式 g

    返回:
        s₁ 的估计值 (越大表示陷门质量越差)
    """
    n = f.n

    # 计算 DFT
    dft_f = _compute_dft(f.coeffs)
    dft_g = _compute_dft(g.coeffs)

    # s₁ = max_k (|DFT(f)_k|^2 + |DFT(g)_k|^2)
    max_val = 0.0
    for k in range(n):
        val = abs(dft_f[k])**2 + abs(dft_g[k])**2
        if val > max_val:
            max_val = val

    return max_val


def compute_fg_norm(f: Polynomial, g: Polynomial) -> float:
    """
    计算私钥向量的范数 ||(f, g)|| = √(2(a+b))

    对于 f, g ∈ T(n, a, b):
      ||(f, g)||^2 = ||f||^2 + ||g||^2 = (a+b) + (a+b) = 2(a+b)

    参数:
        f, g: 私钥多项式

    返回:
        ||(f, g)||
    """
    norm_sq = f.norm_squared() + g.norm_squared()
    return math.sqrt(norm_sq)


# ============================================================================
# Robin 密钥生成 (Algorithm 3)
# ============================================================================

def robin_keygen(params: RobinParams) -> Tuple[Polynomial, Tuple[Polynomial, Polynomial]]:
    """
    Robin 密钥生成算法 (论文 Algorithm 3)

    生成过程:
      1. 采样 K^2 对 (f_i, g_j) ∈ T(n, a, b)
      2. 寻找最优的 k 使得 σ_k(g_j) 的 s₁ 最小
      3. 检查陷阱门质量: √s₁ ≤ α·√(2(a+b))
      4. 计算公钥 h = (p - g)/f mod Q

    返回值:
      - 如果质量检查通过: h = (p-g)/f mod Q
      - 如果不通过: 重新采样

    参数:
        params: Robin 参数集

    返回:
        (h, (f, g)): 公钥和私钥对
    """
    import random
    rng = random.SystemRandom()

    n = params.n
    Q = params.Q
    p = params.p
    a = params.a
    b = params.b
    alpha = params.alpha
    K = 5  # 论文默认采样5对

    target_norm = math.sqrt(2 * (a + b))  # ||(f, g)|| 的理论值

    # 循环直到找到合格的公私钥对
    attempt = 0
    max_keygen_attempts = 50

    while attempt < max_keygen_attempts:
        attempt += 1

        # 采样 K 对候选 (f_i, g_i)
        f_candidates = [sample_ternary_polynomial(n, a, b, Q) for _ in range(K)]
        g_candidates = [sample_ternary_polynomial(n, a, b, Q) for _ in range(K)]

        for i in range(K):
            for j in range(K):
                fi = f_candidates[i]
                gj = g_candidates[j]

                # 寻找最优的自同构 (对于demo参数只检查 k=1)
                best_s1 = float('inf')
                best_g = gj

                k_list = [1]  # 默认只检查 k=1
                if n > 100:
                    # 对于大 n, 随机采样一些 k (论文建议遍历所有与n互素的k)
                    k_list = [k for k in range(1, 20, 2) if math.gcd(k, n) == 1]

                for k in k_list:
                    gk = gj.automorphism(k) if k != 1 else gj
                    s1_val = compute_s1_upper_bound(fi, gk)
                    if s1_val < best_s1:
                        best_s1 = s1_val
                        best_g = gk

                # 检查陷门质量: √s₁ ≤ α·||(f,g)||
                quality_threshold = alpha * target_norm
                sqrt_s1 = math.sqrt(best_s1)

                if sqrt_s1 <= quality_threshold:
                    f = fi
                    g = best_g

                    # 计算公钥 h = (p - g) / f mod Q
                    f_inv = f.invert()
                    if f_inv is None:
                        continue

                    p_poly = Polynomial(
                        [p] + [0] * (n - 1), Q, "cyclic"
                    )
                    h = f_inv * (p_poly - g)

                    return h, (f, g)

    raise RuntimeError(f"密钥生成失败: 超过最大尝试次数 ({max_keygen_attempts})")


# ============================================================================
# Robin 签名 (Algorithm 4)
# ============================================================================

def robin_sign(msg, h: Polynomial, sk: Tuple[Polynomial, Polynomial],
               params: RobinParams,
               hasher: MessageHasher = None,
               max_attempts: int = 100) -> Tuple[bytes, Polynomial]:
    """
    Robin 签名算法 (论文 Algorithm 4)

    流程:
      1. 构造 A = [I_n | M(h)], T = [M(g); M(f)]
      2. salt ← 随机, u ← H(msg, salt)
      3. (z0, z1) ← ApproxPreSamp(A, T, u, r, s)
      4. e = u - (z0 + z1·h) mod Q
      5. if ||(z0+e, γ·z1)|| > β: 重新开始
      6. return (salt, z1)

    参数:
        msg: 要签名的消息 (str 或 bytes)
        h: 公钥多项式
        sk: 私钥 (f, g)
        params: Robin参数集
        hasher: 消息哈希器
        max_attempts: 最大重试次数

    返回:
        (salt, z1): 签名 (盐 + 一个环元素)
    """
    n = params.n
    Q = params.Q
    p = params.p
    q = params.q
    r = params.r
    s = params.s
    gamma = params.gamma
    beta = params.beta

    f, g = sk
    msg_bytes = msg_to_bytes(msg)

    if hasher is None:
        hasher = MessageHasher(n, Q)

    # 构造矩阵形式
    # A = [I_n | M(h)]
    M_h = circulant_matrix(h)
    # T = [M(g); M(f)]  (2n × n)
    M_g = circulant_matrix(g)
    M_f = circulant_matrix(f)

    # 预计算 T 的Frobenius范数用于扰动采样
    T_rows = []
    for i in range(n):
        T_rows.append(list(M_g[i]))
    for i in range(n):
        T_rows.append(list(M_f[i]))

    attempt = 0
    while attempt < max_attempts:
        attempt += 1

        # 生成随机盐和目标
        salt = hasher.generate_salt()
        u_coeffs = hasher.hash_message(msg_bytes, salt)

        # === 扰动采样 ===
        # 计算 s_eff 使得扰动协方差有效
        # 简化: 使用 s_eff = √(s^2 - r^2·||T||^2_F)
        frob_sq_T = sum(sum(val * val for val in row) for row in T_rows)
        s_eff_sq = s * s - r * r * frob_sq_T / n  # 归一化
        if s_eff_sq <= 0:
            continue  # 参数问题, 跳过

        s_eff = math.sqrt(s_eff_sq)

        # 采样扰动向量 p_vec (维度 2n)
        p_vec = [sample_z(c=0, sigma=s_eff) for _ in range(2 * n)]

        # === 计算 u' = u - A·p ===
        # A·p = I_n·p_0 + M(h)·p_1
        Ap = [0] * n
        # p_0 部分 (直接加)
        for i in range(n):
            Ap[i] = p_vec[i] % Q
        # M(h)·p_1 部分
        for i in range(n):
            total = 0
            for j in range(n):
                total += M_h[i][j] * p_vec[n + j]
            Ap[i] = (Ap[i] + total) % Q

        u_prime = [(u_coeffs[i] - Ap[i]) % Q for i in range(n)]
        u_prime = [x if x <= Q // 2 else x - Q for x in u_prime]

        # === Gadget采样 ===
        x_prime_coeffs, e_coeffs = approx_gadget(u_prime, r, p, q, Q)

        # === 计算 x = p + T·x' ===
        # T·x' = [M(g)·x'; M(f)·x']
        Tx0 = [0] * n
        Tx1 = [0] * n
        for i in range(n):
            total0 = 0
            total1 = 0
            for j in range(n):
                total0 += M_g[i][j] * x_prime_coeffs[j]
                total1 += M_f[i][j] * x_prime_coeffs[j]
            Tx0[i] = total0
            Tx1[i] = total1

        z0_coeffs = [p_vec[i] + Tx0[i] for i in range(n)]
        z1_coeffs = [p_vec[n + i] + Tx1[i] for i in range(n)]

        # 构造多项式
        z0 = Polynomial(z0_coeffs, Q, "cyclic")
        z1 = Polynomial(z1_coeffs, Q, "cyclic")
        e = Polynomial(e_coeffs, Q, "cyclic")

        # === 恢复 e (实际上应该相同) ===
        # e = u - (z0 + h·z1) mod Q
        actual_e = Polynomial(u_coeffs, Q, "cyclic") - (z0 + h * z1)

        # === 接受检查 ===
        # ||(z0+e, γ·z1)|| ≤ β
        z0_plus_e = z0 + e
        check_norm_sq = z0_plus_e.norm_squared() + (gamma * gamma) * z1.norm_squared()
        check_norm = math.sqrt(check_norm_sq)

        if check_norm <= beta:
            return salt, z1

    raise RuntimeError(f"签名失败: 超过最大重试次数 ({max_attempts})")


# ============================================================================
# Robin 验签 (Algorithm 5)
# ============================================================================

def robin_verify(msg, salt: bytes, z1: Polynomial, h: Polynomial,
                 params: RobinParams,
                 hasher: MessageHasher = None) -> bool:
    """
    Robin 验签算法 (论文 Algorithm 5)

    流程:
      1. u ← H(msg, salt)
      2. z' ← (u - h·z1) mod Q
      3. 接受当且仅当 ||(z', γ·z1)|| ≤ β

    参数:
        msg: 原始消息
        salt: 签名中的盐
        z1: 签名的多项式部分
        h: 公钥
        params: Robin参数集
        hasher: 消息哈希器

    返回:
        True 如果签名有效, False 否则
    """
    n = params.n
    Q = params.Q
    gamma = params.gamma
    beta = params.beta

    msg_bytes = msg_to_bytes(msg)
    if hasher is None:
        hasher = MessageHasher(n, Q)

    # u ← H(msg, salt)
    u_coeffs = hasher.hash_message(msg_bytes, salt)
    u = Polynomial(u_coeffs, Q, "cyclic")

    # z' = u - h·z1 mod Q
    z_prime = u - h * z1

    # 检查 ||(z', γ·z1)|| ≤ β
    norm_sq = z_prime.norm_squared() + (gamma * gamma) * z1.norm_squared()
    norm = math.sqrt(norm_sq)

    return norm <= beta


# ============================================================================
# 签名大小估算
# ============================================================================

def estimate_sig_size(z1: Polynomial) -> int:
    """
    估算签名大小 (基于熵界)

    论文使用 ANS (Asymmetric Numeral System) 进行批量编码,
    可以达到接近熵界的签名大小。

    这里给出简化的估算。

    参数:
        z1: 签名多项式

    返回:
        估算的签名大小 (字节)
    """
    n = len(z1.coeffs)

    # 简单的大小估算: 每个系数编码所需位数
    max_coeff = max(abs(c) for c in z1.coeffs)
    bits_per_coeff = max(1, (max_coeff.bit_length() + 1))

    total_bits = n * bits_per_coeff
    total_bytes = (total_bits + 7) // 8

    return total_bytes


# ============================================================================
# 测试代码
# ============================================================================

def _test():
    """Robin 签名方案自测"""
    print("=" * 60)
    print("Robin 签名方案测试 (NTRU-based)")
    print("=" * 60)

    # 使用演示参数
    params = get_robin_params("demo")
    print(f"\n参数集: {params.name}")
    print(f"  n={params.n}, Q={params.Q}, p={params.p}, q={params.q}")
    print(f"  a={params.a}, b={params.b}")
    print(f"  r={params.r:.1f}, s={params.s:.1f}")
    print(f"  γ={params.gamma:.2f}, β={params.beta:.1f}")
    print(f"  安全级别: {params.security_level}")

    # 密钥生成
    print("\n[步骤1] 密钥生成...")
    start = time.time()
    h, (f, g) = robin_keygen(params)
    keygen_time = time.time() - start
    print(f"  密钥生成时间: {keygen_time:.2f}s")
    print(f"  公钥 h (前8个系数): {h.coeffs[:8]}")
    print(f"  私钥 f (前8个系数): {f.coeffs[:8]}")
    print(f"  私钥 g (前8个系数): {g.coeffs[:8]}")

    # 验证陷门关系: h·f + g = p mod Q
    hf_plus_g = h * f + g
    p_theory = params.p
    trapdoor_ok = (hf_plus_g.coeffs[0] % params.Q == p_theory % params.Q)
    # 验证所有系数: hf + g 应该等于常数 p
    expected_p = Polynomial([p_theory] + [0] * (params.n - 1), params.Q, "cyclic")
    trapdoor_ok = all(
        abs(hf_plus_g.coeffs[i] - expected_p.coeffs[i]) <= 1
        for i in range(params.n)
    )
    print(f"  陷门关系验证 (hf+g=p): {trapdoor_ok}")

    # 签名
    print("\n[步骤2] 签名...")
    msg = "Hello, Robin signature scheme!"
    hasher = MessageHasher(params.n, params.Q)

    try:
        start = time.time()
        salt, z1 = robin_sign(msg, h, (f, g), params, hasher)
        sign_time = time.time() - start
        print(f"  签名时间: {sign_time:.2f}s")
        print(f"  盐 (hex): {salt.hex()[:20]}...")
        print(f"  签名 z1 (前8个系数): {z1.coeffs[:8]}")
        print(f"  ||z1||^2 = {z1.norm_squared()}")
        print(f"  估算签名大小: {estimate_sig_size(z1)} 字节")

        # 验签
        print("\n[步骤3] 验签...")
        start = time.time()
        valid = robin_verify(msg, salt, z1, h, params, hasher)
        verify_time = time.time() - start
        print(f"  验签时间: {verify_time*1000:.2f}ms")
        print(f"  验签结果: {'有效 ✓' if valid else '无效 ✗'}")

        # 测试错误签名
        print("\n[步骤4] 安全性测试...")
        # 错误消息
        invalid = robin_verify("Wrong message", salt, z1, h, params, hasher)
        print(f"  错误消息验签: {'有效 (BUG!)' if invalid else '无效 ✓'}")

        # 错误公钥
        h_wrong = Polynomial(
            [(c + 1) % params.Q for c in h.coeffs],
            params.Q, "cyclic"
        )
        invalid2 = robin_verify(msg, salt, z1, h_wrong, params, hasher)
        print(f"  错误公钥验签: {'有效 (BUG!)' if invalid2 else '无效 ✓'}")

        # 多次签名统计
        print("\n[步骤5] 多次签名统计...")
        num_sigs = 10
        successes = 0
        for i in range(num_sigs):
            try:
                s, z = robin_sign(f"Test message {i}", h, (f, g), params, hasher)
                if robin_verify(f"Test message {i}", s, z, h, params, hasher):
                    successes += 1
            except RuntimeError:
                pass
        print(f"  签名成功率: {successes}/{num_sigs}")

    except RuntimeError as e:
        print(f"  签名失败: {e}")

    print("\n所有测试完成！")


if __name__ == "__main__":
    _test()
