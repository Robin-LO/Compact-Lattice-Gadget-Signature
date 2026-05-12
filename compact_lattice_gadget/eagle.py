# -*- coding: utf-8 -*-
"""
Eagle 签名方案 (Ring-LWE-based Hash-and-Sign)

基于紧凑Gadget和Ring-LWE假设的Hash-and-Sign签名方案。

方案概述:
  - 环: R = Z[x]/(x^n + 1), n 为 2 的幂
  - 公钥: (seed_a, b) 其中 b = p - (a·f + g) mod Q
  - 私钥: (f, g) (两个短的三值多项式)
  - 公钥矩阵: A = [I_n | M(a) | M(b)]  (n × 3n)
  - 陷门矩阵: T = [M(g); M(f); I_n]  (3n × n)
  - 陷门关系: A·T = I_n·M(g) + M(a)·M(f) + M(b)·I_n
                     = M(g + a·f + b) = M(p) = p·I_n mod Q

签名:
  给定消息 msg, 计算 u = H(msg||salt),
  使用私钥采样短向量 (z0, z1, z2) 满足 z0 + a·z1 + b·z2 = u - e mod Q,
  输出 (salt, z1, z2). 验签时恢复 z0+e = u - a·z1 - b·z2 mod Q,
  验证 (z0+e, z1, z2) 的范数是否在界内。

与 Robin 的区别:
  - Eagle 基于 Ring-LWE, Robin 基于 NTRU
  - Eagle 的公钥更大 (需要存储 a 和 b)
  - Eagle 的签名更大 (需要 z1 和 z2, 而非仅 z1)
  - Eagle 的密钥生成更简单 (a 是公开的随机多项式)

算法参考 (论文 Section 6):
  - Algorithm 6: Eagle.KeyGen
  - Algorithm 7: Eagle.Sign
  - Algorithm 8: Eagle.Verify

参考文献: CRYPTO 2023, Section 6
"""

import math
import time
from typing import Tuple, Optional, List
from .polynomials import Polynomial, circulant_matrix
from .compact_gadget import approx_gadget
from .gaussian_sampler import sample_z
from .hash_utils import MessageHasher, msg_to_bytes, expand_seed, random_seed_32
from .parameters import EagleParams, get_eagle_params


# ============================================================================
# 辅助函数: 生成三值多项式 (用于 R^+_n)
# ============================================================================

def sample_ternary_polynomial_plus(n: int, a: int, b: int,
                                     Q: int) -> Polynomial:
    """
    从 T(n, a, b) 中均匀采样三值多项式 (用于 R^+_n 环)

    与 Robin 的采样完全相同，但环类型为 "negacyclic"。

    参数:
        n: 多项式维度 (2的幂)
        a: +1 系数的个数
        b: -1 系数的个数
        Q: 系数模数

    返回:
        采样到的三值多项式
    """
    import random
    rng = random.SystemRandom()

    coeffs = [1] * a + [-1] * b + [0] * (n - a - b)
    rng.shuffle(coeffs)

    return Polynomial(coeffs, Q, "negacyclic")


# ============================================================================
# 计算最大奇异值上界 (R^+_n 版本)
# ============================================================================

def _compute_negacyclic_dft(coeffs: List[int]) -> List[complex]:
    """
    计算 R^+_n 中多项式的 "DFT" (用于估计反循环矩阵的奇异值)

    对于 R^+_n = Z[x]/(x^n+1) 中的 a, 反循环矩阵 M(a) 的奇异值为:
      s_k = |Σ_{j=0}^{n-1} a_j · ω^{(2k+1)j}|
    其中 ω = exp(πi/n) 是本原 2n 次单位根

    这等价于: 将 a 补零到长度 2n, 做 2n 点 FFT, 取奇数索引
    优先使用 numpy.fft (O(n log n)), 否则回退到纯 Python (O(n²))。
    """
    n = len(coeffs)
    try:
        import numpy as np
        # 补零到 2n
        padded = np.zeros(2 * n, dtype=np.float64)
        padded[:n] = np.array(coeffs, dtype=np.float64)
        # 2n 点 FFT
        fft = np.fft.fft(padded)
        # 取奇数索引 (k=1,3,5,...,2n-1)
        result = [complex(fft[2 * k + 1]) for k in range(n)]
        return result
    except ImportError:
        pass
    # 纯 Python 回退
    import math as m
    result = []
    for k in range(n):
        total = complex(0, 0)
        for j in range(n):
            angle = m.pi * (2 * k + 1) * j / n
            total += coeffs[j] * complex(m.cos(angle), m.sin(angle))
        result.append(total)
    return result


def compute_s1_bound_eagle(f: Polynomial, g: Polynomial) -> float:
    """
    计算 s₁(M(f·f̄ + g·ḡ)) 的估计值 (用于 Eagle)

    对于 R^+_n 中的多项式 a, 反循环矩阵 M(a) 的奇异值是
    |NTT(a)_k| (negacyclic DFT)。
    s₁(M(a)) = max_k |NTT(a)_k|。

    DFT(f·f̄)_k = |NTT(f)_k|^2
    所以 s₁ = max_k (|NTT(f)_k|^2 + |NTT(g)_k|^2)

    参数:
        f: 私钥多项式 f
        g: 私钥多项式 g

    返回:
        s₁ 的估计值
    """
    n = f.n

    dft_f = _compute_negacyclic_dft(f.coeffs)
    dft_g = _compute_negacyclic_dft(g.coeffs)

    max_val = 0.0
    for k in range(n):
        val = abs(dft_f[k])**2 + abs(dft_g[k])**2
        if val > max_val:
            max_val = val

    return max_val


# ============================================================================
# Eagle 密钥生成 (Algorithm 6)
# ============================================================================

def eagle_keygen(params: EagleParams) -> Tuple[
    Tuple[bytes, Polynomial],  # 公钥 (seed_a, b)
    Tuple[Polynomial, Polynomial]  # 私钥 (f, g)
]:
    """
    Eagle 密钥生成算法 (论文 Algorithm 6)

    过程:
      1. 生成随机种子 seed_a, 扩展为多项式 a ∈ R_Q
      2. 采样 f,g ∈ T(n, a, b)
      3. 寻找最优自同构 k 使 s₁ 最小
      4. 计算 b = p - (a·f + g) mod Q
      5. 公钥: (seed_a, b), 私钥: (f, g)

    参数:
        params: Eagle 参数集

    返回:
        ((seed_a, b), (f, g)): 公钥和私钥对
    """
    import random
    rng = random.SystemRandom()

    n = params.n
    Q = params.Q
    p = params.p
    a_count = params.a
    b_count = params.b
    alpha = params.alpha
    K = 5

    target_norm = math.sqrt(2 * (a_count + b_count))

    # 循环直到找到合格的公私钥对
    attempt = 0
    max_keygen_attempts = 50

    while attempt < max_keygen_attempts:
        attempt += 1

        # 生成随机多项式 a
        seed_a = random_seed_32()
        a_coeffs = expand_seed(seed_a, n, Q)
        a = Polynomial(a_coeffs, Q, "negacyclic")

        # 采样候选 f,g
        f_candidates = [sample_ternary_polynomial_plus(n, a_count, b_count, Q)
                        for _ in range(K)]
        g_candidates = [sample_ternary_polynomial_plus(n, a_count, b_count, Q)
                        for _ in range(K)]

        for i in range(K):
            for j in range(K):
                fi = f_candidates[i]
                gj = g_candidates[j]

                # 寻找最优自同构 (对于demo只检查k=1)
                best_s1 = float('inf')
                best_g = gj

                k_list = [1]
                if n > 100:
                    k_list = [k for k in range(1, 20, 2) if math.gcd(k, 2*n) == 1]

                for k in k_list:
                    gk = gj.automorphism(k) if k != 1 else gj
                    s1_val = compute_s1_bound_eagle(fi, gk)
                    if s1_val < best_s1:
                        best_s1 = s1_val
                        best_g = gk

                # 检查陷门质量
                quality_threshold = alpha * target_norm
                sqrt_s1 = math.sqrt(best_s1)

                if sqrt_s1 <= quality_threshold:
                    f = fi
                    g = best_g

                    # 计算 b = p - (a·f + g) mod Q
                    af = a * f
                    af_plus_g = af + g
                    p_poly = Polynomial(
                        [p] + [0] * (n - 1), Q, "negacyclic"
                    )
                    b_poly = p_poly - af_plus_g

                    return (seed_a, b_poly), (f, g)

    raise RuntimeError(f"密钥生成失败: 超过最大尝试次数 ({max_keygen_attempts})")


# ============================================================================
# Eagle 签名 (Algorithm 7)
# ============================================================================

def eagle_sign(msg, pk: Tuple[bytes, Polynomial],
               sk: Tuple[Polynomial, Polynomial],
               params: EagleParams,
               hasher: MessageHasher = None,
               max_attempts: int = 100) -> Tuple[bytes, Polynomial, Polynomial]:
    """
    Eagle 签名算法 (论文 Algorithm 7)

    流程:
      1. a ← Expand(seed_a)
      2. salt ← 随机
      3. u ← H(msg, salt)
      4. 构造 A = [I_n | M(a) | M(b)], T = [M(g); M(f); I_n]
      5. (z0, z1, z2) ← ApproxPreSamp(A, T, u, r, s)
      6. e = u - (z0 + a·z1 + b·z2) mod Q
      7. if ||(z0+e, γ·z1, γ·z2)|| > β: 重新开始
      8. return (salt, (z1, z2))

    参数:
        msg: 要签名的消息
        pk: 公钥 (seed_a, b)
        sk: 私钥 (f, g)
        params: Eagle参数集
        hasher: 消息哈希器
        max_attempts: 最大重试次数

    返回:
        (salt, z1, z2): 签名
    """
    n = params.n
    Q = params.Q
    p = params.p
    q = params.q
    r = params.r
    s = params.s
    gamma = params.gamma
    beta = params.beta

    seed_a, b = pk
    f, g = sk
    msg_bytes = msg_to_bytes(msg)

    if hasher is None:
        hasher = MessageHasher(n, Q)

    # 从种子恢复 a
    a_coeffs = expand_seed(seed_a, n, Q)
    a = Polynomial(a_coeffs, Q, "negacyclic")

    # 构造矩阵形式
    M_a = circulant_matrix(a)
    M_b = circulant_matrix(b)
    M_g = circulant_matrix(g)
    M_f = circulant_matrix(f)

    # I_n 的矩阵形式
    I_n_mat = [[0]*n for _ in range(n)]
    for i in range(n):
        I_n_mat[i][i] = 1

    # T_rows = [M(g); M(f); I_n]  (3n × n)
    T_rows = []
    for i in range(n):
        T_rows.append(list(M_g[i]))
    for i in range(n):
        T_rows.append(list(M_f[i]))
    for i in range(n):
        T_rows.append(list(I_n_mat[i]))

    # 计算 T 的 Frobenius 范数
    frob_sq_T = sum(sum(val * val for val in row) for row in T_rows)

    attempt = 0
    while attempt < max_attempts:
        attempt += 1

        # 生成盐和目标
        salt = hasher.generate_salt()
        u_coeffs = hasher.hash_message(msg_bytes, salt)

        # === 扰动采样 ===
        s_eff_sq = s * s - r * r * frob_sq_T / n
        if s_eff_sq <= 0:
            continue
        s_eff = math.sqrt(s_eff_sq)

        # 采样 p (维度 3n)
        p_vec = [sample_z(c=0, sigma=s_eff) for _ in range(3 * n)]

        # === 计算 u' = u - A·p ===
        # A·p = p_0 + M(a)·p_1 + M(b)·p_2
        Ap = [0] * n
        # p_0 部分
        for i in range(n):
            Ap[i] = p_vec[i] % Q
        # M(a)·p_1 部分
        for i in range(n):
            total = 0
            for j in range(n):
                total += M_a[i][j] * p_vec[n + j]
            Ap[i] = (Ap[i] + total) % Q
        # M(b)·p_2 部分
        for i in range(n):
            total = 0
            for j in range(n):
                total += M_b[i][j] * p_vec[2 * n + j]
            Ap[i] = (Ap[i] + total) % Q

        u_prime = [(u_coeffs[i] - Ap[i]) % Q for i in range(n)]
        u_prime = [x if x <= Q // 2 else x - Q for x in u_prime]

        # === Gadget采样 ===
        x_prime_coeffs, e_coeffs = approx_gadget(u_prime, r, p, q, Q)

        # === 计算 x = p + T·x' ===
        # T·x' = [M(g)·x'; M(f)·x'; x']
        Tx0 = [0] * n
        Tx1 = [0] * n
        Tx2 = x_prime_coeffs[:]  # I_n·x' = x'

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
        z2_coeffs = [p_vec[2*n + i] + Tx2[i] for i in range(n)]

        # 构造多项式
        z0 = Polynomial(z0_coeffs, Q, "negacyclic")
        z1 = Polynomial(z1_coeffs, Q, "negacyclic")
        z2 = Polynomial(z2_coeffs, Q, "negacyclic")
        e = Polynomial(e_coeffs, Q, "negacyclic")

        # === 接受检查: ||(z0+e, γ·z1, γ·z2)|| ≤ β ===
        z0_plus_e = z0 + e
        norm_sq = (z0_plus_e.norm_squared() +
                   gamma * gamma * z1.norm_squared() +
                   gamma * gamma * z2.norm_squared())
        check_norm = math.sqrt(norm_sq)

        if check_norm <= beta:
            return salt, z1, z2

    raise RuntimeError(f"签名失败: 超过最大重试次数 ({max_attempts})")


# ============================================================================
# Eagle 验签 (Algorithm 8)
# ============================================================================

def eagle_verify(msg, salt: bytes, z1: Polynomial, z2: Polynomial,
                 pk: Tuple[bytes, Polynomial],
                 params: EagleParams,
                 hasher: MessageHasher = None) -> bool:
    """
    Eagle 验签算法 (论文 Algorithm 8)

    流程:
      1. a ← Expand(seed_a)
      2. u ← H(msg, salt)
      3. z' ← u - a·z1 - b·z2 mod Q
      4. 接受当且仅当 ||(z', γ·z1, γ·z2)|| ≤ β

    参数:
        msg: 原始消息
        salt: 盐
        z1, z2: 签名的两个多项式部分
        pk: 公钥 (seed_a, b)
        params: Eagle参数集
        hasher: 消息哈希器

    返回:
        True 如果签名有效
    """
    n = params.n
    Q = params.Q
    gamma = params.gamma
    beta = params.beta

    seed_a, b = pk
    msg_bytes = msg_to_bytes(msg)

    if hasher is None:
        hasher = MessageHasher(n, Q)

    # 恢复 a
    a_coeffs = expand_seed(seed_a, n, Q)
    a = Polynomial(a_coeffs, Q, "negacyclic")

    # u ← H(msg, salt)
    u_coeffs = hasher.hash_message(msg_bytes, salt)
    u = Polynomial(u_coeffs, Q, "negacyclic")

    # z' = u - a·z1 - b·z2 mod Q
    z_prime = u - a * z1 - b * z2

    # 检查 ||(z', γ·z1, γ·z2)|| ≤ β
    norm_sq = (z_prime.norm_squared() +
               gamma * gamma * z1.norm_squared() +
               gamma * gamma * z2.norm_squared())
    norm = math.sqrt(norm_sq)

    return norm <= beta


# ============================================================================
# 签名大小估算
# ============================================================================

def estimate_eagle_sig_size(z1: Polynomial, z2: Polynomial) -> int:
    """估算 Eagle 签名大小 (字节)"""
    total = 0
    for z in [z1, z2]:
        max_coeff = max(abs(c) for c in z.coeffs)
        bits_per_coeff = max(1, (max_coeff.bit_length() + 1))
        total += (z.n * bits_per_coeff + 7) // 8
    return total


# ============================================================================
# 测试代码
# ============================================================================

def _test():
    """Eagle 签名方案自测"""
    print("=" * 60)
    print("Eagle 签名方案测试 (Ring-LWE-based)")
    print("=" * 60)

    # 使用演示参数
    params = get_eagle_params("demo")
    print(f"\n参数集: {params.name}")
    print(f"  n={params.n}, Q={params.Q}, p={params.p}, q={params.q}")
    print(f"  a={params.a}, b={params.b}")
    print(f"  r={params.r:.1f}, s={params.s:.1f}")
    print(f"  γ={params.gamma:.2f}, β={params.beta:.1f}")
    print(f"  安全级别: {params.security_level}")

    # 密钥生成
    print("\n[步骤1] 密钥生成...")
    start = time.time()
    (seed_a, b), (f, g) = eagle_keygen(params)
    keygen_time = time.time() - start
    print(f"  密钥生成时间: {keygen_time:.2f}s")
    print(f"  seed_a (hex): {seed_a.hex()[:20]}...")
    print(f"  b (前8个系数): {b.coeffs[:8]}")
    print(f"  私钥 f (前8个系数): {f.coeffs[:8]}")
    print(f"  私钥 g (前8个系数): {g.coeffs[:8]}")

    # 验证陷门关系: a·f + g + b = p mod Q
    af = Polynomial(expand_seed(seed_a, params.n, params.Q),
                    params.Q, "negacyclic") * f
    af_plus_g_plus_b = af + g + b
    p_theory = params.p
    trapdoor_ok = (af_plus_g_plus_b.coeffs[0] % params.Q == p_theory % params.Q)
    print(f"  陷门关系验证 (af + g + b = p): {trapdoor_ok}")

    # 签名
    print("\n[步骤2] 签名...")
    msg = "Hello, Eagle signature scheme!"
    hasher = MessageHasher(params.n, params.Q)

    try:
        start = time.time()
        salt, z1, z2 = eagle_sign(msg, (seed_a, b), (f, g), params, hasher)
        sign_time = time.time() - start
        print(f"  签名时间: {sign_time:.2f}s")
        print(f"  盐 (hex): {salt.hex()[:20]}...")
        print(f"  z1 (前4个系数): {z1.coeffs[:4]}")
        print(f"  z2 (前4个系数): {z2.coeffs[:4]}")
        print(f"  ||z1||^2 = {z1.norm_squared()}")
        print(f"  ||z2||^2 = {z2.norm_squared()}")
        print(f"  估算签名大小: {estimate_eagle_sig_size(z1, z2)} 字节")

        # 验签
        print("\n[步骤3] 验签...")
        start = time.time()
        valid = eagle_verify(msg, salt, z1, z2, (seed_a, b), params, hasher)
        verify_time = time.time() - start
        print(f"  验签时间: {verify_time*1000:.2f}ms")
        print(f"  验签结果: {'有效 ✓' if valid else '无效 ✗'}")

        # 安全性测试
        print("\n[步骤4] 安全性测试...")
        tests = [
            ("错误消息", "Wrong", salt, z1, z2, (seed_a, b)),
            ("错误公钥b",
             msg, salt, z1, z2,
             (seed_a, Polynomial([(c + 1) % params.Q for c in b.coeffs],
                                 params.Q, "negacyclic"))),
        ]
        for test_name, *args in tests:
            result = eagle_verify(*args, params, hasher)
            print(f"  {test_name}: {'有效 (BUG!)' if result else '无效 ✓'}")

        # 多次签名
        print("\n[步骤5] 多次签名统计...")
        successes = 0
        for i in range(10):
            try:
                s, zz1, zz2 = eagle_sign(f"Test {i}", (seed_a, b), (f, g), params, hasher)
                if eagle_verify(f"Test {i}", s, zz1, zz2, (seed_a, b), params, hasher):
                    successes += 1
            except RuntimeError:
                pass
        print(f"  签名成功率: {successes}/10")

    except RuntimeError as e:
        print(f"  签名失败: {e}")

    print("\n所有测试完成！")


if __name__ == "__main__":
    _test()
