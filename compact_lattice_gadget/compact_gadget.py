# -*- coding: utf-8 -*-
"""
紧凑Gadget与半随机采样器模块

本模块实现论文的核心创新——紧凑Gadget框架和半随机采样器。

核心思想:
  - 传统Gadget: G = I_n ⊗ g^t (n行nk列), g = (1, b, ..., b^{k-1})
  - 本工作: P = p·I_n (n行n列，紧凑!), Q = q·I_n, 满足 PQ = Q·I_n (即 p·q = Q)

半随机采样器 (Semi-random Sampler) 分为两步:
  1. 确定性错误解码 (Deterministic Error Decoding):
     给定目标 u' ∈ Z^n_Q, 计算 e = u' mod p, c = (u' - e) / p
     使得 u' - e = P·c ∈ Λ(P)

  2. 随机原像采样 (Random Preimage Sampling):
     从陪集分布 D_{Λ(Q) + c, r} 中采样 x'
     满足 P·x' = u' - e mod Q

关键性质 (Lemma 6):
  对于均匀随机的目标 u', (x', e) 的分布可以在不知道陷门的情况下模拟。
  这保证了签名方案的安全性。

参考文献: CRYPTO 2023, Section 4
"""

from typing import Tuple, List
from .gaussian_sampler import sample_polynomial_gadget


# ============================================================================
# 格解码器 (Lattice Decoder)
# ============================================================================

def lattice_decoder(u_coeffs: List[int], p: int, Q: int) -> Tuple[List[int], List[int]]:
    """
    确定性格解码: 将目标向量解码到格 p·Z^n 上

    算法 (Section 4.1, Algorithm 1 第一步):
      输入: u' ∈ Z^n_Q
      输出: (c ∈ Z^n_Q, e ∈ Z^n_p) 满足 u' = p·c + e mod Q

    具体步骤 (逐坐标运算):
      e_i = u_i mod p     (对称取模, 使得 e_i ∈ [-p/2, p/2))
      c_i = (u_i - e_i) / p   (在整数上进行精确除法)

    参数:
        u_coeffs: 目标向量 u' 的 n 个系数 (模 Q)
        p: 小模数 (gadget参数)
        Q: 大模数 (满足 Q = p·q)

    返回:
        (c_coeffs, e_coeffs): 解码结果
          - c_coeffs[i]: c 的第 i 个系数
          - e_coeffs[i]: 错误 e 的第 i 个系数 (在 [-p/2, p/2) 范围内)

    注意:
      由于 Q = p·q, 等式 u' = p·c + e 在整数上也成立 (不仅仅是模 Q)。
    """

    c_coeffs = []
    e_coeffs = []

    for u_i in u_coeffs:
        # 对称取模 [-(p-1)/2, p/2]
        e_i = u_i % p
        if e_i > p // 2:
            e_i -= p

        # 精确除法
        c_i = (u_i - e_i) // p
        c_coeffs.append(c_i)
        e_coeffs.append(e_i)

    return c_coeffs, e_coeffs


# ============================================================================
# 半随机 Gadget 采样器 (Semi-random Sampler)
# ============================================================================

def approx_gadget(u_prime_coeffs: List[int], r: float, p: int, q: int,
                  Q: int) -> Tuple[List[int], List[int]]:
    """
    近似Gadget采样器 (论文 Algorithm 1: ApproxGadget)

    这是半随机采样器的完整实现。

    输入: 目标 u' ∈ Z^n_Q, 参数 r ≥ η_ε(Λ(Q))
    输出: (x'_coeffs, e_coeffs) 满足 P·x' = u' - e mod Q

    算法步骤:
      1. (c, e) ← LatticeDecoder(u', p)
         确定性计算错误 e 和中间值 c
      2. x' ← D_{qZ^n + c, r}
         从陪集高斯分布中采样原像

    对于简单实例化 P = p·I_n, Q = q·I_n:
      - Λ(P) = p·Z^n, 解码使用 mod p 操作
      - E(P) = Z^n_p, 错误均匀分布在 Z^n_p 上
      - Λ(Q) = q·Z^n, 采样使用 D_{qZ + c_i, r}

    参数:
        u_prime_coeffs: 目标向量 u' 的 n 个系数 (在 Z_Q 中)
        r: Gadget高斯宽度 (满足 r ≥ q·η_ε(Z^n))
        p: 小模数
        q: 中等模数 (Q = p·q)
        Q: 大模数

    返回:
        (x_prime_coeffs, e_coeffs):
          - x_prime_coeffs: 原像 x' 的系数 (在 Z^n 上)
          - e_coeffs: 错误 e 的系数 (在 [-p/2, p/2) 范围内, 由LatticeDecoder确定)
    """
    n = len(u_prime_coeffs)

    # === 步骤1: 确定性错误解码 ===
    c_coeffs, e_coeffs = lattice_decoder(u_prime_coeffs, p, Q)

    # === 步骤2: 随机原像采样 ===
    # 从 D_{Λ(Q) + c, r} = D_{qZ^n + c, r} 中采样
    # 等价于逐坐标从 D_{qZ + c_i, r} 中采样
    x_prime_coeffs = sample_polynomial_gadget(n, r, q, c_coeffs)

    return x_prime_coeffs, e_coeffs


# ============================================================================
# 统计信息
# ============================================================================

def compute_error_stats(e_coeffs: List[int], p: int) -> dict:
    """
    计算错误向量的统计信息

    根据论文 Section 4.2, 对于均匀随机目标:
      - e 在 Z^n_p 上均匀分布
      - e 的系数的标准差 ≈ √(p^2-1)/12

    参数:
        e_coeffs: 错误向量系数
        p: 小模数

    返回:
        包含统计信息的字典
    """
    import math
    n = len(e_coeffs)
    mean = sum(e_coeffs) / n
    variance = sum((x - mean) ** 2 for x in e_coeffs) / n
    std = math.sqrt(variance)
    expected_std = math.sqrt((p * p - 1) / 12)

    norm_sq = sum(x * x for x in e_coeffs)

    return {
        'n': n,
        'p': p,
        'mean': mean,
        'std': std,
        'expected_std': expected_std,
        'norm': math.sqrt(norm_sq),
        'norm_sq': norm_sq
    }


# ============================================================================
# 测试代码
# ============================================================================

def _test():
    """紧凑Gadget模块自测"""
    print("=" * 60)
    print("紧凑Gadget模块测试")
    print("=" * 60)

    # 使用小参数测试
    p = 7
    q = 8
    Q = p * q  # 56
    n = 16

    print(f"\n参数: p={p}, q={q}, Q={Q}, n={n}")

    # 测试1: 格解码器
    print("\n[测试1] 格解码器 (Lattice Decoder)")
    import random
    rng = random.Random(42)
    # 生成随机目标
    u_coeffs = [rng.randint(-Q // 2, Q // 2 - 1) for _ in range(n)]
    c_coeffs, e_coeffs = lattice_decoder(u_coeffs, p, Q)

    # 验证: u = p*c + e
    reconstruction_ok = all(
        (p * c_coeffs[i] + e_coeffs[i]) == u_coeffs[i]
        for i in range(n)
    )
    print(f"  解码重构验证: {reconstruction_ok}")

    # 验证 e 在范围内
    e_in_range = all(-p // 2 <= e_i <= p // 2 for e_i in e_coeffs)
    print(f"  错误范围验证 (|e_i| <= {p//2}): {e_in_range}")

    print(f"  u 的前5个系数: {u_coeffs[:5]}")
    print(f"  e 的前5个系数: {e_coeffs[:5]}")
    print(f"  c 的前5个系数: {c_coeffs[:5]}")

    # 测试2: 半随机采样器
    print("\n[测试2] 半随机采样器 (ApproxGadget)")
    r = 10.0  # 高斯宽度
    x_prime_coeffs, e_coeffs = approx_gadget(u_coeffs, r, p, q, Q)

    # 验证: p * x' ≡ u' - e (mod Q)
    verify_ok = all(
        (p * x_prime_coeffs[i]) % Q == (u_coeffs[i] - e_coeffs[i]) % Q
        for i in range(n)
    )
    print(f"  采样正确性验证 (Px' = u' - e mod Q): {verify_ok}")

    # 检查 x' 是否在正确的陪集中
    # x' = q*k + c, 所以 (x' - c) % q = 0
    x_in_coset = all(
        (x_prime_coeffs[i] - c_coeffs[i]) % q == 0
        for i in range(n)
    )
    print(f"  x' 陪集验证: {x_in_coset}")

    print(f"  x' 的前5个系数: {x_prime_coeffs[:5]}")
    print(f"  ||x'||^2 ≈ {sum(c*c for c in x_prime_coeffs)}")

    # 测试3: 统计性质
    print("\n[测试3] 错误分布统计")
    # 生成多个随机目标,观察错误分布
    num_tests = 1000
    all_errors = []
    for _ in range(num_tests):
        u = [rng.randint(0, Q - 1) - Q // 2 for _ in range(n)]
        _, e = lattice_decoder(u, p, Q)
        all_errors.extend(e)

    import math
    mean = sum(all_errors) / len(all_errors)
    var = sum((x - mean) ** 2 for x in all_errors) / len(all_errors)
    std = math.sqrt(var)
    expected_std = math.sqrt((p * p - 1) / 12)

    print(f"  样本数: {len(all_errors)}")
    print(f"  均值: {mean:.4f} (理论: 0)")
    print(f"  标准差: {std:.4f} (理论: {expected_std:.4f})")

    # 测试4: 不同的 r 对 x' 大小的影响
    print("\n[测试4] 不同高斯宽度r对原像大小的影响")
    for r_test in [5.0, 10.0, 20.0, 50.0]:
        u_test = [rng.randint(0, Q - 1) - Q // 2 for _ in range(n)]
        xp, _ = approx_gadget(u_test, r_test, p, q, Q)
        norm_sq = sum(c * c for c in xp)
        print(f"  r={r_test:.1f}: ||x'||^2 ≈ {norm_sq}")

    print("\n所有测试完成！")


if __name__ == "__main__":
    _test()
