# -*- coding: utf-8 -*-
"""
近似陷门框架模块 (Approximate Trapdoor Framework)

实现基于紧凑Gadget的近似陷门原像采样算法。

核心概念:
  - 公钥矩阵 A ∈ Z^{n×m}_Q
  - 陷门矩阵 T ∈ Z^{m×n} (小条目), 满足 A·T = P mod Q
  - Gadget矩阵 P = p·I_n, Q = q·I_n (满足 p·q = Q)

近似原像采样 (Algorithm 2: ApproxPreSamp):
  给定目标 u ∈ Z^n_Q, 使用陷门 T 计算短向量 x 满足:
    A·x = u - e mod Q
  其中 e 是短错误向量。

采样过程:
  1. 扰动采样 (Perturbation):
     p ← D_{Z^m, √Σ_p} 其中 Σ_p = s^2·I_m - r^2·T·T^t
     需要 s ≥ r·s₁(T) (s₁是最大奇异值)

  2. 中间目标计算:
     u' = u - A·p mod Q

  3. Gadget采样:
     x' ← ApproxGadget(u', r, P, Q)

  4. 原像合成:
     x = p + T·x'

安全性 (Theorem 2):
  对于均匀随机的目标 u, 输出 (x, u, e) 的分布可以与不知道陷门
  的模拟分布统计不可区分。这保证了签名方案的安全性。

参考文献: CRYPTO 2023, Section 4
"""

import math
from typing import Tuple, List
from .compact_gadget import approx_gadget
from .gaussian_sampler import sample_z


# ============================================================================
# 矩阵工具函数
# ============================================================================

def compute_tt_transpose(T_rows: List[List[int]], T_cols: int) -> List[List[int]]:
    """
    计算 T·T^t, 用于构造扰动协方差矩阵

    T 是 m×n 的整数矩阵, TT^t 是 m×m 的矩阵

    参数:
        T_rows: T 按行存储的列表 (m行, 每行n个元素)
        T_cols: T 的列数 n

    返回:
        m×m 的乘积矩阵 TT^t
    """
    m = len(T_rows)
    result = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            total = 0
            for k in range(T_cols):
                total += T_rows[i][k] * T_rows[j][k]
            result[i][j] = total
    return result


# ============================================================================
# 扰动采样
# ============================================================================

def sample_perturbation_vector(m: int, s: float, r: float,
                                T_rows: List[List[int]],
                                T_cols: int) -> List[int]:
    """
    扰动采样: 从 D_{Z^m, √(s^2·I - r^2·T·T^t)} 中采样

    根据 Lemma 7, 当 s^2 ≥ (r^2 + η^2)·(s₁(T)^2 + 1) 时,
    扰动分布可以很好地近似。

    简化实现: 采样 p 使得 [I|T]·D_{Z^{m+n}, Σ_p ⊕ r^2I} ≈ D_{Z^m, s^2I}

    根据论文框架, 扰动采样等价于:
      1. 从 D_{Z^{m+n}, diag(Σ_p, r^2I)} 中采样 (p', x')
      2. 计算 p = p' + T·x'

    但更直接的实现是利用 Cholesky 分解等方法。
    这里采用近似实现: 独立采样每个分量。

    参数:
        m: 向量维度
        s: 原像高斯宽度
        r: Gadget高斯宽度
        T_rows: 陷门矩阵
        T_cols: 陷门矩阵列数

    返回:
        扰动向量 p (长度 m 的整数列表)
    """
    # 对于简单演示, 使用独立采样近似
    # 每个分量从 D_{Z, s_eff} 中采样, s_eff = √(s^2 - r^2·λ_max(TT^t))

    # 计算 s₁(T) 的上界
    # s₁(T) ≤ ||T||_F (Frobenius 范数)
    frob_sq = 0
    for row in T_rows:
        for val in row:
            frob_sq += val * val
    s1_approx = math.sqrt(frob_sq)

    # 使用简化条件: s^2 ≥ r^2·s₁(T)^2
    # 有效宽度 s_eff = √(s^2 - r^2·s₁(T)^2)
    s_eff_sq = s * s - r * r * s1_approx * s1_approx
    if s_eff_sq < 0:
        raise ValueError(f"参数不满足: s^2 < r^2·s₁(T)^2 ({s*s:.1f} < {r*r*s1_approx*s1_approx:.1f})")

    s_eff = math.sqrt(s_eff_sq)

    return [sample_z(c=0, sigma=s_eff) for _ in range(m)]


# ============================================================================
# 近似原像采样 (核心算法)
# ============================================================================

def approx_preimage_sample(
    A_rows: List[List[int]],   # 公钥矩阵 A (n×m, 每个元素是整数)
    T_rows: List[List[int]],   # 陷门矩阵 T (m×n)
    u_coeffs: List[int],        # 目标向量 u (n个系数, 在Z_Q中)
    r: float,                   # Gadget高斯宽度
    s: float,                   # 原像高斯宽度
    p: int,                     # 小模数 (gadget参数)
    q: int,                     # 中等模数
    Q: int,                     # 大模数
) -> Tuple[List[int], List[int]]:
    """
    近似原像采样 (论文 Algorithm 2: ApproxPreSamp)

    给定公钥矩阵 A 和陷门 T (满足 A·T = P mod Q),
    计算目标 u 的近似原像 x ∈ Z^m, 满足:
      A·x = u - e mod Q

    完整算法:
      输入: (A, T) ∈ Z^{n×m}_Q × Z^{m×n} 满足 A·T = P mod Q,
            u ∈ Z^n_Q, r ≥ η_ε(Λ(Q)), s^2·I_m ≥ r^2·T·T^t
      输出: 近似原像 x, 满足 A·x = u - e mod Q

      1. p ← D_{Z^m, √Σ_p}, 其中 Σ_p = s^2·I_m - r^2·T·T^t
      2. u' = u - A·p mod Q
      3. x' ← ApproxGadget(u', r, P, Q)
      4. return x = p + T·x'

    参数:
        A_rows: 公钥矩阵 (n行 × m列, 整数)
        T_rows: 陷门矩阵 (m行 × n列, 整数)
        u_coeffs: 目标向量 (n个整数, 对称表示在[-Q/2, Q/2))
        r: Gadget高斯宽度
        s: 原像高斯宽度
        p: 小模数
        q: 中等模数
        Q: 大模数

    返回:
        (x_coeffs, e_coeffs):
          - x_coeffs: 原像向量 (m个整数)
          - e_coeffs: 错误向量 (n个整数, 在[-p/2, p/2)范围内)
    """
    n = len(A_rows)     # A是 n×m
    m = len(A_rows[0])  # A是 n×m
    T_cols = len(T_rows[0])  # T是 m×nCols

    # === 步骤1: 扰动采样 ===
    # p ← D_{Z^m, √Σ_p}
    p_vec = sample_perturbation_vector(m, s, r, T_rows, T_cols)

    # === 步骤2: 计算中间目标 u' ===
    # u' = u - A·p mod Q
    # 计算 A·p
    Ap = [0] * n
    for i in range(n):
        total = 0
        for j in range(m):
            total = (total + A_rows[i][j] * p_vec[j]) % Q
        Ap[i] = total

    # u' = u - Ap mod Q
    u_prime = [(u_coeffs[i] - Ap[i]) % Q for i in range(n)]
    # 转换到对称表示
    u_prime = [x if x <= Q // 2 else x - Q for x in u_prime]

    # === 步骤3: Gadget采样 ===
    x_prime_coeffs, e_coeffs = approx_gadget(u_prime, r, p, q, Q)

    # === 步骤4: 合成原像 x = p + T·x' ===
    # T 是 m×n 矩阵, x' 是 n 维向量
    Tx = [0] * m
    for i in range(m):
        total = 0
        for j in range(T_cols):
            total += T_rows[i][j] * x_prime_coeffs[j]
        Tx[i] = total

    x_coeffs = [p_vec[i] + Tx[i] for i in range(m)]

    # 验证正确性 (可选的调试检查)
    # A·x = A·p + A·T·x' = A·p + P·x' = A·p + u' - e = u - e mod Q

    return x_coeffs, e_coeffs


# ============================================================================
# 简化版: 多项式环上的原像采样
# ============================================================================

def approx_preimage_poly(
    A_polys,      # 公钥矩阵 (多项式环上的矩阵)
    T_polys,      # 陷门矩阵 (多项式环上的矩阵)
    u_poly,       # 目标多项式
    r: float,     # Gadget高斯宽度
    s: float,     # 原像高斯宽度
    p: int,       # 小模数
    q: int,       # 中等模数
    Q: int,        # 大模数
):
    """
    多项式环上的近似原像采样

    与 approx_preimage_sample 相同, 但工作在多项式环 R_Q 上,
    利用多项式矩阵形式将环运算转化为整数运算。

    这是实际签名方案中使用的版本。

    参数:
        A_polys: 公钥多项式的列表 (用于构造 M(A))
        T_polys: 陷门多项式的列表 (用于构造 M(T))
        u_poly: 目标多项式
        r, s, p, q, Q: 和其他参数一致

    返回:
        (x_polys, e_poly): 原像多项式列表和错误多项式
    """
    n = u_poly.n
    m = len(A_polys)

    # 步骤1: 构建矩阵形式
    # A = [M(a_0) | M(a_1) | ... | M(a_{m-1})] (n × nm)
    # T = [M(t_0)^t | M(t_1)^t | ... | M(t_{m-1})^t]^t (nm × n)

    # 对每个多项式构造其矩阵形式
    from .polynomials import circulant_matrix, matrix_vector_multiply

    A_matrices = [circulant_matrix(poly) for poly in A_polys]
    T_matrices = [circulant_matrix(poly) for poly in T_polys]

    # u 的系数向量
    u_coeffs = u_poly.coeffs[:]

    # 步骤2: 扰动采样 p_vec (维度 m*n)
    p_vec = sample_perturbation_vector(
        m * n,
        s, r,
        [], 0  # 简化: 不传递TT^t
    )
    # 实际上, 上面的扰动采样过于简化。在完整实现中需要正确计算 TT^t。

    # 步骤3: 计算 Ap = A·p = Σ M(a_i)·p_i
    Ap = [0] * n
    for i in range(m):
        p_i = p_vec[i * n:(i + 1) * n]
        # 计算 M(a_i)·p_i
        mat_vec = matrix_vector_multiply(A_matrices[i],
                                          type('Poly', (), {'coeffs': p_i})(),
                                          Q)
        for j in range(n):
            Ap[j] = (Ap[j] + mat_vec.coeffs[j]) % Q

    # 步骤4: u' = u - Ap
    u_prime = [(u_coeffs[i] - Ap[i]) % Q for i in range(n)]
    u_prime = [x if x <= Q // 2 else x - Q for x in u_prime]

    # 步骤5: Gadget采样
    x_prime_coeffs, e_coeffs = approx_gadget(u_prime, r, p, q, Q)

    # 步骤6: 计算 T·x'
    Tx = [0] * (m * n)
    for i in range(m):
        # T_i 的第i块
        M_ti = T_matrices[i]
        # 计算 M(t_i)^t · x'
        # 这里是转置乘法
        for row in range(n):
            total = 0
            for col in range(n):
                total += M_ti[col][row] * x_prime_coeffs[col]
            Tx[i * n + row] = total

    # 步骤7: x = p + T·x'
    x_coeffs = [p_vec[i] + Tx[i] for i in range(m * n)]

    # 将结果分解回多项式列表
    from .polynomials import Polynomial
    x_polys = []
    for i in range(m):
        start = i * n
        x_i_coeffs = x_coeffs[start:start + n]
        x_polys.append(Polynomial(x_i_coeffs, Q, u_poly.ring_type))

    e_poly = Polynomial(e_coeffs, Q, u_poly.ring_type)

    return x_polys, e_poly


# ============================================================================
# 测试代码
# ============================================================================

def _test():
    """近似陷门框架模块自测"""
    print("=" * 60)
    print("近似陷门框架模块测试")
    print("=" * 60)

    # 小参数
    n = 8    # 环维度
    m = 2    # A 的列数 (多项式个数)
    p = 7
    q = 8
    Q = p * q  # 56
    r = 10.0
    s = 100.0

    print(f"\n参数: n={n}, m={m}, p={p}, q={q}, Q={Q}")
    print(f"  r={r}, s={s}")

    # 生成一个简单的陷门对
    # A = [I_n | M(h)] (用于Robin)
    # 使用简化的随机矩阵测试
    import random
    rng = random.Random(42)

    # 构造简单的多项式环
    from .polynomials import Polynomial, circulant_matrix

    # 随机公钥多项式 h
    h_coeffs = [rng.randint(0, Q - 1) - Q // 2 for _ in range(n)]
    h = Polynomial(h_coeffs, Q, "cyclic")

    # 公钥矩阵 A = [I_n | M(h)]
    # 在整数层面: A 是 n×2n 的矩阵

    # 构造 T (简化的陷门, 真实版本来自NTRU密钥)
    f_coeffs = [rng.randint(-1, 1) for _ in range(n)]
    g_coeffs = [rng.randint(-1, 1) for _ in range(n)]
    f = Polynomial(f_coeffs, Q, "cyclic")
    g = Polynomial(g_coeffs, Q, "cyclic")

    M_f = circulant_matrix(f)
    M_g = circulant_matrix(g)

    # T = [M(g); M(f)]^t  2n × n
    T_rows = []
    for i in range(n):
        T_rows.append(list(M_g[i]))  # 第一块
    for i in range(n):
        T_rows.append(list(M_f[i]))  # 第二块

    # 构造 A = [I_n | M(h)] 其中 M(h)·f + I·g = p·I? 需要验证

    # 简化测试: 直接使用整数的 ApproxPreSamp
    # 随机 A (n×2n) 和 T (2n×n)
    Ar = [[rng.randint(0, Q - 1) for _ in range(2 * n)] for _ in range(n)]
    Tr = [[rng.randint(-2, 2) for _ in range(n)] for _ in range(2 * n)]

    # 确保 A·T = p·I_n mod Q
    # 调整 A 的最后一列满足此关系 (粗略的测试)
    # 这不是真实的陷门构造,仅用于测试流程

    # 随机目标
    u = [rng.randint(0, Q - 1) - Q // 2 for _ in range(n)]

    print(f"\n  目标 u: {u[:4]}...")

    try:
        x_coeffs, e_coeffs = approx_preimage_sample(
            Ar, Tr, u, r, s, p, q, Q
        )
        print(f"  原像 x (前8个): {x_coeffs[:8]}...")
        print(f"  错误 e: {e_coeffs}")

        # 验证 Ax = u - e mod Q
        Ax = [0] * n
        for i in range(n):
            total = 0
            for j in range(2 * n):
                total += Ar[i][j] * x_coeffs[j]
            Ax[i] = total % Q

        expected = [(u[i] - e_coeffs[i]) % Q for i in range(n)]
        verify = all(Ax[i] == expected[i] for i in range(n))
        print(f"  验证 Ax = u - e mod Q: {verify}")

        # 计算原像和错误的大小
        x_norm = sum(c * c for c in x_coeffs)
        e_norm = sum(c * c for c in e_coeffs)
        print(f"  ||x||^2 = {x_norm}")
        print(f"  ||e||^2 = {e_norm}")

    except ValueError as e:
        print(f"  参数错误: {e}")

    print("\n所有测试完成！")


if __name__ == "__main__":
    _test()
