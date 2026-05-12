# -*- coding: utf-8 -*-
"""
多项式环运算模块

支持两种多项式环：
  - R^-_n = Z[x]/(x^n - 1)  : 用于 Robin (NTRU-based) 签名方案，n为素数
  - R^+_n = Z[x]/(x^n + 1)  : 用于 Eagle (Ring-LWE-based) 签名方案，n为2的幂

提供多项式在模Q下的基本运算：加法、减法、乘法、求逆，
以及矩阵形式 M(a)、自同构 σ_k(a)、共轭 ā 等操作。

参考文献: CRYPTO 2023, "Compact Lattice Gadget and Its Applications to Hash-and-Sign Signatures"
"""

import random
from typing import List, Tuple, Optional

# ============================================================================
# 辅助函数
# ============================================================================

def mod_symmetric(x: int, q: int) -> int:
    """
    对称取模: 将整数x映射到 [-q/2, q/2) 范围内
    这是格密码中常用的取模方式
    """
    r = x % q
    if r > q // 2:
        r -= q
    return r


def mod_positive(x: int, q: int) -> int:
    """正向取模: 将整数x映射到 [0, q) 范围内"""
    return x % q


# ============================================================================
# 多项式类
# ============================================================================

class Polynomial:
    """
    多项式类，表示为系数向量的形式 a_0 + a_1·x + ... + a_{n-1}·x^{n-1}

    参数:
        coeffs: 系数列表 [a_0, a_1, ..., a_{n-1}]
        modulus: 系数模数 Q
        ring_type: 环类型
            - "negacyclic": R^+_n = Z[x]/(x^n + 1) (默认, 用于Eagle)
            - "cyclic":     R^-_n = Z[x]/(x^n - 1) (用于Robin)
    """

    def __init__(self, coeffs: List[int], modulus: int,
                 ring_type: str = "negacyclic"):
        self.n = len(coeffs)
        self.modulus = modulus  # 通常记作 Q
        self.ring_type = ring_type

        # 将系数存储在对称表示中 [-Q/2, Q/2)
        self.coeffs = [mod_symmetric(c, modulus) for c in coeffs]

    # ----- 基本操作 -----

    def __add__(self, other: 'Polynomial') -> 'Polynomial':
        """多项式加法 (逐系数)"""
        if self.n != other.n or self.modulus != other.modulus:
            raise ValueError("多项式维度或模数不匹配")
        return Polynomial(
            [(a + b) % self.modulus for (a, b) in zip(self.coeffs, other.coeffs)],
            self.modulus, self.ring_type
        )

    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        """多项式减法 (逐系数)"""
        if self.n != other.n or self.modulus != other.modulus:
            raise ValueError("多项式维度或模数不匹配")
        return Polynomial(
            [(a - b) % self.modulus for (a, b) in zip(self.coeffs, other.coeffs)],
            self.modulus, self.ring_type
        )

    def __neg__(self) -> 'Polynomial':
        """多项式取负"""
        return Polynomial(
            [(-c) % self.modulus for c in self.coeffs],
            self.modulus, self.ring_type
        )

    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        """
        多项式乘法 (模 x^n ± 1 和模 Q)

        使用直接卷积方法。对于大规模参数，可用NTT加速。
        这里先用朴素实现以保证正确性。
        """
        if self.n != other.n or self.modulus != other.modulus:
            raise ValueError("多项式维度或模数不匹配")
        n = self.n
        q = self.modulus
        result = [0] * n

        for i in range(n):
            if self.coeffs[i] == 0:
                continue
            for j in range(n):
                if other.coeffs[j] == 0:
                    continue
                k = i + j
                if k >= n:
                    k -= n
                    # 根据环类型决定是加还是减
                    if self.ring_type == "negacyclic":
                        # R^+: x^n ≡ -1, 所以 x^k ≡ -x^{k-n}
                        result[k] = (result[k] - self.coeffs[i] * other.coeffs[j]) % q
                    else:
                        # R^-: x^n ≡ 1, 所以 x^k ≡ x^{k-n}
                        result[k] = (result[k] + self.coeffs[i] * other.coeffs[j]) % q
                else:
                    result[k] = (result[k] + self.coeffs[i] * other.coeffs[j]) % q

        return Polynomial(result, q, self.ring_type)

    def __eq__(self, other: 'Polynomial') -> bool:
        if not isinstance(other, Polynomial):
            return False
        return (self.n == other.n and
                self.modulus == other.modulus and
                self.ring_type == other.ring_type and
                self.coeffs == other.coeffs)

    def __repr__(self) -> str:
        return f"Polynomial(n={self.n}, Q={self.modulus}, type={self.ring_type})"

    def copy(self) -> 'Polynomial':
        """深拷贝"""
        return Polynomial(self.coeffs[:], self.modulus, self.ring_type)

    # ----- 标量运算 -----

    def scale(self, c: int) -> 'Polynomial':
        """标量乘法: c * poly mod Q"""
        return Polynomial(
            [(c * a) % self.modulus for a in self.coeffs],
            self.modulus, self.ring_type
        )

    def add_constant(self, c: int) -> 'Polynomial':
        """常数项加上常量: poly + c"""
        new_coeffs = self.coeffs[:]
        new_coeffs[0] = (new_coeffs[0] + c) % self.modulus
        return Polynomial(new_coeffs, self.modulus, self.ring_type)

    # ----- 范数计算 -----

    def norm_squared(self) -> int:
        """计算 ||a||^2 = Σ a_i^2 (系数的平方和)"""
        return sum(c * c for c in self.coeffs)

    # ----- 多项式特殊运算 -----

    def conjugate(self) -> 'Polynomial':
        """
        共轭 ā = a(x^{-1})
        对于 R^-_n: ā = a_0 + Σ_{i=1}^{n-1} a_{n-i}·x^i
        对于 R^+_n: ā = a_0 - Σ_{i=1}^{n-1} a_{n-i}·x^i
        """
        new_coeffs = [self.coeffs[0]]
        for i in range(1, self.n):
            val = self.coeffs[self.n - i]
            if self.ring_type == "negacyclic":
                val = -val
            new_coeffs.append(val)
        return Polynomial(new_coeffs, self.modulus, self.ring_type)

    def automorphism(self, k: int) -> 'Polynomial':
        """
        自同构 σ_k(a) = a(x^k) mod (x^n ± 1)
        k 必须与 2n 互素 (对于 R^+) 或与 n 互素 (对于 R^-)
        """
        new_coeffs = [0] * self.n
        for i in range(self.n):
            new_idx = (i * k) % self.n
            val = self.coeffs[i]
            # 对于 R^+，需要考虑符号翻转
            if self.ring_type == "negacyclic":
                # 当 (i * k) // n 为奇数时取负
                if ((i * k) // self.n) % 2 == 1:
                    val = -val
            new_coeffs[new_idx] = mod_symmetric(val, self.modulus)
        return Polynomial(new_coeffs, self.modulus, self.ring_type)

    # ----- 多项式求逆 (模 Q) -----

    def invert(self) -> Optional['Polynomial']:
        """
        计算多项式在模 x^n ± 1 和模 Q 下的乘法逆元
        使用扩展欧几里得算法 (仅适用于 R^+_n 和 Q 为素数的幂的情况)

        返回 None 如果逆元不存在
        """
        # 目前使用简单实现：转换为幂级数求逆
        # 对于实际应用，应使用 NTT-based 方法
        return self._invert_by_newton()

    # ===== 多项式求逆 =====

    def invert(self) -> Optional['Polynomial']:
        """
        计算多项式在 Z[x]/(x^n ± 1) 模 Q 下的乘法逆元

        根据 Q 的类型选择不同方法:
          - Q 为素数 → 使用扩展欧几里得算法 (EEA)
          - Q 为 2 的幂 → 使用 Hensel 提升 (从模2开始逐步提升)
        """
        if self.n <= 128 and self.modulus > 2:
            # 检查 Q 是否为 2 的幂 (论文中的典型情况)
            q_val = self.modulus
            if q_val & (q_val - 1) == 0 and q_val > 2:
                # Q 是 2 的幂, 使用 Hensel 提升
                return self._invert_by_hensel()
            else:
                return self._invert_by_gcd()
        else:
            return self._invert_by_newton()

    def _invert_by_hensel(self) -> Optional['Polynomial']:
        """
        使用 Hensel 提升法求多项式在 Z_{2^k}[x]/(x^n ± 1) 中的逆元

        算法:
          1. 在 Z_2[x]/(x^n ± 1) 中用 EEA 求逆 g_1
          2. 逐步提升: g_{k+1} = g_k * (2 - f * g_k) mod 2^{k+1}
        其中乘法使用普通截断(非环乘法), 最后进行环修正。

        这是 NTRU 实现中的标准方法。
        """
        q = self.modulus
        n = self.n
        ring = self.ring_type

        # 步骤1: 在 Z_2 中求逆
        try:
            f_mod2 = Polynomial([c % 2 for c in self.coeffs], 2, ring)
            g_mod2 = f_mod2._invert_by_gcd()
            if g_mod2 is None:
                return None
        except Exception:
            return None

        # 步骤2: Hensel 提升
        # g_{k+1} = g_k * (2 - f * g_k) mod 2^{k+1}
        # 在普通截断多项式环中计算 (非环乘法), 然后做环修正

        g_coeffs = g_mod2.coeffs[:]
        current_q = 2  # 当前模数

        while current_q < q:
            next_q = min(current_q * current_q, q)

            # 构造中间多项式 (模 next_q)
            f_cur = Polynomial(
                [c % next_q for c in self.coeffs], next_q, ring
            )
            g_cur = Polynomial(
                [c % next_q for c in g_coeffs], next_q, ring
            )
            two_poly = Polynomial([2] + [0]*(n-1), next_q, ring)

            # 在环中计算 f*g
            fg = f_cur * g_cur

            # 2 - f*g
            diff = two_poly - fg

            # g * (2 - f*g)
            g_next = g_cur * diff

            # 更新系数 (保持模 next_q)
            g_coeffs = g_next.coeffs[:]
            current_q = next_q

        # 最终结果 (模 Q)
        result = Polynomial(
            [c % q for c in g_coeffs], q, ring
        )

        # 验证
        prod = self * result
        if prod.coeffs[0] == 1 and all(c == 0 for c in prod.coeffs[1:]):
            return result
        return None

    def _invert_by_gcd(self) -> Optional['Polynomial']:
        """
        使用多项式扩展欧几里得算法求逆 (适用于素数模数)

        原理: 在 Z_Q[x] 中, 如果 gcd(f, x^n ± 1) = 1,
        则 ∃ g, k 使得 f·g + k·(x^n ± 1) = 1
        那么 g = f^{-1} mod (x^n ± 1)
        """
        q = self.modulus
        n = self.n

        # 取模环多项式: x^n - 1 或 x^n + 1
        if self.ring_type == "negacyclic":
            modulus_poly = [1] + [0] * (n - 1) + [1]  # x^n + 1
        else:
            modulus_poly = [-1] + [0] * (n - 1) + [1]  # x^n - 1

        a = self.coeffs[:]  # f(x)
        b = modulus_poly[:]  # x^n ± 1

        # 扩展欧几里得算法
        old_r, r = a[:], b[:]  # 余数
        old_u, u = [1], [0]    # Bezout 系数
        old_v, v = [0], [1]

        while not _poly_is_zero(r):
            quotient = _poly_div(old_r, r, q)
            if quotient is None:
                return None  # 除法的首项不可逆

            new_r = _poly_sub(old_r, _poly_mul(quotient, r, q), q)
            old_r, r = _poly_normalize(r), _poly_normalize(new_r)

            new_u = _poly_sub(old_u, _poly_mul(quotient, u, q), q)
            old_u, u = _poly_normalize(u), _poly_normalize(new_u)

            new_v = _poly_sub(old_v, _poly_mul(quotient, v, q), q)
            old_v, v = _poly_normalize(v), _poly_normalize(new_v)

        # old_r 是 gcd, 应该是一个非零常数
        if len(old_r) != 1 or old_r[0] == 0:
            return None

        gcd = old_r[0]
        try:
            inv_gcd = pow(gcd, -1, q)
        except ValueError:
            return None

        inv_coeffs = [(c * inv_gcd) % q for c in old_u]
        if len(inv_coeffs) > n:
            inv_coeffs = inv_coeffs[:n]
        elif len(inv_coeffs) < n:
            inv_coeffs += [0] * (n - len(inv_coeffs))

        return Polynomial(inv_coeffs, q, self.ring_type)

    @staticmethod
    def _mul_plain_trunc(a_coeffs: List[int], b_coeffs: List[int],
                         m: int, q: int) -> List[int]:
        """普通多项式乘法 (截断到第 m 项), 不进行 x^n ± 1 的模约化"""
        result = [0] * m
        for i in range(min(m, len(a_coeffs))):
            if a_coeffs[i] == 0:
                continue
            for j in range(min(m - i, len(b_coeffs))):
                if b_coeffs[j] == 0:
                    continue
                result[i + j] = (result[i + j] + a_coeffs[i] * b_coeffs[j]) % q
        return result

    def _invert_by_newton(self) -> Optional['Polynomial']:
        """
        使用牛顿迭代法求逆 (适用于较大的 n)

        对于 R^-_n (cyclic): 在普通截断环中迭代, 最后做环修正
        对于 R^+_n (negacyclic): 同样需要修正步骤
        """
        q = self.modulus
        n = self.n

        try:
            g0 = pow(self.coeffs[0], -1, q)
        except ValueError:
            return None

        g_coeffs = [g0]
        m = 1

        while m < n:
            next_m = min(2 * m, n)
            g_ext = g_coeffs + [0] * (next_m - len(g_coeffs))
            f_trunc = self.coeffs[:next_m]

            fg = self._mul_plain_trunc(f_trunc, g_ext, next_m, q)
            two_minus_fg = [(2 - fg[0]) % q]
            two_minus_fg += [(-fg[k]) % q for k in range(1, next_m)]

            g_coeffs = self._mul_plain_trunc(g_ext, two_minus_fg, next_m, q)
            m = next_m

        # 环修正: 在环中计算 f*g, 如果 ≠ 1, 做整体缩放
        tmp = Polynomial(self.coeffs, q, self.ring_type) * Polynomial(g_coeffs, q, self.ring_type)
        if tmp.coeffs[0] != 1 or any(c != 0 for c in tmp.coeffs[1:]):
            try:
                k_inv = pow(tmp.coeffs[0], -1, q)
            except ValueError:
                return None
            g_coeffs = [(c * k_inv) % q for c in g_coeffs]

        return Polynomial(g_coeffs, q, self.ring_type)


# ============================================================================
# 多项式扩展欧几里得算法辅助函数
# ============================================================================

def _poly_normalize(coeffs: List[int]) -> List[int]:
    """去除多项式尾部零系数"""
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs = coeffs[:-1]
    return coeffs


def _poly_is_zero(coeffs: List[int]) -> bool:
    """检查多项式是否为零"""
    return len(coeffs) == 1 and coeffs[0] == 0


def _poly_degree(coeffs: List[int]) -> int:
    """多项式的次数"""
    return len(coeffs) - 1


def _poly_div(a: List[int], b: List[int], q: int) -> Optional[List[int]]:
    """
    多项式长除法: a / b 在 Z_q 上

    返回商多项式, 如果除法的首项不可逆则返回 None
    """
    a = _poly_normalize(a)
    b = _poly_normalize(b)

    deg_a = _poly_degree(a)
    deg_b = _poly_degree(b)

    if deg_a < deg_b:
        return [0]

    # 除数的首项系数
    lead_b = b[-1]
    try:
        inv_lead_b = pow(lead_b, -1, q)
    except ValueError:
        return None

    quotient = [0] * (deg_a - deg_b + 1)
    remainder = a[:]

    for i in range(deg_a - deg_b, -1, -1):
        if len(remainder) <= deg_b + i - 1:
            break
        # 商系数 = 余数首项 / 除数首项
        coef = (remainder[deg_b + i] * inv_lead_b) % q
        quotient[i] = coef

        # 余数 -= coef * x^i * 除数
        for j in range(deg_b + 1):
            remainder[i + j] = (remainder[i + j] - coef * b[j]) % q

    return _poly_normalize(quotient)


def _poly_mul(a: List[int], b: List[int], q: int) -> List[int]:
    """普通多项式乘法 (不截断, 不进行模约化)"""
    deg_a = _poly_degree(a)
    deg_b = _poly_degree(b)
    result = [0] * (deg_a + deg_b + 1)
    for i in range(deg_a + 1):
        if a[i] == 0:
            continue
        for j in range(deg_b + 1):
            if b[j] == 0:
                continue
            result[i + j] = (result[i + j] + a[i] * b[j]) % q
    return result


def _poly_sub(a: List[int], b: List[int], q: int) -> List[int]:
    """多项式减法"""
    max_len = max(len(a), len(b))
    result = [0] * max_len
    for i in range(max_len):
        va = a[i] if i < len(a) else 0
        vb = b[i] if i < len(b) else 0
        result[i] = (va - vb) % q
    return _poly_normalize(result)


# ============================================================================
# 多项式矩阵形式
# ============================================================================

def circulant_matrix(poly: Polynomial) -> List[List[int]]:
    """
    计算多项式的循环矩阵 M(a)

    对于 R^-_n (cyclic):
        M(a) = [v(a), v(a·x), ..., v(a·x^{n-1})]
    其中每列是 v(a·x^k), 因此 M(a)[i][j] = a_{(i-j) mod n}

    此矩阵满足性质: M(a) * v(b) = v(a * b)

    对于 R^+_n (negacyclic):
        M(a) 是反循环矩阵 (anti-circulant)
    """
    n = poly.n
    q = poly.modulus
    matrix = []

    for i in range(n):
        row = []
        for j in range(n):
            # M(a)[i][j] = 系数 a_{k}, 其中 k = (i-j) mod n
            k = (i - j) % n
            val = poly.coeffs[k]
            if poly.ring_type == "negacyclic":
                # 对于反循环矩阵: 当 i < j 时取负
                if i < j:
                    val = -val
            row.append(val)
        matrix.append(row)

    return matrix


def matrix_vector_multiply(matrix: List[List[int]],
                           vec: 'Polynomial',
                           modulus: int) -> 'Polynomial':
    """
    多项式矩阵 × 系数向量的乘法
    M(a) · v (在Z_modulus下运算)

    注意：这里 vec 被视为长度为 n 的整数向量 (不是多项式)
    """
    n = len(matrix)
    q = modulus
    result = [0] * n

    for i in range(n):
        row = matrix[i]
        total = 0
        for j in range(n):
            total = (total + row[j] * vec.coeffs[j]) % q
        result[i] = total

    return Polynomial(result, q, "cyclic")  # ring_type 此处不关键


# ============================================================================
# 多项式环上的向量和矩阵
# ============================================================================

class PolyVector:
    """
    多项式向量: (p_0, p_1, ..., p_{m-1}), 每个 p_i ∈ R_Q
    """

    def __init__(self, polys: List[Polynomial]):
        if not polys:
            raise ValueError("多项式向量不能为空")
        self.polys = polys
        self.m = len(polys)
        self.n = polys[0].n
        self.modulus = polys[0].modulus
        self.ring_type = polys[0].ring_type

    def __add__(self, other: 'PolyVector') -> 'PolyVector':
        return PolyVector([a + b for a, b in zip(self.polys, other.polys)])

    def __sub__(self, other: 'PolyVector') -> 'PolyVector':
        return PolyVector([a - b for a, b in zip(self.polys, other.polys)])

    def to_coeff_vector(self) -> List[int]:
        """将所有多项式串联为一个系数向量 (用于矩阵运算)"""
        result = []
        for p in self.polys:
            result.extend(p.coeffs)
        return result

    @staticmethod
    def from_coeff_vector(coeffs: List[int], m: int, n: int,
                          modulus: int, ring_type: str) -> 'PolyVector':
        """从系数向量恢复多项式向量"""
        polys = []
        for i in range(m):
            start = i * n
            p_coeffs = coeffs[start:start + n]
            polys.append(Polynomial(p_coeffs, modulus, ring_type))
        return PolyVector(polys)

    def norm_squared(self) -> int:
        """计算向量范数的平方 Σ||p_i||^2"""
        return sum(p.norm_squared() for p in self.polys)

    def __repr__(self) -> str:
        return f"PolyVector(m={self.m}, n={self.n}, Q={self.modulus})"


# ============================================================================
# 测试代码
# ============================================================================

def _test():
    """多项式模块自测"""
    print("=" * 60)
    print("多项式模块测试")
    print("=" * 60)

    # 测试 R^-_n (cyclic)
    print("\n[测试1] R^-_n = Z[x]/(x^n-1) 基本运算")
    n, Q = 8, 65536
    a = Polynomial([1, 2, 3, 4, 5, 6, 7, 8], Q, "cyclic")
    b = Polynomial([8, 7, 6, 5, 4, 3, 2, 1], Q, "cyclic")
    c = a + b
    print(f"  a + b = {c.coeffs}")
    d = a * b
    print(f"  a * b = {d.coeffs[:4]}... (截断显示)")

    # 验证 (x^n - 1) 关系
    x = Polynomial([0, 1] + [0]*6, Q, "cyclic")
    x_pow_n = Polynomial([1] + [0]*7, Q, "cyclic")
    for _ in range(n):
        x_pow_n = x_pow_n * x
    print(f"  x^n mod (x^n - 1) = {x_pow_n.coeffs[:4]}... (应为 [1,0,0,0]...)")

    # 测试 R^+_n (negacyclic)
    print("\n[测试2] R^+_n = Z[x]/(x^n+1) 基本运算")
    a2 = Polynomial([1, 2, 3, 4, 5, 6, 7, 8], Q, "negacyclic")
    b2 = Polynomial([8, 7, 6, 5, 4, 3, 2, 1], Q, "negacyclic")
    c2 = a2 + b2
    print(f"  a + b = {c2.coeffs}")
    d2 = a2 * b2
    print(f"  a * b = {d2.coeffs[:4]}... (截断显示)")

    # 验证 (x^n + 1) 关系
    x2 = Polynomial([0, 1] + [0]*6, Q, "negacyclic")
    x_pow_n2 = Polynomial([1] + [0]*7, Q, "negacyclic")
    for _ in range(n):
        x_pow_n2 = x_pow_n2 * x2
    print(f"  x^n mod (x^n + 1) = {x_pow_n2.coeffs[:4]}... (应为 [-1,0,0,0]...)")

    # 测试共轭
    print("\n[测试3] 共轭 a(x^{-1})")
    a_bar = a.conjugate()
    print(f"  a = {a.coeffs}")
    print(f"  ā = {a_bar.coeffs}")

    # 测试自同构
    print("\n[测试4] 自同构 σ_k(a) = a(x^k)")
    sigma = a.automorphism(3)
    print(f"  σ_3(a) = {sigma.coeffs}")

    # 测试多项式求逆
    print("\n[测试5] 多项式求逆 (在 Z[x]/(x^n+1) 中)")
    small_q = 97  # 小素数方便测试
    f = Polynomial([3, 0, 1, 0, 0, 0, 0, 0], small_q, "negacyclic")
    f_inv = f.invert()
    if f_inv:
        prod = f * f_inv
        print(f"  f = {f.coeffs}")
        print(f"  f^{-1} = {f_inv.coeffs[:4]}...")
        print(f"  f * f^{-1} mod Q = {prod.coeffs}")
        print(f"  验证通过: {prod.coeffs[0] == 1 and all(c == 0 for c in prod.coeffs[1:])}")
    else:
        print("  求逆失败")

    # 测试矩阵形式
    print("\n[测试6] 循环矩阵 M(a)")
    mat = circulant_matrix(a)
    print(f"  M(a) 第1行: {mat[0]}")
    print(f"  M(a) 第2行: {mat[1]}")

    print("\n所有测试完成！")


if __name__ == "__main__":
    _test()
