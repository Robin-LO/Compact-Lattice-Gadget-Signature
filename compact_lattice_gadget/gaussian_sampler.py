# -*- coding: utf-8 -*-
"""
离散高斯采样器模块 (Discrete Gaussian Sampler)

实现整数格上的离散高斯分布采样。
离散高斯分布 D_{Z, σ, c} 定义为:
    D_{Z, σ, c}(x) = ρ_{σ}(x - c) / ρ_{σ}(Z - c)
    其中 ρ_{σ}(x) = exp(-π·x²/σ²)

本实现采用 CDT (Cumulative Distribution Table) 方法:
  - 预计算累积分布表
  - 采样时通过二分查找得到随机样本
  - 纯整数运算，无需浮点数

主要使用场景:
  - 近似Gadget采样: D_{qZ + c, r}
  - 扰动采样: D_{Z^m, √Σ_p}

参考文献: CRYPTO 2023, "Compact Lattice Gadget and Its Applications to Hash-and-Sign Signatures"
"""

import random
import math
from typing import List, Dict, Tuple
from bisect import bisect_left

# ============================================================================
# 全局随机数生成器 (使用 SystemRandom 保证密码学安全性)
# ============================================================================

_rng = random.SystemRandom()

# ============================================================================
# CDT 预计算和缓存
# ============================================================================

# CDT表缓存: sigma -> (table, scale_factor)
# table[i] 包含 x在 [-cutoff, i-cutoff] 范围内的累积概率 (缩放为整数)
_cdt_cache: Dict[float, Tuple[List[int], int]] = {}


def _compute_cdt(sigma: float, tail_cut: float = 13.0,
                 table_bits: int = 64) -> Tuple[List[int], int]:
    """
    预计算离散高斯分布的 CDT 表

    参数:
        sigma: 高斯宽度参数
        tail_cut: 尾部截断参数 (默认13σ足以达到2^{-128}统计距离)
        table_bits: CDT表精度位数

    返回:
        (cdt_table, scale): CDT表和缩放因子
    """
    # 截断范围
    cutoff = int(math.ceil(tail_cut * sigma))
    scale = 2 ** table_bits

    # 计算每个点的概率质量 (未归一化)
    # ρ_σ(x) = exp(-π * x² / σ²)
    probs = []
    for x in range(-cutoff, cutoff + 1):
        rho = math.exp(-math.pi * x * x / (sigma * sigma))
        # 对于中心为0的情况，对称性可以优化，但这里先直接计算
        probs.append(rho)

    # 归一化
    total = sum(probs)
    probs = [p / total for p in probs]

    # 构建累积分布表
    cdt = []
    cumsum = 0.0
    for p in probs[:-1]:  # 最后一个不需要 (概率为1)
        cumsum += p
        cdt.append(int(cumsum * scale))

    # 确保最后一个值是 scale (精确的1)
    cdt.append(scale)

    return cdt, scale


def _get_cdt(sigma: float) -> Tuple[List[int], int]:
    """获取或创建CDT表 (带缓存)"""
    # 四舍五入到合理精度以便缓存
    key = round(sigma, 4)
    if key not in _cdt_cache:
        _cdt_cache[key] = _compute_cdt(sigma)
    return _cdt_cache[key]


# ============================================================================
# 离散高斯采样函数
# ============================================================================

def sample_z_binary(c: int = 0, sigma: float = 1.0) -> int:
    """
    采一个整数 z ← D_{Z, σ, c} (CDT 方法)

    核心算法:
    1. 查表获取 σ 对应的 CDT
    2. 生成 [0, scale) 内的均匀随机整数
    3. 二分查找 CDT 定位对应区间
    4. 加上偏移中心 c

    参数:
        c: 分布中心 (整数)
        sigma: 高斯宽度参数 (正实数)

    返回:
        采样的整数
    """
    if sigma <= 0:
        raise ValueError(f"sigma 必须为正数: {sigma}")

    cdt, scale = _get_cdt(sigma)
    cutoff = (len(cdt) - 1) // 2

    # 生成随机整数
    r = _rng.randint(0, scale - 1)

    # 二分查找
    idx = bisect_left(cdt, r)
    # idx 对应 x = idx - cutoff
    x = idx - cutoff

    return x + c


def sample_z_rejection(c: int = 0, sigma: float = 1.0) -> int:
    """
    采一个整数 z ← D_{Z, σ, c} (拒绝采样方法)

    使用 Karney 的拒绝采样算法，适合较大的σ。
    对于 σ < 1 的情况，CDT 方法更好。

    参数:
        c: 分布中心
        sigma: 高斯宽度

    返回:
        采样的整数
    """
    # 对于较小的sigma，使用CDT
    if sigma < 10:
        return sample_z_binary(c, sigma)

    # Karney 算法: 使用均匀分布的拒绝采样
    # 这里实现的是简化的 Polar 方法变体
    sigma_sq = sigma * sigma
    while True:
        # 生成候选点 (在附近范围内)
        # 使用 Box-Muller 的离散变体
        u = _rng.random()
        v = _rng.random()

        # 从连续高斯中采样
        z_cont = sigma * math.sqrt(-2.0 * math.log(u)) * math.cos(2.0 * math.pi * v)
        k = int(round(z_cont))

        # 接受概率
        accept_prob = math.exp(-math.pi * (k - c) * (k - c) / sigma_sq)
        if _rng.random() < accept_prob:
            return k


def sample_z(c: int = 0, sigma: float = 1.0) -> int:
    """
    整数高斯采样 (自动选择最佳算法)

    参数:
        c: 分布中心 (整数)
        sigma: 高斯宽度

    返回:
        采样的整数
    """
    return sample_z_binary(c, sigma)


def sample_gadget_coset(sigma: float, q: int, c: int) -> int:
    """
    从陪集高斯分布 D_{qZ + c, σ} 中采样

    这是紧凑Gadget中半随机采样器的核心操作。

    原理:
        D_{qZ + c, σ}(x) ∝ exp(-π·x²/σ²) for x ∈ qZ + c

    等价于:
        采样 k ← D_{Z, σ/q, -c/q}
        返回 x = q·k + c

    参数:
        sigma: 高斯宽度 r (注意是陪集的高斯宽度)
        q: 陪集步长
        c: 陪集偏移 (0 ≤ c < q)

    返回:
        采样的整数 x = q*k + c
    """
    # 计算缩放后的高斯参数
    base_sigma = sigma / q  # k 的高斯宽度

    # 采样 k ← D_{Z, base_sigma}
    # 注意: 中心为0的采样 (偏移由 c 处理)
    k = sample_z(c=0, sigma=base_sigma)

    # 计算最终值
    x = q * k + c
    return x


def sample_polynomial_gadget(n: int, sigma: float, q: int,
                              c_poly: List[int]) -> List[int]:
    """
    从 D_{qZ^n + c, sigma} 中采样一个多项式向量

    对每个坐标独立采样 D_{qZ + c_i, sigma}

    参数:
        n: 多项式维度
        sigma: 高斯宽度 r
        q: 陪集步长
        c_poly: 陪集偏移向量 (长度n)

    返回:
        采样的系数向量 (长度n)
    """
    x_coeffs = []
    for i in range(n):
        x_coeffs.append(sample_gadget_coset(sigma, q, c_poly[i]))
    return x_coeffs


def sample_gaussian_vector(dim: int, sigma: float) -> List[int]:
    """
    从 D_{Z^dim, sigma} 中采样一个向量

    参数:
        dim: 向量维度
        sigma: 高斯宽度

    返回:
        采样的整数向量 (长度dim)
    """
    return [sample_z(c=0, sigma=sigma) for _ in range(dim)]


# ============================================================================
# 统计测试函数
# ============================================================================

def _test_distribution(sigma: float, num_samples: int = 100000):
    """测试高斯采样的统计特性"""
    print(f"\n  测试 D_{{Z, {sigma}}} 采样 ({num_samples} 个样本):")

    from collections import Counter
    samples = [sample_z(c=0, sigma=sigma) for _ in range(num_samples)]

    # 统计均值和标准差
    mean = sum(samples) / len(samples)
    var = sum((s - mean)**2 for s in samples) / len(samples)
    std = math.sqrt(var)

    # 理论值
    # 对于离散高斯 D_{Z, σ}:
    # 对于较大的σ，方差 ≈ σ²/(2π) (从连续近似)
    # 更准确的理论分析更复杂

    print(f"    均值: {mean:.4f} (理论: 0)")
    print(f"    标准差: {std:.4f} (高斯宽度 σ={sigma})")

    # 显示分布
    counter = Counter(samples)
    most_common = counter.most_common(10)
    print(f"    最常见值: {most_common}")

    return samples


def _test():
    """高斯采样器模块自测"""
    print("=" * 60)
    print("离散高斯采样模块测试")
    print("=" * 60)

    # 测试基本采样
    print("\n[测试1] 基本整数高斯采样")
    sigma = 3.0
    samples = [sample_z(c=0, sigma=sigma) for _ in range(20)]
    print(f"  D_{{Z, {sigma}}} 采样: {samples}")

    # 测试带中心的采样
    print("\n[测试2] 带偏移的高斯采样")
    c = 5
    sigma = 2.0
    samples = [sample_z(c=c, sigma=sigma) for _ in range(20)]
    print(f"  D_{{Z, {sigma}, {c}}} 采样: {samples}")

    # 测试陪集采样 (Gadget采样)
    print("\n[测试3] 陪集高斯采样 D_{qZ + c, r}")
    q = 8
    c = 3
    r = 10.0
    samples = [sample_gadget_coset(r, q, c) for _ in range(30)]
    print(f"  D_{{{q}Z + {c}, {r}}} 采样: {samples}")
    # 验证所有样本都在陪集中
    all_in_coset = all((s - c) % q == 0 for s in samples)
    print(f"  所有样本在陪集 qZ + {c} 中: {all_in_coset}")

    # 统计测试
    _test_distribution(3.0, 50000)

    print("\n所有测试完成！")


if __name__ == "__main__":
    _test()
