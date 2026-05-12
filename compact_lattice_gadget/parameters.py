# -*- coding: utf-8 -*-
"""
参数集定义模块

包含论文中不同安全级别的参数集。

Robin (NTRU-based) 参数集 (论文 Table 5):
  - Robin-701:  NIST-I 安全级别
  - Robin-1061: NIST-III 安全级别
  - Robin-1279: NIST-V 安全级别

Eagle (Ring-LWE-based) 参数集 (论文 Table 7):
  - Eagle-512:  80-bit 安全级别
  - Eagle-1024: NIST-III 安全级别

此外，还包含用于演示和测试的小型参数集 (Demo)。

参考文献:
  CRYPTO 2023, Tables 5 & 7
  ePrint: https://eprint.iacr.org/2023/729
"""

from dataclasses import dataclass
from typing import Tuple


@dataclass
class RobinParams:
    """
    Robin 签名方案的参数

    属性:
        name: 参数集名称
        n: 环维度 (素数, Z[x]/(x^n-1))
        Q: 大模数
        p: 小模数 (gadget参数)
        q: 中等模数 (Q = p*q)
        a: 系数中 +1 的个数
        b: 系数中 -1 的个数
        alpha: 陷门质量控制参数
        r: Gadget高斯宽度
        s: 原像高斯宽度
        gamma: 扭转因子 (平衡 z0+e 和 z1 的大小)
        beta: 接受界
        security_level: 目标安全级别描述
    """
    name: str
    n: int
    Q: int
    p: int
    q: int
    a: int
    b: int
    alpha: float
    r: float
    s: float
    gamma: float
    beta: float
    security_level: str

    @property
    def ring_type(self) -> str:
        """Robin 使用 cyclic 环 R^-_n = Z[x]/(x^n-1)"""
        return "cyclic"


@dataclass
class EagleParams:
    """
    Eagle 签名方案的参数

    属性:
        name: 参数集名称
        n: 环维度 (2的幂, Z[x]/(x^n+1))
        Q: 大模数
        p: 小模数 (gadget参数)
        q: 中等模数 (Q = p*q)
        a: 系数中 +1 的个数
        b: 系数中 -1 的个数
        alpha: 陷门质量控制参数
        r: Gadget高斯宽度
        s: 原像高斯宽度
        gamma: 扭转因子
        beta: 接受界
        security_level: 目标安全级别描述
    """
    name: str
    n: int
    Q: int
    p: int
    q: int
    a: int
    b: int
    alpha: float
    r: float
    s: float
    gamma: float
    beta: float
    security_level: str

    @property
    def ring_type(self) -> str:
        """Eagle 使用 negacyclic 环 R^+_n = Z[x]/(x^n+1)"""
        return "negacyclic"


# ============================================================================
# Robin 参数集 (论文 Table 5)
# ============================================================================

ROBIN_PARAMS = {
    "Robin-701": RobinParams(
        name="Robin-701",
        n=701,
        Q=16384,
        p=2048,
        q=8,
        a=176,
        b=175,
        alpha=1.65,
        r=10.22,
        s=449.8,
        gamma=1.65,
        beta=28928.7,
        security_level="NIST-I (C/Q: 116/105)"
    ),
    "Robin-1061": RobinParams(
        name="Robin-1061",
        n=1061,
        Q=32768,
        p=4096,
        q=8,
        a=266,
        b=265,
        alpha=1.70,
        r=10.28,
        s=573.8,
        gamma=2.29,
        beta=62965.5,
        security_level="NIST-III (C/Q: 181/165)"
    ),
    "Robin-1279": RobinParams(
        name="Robin-1279",
        n=1279,
        Q=32768,
        p=4096,
        q=8,
        a=320,
        b=319,
        alpha=1.75,
        r=10.31,
        s=650.4,
        gamma=2.07,
        beta=70983.7,
        security_level="NIST-V (C/Q: 228/207)"
    ),
}

# ============================================================================
# Eagle 参数集 (论文 Table 7)
# ============================================================================

EAGLE_PARAMS = {
    "Eagle-512": EagleParams(
        name="Eagle-512",
        n=512,
        Q=16000,
        p=2000,
        q=8,
        a=128,
        b=128,
        alpha=1.70,
        r=10.17,
        s=394.2,
        gamma=1.36,
        beta=28493.5,
        security_level="80-bit Classic (C/Q: 79/71)"
    ),
    "Eagle-1024": EagleParams(
        name="Eagle-1024",
        n=1024,
        Q=32400,
        p=2700,
        q=12,
        a=256,
        b=256,
        alpha=1.70,
        r=15.42,
        s=841.5,
        gamma=1.19,
        beta=66118.5,
        security_level="NIST-III (C/Q: 176/160)"
    ),
}

# ============================================================================
# 演示用参数 (小规模, 用于快速测试)
# ============================================================================

DEMO_ROBIN_PARAMS = RobinParams(
    name="Robin-Demo",
    n=67,           # 小素数
    Q=16384,        # 2^14 (与论文一致, 虽然大但可工作)
    p=2048,
    q=8,
    a=30,           # 较高权重的三值多项式, 提高可逆概率
    b=29,           # b = a-1 确保 f(1) ≠ 0
    alpha=2.5,      # 宽松的质量要求
    r=10.0,         # Gadget 宽度
    s=300.0,        # 原像宽度
    gamma=1.5,
    beta=12000.0,   # 接受界
    security_level="Demo (NOT SECURE)"
)

DEMO_EAGLE_PARAMS = EagleParams(
    name="Eagle-Demo",
    n=64,           # 2的幂 (较小)
    Q=16384,        # 2^14
    p=2048,
    q=8,
    a=30,           # 较高权重
    b=29,
    alpha=2.5,
    r=10.0,
    s=300.0,
    gamma=1.3,
    beta=12000.0,
    security_level="Demo (NOT SECURE)"
)


# ============================================================================
# 参数查询函数
# ============================================================================

def get_robin_params(name: str = "Robin-701") -> RobinParams:
    """获取指定名称的 Robin 参数集"""
    if name in ROBIN_PARAMS:
        return ROBIN_PARAMS[name]
    if name == "demo":
        return DEMO_ROBIN_PARAMS
    raise ValueError(f"未知的 Robin 参数集: {name}")


def get_eagle_params(name: str = "Eagle-512") -> EagleParams:
    """获取指定名称的 Eagle 参数集"""
    if name in EAGLE_PARAMS:
        return EAGLE_PARAMS[name]
    if name == "demo":
        return DEMO_EAGLE_PARAMS
    raise ValueError(f"未知的 Eagle 参数集: {name}")


def list_params():
    """列出所有可用的参数集"""
    print("Robin 参数集 (NTRU-based):")
    for name, p in ROBIN_PARAMS.items():
        print(f"  {name}: n={p.n}, Q={p.Q}, p={p.p}, q={p.q}")
        print(f"         a={p.a}, b={p.b}, α={p.alpha}")
        print(f"         r={p.r}, s={p.s}, γ={p.gamma}, β={p.beta}")
        print(f"         安全级别: {p.security_level}")
    print()
    print("Eagle 参数集 (Ring-LWE-based):")
    for name, p in EAGLE_PARAMS.items():
        print(f"  {name}: n={p.n}, Q={p.Q}, p={p.p}, q={p.q}")
        print(f"         a={p.a}, b={p.b}, α={p.alpha}")
        print(f"         r={p.r}, s={p.s}, γ={p.gamma}, β={p.beta}")
        print(f"         安全级别: {p.security_level}")
    print()
    print("Demo 参数集 (仅用于测试):")
    dp = DEMO_ROBIN_PARAMS
    print(f"  {dp.name}: n={dp.n}, Q={dp.Q}, p={dp.p}, q={dp.q}")


if __name__ == "__main__":
    list_params()
