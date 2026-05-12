# -*- coding: utf-8 -*-
"""
哈希工具模块

提供签名方案所需的哈希函数和消息处理功能:
  - SHAKE-256 哈希函数 (用于随机预言机)
  - 消息到多项式环元素的映射
  - 随机盐 (salt) 生成

参考文献: CRYPTO 2023, Sections 5-6
"""

import hashlib
import os
from typing import List, Tuple


# ============================================================================
# 哈希函数
# ============================================================================

def shake256(data: bytes, output_len: int) -> bytes:
    """
    SHAKE-256 扩展输出函数 (XOF)

    这是论文中使用的哈希函数。SHAKE-256 是 SHA-3 的可扩展输出函数,
    可以提供任意长度的输出。

    参数:
        data: 输入数据
        output_len: 输出字节数

    返回:
        output_len 字节的哈希输出
    """
    return hashlib.shake_256(data).digest(output_len)


def hash_to_bytes(msg: bytes, salt: bytes, output_bytes: int) -> bytes:
    """
    对消息进行哈希, 输出指定长度的字节串

    H(msg, salt) = SHAKE-256(salt || msg, output_bytes)

    参数:
        msg: 消息 (字节串)
        salt: 盐 (40字节 = 320位, 论文中的参数)
        output_bytes: 输出字节数

    返回:
        哈希值
    """
    return shake256(salt + msg, output_bytes)


def bytes_to_poly_coeffs(data: bytes, n: int, Q: int) -> List[int]:
    """
    将字节串转换为多项式系数 (在 Z_Q 中)

    使用简单的解析方法:
      - 将字节视为大整数
      - 逐系数取模 Q

    参数:
        data: 输入字节串
        n: 多项式维度
        Q: 系数模数

    返回:
        n 个在 [-Q/2, Q/2) 范围内的系数
    """
    # 将字节串转为大整数
    big_int = int.from_bytes(data, byteorder='little')

    coeffs = []
    for _ in range(n):
        coeff = big_int % Q
        # 转换到对称表示
        if coeff > Q // 2:
            coeff -= Q
        coeffs.append(coeff)
        big_int //= Q

    return coeffs


# ============================================================================
# 消息哈希 (到多项式环)
# ============================================================================

class MessageHasher:
    """
    消息哈希器: 将任意消息映射到多项式环 R_Q 中的目标 u

    用于签名方案的随机预言机:
      u = H(msg, salt) ∈ Z^n_Q
    """

    def __init__(self, n: int, Q: int):
        """
        参数:
            n: 多项式维度
            Q: 系数模数
        """
        self.n = n
        self.Q = Q
        # 计算需要的字节数 (每个系数最多需要 ceil(log2(Q)) 位)
        self.bits_per_coeff = (Q.bit_length() + 7) // 8
        self.output_bytes = n * self.bits_per_coeff

    def generate_salt(self) -> bytes:
        """
        生成随机盐 (40字节 = 320位, 论文默认)

        使用 os.urandom 生成密码学安全的随机盐
        """
        return os.urandom(40)

    def hash_message(self, msg: bytes, salt: bytes) -> List[int]:
        """
        将消息哈希到多项式系数

        返回 n 个在 [-Q/2, Q/2) 范围内的系数

        参数:
            msg: 消息 (任意字节串)
            salt: 随机盐

        返回:
            目标多项式 u 的系数列表
        """
        # SHAKE-256 输出
        hash_output = hash_to_bytes(msg, salt, self.output_bytes)

        # 转换为多项式系数
        return bytes_to_poly_coeffs(hash_output, self.n, self.Q)


# ============================================================================
# 种子扩展 (用于Eagle的公钥)
# ============================================================================

def expand_seed(seed: bytes, n: int, Q: int) -> List[int]:
    """
    从种子扩展生成均匀随机的多项式 a ∈ R_Q

    用于 Eagle 方案中从 32 字节种子生成公钥多项式。
    使用 SHAKE-256 作为伪随机生成器。

    参数:
        seed: 32字节种子
        n: 多项式维度
        Q: 系数模数

    返回:
        多项式 a 的 n 个系数 (在 Z_Q 中均匀分布)
    """
    bits_per_coeff = (Q.bit_length() + 7) // 8
    output_bytes = n * bits_per_coeff
    random_bytes = shake256(seed, output_bytes)
    return bytes_to_poly_coeffs(random_bytes, n, Q)


# ============================================================================
# 辅助函数
# ============================================================================

def random_seed_32() -> bytes:
    """生成32字节的随机种子"""
    return os.urandom(32)


def msg_to_bytes(msg) -> bytes:
    """将消息转为字节串 (支持 str 和 bytes)"""
    if isinstance(msg, str):
        return msg.encode('utf-8')
    elif isinstance(msg, bytes):
        return msg
    else:
        return str(msg).encode('utf-8')


# ============================================================================
# 测试代码
# ============================================================================

def _test():
    """哈希工具模块自测"""
    print("=" * 60)
    print("哈希工具模块测试")
    print("=" * 60)

    n, Q = 16, 65536

    # 测试1: 基本SHAKE-256
    print("\n[测试1] SHAKE-256 基本测试")
    output = shake256(b"Hello, Lattice!", 32)
    print(f"  SHAKE-256('Hello, Lattice!', 32B) = {output.hex()[:32]}...")

    # 测试2: 消息哈希器
    print("\n[测试2] 消息到多项式的哈希")
    hasher = MessageHasher(n, Q)
    msg = b"Test message for signature"
    salt = hasher.generate_salt()
    print(f"  盐 (hex): {salt.hex()[:20]}...")

    u_coeffs = hasher.hash_message(msg, salt)
    print(f"  目标 u (前8个系数): {u_coeffs[:8]}")

    # 验证系数在范围内
    in_range = all(-Q // 2 <= c < Q // 2 for c in u_coeffs)
    print(f"  系数范围验证: {in_range}")

    # 测试3: 同一消息同一盐应产生相同输出
    print("\n[测试3] 哈希确定性")
    u2 = hasher.hash_message(msg, salt)
    consistent = u_coeffs == u2
    print(f"  相同输入产生相同输出: {consistent}")

    # 测试4: 种子扩展
    print("\n[测试4] 种子扩展 (用于Eagle)")
    seed = random_seed_32()
    a_coeffs = expand_seed(seed, 8, 65536)
    print(f"  种子 (hex): {seed.hex()}")
    print(f"  扩展的多项式系数: {a_coeffs[:4]}...")

    print("\n所有测试完成！")


if __name__ == "__main__":
    _test()
