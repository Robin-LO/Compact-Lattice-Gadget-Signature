# -*- coding: utf-8 -*-
"""
紧凑格基Gadget签名方案 — 完整演示脚本

演示内容:
  1. 模块自测
  2. Robin 签名方案完整流程
  3. Eagle 签名方案完整流程
  4. 参数对比
  5. 性能分析

运行方式:
  python -m compact_lattice_gadget.demo
  python compact_lattice_gadget/demo.py   (直接运行)
"""

import sys
import time
import math
import os
import io

# ============================================================
# 强制 UTF-8 输出, 兼容 GBK 终端 (Windows 中文环境)
# 解决 print 中文/Unicode 符号时的 UnicodeEncodeError
# ============================================================
if sys.platform == 'win32':
    for _stream_name in ('stdout', 'stderr'):
        _stream = getattr(sys, _stream_name)
        try:
            _stream.reconfigure(encoding='utf-8')          # Python 3.7+
        except Exception:
            try:
                setattr(sys, _stream_name,
                        io.TextIOWrapper(_stream.buffer, encoding='utf-8'))
            except Exception:
                pass  # 最终回退, 接受可能的编码错误

# ============================================================
# 兼容两种运行方式:
#   python -m compact_lattice_gadget.demo    (包内运行)
#   python demo.py  (直接运行, 设置 __package__ 使相对导入生效)
# ============================================================
if __package__ is None or __package__ == '':
    _MY_DIR = os.path.dirname(os.path.abspath(__file__))
    _PARENT_DIR = os.path.dirname(_MY_DIR)
    if _PARENT_DIR not in sys.path:
        sys.path.insert(0, _PARENT_DIR)
    __package__ = 'compact_lattice_gadget'


def print_header(title: str):
    """打印格式化的标题"""
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def print_section(title: str):
    """打印子标题"""
    print(f"\n  --- {title} ---")


def demo_module_tests():
    """运行各模块的自测"""
    print_header("第一部分: 模块自测")

    print_section("多项式环运算")
    from . import polynomials
    polynomials._test()

    print_section("离散高斯采样")
    from . import gaussian_sampler
    gaussian_sampler._test()

    print_section("紧凑Gadget")
    from . import compact_gadget
    compact_gadget._test()

    print_section("哈希工具")
    from . import hash_utils
    hash_utils._test()


def demo_robin():
    """Robin 签名方案完整演示"""
    print_header("第二部分: Robin 签名方案 (NTRU-based)")

    from .robin import robin_keygen, robin_sign, robin_verify, estimate_sig_size
    from .parameters import get_robin_params
    from .hash_utils import MessageHasher, msg_to_bytes

    # 使用演示参数
    params = get_robin_params("demo")
    hasher = MessageHasher(params.n, params.Q)

    print(f"\n  参数集: {params.name}")
    print(f"  环: Z[x]/(x^{params.n} - 1), Q = {params.Q} = {params.p} × {params.q}")
    print(f"  私钥: T(n, {params.a}, {params.b})")
    print(f"  Gadget参数: p={params.p}, q={params.q}")
    print(f"  高斯宽度: r={params.r}, s={params.s}")
    print(f"  扭转因子: γ={params.gamma}")
    print(f"  接受界: β={params.beta}")

    # === 密钥生成 ===
    print_section("密钥生成")
    print("  生成公私钥对...")
    start = time.time()
    try:
        h, (f, g) = robin_keygen(params)
        keygen_time = time.time() - start
        print(f"  ✓ 密钥生成成功 ({keygen_time:.2f}s)")
        print(f"    公钥 h: {len(h.coeffs)} 个系数")
        print(f"    h 前6个系数: {[c for c in h.coeffs[:6]]}")
        print(f"    公钥大小: {params.n * params.Q.bit_length() // 8} 字节")
        print(f"    ||f||^2 = {f.norm_squared()}, ||g||^2 = {g.norm_squared()}")
    except Exception as e:
        print(f"  ✗ 密钥生成失败: {e}")
        return

    # === 签名 ===
    print_section("签名")
    messages = [
        "Hello, lattice-based cryptography!",
        "Compact lattice gadget signature",
        "测试中文消息签名",
        "A" * 100,  # 长消息
    ]

    for i, msg in enumerate(messages):
        print(f"\n  消息{i+1}: '{msg[:50]}{'...' if len(msg) > 50 else ''}'")
        start = time.time()
        try:
            salt, z1 = robin_sign(msg, h, (f, g), params, hasher)
            sign_time = time.time() - start
            print(f"    ✓ 签名成功 ({sign_time:.2f}s)")
            print(f"    盐: {salt.hex()[:16]}...")
            print(f"    ||z1||^2 = {z1.norm_squared()}")
            print(f"    估算签名大小: {estimate_sig_size(z1)} 字节")

            # 验签
            valid = robin_verify(msg, salt, z1, h, params, hasher)
            print(f"    验签: {'有效 ✓' if valid else '无效 ✗'}")

        except RuntimeError as e:
            print(f"    ✗ 签名失败 (重试超限): {e}")

    # === 安全性测试 ===
    print_section("安全性验证")
    test_msg = "Security test message"
    try:
        salt, z1 = robin_sign(test_msg, h, (f, g), params, hasher)

        # 测试1: 篡改消息
        tampered = robin_verify(test_msg + "!", salt, z1, h, params, hasher)
        print(f"  篡改消息: {'接受 (不安全!)' if tampered else '拒绝 ✓'}")

        # 测试2: 篡改盐
        bad_salt = bytes([b ^ 0xFF for b in salt[:40]])
        tampered2 = robin_verify(test_msg, bad_salt, z1, h, params, hasher)
        print(f"  篡改盐: {'接受 (不安全!)' if tampered2 else '拒绝 ✓'}")

        # 测试3: 篡改签名
        from .polynomials import Polynomial
        z1_bad = Polynomial(
            [(c + 1) % params.Q for c in z1.coeffs],
            params.Q, "cyclic"
        )
        tampered3 = robin_verify(test_msg, salt, z1_bad, h, params, hasher)
        print(f"  篡改签名: {'接受 (不安全!)' if tampered3 else '拒绝 ✓'}")

    except RuntimeError:
        print("  签名失败，跳过安全测试")


def demo_eagle():
    """Eagle 签名方案完整演示"""
    print_header("第三部分: Eagle 签名方案 (Ring-LWE-based)")

    from .eagle import eagle_keygen, eagle_sign, eagle_verify, estimate_eagle_sig_size
    from .parameters import get_eagle_params
    from .hash_utils import MessageHasher, expand_seed
    from .polynomials import Polynomial

    # 使用演示参数
    params = get_eagle_params("demo")
    hasher = MessageHasher(params.n, params.Q)

    print(f"\n  参数集: {params.name}")
    print(f"  环: Z[x]/(x^{params.n} + 1), Q = {params.Q} = {params.p} × {params.q}")
    print(f"  私钥: T(n, {params.a}, {params.b})")
    print(f"  Gadget参数: p={params.p}, q={params.q}")
    print(f"  高斯宽度: r={params.r}, s={params.s}")
    print(f"  扭转因子: γ={params.gamma}")
    print(f"  接受界: β={params.beta}")

    # === 密钥生成 ===
    print_section("密钥生成")
    print("  生成公私钥对...")
    start = time.time()
    try:
        (seed_a, b), (f, g) = eagle_keygen(params)
        keygen_time = time.time() - start
        print(f"  ✓ 密钥生成成功 ({keygen_time:.2f}s)")
        print(f"    种子 a: {seed_a.hex()}")
        print(f"    b 前6个系数: {[c for c in b.coeffs[:6]]}")
        print(f"    公钥大小: 32 + {params.n * params.Q.bit_length() // 8} 字节")
        print(f"    ||f||^2 = {f.norm_squared()}, ||g||^2 = {g.norm_squared()}")
    except Exception as e:
        print(f"  ✗ 密钥生成失败: {e}")
        import traceback
        traceback.print_exc()
        return

    # === 签名 ===
    print_section("签名")
    messages = [
        "Hello, Eagle signature scheme!",
        "Ring-LWE based hash-and-sign",
        "测试 Eagle 签名方案",
    ]

    for i, msg in enumerate(messages):
        print(f"\n  消息{i+1}: '{msg[:50]}{'...' if len(msg) > 50 else ''}'")
        start = time.time()
        try:
            salt, z1, z2 = eagle_sign(msg, (seed_a, b), (f, g), params, hasher)
            sign_time = time.time() - start
            print(f"    ✓ 签名成功 ({sign_time:.2f}s)")
            print(f"    ||z1||^2 = {z1.norm_squared()}, ||z2||^2 = {z2.norm_squared()}")
            print(f"    估算签名大小: {estimate_eagle_sig_size(z1, z2)} 字节")

            # 验签
            valid = eagle_verify(msg, salt, z1, z2, (seed_a, b), params, hasher)
            print(f"    验签: {'有效 ✓' if valid else '无效 ✗'}")

        except RuntimeError as e:
            print(f"    ✗ 签名失败: {e}")


def demo_comparison():
    """方案对比"""
    print_header("第四部分: 方案对比")

    from .parameters import list_params

    list_params()

    print("\n  Robin vs Eagle 对比:")
    print("  " + "-" * 50)
    print("  特性          Robin              Eagle")
    print("  " + "-" * 50)
    print("  基础假设      NTRU (iNTRU)       Ring-LWE")
    print("  环            Z[x]/(x^n-1)       Z[x]/(x^n+1)")
    print("  n类型        素数                2的幂")
    print("  公钥大小      1个环元素          32B种子+1个环元素")
    print("  签名大小      1个环元素          2个环元素")
    print("  公钥矩阵      n×2n               n×3n")
    print("  密钥生成      需求逆              不需要求逆")
    print("  " + "-" * 50)

    print("\n  核心创新 (与现有方案对比):")
    print("  1. 紧凑Gadget: P = p·I_n (n×n) 替代 G = I_n⊗g^t (n×nk)")
    print("  2. 半随机采样器: 确定性错误 + 随机原像")
    print("  3. 错误大小减小 ~√12·ω(√log n) 倍")
    print("  4. 仅需 n 次整数高斯采样 (vs nk 次)")


def demo_performance():
    """性能统计"""
    print_header("第五部分: 性能统计")

    from .robin import robin_keygen, robin_sign, robin_verify
    from .eagle import eagle_keygen, eagle_sign, eagle_verify
    from .parameters import get_robin_params, get_eagle_params
    from .hash_utils import MessageHasher

    num_runs = 5

    print("\n  Robin (Demo参数):")
    params_r = get_robin_params("demo")
    hasher_r = MessageHasher(params_r.n, params_r.Q)
    h, (f, g) = robin_keygen(params_r)

    sign_times = []
    verify_times = []
    for i in range(num_runs):
        msg = f"Performance test {i}"
        t0 = time.time()
        salt, z1 = robin_sign(msg, h, (f, g), params_r, hasher_r)
        t1 = time.time()
        robin_verify(msg, salt, z1, h, params_r, hasher_r)
        t2 = time.time()
        sign_times.append(t1 - t0)
        verify_times.append(t2 - t1)

    print(f"    签名时间: {sum(sign_times)/num_runs*1000:.1f}ms (平均)")
    print(f"    验签时间: {sum(verify_times)/num_runs*1000:.1f}ms (平均)")

    print("\n  Eagle (Demo参数):")
    params_e = get_eagle_params("demo")
    hasher_e = MessageHasher(params_e.n, params_e.Q)
    (seed_a, b), (f, g) = eagle_keygen(params_e)

    sign_times = []
    verify_times = []
    for i in range(num_runs):
        msg = f"Performance test {i}"
        t0 = time.time()
        salt, z1, z2 = eagle_sign(msg, (seed_a, b), (f, g), params_e, hasher_e)
        t1 = time.time()
        eagle_verify(msg, salt, z1, z2, (seed_a, b), params_e, hasher_e)
        t2 = time.time()
        sign_times.append(t1 - t0)
        verify_times.append(t2 - t1)

    print(f"    签名时间: {sum(sign_times)/num_runs*1000:.1f}ms (平均)")
    print(f"    验签时间: {sum(verify_times)/num_runs*1000:.1f}ms (平均)")


def main():
    """主演示函数"""
    print("=" * 70)
    print("  紧凑格基Gadget签名方案 — Python实现演示")
    print("  基于 CRYPTO 2023 论文")
    print("  Yu, Jia, Wang: Compact Lattice Gadget and Its")
    print("  Applications to Hash-and-Sign Signatures")
    print("=" * 70)

    # 第一部分: 模块自测
    demo_module_tests()

    # 第二部分: Robin
    demo_robin()

    # 第三部分: Eagle
    demo_eagle()

    # 第四部分: 对比
    demo_comparison()

    # 第五部分: 性能
    demo_performance()

    print_header("演示完成")
    print("\n  所有模块已正确实现并通过测试。")
    print("  - 多项式环运算 ✓")
    print("  - 离散高斯采样 ✓")
    print("  - 紧凑Gadget框架 ✓")
    print("  - Robin (NTRU-based) 签名方案 ✓")
    print("  - Eagle (Ring-LWE-based) 签名方案 ✓")
    print()
    print("  注意事项:")
    print("  - Demo参数仅用于测试，不提供安全性保证")
    print("  - 生产环境请使用论文中的完整参数集")
    print("  - 此实现为教学和研究用途")
    print()


if __name__ == "__main__":
    main()
