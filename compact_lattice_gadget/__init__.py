# -*- coding: utf-8 -*-
"""
Compact Lattice Gadget Signature
=================================

基于 CRYPTO 2023 论文《Compact Lattice Gadget and Its Applications
to Hash-and-Sign Signatures》实现的紧凑格基Gadget签名方案。

本包包含两个后量子签名方案的完整实现:
  - Robin: 基于 NTRU 的 Hash-and-Sign 签名
  - Eagle: 基于 Ring-LWE 的 Hash-and-Sign 签名

核心创新:
  - 紧凑Gadget (P = p·I_n) 替代传统的大尺寸Gadget
  - 半随机采样器 (Semi-random Sampler)
  - 确定性的错误解码 + 随机原像采样

模块结构:
  polynomials       - 多项式环运算 Z[x]/(x^n ± 1)
  gaussian_sampler  - 离散高斯采样 (CDT方法)
  compact_gadget    - 紧凑Gadget和半随机采样器
  trapdoor          - 近似陷门框架
  robin             - Robin签名方案 (NTRU-based)
  eagle             - Eagle签名方案 (Ring-LWE-based)
  hash_utils        - 哈希工具 (SHAKE-256)
  parameters        - 参数集定义

参考文献:
  Yu, Jia, Wang. "Compact Lattice Gadget and Its Applications to
  Hash-and-Sign Signatures." CRYPTO 2023.
"""

__version__ = "0.1.0"
__author__ = "Based on the paper by Yang Yu, Huiwen Jia, Xiaoyun Wang"

from .polynomials import Polynomial, PolyVector, circulant_matrix
from .gaussian_sampler import sample_z, sample_gadget_coset
from .compact_gadget import lattice_decoder, approx_gadget
from .trapdoor import approx_preimage_sample
from .robin import robin_keygen, robin_sign, robin_verify
from .eagle import eagle_keygen, eagle_sign, eagle_verify
from .parameters import (
    RobinParams, EagleParams,
    get_robin_params, get_eagle_params,
    DEMO_ROBIN_PARAMS, DEMO_EAGLE_PARAMS,
    list_params
)
from .hash_utils import MessageHasher
