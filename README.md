# 紧凑格基Gadget签名方案 — Python实现

## 论文信息

- **标题**: Compact Lattice Gadget and Its Applications to Hash-and-Sign Signatures
- **作者**: Yang Yu (Tsinghua), Huiwen Jia (Guangzhou Univ.), Xiaoyun Wang (Tsinghua)
- **发表**: CRYPTO 2023, LNCS 14085, pp. 390–420
- **ePrint**: https://eprint.iacr.org/2023/729

## 概述

本实现基于 CRYPTO 2023 论文, 完整实现了两个后量子数字签名方案:

| 方案 | 基础假设 | 环 | 特点 |
|------|----------|-----|------|
| **Robin** | NTRU (iNTRU) | $\mathbb{Z}[x]/(x^n-1)$, $n$ 为素数 | 签名仅需 1 个环元素 |
| **Eagle** | Ring-LWE | $\mathbb{Z}[x]/(x^n+1)$, $n$ 为 2 的幂 | 安全假设更广泛 |

### 核心创新

论文的核心贡献是**紧凑Gadget**框架和**半随机采样器**:

1. **紧凑Gadget**: 使用 $\mathbf{P} = p \cdot \mathbf{I}_n$ ($n \times n$ 方阵) 替代传统方案中 $\mathbf{I}_n \otimes \mathbf{g}^t$ ($n \times nk$ 长矩阵), 大幅减小公钥和签名尺寸。

2. **半随机采样器** (Semi-random Sampler):
   - **确定性错误解码**: $e = u \bmod p$ (逐系数模运算)
   - **随机原像采样**: $\mathbf{x}' \leftarrow D_{\Lambda(\mathbf{Q}) + \mathbf{c},\, r}$ (陪集高斯采样)

3. **错误尺寸减小**: 相比 [CGM19] 的方案, 错误标准差减小约 $\sqrt{12} \cdot \omega(\sqrt{\log n})$ 倍。

## 项目结构

```
compact_lattice_gadget/
├── __init__.py              # 包入口
├── polynomials.py           # 多项式环运算 Z[x]/(x^n±1)
├── gaussian_sampler.py      # 离散高斯采样 (CDT方法)
├── compact_gadget.py        # 紧凑Gadget与半随机采样器
├── trapdoor.py              # 近似陷门框架 (ApproxPreSamp)
├── robin.py                 # Robin签名方案 (NTRU-based)
├── eagle.py                 # Eagle签名方案 (Ring-LWE-based)
├── hash_utils.py            # 哈希工具 (SHAKE-256)
├── parameters.py            # 参数集定义
├── demo.py                  # 演示与测试脚本
└── README.md                # 本文件
```

## 快速开始

### 环境要求

- Python 3.8+ (仅需标准库, 可选 numpy 加速)

### 运行演示

```bash
# 完整演示 (包含所有模块)
python -m compact_lattice_gadget.demo

# 单独测试各模块
python -m compact_lattice_gadget.polynomials
python -m compact_lattice_gadget.gaussian_sampler
python -m compact_lattice_gadget.compact_gadget

# 测试 Robin 签名方案
python -m compact_lattice_gadget.robin

# 测试 Eagle 签名方案
python -m compact_lattice_gadget.eagle
```

### 基本使用

```python
from compact_lattice_gadget import (
    robin_keygen, robin_sign, robin_verify,
    get_robin_params, MessageHasher
)

# 1. 获取参数 (demo参数仅用于测试)
params = get_robin_params("demo")
hasher = MessageHasher(params.n, params.Q)

# 2. 密钥生成
h, (f, g) = robin_keygen(params)

# 3. 签名
msg = "Hello, lattice-based cryptography!"
salt, signature = robin_sign(msg, h, (f, g), params, hasher)

# 4. 验签
is_valid = robin_verify(msg, salt, signature, h, params, hasher)
print(f"签名有效: {is_valid}")
```

## 算法详解

### 1. 多项式环

支持两种多项式环:

- **$R_n^- = \mathbb{Z}[x]/(x^n-1)$** (Robin): $n$ 为素数, 循环卷积
- **$R_n^+ = \mathbb{Z}[x]/(x^n+1)$** (Eagle): $n$ 为 2 的幂, 反循环卷积

关键操作:
- 加法 / 减法 / 乘法 (模 $Q$ 和模 $x^n \pm 1$)
- 共轭: $\bar{a} = a(x^{-1})$
- 自同构: $\sigma_k(a) = a(x^k)$
- 矩阵形式 $M(a)$: 循环矩阵 ($R_n^-$) 或反循环矩阵 ($R_n^+$)
- 多项式求逆 (Hensel 提升法, 适用于 $Q$ 为 2 的幂)

### 2. 离散高斯采样

采用 CDT (Cumulative Distribution Table) 方法进行离散高斯采样:

$$D_{\mathbb{Z},\, \sigma,\, c}(x) \propto \exp\!\left(-\pi \cdot \frac{(x-c)^2}{\sigma^2}\right)$$

- 预计算 CDT 表并缓存
- 支持陪集高斯采样 $D_{q\mathbb{Z} + c,\, r}$ (用于 Gadget 采样)
- 纯整数运算, 无需浮点数

### 3. 紧凑Gadget 与半随机采样器

> 以下为论文 Algorithm 1 (ApproxGadget) 的核心流程。

**格解码器** (Lattice Decoder):

给定目标 $\mathbf{u}' \in \mathbb{Z}_Q^n$, 逐坐标计算:

$$\begin{aligned}
e_i &= u'_i \bmod p \quad \in [-\tfrac{p}{2}, \tfrac{p}{2}) \\[4pt]
c_i &= \frac{u'_i - e_i}{p}
\end{aligned}$$

满足 $\mathbf{u}' = p \cdot \mathbf{c} + \mathbf{e}$ (在整数上精确成立, 而不仅是模 $Q$)。

**半随机采样器** (ApproxGadget):

| 步骤 | 操作 | 说明 |
|------|------|------|
| 1 | $(\mathbf{c}, \mathbf{e}) \leftarrow \text{LatticeDecoder}(\mathbf{u}', p)$ | 确定性错误解码 |
| 2 | $\mathbf{x}' \leftarrow D_{q\mathbb{Z}^n + \mathbf{c},\, r}$ | 随机原像采样 (逐坐标陪集高斯) |

输出 $(\mathbf{x}', \mathbf{e})$ 满足:

$$\mathbf{P} \cdot \mathbf{x}' = \mathbf{u}' - \mathbf{e} \pmod Q$$

### 4. Robin 签名方案 (NTRU-based)

**密钥生成** (Algorithm 3):

1. 采样三值多项式 $f, g \in \mathcal{T}(n, a, b)$ — 含 $a$ 个 $+1$, $b$ 个 $-1$, 其余为 $0$
2. 优化 $s_1\!\left(M(f\bar{f} + g\bar{g})\right)$ 使得陷门质量足够好
3. 计算公钥: $h = (p - g) \cdot f^{-1} \bmod Q$
4. 私钥: $(f, g)$, 满足 $h \cdot f + g = p \pmod Q$

其中 $s_1(\cdot)$ 表示矩阵的最大奇异值, 通过 DFT 谱范数估计。

**签名** (Algorithm 4):

1. 构造公钥矩阵 $\mathbf{A} = [\,\mathbf{I}_n \mid M(h)\,]$ 和陷门 $\mathbf{T} = \begin{bmatrix}M(g) \\ M(f)\end{bmatrix}$
2. $\text{salt} \leftarrow \{0,1\}^{320}$, $\mathbf{u} \leftarrow H(\text{msg}, \text{salt})$
3. $(\mathbf{z}_0, \mathbf{z}_1) \leftarrow \text{ApproxPreSamp}(\mathbf{A}, \mathbf{T}, \mathbf{u}, r, s)$
4. 检查 $\|(\mathbf{z}_0 + \mathbf{e},\, \gamma \cdot \mathbf{z}_1)\| \leq \beta$
5. 输出 $(\text{salt}, \mathbf{z}_1)$

**验签** (Algorithm 5):

1. $\mathbf{u} \leftarrow H(\text{msg}, \text{salt})$
2. $\mathbf{z}' = \mathbf{u} - h \cdot \mathbf{z}_1 \bmod Q$
3. 接受当且仅当 $\|(\mathbf{z}',\, \gamma \cdot \mathbf{z}_1)\| \leq \beta$

**公钥矩阵关系的正确性验证**:

$$\mathbf{A} \cdot \mathbf{T}
= M(g) + M(h) \cdot M(f)
= M(g + h \cdot f)
= M(p)
= p \cdot \mathbf{I}_n$$

### 5. Eagle 签名方案 (Ring-LWE-based)

**密钥生成** (Algorithm 6):

1. $\text{seed}_a \leftarrow \{0,1\}^{256}$, $a \leftarrow \text{Expand}(\text{seed}_a) \in R_Q$
2. 采样 $f, g \in \mathcal{T}(n, a, b)$
3. $b = p - (a \cdot f + g) \bmod Q$
4. 公钥: $(\text{seed}_a, b)$, 私钥: $(f, g)$, 满足 $a \cdot f + g + b = p \pmod Q$

**签名** (Algorithm 7):

1. 构造 $\mathbf{A} = [\,\mathbf{I}_n \mid M(a) \mid M(b)\,]$, $\mathbf{T} = \begin{bmatrix}M(g) \\ M(f) \\ \mathbf{I}_n\end{bmatrix}$
2. $(\mathbf{z}_0, \mathbf{z}_1, \mathbf{z}_2) \leftarrow \text{ApproxPreSamp}(\mathbf{A}, \mathbf{T}, \mathbf{u}, r, s)$
3. 输出 $(\text{salt}, \mathbf{z}_1, \mathbf{z}_2)$

**验签** (Algorithm 8):

1. $\mathbf{z}' = \mathbf{u} - a \cdot \mathbf{z}_1 - b \cdot \mathbf{z}_2 \bmod Q$
2. 接受当且仅当 $\|(\mathbf{z}',\, \gamma \cdot \mathbf{z}_1,\, \gamma \cdot \mathbf{z}_2)\| \leq \beta$

**公钥矩阵关系的正确性验证**:

$$\mathbf{A} \cdot \mathbf{T}
= M(g) + M(a) \cdot M(f) + M(b)
= M(g + a \cdot f + b)
= M(p)
= p \cdot \mathbf{I}_n$$

## 参数集

### Robin 参数 (论文 Table 5)

| 参数集 | $n$ | $Q$ | $p$ | $q$ | $a$ | $b$ | 安全级别 |
|--------|-----|-----|-----|-----|-----|-----|---------|
| Robin-701 | 701 | 16384 | 2048 | 8 | 176 | 175 | NIST-I |
| Robin-1061 | 1061 | 32768 | 4096 | 8 | 266 | 265 | NIST-III |
| Robin-1279 | 1279 | 32768 | 4096 | 8 | 320 | 319 | NIST-V |

### Eagle 参数 (论文 Table 7)

| 参数集 | $n$ | $Q$ | $p$ | $q$ | $a$ | $b$ | 安全级别 |
|--------|-----|-----|-----|-----|-----|-----|---------|
| Eagle-512 | 512 | 16000 | 2000 | 8 | 128 | 128 | 80-bit Classic |
| Eagle-1024 | 1024 | 32400 | 2700 | 12 | 256 | 256 | NIST-III |

### Demo 参数 (仅用于测试, 不提供安全保证)

| 参数集 | $n$ | $Q$ | $p$ | $q$ | $a$ | $b$ |
|--------|-----|-----|-----|-----|-----|-----|
| Robin-Demo | 67 | 16384 | 2048 | 8 | 30 | 29 |
| Eagle-Demo | 64 | 16384 | 2048 | 8 | 30 | 29 |

## 安全性

两个方案的 EU-CMA (存在不可伪造性) 安全性在随机预言机模型下归约到:

- **Robin**: NTRU-SIS 问题 (扭曲范数版本)
- **Eagle**: Ring-SIS 问题 + Ring-LWE 问题

核心安全性质 (Theorem 2):

> 对均匀随机的目标 $\mathbf{u}$, 近似原像 $(\mathbf{x}, \mathbf{e})$ 的分布在统计上与不使用陷门的模拟分布不可区分。

## 性能对比

根据论文数据, 在 NIST-III 安全级别:

| 方案 | 公钥大小 | 签名大小 |
|------|---------|---------|
| Falcon-1024 | 1792 B | 1249 B |
| Mitaka-1024 | 1792 B | 1376 B |
| **Robin-1061** | 1990 B | 1527 B |
| Dilithium-3 | 1952 B | 3293 B |
| [CGM19] | 7712 B | 7172 B |
| **Eagle-1024** | 1952 B | 3052 B |

## 与相关工作的比较

### 与 [CGM19] (Approximate Trapdoor) 的对比

| 特性 | [CGM19] | 本工作 |
|------|---------|--------|
| Gadget 矩阵 | $\mathbf{I}_n \otimes (b^\ell, \dots, b^{k-1})^t$ | $p \cdot \mathbf{I}_n$ |
| Gadget 尺寸 | $n \times n(k-\ell)$ | $n \times n$ |
| 错误分布 | 高斯分布 | 均匀分布 (确定性) |
| 错误尺寸 | $\approx p \cdot \omega(\sqrt{n \log n})$ | $\approx p \cdot \sqrt{n/12}$ |
| 整数采样次数 | $n \cdot k$ | $n$ |

### 与 Falcon / Mitaka 的对比

- **优势**: 更简单的密钥生成 (无需生成 NTRU 完整陷门基)
- **优势**: 在线/离线结构便于侧信道保护
- **优势**: 纯整数运算实现, 无需浮点数
- **劣势**: 签名和公钥略大 (约 20%–40%)

## 实现说明

### 多项式求逆

对于 $Q$ 为 2 的幂的情况 (论文的标准参数), 使用 **Hensel 提升法**:

1. 在 $\mathbb{Z}_2[x]/(x^n \pm 1)$ 中用扩展欧几里得算法求逆 (模 2)
2. 迭代提升:

$$g_{k+1} = g_k \cdot (2 - f \cdot g_k) \bmod 2^{k+1}$$

### 高斯采样

使用 CDT (累积分布表) 方法, 预计算不同 $\sigma$ 值的 CDT 表。支持:

- 整数格采样: $D_{\mathbb{Z},\, \sigma}$
- 陪集采样: $D_{q\mathbb{Z} + c,\, r} = q \cdot D_{\mathbb{Z},\, r/q} + c$

### DFT 谱范数估计

利用循环矩阵的谱性质, 通过 DFT 高效计算 $M(f\bar{f} + g\bar{g})$ 的最大奇异值:

$$s_1\!\left(M(a)\right) = \max_{j \in \mathbb{Z}_n} \left|\hat{a}(\omega_n^j)\right|$$

其中 $\hat{a}(\omega)$ 为多项式 $a$ 在单位根处的 DFT 值。当 numpy 可用时自动使用 FFT 加速 ($O(n \log n)$), 否则回退到纯 Python 实现 ($O(n^2)$)。

## 扩展方向

基于本实现可以进一步研究:

1. **参数优化**: 使用论文中的完整参数集实现安全级别的签名方案
2. **NTT 加速**: 对大规模 $n$ 使用 NTT 加速多项式乘法
3. **侧信道保护**: 实现常数时间的 CDT 采样
4. **熵编码**: 使用 ANS 编码减小签名尺寸
5. **高级原语**: 基于紧凑Gadget构建群签名、属性基加密等高级方案
6. **其他格**: 探索其他 $\Lambda(\mathbf{P})$ 和 $\Lambda(\mathbf{Q})$ 的实例化

## 参考文献

1. Yu, Jia, Wang. "Compact Lattice Gadget and Its Applications to Hash-and-Sign Signatures." CRYPTO 2023.
2. Micciancio, Peikert. "Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller." EUROCRYPT 2012.
3. Chen, Genise, Mukherjee. "Approximate Trapdoors for Lattices and Smaller Hash-and-Sign Signatures." ASIACRYPT 2019.
4. Gentry, Peikert, Vaikuntanathan. "Trapdoors for Hard Lattices and New Cryptographic Constructions." STOC 2008.
5. Falcon: https://falcon-sign.info/
6. Dilithium: https://pq-crystals.org/dilithium/

## 免责声明

此实现仅用于教学和研究目的, 未经过正式安全审计, 不应用于生产环境。Demo 参数集不提供任何安全保证。
