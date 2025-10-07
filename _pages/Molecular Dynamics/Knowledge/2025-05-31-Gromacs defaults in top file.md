---
title: "GROMACS Defaults in Topology Files: Understanding comb-rule and nonbond_params"
date: "2025-05-31"
description: "详细解析 GROMACS 拓扑文件中 defaults 指令下的组合规则和非键参数。深入理解分子动力学力场参数，为模拟配置和参数优化提供完整指南。"
tags: [gromacs, topology, nonbond-params, comb-rule, force-field, molecular-dynamics, parameters]
thumbnail: "/assets/img/thumbnail/dsygx.png"
image: "/assets/img/thumbnail/dsygx.png"
---
# GROMACS 中 comb-rule 与 [nonbond_params] 参数解析

本文档旨在详细解释 GROMACS 拓扑文件中 `[defaults]` 指令下的 `comb-rule`（组合规则）以及 `[atomtypes]` 和 `[nonbond_params]` 部分中非键参数（特别是 Lennard-Jones 参数）的含义和解释方式。

## 一、`[defaults]` 指令详解

在 GROMACS 的拓扑文件（通常是 `.top` 文件或力场主 `.itp` 文件）中，`[defaults]` 指令用于设定非键相互作用的全局默认行为。

### 示例

```
[ defaults ]
; nbfunc      comb-rule       gen-pairs       fudgeLJ fudgeQQ
    1           2               no              1.0     1.0
```

### 参数解释

* **`nbfunc`** (Non-bonded function type)：定义非键势函数类型。
  * `1`：Lennard-Jones 势。这是**绝大多数经典力场**（如 AMBER, CHARMM, OPLS, Martini）使用的形式。
  
  * `2`：Buckingham 势。**注意**：根据 GROMACS 文档和社区讨论，Buckingham 势 (`nbfunc = 2`) 自 GROMACS 2019 版本后可能已被弃用或不再完全支持。
    
    参考链接：https://gromacs.bioexcel.eu/t/how-use-desired-mixing-rule-in-gromacs/10409/3
  
* **`comb-rule`** (Combination rule)：定义当 `[nonbond_params]` 部分**没有显式给出**不同原子类型 `i` 和 `j` 之间的非键参数时，如何从各自的原子类型参数（`[atomtypes]` 部分的参数）计算出交叉项参数。

* **`gen-pairs`** (Generate 1-4 pairs)：决定是否自动生成1-4相互作用对（即通过3个键连接的原子对）。
  * `yes`：根据成键信息自动生成，并通常与 `fudgeLJ` 和 `fudgeQQ` 联用。
  * `no`：不自动生成，1-4相互作用需要在 `[pairs]` 或 `[pairtypes]` 部分显式定义，或者由力场设计本身通过其他方式处理（如Martini）。

* **`fudgeLJ`**：如果 `gen-pairs = yes`，此参数定义了1-4相互作用中 Lennard-Jones 部分的缩放因子。

* **`fudgeQQ`**：如果 `gen-pairs = yes`，此参数定义了1-4相互作用中静电部分的缩放因子。

### GROMACS `comb-rule`：对 `[atomtypes]` 参数的解释及交叉项的计算

`comb-rule` 的设置**直接影响** GROMACS 如何解释 `[atomtypes]` 部分中的 `V` 和 `W` 列参数，以及在 `[nonbond_params]` 中没有显式定义一对原子类型间的非键参数时，如何计算这些交叉项参数。

https://manual.gromacs.org/current/reference-manual/topologies/parameter-files.html#non-bonded-parameters

#### 1. `[atomtypes]` 中 `V` 和 `W` 参数的解释

根据 GROMACS 手册:

* **如果 `comb-rule = 1`**:
  * $V_{ii}$ 代表 $C_{6,ii} = 4 \epsilon_{ii} \sigma_{ii}^6$ (单位：kJ mol⁻¹ nm⁶)
  * $W_{ii}$ 代表 $C_{12,ii} = 4 \epsilon_{ii} \sigma_{ii}^{12}$ (单位：kJ mol⁻¹ nm¹²)
  * 此时 Lennard-Jones 势能通常写作：
  
    $$
    V_{LJ}(r) = \frac{C_{12,ij}}{r^{12}} - \frac{C_{6,ij}}{r^6}
    $$
  
* **如果 `comb-rule = 2` 或 `3`**:
  * $V_{ii}$ 直接代表 $\sigma_{ii}$ (单位：nm)
  * $W_{ii}$ 直接代表 $\epsilon_{ii}$ (单位：kJ mol⁻¹)
  * 此时 Lennard-Jones 势能通常写作：
  
    $$
    V_{LJ}(r) = 4 \epsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6\right]
    $$

#### 2. 交叉项参数的计算 (如果未在 `[nonbond_params]` 中显式定义)

* **对于 `comb-rule = 1` 和 `3`**:
  * GROMACS 使用**几何平均**来组合 $C_6$ 和 $C_{12}$ 参数：
    
    $$
    C_{6,ij} = \sqrt{C_{6,ii} \times C_{6,jj}}
    $$
    
    $$
    C_{12,ij} = \sqrt{C_{12,ii} \times C_{12,jj}}
    $$
    
  * 注意：如果 `comb-rule = 3`，`[atomtypes]` 中的 $V_{ii}$ 和 $W_{ii}$ 被解释为 $\sigma_{ii}$ 和 $\epsilon_{ii}$。GROMACS 内部会先将它们转换为 $C_{6,ii}$ 和 $C_{12,ii}$，然后再应用上述几何平均规则。
  
* **对于 `comb-rule = 2` (Lorentz-Berthelot 规则)**:
  
  * GROMACS 使用**算术平均**组合 $\sigma$ 参数，使用**几何平均**组合 $\epsilon$ 参数：
    
    $$
    \sigma_{ij} = \frac{\sigma_{ii} + \sigma_{jj}}{2}
    $$
    
    $$
    \epsilon_{ij} = \sqrt{\epsilon_{ii} \times \epsilon_{jj}}
    $$
    
    

### 关于常见力场的组合规则说明

**注意**：常见力场（CHARMM、AMBER、OPLS等）与 GROMACS 中 `comb-rule` 参数的对应关系在文献中并不十分明晰，以下信息基于有限的资料整理推测：

| 力场       | σ 组合规则 | ε 组合规则 | 可能的 GROMACS 设置 | 备注                                                         |
| ---------- | ---------- | ---------- | ------------------- | ------------------------------------------------------------ |
| **CHARMM** | 算术平均   | 几何平均   | `comb-rule = 2`     | 如果 `[atomtypes]` 中提供的是 $\sigma_{ii}$ 和 $\epsilon_{ii}$ |
| **AMBER**  | 算术平均   | 几何平均   | `comb-rule = 2`     | 明确使用 Lorentz-Berthelot 规则                              |
| **OPLS**   | 几何平均   | 几何平均   | `comb-rule = 3`     | 通常在 `[nonbond_params]` 中显式定义所有交叉项               |

> 算术平均是Lorentz提出的，几何平均是Berthelot提出的

也就是说，`comb-rule = 1`当然是万能的，但全原子一般是给出 $\sigma$ 和 $\epsilon$，其中`comb-rule = 2` 即Lorentz-Berthelot 规则，`comb-rule = 3` 即均为几何平均。

* **CHARMM**：使用 Lorentz-Berthelot 规则。对 $\sigma$ (或NAMD里面，等效的 $R_{min}$) 使用算术平均，对 $\epsilon$ 使用几何平均。
  $$
  R_{min,ij} = \frac{R_{min,ii} + R_{min,jj}}{2}  \text{(等效于 $\sigma$ 的算术平均)}
  $$

  $$
  \epsilon_{ij} = \sqrt{\epsilon_{ii} \times \epsilon_{jj}}
  $$
  
  * 参考：NAMD Mailing List - https://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2009-2010/3885.html
  
* **AMBER**：明确使用 Lorentz-Berthelot 规则。根据 AMBER 手册节选："For Amber force fields, cross terms involving different atom types i and j are evaluated according to the Lorentz/Berthelot mixing rules..."，可以自行查找

* **OPLS**：OPLS 力场通常对 Lennard-Jones 参数 $\sigma$ 和 $\epsilon$ **都使用几何平均**。
  
  * OPLS 力场的 GROMACS 实现通常没有 `[nonbond_params]` 。
  * NAMD参考文末
  

## 二、`[atomtypes]` 和 `[nonbond_params]` 中的参数解释

GROMACS 通过 `[atomtypes]` 和 `[nonbond_params]` (或 `[pairtypes]`) 这两个主要部分来定义非键相互作用参数。

参考：GROMACS Manual - Non-bonded parameters - https://manual.gromacs.org/current/reference-manual/topologies/parameter-files.html#non-bonded-parameters

### `[atomtypes]` 部分

此部分定义了每种原子类型**自身** (`ii`) 的基本非键参数。这些参数的解释（是 $\sigma$, $\epsilon$ 还是 $C_6$, $C_{12}$）**取决于** `[defaults]` 中设置的 `comb-rule`。

#### 示例 (OPLS-AA 风格，通常 `comb-rule = 1`，意味着 V, W 是 $C_6$, $C_{12}$)

```
[ atomtypes ]
;name   at.num    mass      charge   ptype          V(c6)          W(c12)      ; V 和 W 的含义取决于 comb-rule
O      8   15.99940      0.000       A     0.22617E-02    0.74158E-06   ; V(c6) = C6_ii, W(c12) = C12_ii
OM      8   15.99940      0.000       A     0.22617E-02    0.74158E-06
...
```

### `[nonbond_params]` 部分

此部分用于**显式定义**特定原子类型对 `i` 和 `j` 之间的非键相互作用参数。这里定义的参数将**覆盖任何通过组合规则计算得到的参数**。

#### 示例1 ( `comb-rule = 1` 配合，参数为直接的 $C_{6,ij}$ 和 $C_{12,ij}$)

```
[ nonbond_params ]
; i    j func          V(c6)          W(c12)    ; 列标题指明了是 C6 和 C12
O    O    1   0.22617E-02    0.74158E-06   ; O-O 相互作用的 C6_ij 和 C12_ij
O   OA    1   0.22617E-02    0.13807E-05   ; O-OA 相互作用的 C6_ij 和 C12_ij
...
```

* `V(c6)`：该原子类型对的 Lennard-Jones $C_{6,ij}$ 参数 (单位：kJ mol⁻¹ nm⁶)。
* `W(c12)`：该原子类型对的 Lennard-Jones $C_{12,ij}$ 参数 (单位：kJ mol⁻¹ nm¹²)。

#### 示例2 (Martini 风格，参数为直接的 $\sigma_{ij}$ 和 $\epsilon_{ij}$)

```
[ nonbond_params ]
; i     j     func   sigma      epsilon   ; 列标题通常会指明是 sigma 和 epsilon
P6    P6    1      0.470      4.990     ; P6-P6 相互作用的 sigma_ij 和 epsilon_ij
P6    P5    1      0.470      4.730     ; P6-P5 相互作用的 sigma_ij 和 epsilon_ij
...
```

* `i, j`：相互作用的原子类型。
* `func`：函数类型，`1` 表示 Lennard-Jones 12-6 势。
* `sigma`：该原子类型对的 Lennard-Jones $\sigma_{ij}$ 参数 (单位：nm)。
* `epsilon`：该原子类型对的 Lennard-Jones $\epsilon_{ij}$ 参数 (单位：kJ/mol)。

> **关键点**：`[nonbond_params]` 中参数的含义（是 $\sigma$/$\epsilon$ 还是 $C_6$/$C_{12}$）直接由该力场文件在该部分的**列定义**（通常通过注释中的列标题）决定。`func=1` 只是表示它是一个12-6型的Lennard-Jones势，但参数的表达形式可以有两种。

## 三、Martini 力场的特殊性

对于 Martini 力场 (例如 `martini_v3.0.0.itp`)：

参考文献：PCT Souza, et al., Nat. Methods, 2021. DOI：10.1038/s41592-021-01098-3 （看SI的表）

![image-20250531010245466](E:\GitHub-repo\mendelevium\_posts\Gromacs defaults in top file.assets\image-20250531010245466.png)

### `[defaults]` 指令

Martini 3 的主 `.itp` 文件通常包含：

```
[ defaults ]
; nbfunc      comb-rule 
  1           2           ; (通常 gen-pairs no, fudgeLJ/QQ 不适用或设为1.0)
```

这里的 `comb-rule = 2` 设定了默认的参数类型。

### `[atomtypes]` 部分（真实示例）

在 Martini 3 中，`[atomtypes]` 部分的 $\sigma$ 和 $\epsilon$ 值都设为 **0.0**，因为 **Martini 的核心在于珠子类型之间的相互作用矩阵**：

```
[ atomtypes ]
; name  mass    charge   ptype          sigma      epsilon
P6      72.0    0.000    A              0.0        0.0
P5      72.0    0.000    A              0.0        0.0
...
```

这里的 `sigma` 和 `epsilon` 都是 0.0，表明它们仅是占位符。

### `[nonbond_params]` 部分（真实示例）

这是 **Martini 力场定义非键相互作用的关键**。Martini **不依赖** GROMACS 的组合规则来生成不同珠子类型之间的相互作用参数。相反，它在 `[nonbond_params]` 部分**显式地定义**每一对珠子类型之间的 $\sigma_{ij}$ 和 $\epsilon_{ij}$：

```
[ nonbond_params ]
    P6    P6  1 4.700000e-01    4.990000e+00
    P6    P5  1 4.700000e-01    4.730000e+00
    ...
```

注意这里：
* 没有列标题注释，但根据 Martini 文档，这些参数是 $\sigma_{ij}$ (第4列) 和 $\epsilon_{ij}$ (第5列)
* 所有珠子对的相互作用都被显式定义

因此，当 `grompp` 处理 Martini 拓扑时，它会**优先使用** `[nonbond_params]` 中为特定珠子对定义的 $\sigma_{ij}$ 和 $\epsilon_{ij}$。只有当某一对珠子类型的相互作用没有在 `[nonbond_params]` 中显式定义时，才会退回到使用 `[defaults]` 中指定的 `comb-rule` 和 `[atomtypes]` 中的参数来尝试计算（但由于 `[atomtypes]` 中的值都是 0.0，实际上不会产生有意义的相互作用）。

详见上一篇：

### 总结

对于标准的 Martini 3 力场文件：

1. **`[atomtypes]` 中的 $\sigma$/$\epsilon$ 都是 0.0**：它们是占位符，不用于计算。

2. **核心的异类珠子对相互作用参数来自 `[nonbond_params]`**：这是Martini设计的核心。

3. **`[nonbond_params]` 中提供的是针对特定珠子对 `ij` 的 $\sigma_{ij}$ 和 $\epsilon_{ij}$**：这些不是 $C_{6,ij}$ 和 $C_{12,ij}$。

4. **`[defaults]` 中的 `comb-rule = 2` 在 Martini 中更多的是一个形式上的设定**：因为所有相关的珠子对相互作用参数都是在 `[nonbond_params]` 中显式提供的。

## 四、总结：如何判断参数类型

判断 `.itp` 文件中非键参数是 ($\sigma$, $\epsilon$) 还是 ($C_6$, $C_{12}$) 的**关键步骤**：

### 1. 查看 `[defaults]` 指令中的 `comb-rule`

* 如果 **`comb-rule = 1`**，那么 `[atomtypes]` 中的 `V` 和 `W` 列倾向于是 **$C_{6,ii}$** 和 **$C_{12,ii}$**。
* 如果 **`comb-rule = 2` 或 `3`**，那么 `[atomtypes]` 中的 `V` 和 `W` 列倾向于是 **$\sigma_{ii}$** 和 **$\epsilon_{ii}$**。

### 2. 仔细阅读 `[atomtypes]` 和 `[nonbond_params]` 部分的列标题注释

* 如果列标题明确写着 **`sigma` 和 `epsilon`**，那么这些值就是 $\sigma$ 和 $\epsilon$。
* 如果列标题明确写着 **`V(c6)` 和 `W(c12)`**，那么这些值就是 $C_6$ 和 $C_{12}$。
* 假定开发者不至于搞错，**这是最直接的判断依据。**

### 3. 查阅相应力场的原始文献和手册

这是**最权威**的判断依据。力场开发者会明确说明其参数的定义和使用方式。

### 实用建议

对于您的脚本而言，如果它需要同时处理可能来自不同力场的 `.itp` 文件，建议：

* 通过一个参数来指定当前处理的ITP文件中的非键参数是哪种类型
* 或者通过智能解析列标题来判断
* 对于 Martini 这样的特殊情况（`[atomtypes]` 中都是 0.0），直接使用 `[nonbond_params]` 中的参数

## 其他参考资料

GROMACS Manual - MDP Options for LJ-PME combination rule:
https://manual.gromacs.org/current/user-guide/mdp-options.html#mdp-lj-pme-comb-rule

### NAMD 的相关设置

* NAMD Mailing List：https://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2009-2010/3885.html
  "Yes, as is standard for the CHARMM force field NAMD uses arithmetic mean for sigma, geometric mean for epsilon by default. You can change this by **adding 'vdwGeometricSigma yes' in the config file to support, e.g., OPLS.**"

* NAMD User Guide：https://www.ks.uiuc.edu/Research/namd/3.0.1/ug/node25.html#7012
  "vdwGeometricSigma：Use geometric mean, as required by OPLS, rather than traditional arithmetic mean when combining Lennard-Jones sigma parameters for different atom types."

