---
title: "理解结构因子S(q)及其在聚电解质相分离研究中的应用"
date: "2025-10-08"
tags: [structure-factor, scattering, phase-separation, polyelectrolyte, physics, theory, computational-physics]
description: "详解结构因子S(q)的物理意义与计算方法，从X射线散射到分子动力学模拟，深入分析其在聚电解质液-液相分离研究中的应用，揭示相分离动力学的特征长度和时间演化"
thumbnail: "/assets/img/thumbnail_mine/wh-dpe6lm.jpg"
image: "/assets/img/thumbnail_mine/wh-dpe6lm.jpg"
---
# 理解结构因子 $S(q)$ 及其在聚电解质相分离研究中的应用

本文的主要参考文献为[^1,2]，内容由AI生成，如有错误恳请指出。


## 一、结构因子 $S(q)$ 理论与计算详解

### 1.1 什么是结构因子 $S(q)$？为什么要计算它？

结构因子，通常表示为 $S(\mathbf{q})$，是一个关键的物理量，用于描述材料内部原子或分子在不同空间尺度上的**密度不均匀性或有序性**。如果材料是完全均匀的，那么各个点的密度都一样；但实际材料总会有涨落。结构因子就是用来量化这种不均匀程度的。

在实验上，结构因子可以通过X射线、电子衍射和中子衍射等得到[^3]（注：我们主要讨论的是 $S(\mathbf{q})$ ）。当一束波（X射线、中子、光）入射到材料上时，会因为材料内部的密度不均匀而发生散射。测量到的**散射强度 $I(\mathbf{q})$ 与结构因子 $S(\mathbf{q})$ 直接相关**（通常是成正比，$I(\mathbf{q}) \propto S(\mathbf{q})$）。

这里的 $\mathbf{q}$ 是**散射波矢**（scattering wavevector）。它联系了实验测量的散射角与我们关心的材料内部结构尺度。
* $\mathbf{q}$ 的**方向**与入射波和散射波方向的差异有关，反映了我们探测的结构在空间中的取向。
* 其**模长** $q = |\mathbf{q}|$ 的大小与我们探测的**空间尺度 $l$ 成反比**，通常可以近似认为 $l \approx 2\pi/q$。这意味着：
    * **小的 $q$ 值对应大的空间尺度**（例如，大的团簇、相畴尺寸）。
    * **大的 $q$ 值对应小的空间尺度**（例如，原子间距、键长）。

通过分析散射强度随 $q$ 的变化，我们就能反推出材料在不同长度尺度上的结构信息。例如，如果 $S(q)$ 在某个 $q_0$ 处出现峰值，则表明体系中存在一个以 $l_0 \approx 2\pi/q_0$ 为特征长度的显著结构。

对于**各向同性的体系**（如液体、无序的聚合物溶液或粉末样品），其内部结构在所有方向上统计平均是相同的。因此，**结构因子仅依赖于波矢的模长** $q$，可以简化记作 $S(q)$。

---
#### 1.1.1 结构因子一阶矩 $\langle q \rangle (t)$ 的物理意义

**特征波数 $\langle q \rangle (t)$ (Characteristic Wavenumber)** 是用来定量表征相分离过程中结构特征尺寸的一个关键物理量。

* **结构因子 $S(q,t)$ 的一阶矩**: 论文中明确指出，$\langle q \rangle$ 是通过含时结构因子 $S(q,t)$ 的一阶矩来定义的。其计算公式为：
    $$
    \langle q \rangle (t) = \frac{\int_0^\infty q S(q,t) \, dq}{\int_0^\infty S(q,t) \, dq}
    $$
    
    这个公式的意义是：对所有波数 $q$ 进行积分（或在离散数据中求和），每个 $q$ 的"权重"是其对应的结构因子强度 $S(q,t)$。分子是加权后的波数总和，分母是总的结构因子强度（起到归一化作用）。
    
* **倒易空间中的平均特征尺度**: $\langle q \rangle$ 作为 $S(q,t)$ 的加权平均波数值，反映了在时刻 $t$ 体系中**占主导地位的结构特征**（如网络的平均线宽或孔洞大小，或者液滴的平均尺寸）所对应的平均波数值。

* **与特征长度成反比**: 波矢 $q$ 与实空间中的特征长度尺度 $l$ 成反比关系，通常可以认为 $l =2\pi/\langle q \rangle$。因此，特征波数 $\langle q \rangle$ 的**减小**直接反映了实空间中相畴（网络或液滴）特征尺寸 $l$ 的**增大**。

* **描述畴粗化过程**: 在相分离后期，小的相畴会逐渐合并变大，这个过程称为**"畴粗化" (Domain Coarsening)**。在这个过程中，特征长度 $l$ 会随时间 $t$ 增大，因此，特征波数 $\langle q \rangle$ 会随时间 $t$ 减小。通过追踪 $\langle q \rangle$ 随时间的变化，可以定量地研究畴粗化的动力学过程及其标度行为。

### 1.2 从密度涨落到结构因子：数学定义

#### 1.2.1 瞬时粒子密度及其傅里叶分量

要从微观层面理解材料或流体的结构，我们首先需要一种描述体系中粒子空间分布的方法。考虑一个包含 $N$ 个粒子的体系，在某个特定时刻 $t$，其**瞬时单粒子密度** $\rho(\mathbf{r},t)$ 可以被精确地表示为体系中所有粒子在各自位置 $\mathbf{r}_i(t)$ 上的贡献之和。数学上，这通常通过**狄拉克 $\delta$ 函数**或**双曲正切函数映射**来实现：

$$
\rho(\mathbf{r},t) = \sum_{i=1}^{N} \delta(\mathbf{r} - \mathbf{r}_i(t))
$$

这里的 $\delta(\mathbf{r} - \mathbf{r}_i(t))$ 是一个三维狄拉克 $\delta$ 函数。它的核心特性是：当 $\mathbf{r} = \mathbf{r}_i(t)$ 时，其值为无穷大；而当 $\mathbf{r} \neq \mathbf{r}_i(t)$ 时，其值为零。然而，它在整个空间的积分为1，即 
$$
\int \delta(\mathbf{r} - \mathbf{r}_i(t)) d\mathbf{r} = 1
$$
因此，这个表达式的物理意义是：在粒子 $i$ 所在的位置 $\mathbf{r}_i(t)$ 处密度是无穷集中的，而在其他任何没有粒子的地方密度为零。

直接在实空间处理 $\rho(\mathbf{r},t)$ 来分析跨越不同空间尺度的结构特征（例如，从单个粒子的大小到宏观聚集体的尺寸）往往非常复杂。为了更有效地揭示这些结构信息，我们通常将其转换到**倒易空间**（reciprocal space），也常被称为**傅里叶空间**或 **q-空间**。这一转换通过**傅里叶变换**完成，它将实空间中的密度函数 $\rho(\mathbf{r},t)$ 分解为一系列具有不同波矢 $\mathbf{q}$ 的平面波（也称为密度波）的线性叠加。每个波矢 $\mathbf{q}$ 对应着一个特定的空间尺度（波长 $\lambda = 2\pi/|\mathbf{q}|$）和方向。

密度场 $\rho(\mathbf{r},t)$ 在特定波矢 $\mathbf{q}$ 上的**傅里叶分量** $\rho_{\mathbf{q}}(t)$（也常被称为密度涨落的傅里叶模式）定义为：

$$
\rho_{\mathbf{q}}(t) = \int_{\text{box}} e^{-i\mathbf{q}\cdot\mathbf{r}} \rho(\mathbf{r},t) \, d\mathbf{r}
$$

其中积分在体系的整个体积（"box"）上进行，$i$ 是虚数单位。将上面瞬时密度的表达式代入此定义，我们可以得到 $\rho_{\mathbf{q}}(t)$ 的一个更直接的计算形式：

$$
\rho_{\mathbf{q}}(t) = \int_{\text{box}} e^{-i\mathbf{q}\cdot\mathbf{r}} \left( \sum_{j=1}^{N} \delta(\mathbf{r} - \mathbf{r}_j(t)) \right) \, d\mathbf{r}
$$

利用狄拉克 $\delta$ 函数的筛选性质（即 $\int f(\mathbf{r})\delta(\mathbf{r}-\mathbf{a})d\mathbf{r} = f(\mathbf{a})$），上式简化为：

$$
\rho_{\mathbf{q}}(t) = \sum_{j=1}^{N} e^{-i\mathbf{q}\cdot\mathbf{r}_j(t)}
$$

$\rho_{\mathbf{q}}(t)$ 是一个复数。它的**模长** $| \rho_{\mathbf{q}}(t) |$ 反映了体系在波矢 $\mathbf{q}$ 所对应的空间尺度和方向上密度起伏的**幅度**，而它的**相位**则给出了这些密度波的相对位置信息。

#### 1.2.2 结构因子 $S(q)$ 的定义

有了密度的傅里叶分量，我们就可以定义**结构因子**，它是表征材料平均结构的关键物理量。

**静态结构因子 $S(\mathbf{q})$** 通常被定义为密度傅里叶分量的均方涨落，并除以粒子总数 $N$ 进行归一化：
$$
S(\mathbf{q}) = \frac{1}{N} \langle \rho_{\mathbf{q}}(t) \rho_{-\mathbf{q}}(t) \rangle
$$

其中 $\langle \dots \rangle$ 表示对系统进行**系综平均**（例如，在平衡态下对所有可能的微观状态进行平均）或在足够长的时间内进行**时间平均**。
我们注意到 $\rho_{-\mathbf{q}}(t)$ 与 $\rho_{\mathbf{q}}(t)$ 的复共轭 $\rho_{\mathbf{q}}^*(t)$ 之间存在一个简单的关系。其推导如下：
$$
\rho_{-\mathbf{q}}(t) = \sum_{j=1}^{N} e^{-i(-\mathbf{q})\cdot\mathbf{r}_j(t)} = \sum_{j=1}^{N} e^{i\mathbf{q}\cdot\mathbf{r}_j(t)}
$$

另一方面，$\rho_{\mathbf{q}}(t)$ 的复共轭是：
$$
\rho_{\mathbf{q}}^*(t) = \left(\sum_{j=1}^{N} e^{-i\mathbf{q}\cdot\mathbf{r}_j(t)}\right)^* = \sum_{j=1}^{N} (e^{-i\mathbf{q}\cdot\mathbf{r}_j(t)})^* = \sum_{j=1}^{N} e^{i\mathbf{q}\cdot\mathbf{r}_j(t)}
$$
因此，我们得到 $\rho_{-\mathbf{q}}(t) = \rho_{\mathbf{q}}^*(t)$。

利用这个关系，静态结构因子可以更直观地写成：
$$
S(\mathbf{q}) = \frac{1}{N} \langle |\rho_{\mathbf{q}}(t)|^2 \rangle
$$
这个定义清晰地表明，$S(\mathbf{q})$ **衡量了在波矢 $\mathbf{q}$ （即特定尺度和方向）上密度涨落的平均强度**。在实验中（如X射线或中子散射），$S(\mathbf{q})$ 与散射强度直接相关。$S(\mathbf{q})$ 的峰值位置揭示了体系中占主导地位的特征长度或周期性结构。

在研究动态过程（例如相分离动力学、玻璃化转变等）时，我们更关心的是结构如何随时间演化。这时，**含时结构因子 $S(\mathbf{q}, t)$** 成为一个重要的分析工具。在 Yuan & Tanaka (2025) 的研究中 [^1]，它被类似地定义（通常，如果体系是各向同性的，或者我们只关心不同尺度上的平均结构演化，结构因子可以只表示为波矢模长 $q = |\mathbf{q}|$ 的函数，此时 $S(q,t)$ 是对所有方向的 $\mathbf{q}$ 进行平均的结果）：
$$
S(q,t) = \frac{\langle \rho_q(t) \rho_{-q}(t) \rangle}{N}
$$
与上式是等价的。

$S(q,t)$ **描述了在特定空间尺度 $q$ 上的结构特征如何随时间 $t$ 演变**。例如，在相分离过程中，特征峰的位置会向更小的 $q$ 值移动（对应更大的结构尺寸），峰高也会增加。

### 1.3 高分子体系中结构因子的计算细节

在分子动力学 (MD) 或粗粒化 (CG) 模拟中计算 $S(q,t)$ 时，通常涉及以下步骤：

1.  **粒子选择**:
    * **原子级模拟**：通常选择重原子或分子的质心 (COM)。
    * **粗粒化模拟**：选择代表性的粗粒化珠子 (bead)。

2.  **密度场计算**:
    * 将选定粒子的坐标映射到三维网格上，得到离散的密度场 $\rho(\mathbf{r},t)$。
    * 可使用**高斯平滑**或其他平滑函数（如论文[1]中提到的双曲正切函数）处理点粒子密度，获得更连续的密度场。

3.  **傅里叶变换**:
    * 使用**快速傅里叶变换 (FFT)** 算法计算离散密度场的傅里叶分量 $\rho_{\mathbf{q}}(t)$。
    * 通常会去除 $\mathbf{q}=0$ 的分量（直流分量），因为它代表体系的平均密度。

4.  **计算 $S(q,t)$**:
    * 根据 $S(\mathbf{q},t) = \frac{1}{N} |\rho_{\mathbf{q}}(t)|^2$ 计算。
    * 对于各向同性体系，进行**球面平均 (spherical averaging)**，得到仅依赖于 $q$ 的标量函数 $S(q,t)$。

5.  **时间平均或演化**:
    * 对于**静态结构因子 $S(q)$**，需对多个时间步或独立轨迹的 $S(q,t)$ 进行平均。
    * 研究**动力学过程**时，观察 $S(q,t)$ 或其导出量随时间 $t$ 的演化。

------

## 二、聚电解质相分离研究中的结构因子与特征波数[^1]

### 2.1 引言概述：从传统认知到新发现的科学突破

在生物体系和材料科学中，相反电荷聚电解质（PEs）通过相分离形成的凝聚层（coacervates）扮演着至关重要的角色。这些凝聚层不仅是理解生物凝聚体（如无膜细胞器）形成机制的关键，也为开发响应性智能材料提供了新思路。

**传统观点的局限性**：长期以来，科学界普遍认为聚电解质凝聚层主要形成**球形液滴**，其生长动力学遵循经典的液-液相分离（LLPS）机制。在这种传统框架下，凝聚层被视为简单的液滴，**通过蒸发-凝结或碰撞-合并等机制长大**，其特征尺寸遵循 $l \propto t^{1/3}$ 的生长规律。

**革命性的新发现**：然而，Yuan & Tanaka (2025) 通过包含流体动力学相互作用（HI）和静电相互作用的流体粒子动力学（FPD）模拟，彻底颠覆了这一传统认知。他们的研究揭示了一个惊人的现象：即使在半稀溶液中（体积分数仅约2.3%），相反电荷的聚电解质也能自发形成**贯通的网络结构**，而非传统认为的孤立液滴。

**独特的生长规律**：更令人瞩目的是，该网络结构在粗化过程中遵循一个**独特的生长规律 $l \propto t^{1/2}$**（其中 $l$ 是特征长度，$t$ 是时间）。这种自相似的生长行为在中性聚合物的不良溶剂体系中通常不存在。其背后的物理机制源于良溶剂中的聚电解质在整体电中性的约束下，由于空间电荷的不均匀性，表现出**更弱但更长程的有效吸引力，导致形成的聚电解质富集相密度较低（约40%），界面张力也显著降低**。

**研究的核心科学问题**：

1. **相形态的决定因素**：在何种条件下会发生液滴状或网络状的相分离？初始状态、体积分数、链长等因素如何影响最终形态？

2. **畴粗化的自相似性**：网络状相分离的畴粗化过程是否存在自相似性？其背后的物理机制是什么？

3. **静电相互作用的独特作用**：静电荷及其对称性对网络状相分离有何影响？与中性聚合物体系有何本质区别？

**研究的重要意义**：这项研究不仅挑战了我们对聚电解质凝聚层形成机制的基本认识，还为理解生物体系中的网络状凝聚体（如中心体组装、蛋白质颗粒等）提供了新的理论基础。同时，通过调控电荷不对称性来稳定网络结构的发现，为设计新型多孔材料和生物响应材料开辟了新途径。

### 2.2 模拟参数说明：$\sigma$、$\tau_{BD}$ 及无量纲化处理

在Yuan & Tanaka (2025) 的粗粒化模拟研究中，为了使结果具有普适性并便于比较，采用了无量纲化的处理方法。

#### 2.2.1 基本长度单位 $\sigma$

**$\sigma$ (sigma)** 代表粗粒化模型中单体（monomer）或离子（ion）的**直径**。论文中设定：

$$
\sigma = 0.72 \text{ nm}
$$

这个尺度对应于典型的水合离子直径，被用作基本的长度单位。在模拟中，所有的长度量都以 $\sigma$ 为单位进行无量纲化处理。

#### 2.2.2 布朗时间 $\tau_{BD}$

**$\tau_{BD}$** 是**布朗时间 (Brownian time)**，代表一个粒子由于热运动扩散其自身直径 $\sigma$ 距离所需的特征时间尺度。根据Stokes-Einstein关系，布朗时间定义为：
$$
\tau_{BD} = \frac{\pi \sigma^3 \eta}{8 k_B T}
$$

其中：
- $\eta$ 是溶剂粘度
- $k_B$ 是玻尔兹曼常数
- $T$ 是绝对温度

对于室温下的水（$\eta \approx 10^{-3}$ Pa·s），计算得到：

$$
\tau_{BD} \approx 0.035 \text{ ns}
$$

这意味着模拟的时间尺度从1微秒到10微秒不等，这在原子级模拟中是极具挑战性的。

#### 2.2.3 无量纲特征波数 $\langle q \rangle \sigma / 2\pi$

在图二中，y轴显示的是无量纲化的特征波数 $\langle q \rangle \sigma / 2\pi$。这个量的物理意义可以通过以下推导理解：

由于特征长度 $l \approx 2\pi / \langle q \rangle$，我们有：

$$
\frac{\langle q \rangle \sigma}{2\pi} = \frac{\sigma}{2\pi / \langle q \rangle} = \frac{\sigma}{l}
$$

因此，**$\langle q \rangle \sigma / 2\pi$ 表示的是单体直径 $\sigma$ 与体系特征长度 $l$ 的比值**。

- 当相畴较小时，$l$ 小，$\langle q \rangle \sigma / 2\pi$ 大
- 随着相畴粗化，$l$ 增大，$\langle q \rangle$ 减小，$\langle q \rangle \sigma / 2\pi$ 随时间减小

这种无量纲化处理使得不同参数条件下的结果可以在同一坐标系中进行比较，揭示普适的标度行为。

### 2.3 双对数坐标图分析与核心发现

#### 2.3.1 为什么使用双对数坐标图？

论文中图二 (Fig. 2) 将 $\langle q \rangle \sigma / 2\pi$ 对 $t/\tau_{BD}$ 绘制在双对数坐标上。这种作图方式的主要目的是**检验数据是否满足幂律关系 (Power Law)**，即形如 $y = A \cdot x^m$ 的关系。若满足幂律关系，在双对数图上数据点会落在一条直线上，其**斜率 (slope) 即为幂指数 $m$**。

#### 2.3.2 图二的关键发现

![image-20250603123852532](E:\GitHub-repo\mendelevium\_posts\structure-factor.assets\image-20250603123852532.png)

图 2 | 网络形成相分离过程中的区域粗化和时间尺度表征。

**a, c, e**: 在流体粒子动力学 (FPD) 模拟中，不同比耶鲁姆长度（Bjerrum length）lB=1.1σ (a)、lB=2σ (c) 和 lB=3σ (e) 条件下，特征波数 ⟨q⟩（定义为结构因子 S(q,t) 的一阶矩）随时间的演化过程 。

**b, d, f**: 通过布朗动力学 (BD) 模拟得到的结果 。误差棒代表了根据四次独立模拟计算得到的标准误差 。在电荷对称条件 (Na=Nc=40) 下，区域粗化过程在 FPD 模拟中遵循 ⟨q⟩∼t−1/2 的规律，而在 BD 模拟中则遵循 ⟨q⟩∼t−1/3 的规律 。

**g**: 二元带电聚电解质（PE）溶液在链长为 (Nc,Na)=(40,40)、比耶鲁姆长度分别为 lB=2σ（体积分数 ϕ≈0.38）和 lB=3σ（体积分数 ϕ≈0.42）时，其致密相的自中间散射函数 Fs(q,t) 。其中，q 选为结构因子 S(q) 第一个峰值对应的波数 。S(q) 和 Fs(q,t) 的定义参见“方法”部分。结构弛豫时间 τα 定义为 Fs(q,t) 衰减到 1/e 时的时间 。我们发现 τα≈70∼100τBD 。

**h, i**: 在比耶鲁姆长度分别为 lB=2σ (h) 和 lB=3σ (i) 条件下，相同聚电解质溶液的整体变形和剪切变形特征时间尺度（应变率的倒数，Δt/∣ϵbulk∣ 和 Δt/∣ϵshear∣）随时间的变化 。估算的区域变形时间尺度 τdef 约为 5∼10τBD 。这些结果表明，以 τα 为特征的粒子重排过程慢于区域变形过程，这说明网络粗化过程是由机械弛豫控制的 。

图二展示了不同条件下特征波数的时间演化：

1. **电荷对称条件 (Na = Nc = 40)**：
   - **FPD模拟（含HI）**：图 2a, c, e 中，数据点在双对数图上呈现良好的线性关系，斜率为 **-1/2**
   - 这表示 $\langle q \rangle \propto t^{-1/2}$，即特征长度 $l \propto t^{1/2}$
   - 这个标度关系在不同的Bjerrum长度（$l_B = 1.1\sigma$ 到 $3\sigma$）下保持一致

   - **BD模拟（不含HI）**：图 2b, d, f 中，斜率为 **-1/3**
   - 表示 $\langle q \rangle \propto t^{-1/3}$，即特征长度 $l \propto t^{1/3}$
   - 虽然也形成网络结构，但粗化动力学不同

2. **流体动力学相互作用的关键作用**：
   - **HI对实现 $t^{1/2}$ 幂律至关重要**
   - 这一发现强调了在模拟聚电解质相分离时包含流体动力学效应的必要性

3. **自相似性的体现**：
   - 幂律关系的存在通常意味着系统在粗化过程中表现出**自相似性 (self-similarity)**
   - 这种自相似性在图3b中得到进一步验证：不同时刻的标度结构因子塌缩到同一主曲线上

### 2.4 物理机制解析：粘弹性相分离

通过对特征波数 $\langle q \rangle (t)$ 演化的分析，结合对弛豫时间尺度的比较（图2g-i），Yuan & Tanaka揭示了网络形成的物理机制：

1. **动力学不对称性**：
   - 结构弛豫时间 $\tau_\alpha \approx 70-100\tau_{BD}$
   - 畴变形时间 $\tau_{def} \approx 5-10\tau_{BD}$
   - 由于 $\tau_\alpha \gg \tau_{def}$，密集相中的粒子重排跟不上快速的畴变形

2. **粘弹性相分离 (VPS)**：
   - 这种动力学不对称性激活了粘弹性效应
   - 导致形成瞬态网络结构，而非传统的液滴
   - 网络粗化由其力学弛豫控制，该弛豫受限于溶剂在网络中的渗透流动（孔隙弹性弛豫）

3. **与中性聚合物的区别**：
   - 聚电解质在良溶剂中的有效吸引力较弱
   - 形成的富集相密度较低（约40%，而中性聚合物约50%）
   - 这种较松散的堆积促进了局部键弛豫，维持了自相似生长

### 2.5 电荷不对称的影响

当引入电荷不对称（如 Nc = 50, Na = 30）时：

1. **粗化动力学显著减慢**：
   - 偏离 $t^{-1/2}$ 幂律
   - 后期出现动力学慢化趋势

2. **物理机制**：
   - 网络表面积累净电荷（图3d, 4b, 4d）
   - 静电排斥阻碍进一步粗化
   - 与中性聚合物的VPS不同，电荷不对称可以稳定网络结构

3. **应用前景**：
   - 通过调节电荷不对称性可控制网络稳定性
   - 为设计稳定的多孔材料提供新途径

------

## 三、实操指南：从模拟数据计算 S(q) 和拟合幂律指数

### 3.1 Python 代码实现与详细解读

以下Python函数展示了如何使用MDAnalysis和SciPy/NumPy从模拟轨迹计算 $S(q)$ 和特征波数 $\langle q \rangle$：

```python
def calculate_structure_factor(
    u: mda.Universe,
    frame_index: int,
    selection: str,
    n_bins: int = 64,
    q_max_factor: float = 0.5,  # 计算到 q_max = q_max_factor * Nyquist频率
    density_method: str = 'histogram'
    ) -> tuple[np.ndarray | None, np.ndarray | None, float | None]:
    """
    计算特定帧和原子选择的静态结构因子 S(q) 和特征波数 <q>
    
    参数:
        u (mda.Universe): MDAnalysis Universe对象，包含轨迹
        frame_index (int): 要分析的帧索引
        selection (str): MDAnalysis选择字符串（如 'resname HA and name A'）
        n_bins (int): 密度网格每个维度的格子数，默认64
        q_max_factor (float): 计算q的最大值相对于Nyquist频率的比例
        density_method (str): 密度计算方法，目前仅支持'histogram'
    
    返回:
        tuple: (q_bin_centers, S_q_radially_averaged, char_q)
               - q_bin_centers: q值的数组
               - S_q_radially_averaged: 对应的S(q)值
               - char_q: 特征波数<q>
    """
    if density_method != 'histogram':
        raise NotImplementedError("Only 'histogram' density method is currently supported.")

    try:
        # 确保轨迹定位到正确的帧
        # 这在重复调用时很关键
        u.trajectory[frame_index]
    except IndexError:
        print(f"Error: Frame index {frame_index} is out of bounds.")
        return None, None, None

    # --- 1. 选择原子并获取盒子尺寸 ---
    ag = u.select_atoms(selection)
    N = len(ag)
    if N == 0:
        return None, None, None

    from scipy.fft import fftn, fftshift, fftfreq
    from scipy import stats  # 用于径向平均的binned_statistic
    
    coords = ag.positions
    # 假设是正交盒子，从dimensions属性获取
    box_dims = u.dimensions[:3]
    if box_dims is None or np.any(box_dims <= 0):
         print(f"Error: Invalid box dimensions {box_dims} at frame {frame_index}.")
         return None, None, None

    # --- 2. 计算密度场 (rho_r) ---
    # 使用3D直方图将粒子坐标转换为密度场
    ranges = [[0, L] for L in box_dims]
    try:
        rho_r, edges = np.histogramdd(
            coords,
            bins=n_bins,
            range=ranges,
            density=False  # 获取计数，而非概率密度
        )
    except ValueError as e:
        print(f"Error during histogramming for frame {frame_index}: {e}")
        return None, None, None

    delta_xyz = box_dims / n_bins  # 每个格子的尺寸

    # --- 3. 计算 S(q) 网格 ---
    # 对密度场进行FFT得到傅里叶分量
    rho_q = fftn(rho_r)
    # S(q) = |rho_q|^2 / N
    S_q_grid = (np.abs(rho_q)**2) / N if N > 0 else np.zeros_like(rho_r, dtype=float)

    # --- 4. 计算 q 向量和模长 ---
    # fftfreq给出归一化的频率，需要乘以2π/d得到波矢
    qx = 2 * np.pi * fftfreq(n_bins, d=delta_xyz[0])
    qy = 2 * np.pi * fftfreq(n_bins, d=delta_xyz[1])
    qz = 2 * np.pi * fftfreq(n_bins, d=delta_xyz[2])
    # 创建3D网格
    qxg, qyg, qzg = np.meshgrid(qx, qy, qz, indexing='ij')
    q_magnitude_grid = np.sqrt(qxg**2 + qyg**2 + qzg**2)

    # --- 5. 径向平均 ---
    # 将FFT结果移动到中心（低频在中心）
    S_q_grid_shifted = fftshift(S_q_grid)
    q_magnitude_grid_shifted = fftshift(q_magnitude_grid)
    
    # 展平为1D数组以便进行统计
    q_values_flat = q_magnitude_grid_shifted.ravel()
    S_q_values_flat = S_q_grid_shifted.ravel()

    # 确定q的范围和分辨率
    q_min_res = np.min([np.min(np.abs(qi[qi!=0])) for qi in [qx, qy, qz] if np.any(qi!=0)]) if np.any([np.any(qi!=0) for qi in [qx, qy, qz]]) else 0.01
    q_nyquist = np.min(np.pi / delta_xyz) if np.all(delta_xyz > 0) else 1.0
    q_max_calc = q_max_factor * q_nyquist
    delta_q = q_min_res / 2.0
    if delta_q <= 0:
        delta_q = q_max_calc / (n_bins // 2) if q_max_calc > 0 else 0.01

    # 创建q的bins用于径向平均
    if q_max_calc <= delta_q:
         q_bins = np.array([0, q_max_calc + delta_q]) if q_max_calc > 0 else np.array([0, 0.1])
    else:
        q_bins = np.arange(0, q_max_calc + delta_q, delta_q)

    # 对每个q区间内的S(q)值求和
    S_q_sum, _, binnumber = stats.binned_statistic(
        q_values_flat, S_q_values_flat, statistic='sum', bins=q_bins
    )
    # 计算每个区间内的点数
    counts, _, _ = stats.binned_statistic(
        q_values_flat, q_values_flat, statistic='count', bins=q_bins
    )
    # 径向平均 = 总和 / 计数
    S_q_radially_averaged = np.divide(S_q_sum, counts, out=np.zeros_like(S_q_sum), where=counts != 0)
    q_bin_centers = (q_bins[:-1] + q_bins[1:]) / 2

    # --- 6. 计算特征波数 <q> ---
    # <q> = ∫q·S(q)dq / ∫S(q)dq
    if len(q_bin_centers) > 1:
        # 排除q=0的点（通常对应均匀背景）
        q_relevant = q_bin_centers[1:]
        S_q_relevant = S_q_radially_averaged[1:]
        
        # 只考虑S(q)显著大于0的点
        valid_indices = np.where(S_q_relevant > 1e-9)[0]
        if len(valid_indices) > 0:
            q_relevant = q_relevant[valid_indices]
            S_q_relevant = S_q_relevant[valid_indices]
            
            # 计算一阶矩
            numerator = np.sum(q_relevant * S_q_relevant)
            denominator = np.sum(S_q_relevant)
            char_q = numerator / denominator if denominator > 0 else np.nan
        else:
            char_q = np.nan
    else:
        char_q = np.nan

    return q_bin_centers, S_q_radially_averaged, char_q


def calculate_sq_trajectory(
    u: mda.Universe,
    selection: str = 'resname HA and name A',
    n_bins: int = 64,
    start_frame: int = 0,
    stop_frame: int | None = None,
    step: int = 1,
    show_progress: bool = True,
    **kwargs  # 传递额外参数给calculate_structure_factor
    ) -> np.ndarray:
    """
    计算整个轨迹的特征波数 <q> 随时间的演化
    
    通过对每个指定帧调用 calculate_structure_factor 并收集特征波数
    
    参数:
        u (mda.Universe): MDAnalysis Universe对象
        selection (str): 原子选择字符串
        n_bins (int): 密度网格的bins数
        start_frame (int): 起始帧索引
        stop_frame (int | None): 结束帧索引（不包含）
        step (int): 帧间隔
        show_progress (bool): 是否显示进度条
        **kwargs: 传递给calculate_structure_factor的额外参数
    
    返回:
        np.ndarray: 包含每帧特征波数<q>的数组
    """
    all_char_q = []
    n_frames_total = len(u.trajectory)
    
    if stop_frame is None:
        stop_frame = n_frames_total
    else:
        stop_frame = min(stop_frame, n_frames_total)  # 确保不超过轨迹长度
    
    frame_indices = range(start_frame, stop_frame, step)

    # 设置进度条
    iterator = frame_indices
    if show_progress:
        try:
            # 尝试自动检测是否在notebook环境
            if 'ipykernel' in str(type(getattr(__builtins__, '__dict__', {}).get('get_ipython'))):
                 from tqdm.notebook import tqdm
            else:
                 from tqdm import tqdm
            iterator = tqdm(frame_indices, desc="Calculating <q> per frame")
        except ImportError:
            print("tqdm library not found. Progress bar disabled.")

    # 遍历指定的帧
    for frame_idx in iterator:
        q_bins, S_q, char_q = calculate_structure_factor(
            u=u,
            frame_index=frame_idx,
            selection=selection,
            n_bins=n_bins,
            **kwargs  # 传递额外参数如q_max_factor
        )
        all_char_q.append(char_q if char_q is not None else np.nan)

    return np.array(all_char_q)
```

### 3.2 代码解读要点

1. **密度场构建**：
   - 使用3D直方图将离散的粒子坐标转换为连续的密度场
   - 格子大小影响q空间的分辨率和最大可探测的q值

2. **FFT计算**：
   - 使用快速傅里叶变换计算密度场的傅里叶分量
   - $S(q) = |\rho_q|^2 / N$ 给出了每个q模式的强度

3. **径向平均**：
   - 对于各向同性体系，将3D的S(q)数据按q的模长进行平均
   - 使用binned_statistic高效实现

4. **特征波数计算**：
   - 排除q=0的贡献（对应均匀背景）
   - 只考虑S(q)显著的区域，避免噪声影响

### 3.3 幂律拟合实操指南

当您追踪 $\langle q \rangle (t)$ 并希望通过线性拟合其双对数图来确定幂律指数 $m$ (即 $\langle q \rangle \propto t^m$) 时，以下是重要的实操考虑：

1. **观察双对数图**：
   ```python
   import numpy as np
   import matplotlib.pyplot as plt
   
   # 假设已经计算得到时间和特征波数数据
   time_values = np.array([...])  # 时间数据
   char_q_values = np.array([...])  # 特征波数数据
   
   # 绘制双对数图
   plt.figure(figsize=(8, 6))
   plt.loglog(time_values, char_q_values, 'o-', label='Data')
   plt.xlabel('Time t')
   plt.ylabel('Characteristic wavenumber <q>')
   plt.grid(True, which="both", ls="-", alpha=0.2)
   plt.legend()
   plt.show()
   ```

2. **选择合适的拟合区域**：
   - **并非所有数据点都适用于拟合**。相分离过程复杂，通常仅在特定阶段表现清晰幂律
   - **后期粗化阶段**：这是通常关注的阶段。当相畴形成并开始粗化时，体系常进入自相似生长
   - **避免早期和极后期**：早期（成核/旋节线分解初期）或极后期（有限尺寸效应/平衡）可能偏离幂律
   - **目视检查**：找出数据点在双对数图上近似排列成直线的时间区间

3. **拟合方法**：
   ```python
   # 选择拟合区间（例如，帧100到帧400）
   fit_start_frame = 100
   fit_end_frame = 400
   
   # 提取拟合区间的数据
   fit_mask = (frame_indices >= fit_start_frame) & (frame_indices <= fit_end_frame)
   t_fit = time_values[fit_mask]
   q_fit = char_q_values[fit_mask]
   
   # 取对数
   log_t = np.log10(t_fit)
   log_q = np.log10(q_fit)
   
   # 线性拟合
   coeffs = np.polyfit(log_t, log_q, 1)
   slope = coeffs[0]  # 这就是幂律指数m
   intercept = coeffs[1]
   
   # 计算拟合线
   fit_line = 10**(slope * log_t + intercept)
   
   print(f"幂律指数 m = {slope:.3f}")
   print(f"特征长度生长指数 ν = {-slope:.3f}")
   ```

4. **结果解释与物理意义**：
   - 拟合得到的**斜率 $m$ 就是幂律关系 $\langle q \rangle \propto t^m$ 中的指数**
   - 特征长度的生长指数 $\nu = -m$，因为 $l \propto 1/\langle q \rangle$
   - 根据Yuan & Tanaka (2025)：
     * 若 $m \approx -0.5$ (即 $\nu=0.5$)：指示由流体动力学和孔隙弹性主导的粘弹性相分离
     * 若 $m \approx -0.33$ (即 $\nu=1/3$)：对应经典扩散控制的粗化或无HI的情况
   - 通过比较拟合斜率与理论/文献值，可推断体系主导的粗化机制

5. **注意事项**：
   - 确保拟合区间有足够的数据点（通常至少跨越一个数量级的时间）
   - 检查拟合的R²值，确保线性关系良好
   - 考虑多次运行的统计误差
   - 对于有噪声的数据，可以先进行适当的平滑处理

------

## 参考文献

[1] Yuan J.; Tanaka H. *Network-forming phase separation of oppositely charged polyelectrolytes forming coacervates in a solvent*. **Nat. Commun.** **2025**, *16*, 1517. (DOI: https://doi.org/10.1038/s41467-025-56583-6)

[2] Hansen, J.-P.; McDonald, I. R. *Theory of Simple Liquids*, 4th ed.; Academic Press, 2013.

[3] https://zh.wikipedia.org/wiki/%E7%BB%93%E6%9E%84%E5%9B%A0%E5%AD%90  or https://en.wikipedia.org/wiki/Structure_factor



> 本文编辑：摸鱼的帆仔
>
> 校对：AIB001