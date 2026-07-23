---
title: "聚合物自组装体系的自由能面构建：从指标到景观的综合综述"
date: "2026-06-11"
last_modified_at: "2026-06-11"
tags: [self-assembly, free-energy-surface, molecular-dynamics, metadynamics, polymer-nanoparticle, umbrella-sampling, coarse-grained, review]
description: "系统综述构建聚合物自组装体系自由能面所需的指标选取、坐标组合、采样方法与可视化策略，涵盖二维FES、PMF计算与多维景观的最新文献进展"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-1kxerg.jpg"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-1kxerg.jpg"
author: Xufan Gao
lang: zh-CN
---

# 聚合物自组装体系的自由能面构建——指标选取、采样方法与景观解读的综合综述

> 由AI调研和总结，请自行甄别信息正确性

## 摘要

> 自组装体系的自由能面（Free Energy Surface，FES）是理解纳米粒子稳定性、动力学路径和可控制备的核心工具。本文系统综述了构建自组装体系自由能面所需的指标体系、坐标组合、采样方法和可视化策略，重点关注**类蛋白折叠漏斗状势能面**在聚合物纳米粒子体系中的应用。**候选指标**涵盖最大团簇占比$f_{LCC}$、异质/同质接触数$C_{AB}/C_{AA}$、混合度指数$\chi_{mix}$、回转半径$R_g$、溶剂可及表面积SASA、结构因子$S(q)$、径向分布函数$g(r)$、配位数分布、网络指标、相互作用能分解、构型熵估计等。对于每项指标，给出数学定义、物理意义、与聚集稳定性的相关性、对噪声和采样的敏感性、计算复杂度，以及是否可由常见软件直接输出。**推荐坐标组合**包括$(f_{LCC}, R_g)$、$(C_{AB}, \chi_{mix})$、$(S(q), f_{LCC})$等。本文还提出漏斗质量评分函数，综合考虑自由能差、陷阱数目、粗糙度、产物生成概率等因素。对于软件可能缺失的功能，给出利用PLUMED、WHAM、MDAnalysis脚本等补充实现要点，并提供可视化推荐与Python/MDAnalysis代码示例。

### 核心结论

- 自组装FES构建的核心挑战是**CV选择**：单一距离型CV会混淆不同机制，需要**”聚集程度坐标+构象紧凑度坐标”的二维设计**
- **CV选择的核心逻辑**：一个CV反映”聚集到什么程度”（如$f_{LCC}$、接触数、cluster size），另一个CV反映”构象是否紧密”（如$R_g$、链端距离、coordination number）
- 推荐主图坐标组合为$(f_{LCC}, R_g)$，备选包括$(C_{AB}, \chi_{mix})$、$(S(q), f_{LCC})$、$(\langle z \rangle, E_{int})$等
- 增强采样（Metadynamics、伞形采样、REST2等）对自组装体系至关重要，特别是路径复杂、能垒高的多步组装过程
- 评估FES质量应使用**漏斗评分函数**，综合考虑全局稳定性、陷阱深度与数量、产物生成概率
- 文献中明确构建自组装二维FES的工作仍不多，Varner等2025年的”距离+链构象”二维FES是最直接的方法学参考，体现了”过程+构象”的CV设计思想

## 背景

自组装作为软物质、纳米材料和生物大分子领域的核心现象，其热力学驱动力和动力学路径的理解是建立**可推广结构-性质关系**的关键。蛋白折叠领域的“漏斗势能面”概念为理解自组装提供了理论框架：**折叠态对应能量漏斗底部，分散态对应漏斗顶部，中间过程需穿越不同深度的能垒**。但与蛋白折叠不同，自组装体系（如聚合物胶束、纳米颗粒聚集体）涉及**多体相互作用、组分多样性和可调参数空间**，其自由能面构建面临独特挑战。

近10-15年，分子模拟在自组装研究中扮演越来越重要的角色。从全原子MD到DPD粗粒化、从直接Boltzmann反演到Metadynamics增强采样，研究者发展了多种构建FES的方法。然而，**自组装FES的系统综述仍相对缺乏**，特别是在指标选择、坐标组合、采样策略和结果解读方面缺乏统一指导。

本综述基于近10-15年文献，系统梳理构建自组装体系自由能面的指标与方法，重点关注**聚合物-聚合物共组装形成的纳米粒子**这一类重要体系。

## 研究内容

### 一、指标清单与评估

构建自组装体系自由能面需要选择合适的**集体变量（Collective Variable, CV）**。以下逐项列出候选指标的定义、物理意义、计算表达式和评价。

#### 最大团簇占比 $f_{LCC}$

定义为最大簇中粒子数$S_{\max}$与系统总粒子数$N$之比：

$$
f_{LCC} = S_{\max} / N
$$

反映体系聚集程度，值越接近1说明大部分粒子聚成一团；在多体聚集时常用于区分集聚与分散态。物理上与相互作用强弱、温度、浓度相关。计算简单，可通过并查集算法或GROMACS中`gmx cluster`模块得到每帧$S_{\max}$。对含噪轨迹较稳健，但需足够采样显著信号。此指标可直接用于二维FES绘制，例如与能量或$R_g$联立分析。

#### 异质/同质接触数 $C_{AB}$, $C_{AA}$

定义为不同种类（A-B）或同种类（A-A）粒子对在给定截断距离$r_c$内的数目：

$$
C_{AB} = \sum_{i \in A, j \in B} \Theta(r_c - r_{ij})
$$

其中$\Theta$为阶跃函数。物理意义为互补组分间或自身间的结合程度。高$C_{AB}$意味着A/B混合良好，$C_{AA}$大则表明A粒子自团聚明显。该指标线性依赖截断参数，需经验选取，一般取第一近邻距离。计算可用MDAnalysis、MDTraj直接累加距离判据。对噪声较敏感，但能揭示组分间亲和性。适合与$f_{LCC}$组合绘制二维FES。

#### 混合度指数 $\chi_{mix}$

定义为混合接触占比：

$$
\chi_{mix} = C_{AB} / (C_{AB} + C_{AA} + C_{BB})
$$

反映A/B异质混合程度。若系统完全混合，$\chi_{mix} \to 1$，若自分相，$\chi_{mix} \to 0$。物理上表达共组装质量；易于从轨迹计算，只需接触数统计。对系统大小和配比敏感，需注意正则化（当某类接触极少时易发散）。通常与全局混合能判断结果一致。

#### 质心距/簇半径 $R_g$

回转半径定义为所有原子到质心距离的均方平均：

$$
R_g = \sqrt{\frac{1}{M} \sum_i m_i |\mathbf{r}_i - \mathbf{r}_{CM}|^2}
$$

其中$M = \sum m_i$为总质量。$R_g$衡量结构整体尺寸和紧凑度，值小表示高度聚集。可直接用GROMACS（`gmx gyrate`）、MDTraj或MDAnalysis计算。受温度和形状变化影响，对于非团簇或高噪声轨迹可能误差较大。通常与接触数或簇大小联合使用。

#### 溶剂可及表面积（SASA）与结合面面积（AP）

计算分子/纳米粒子聚集物的溶剂可及表面积，常用GROMACS（`gmx sasa`）或MDTraj方法。SASA减少通常意味着疏水驱动的聚集增强。AP一般指两组分接触界面的面积，也可用切换函数计算。公式较复杂（求球面网格交点或解析式），一般通过程序得出。受粒子形状、定义截断面影响；对反映暴露/接触面有用，但计算量较大。

#### 结构因子 $S(q$) 与峰位 $q^*$

结构因子度量体系在动量空间的有序性：

$$
S(\mathbf{q}) = \frac{1}{N} \left\langle \left| \sum_j e^{-i \mathbf{q} \cdot \mathbf{r}_j} \right|^2 \right\rangle
$$

峰位$q^*$对应主要结构周期（如晶格常数约$2\pi/q^*$）。在无序系统$S(q) \approx 1$平坦，有序聚集时出现峰值。计算可用FFT方法或PLUMED的`STRUCTURE_FACTOR`功能得到角平均$S(q$)谱。$S(q$)对体系有序性敏感，对噪声和有限尺寸影响较大。适合揭示长程有序结构，但不适合小团簇局部自由能面。峰位$q^*$可作为聚集间距指标。

#### 径向分布函数 $g(r$)

定义为单位密度下在距离$r$处的粒子对分布概率：

$$
g_{AB}(r) = \frac{1}{4\pi r^2 \langle \rho_B \rangle} \frac{1}{N_A} \sum_{i \in A, j \in B} \delta(r_{ij} - r)
$$

$g(r$)刻画短程结构，如第一个峰对应近邻距。可用GROMACS（`gmx rdf`）或MDAnalysis计算。$g(r$)有助于确定配位数（通过积分至第一个谷）并作为坐标之一用于FES投影（例如$(R_g, q$)或$(R_g, g(r^*)$)）。对统计采样要求较高，曲线平滑性决定计算精度。

#### 配位数分布

统计每个粒子在某截断半径$r_c$内邻居数的分布（histogram of coordination number）。可定义粒子$i$的配位数$z_i = \sum_{j \neq i} H(r_c - r_{ij}$)，然后统计$P(z$)。反映局部结构多样性，分布宽度增大意味着结构无序度高。通常采用第一近邻截断。易用MDAnalysis或numpy histogram快速得到。对噪声敏感，但适用于表征局部稳定性和链状聚集。

#### 邻接矩阵/网络指标

将体系视为图，节点为粒子，边存在于两粒子$< r_c$。常用指标：节点度（平均度$\langle k \rangle$）、连通性（是否连通）、社区结构（模块度modularity）等。模块度定义为社区内连边密度与随机模型差异。高模块度表明粒子可分为几个紧密子簇。度分布、聚类系数等也可反映聚集组织。计算可借助NetworkX处理邻接矩阵。此类指标能够捕捉复杂结构性特征，但直观关联自由能不易量化，多作为定性辅助。

#### 相互作用能分解

将系统总能分解为部件间相互作用，如A–B和A–A、B–B能；或范德华/静电分量。常定义A–B相互作用能为：

$$
E_{int} = \sum_{i \in A, j \in B} [V_{LJ}(r_{ij}) + V_{elec}(r_{ij})]
$$

GROMACS可用`gmx energy`或自写脚本分析能量输出。物理上直接反映组分间吸引力或排斥，适合评估稳定性。需要平衡剪切截断误差和静电处理方式。对小改动灵敏；作为坐标使用时一般与结构指标联合。

#### 熵估计方法

构型熵可从采样的状态分布估算：如以简化状态簇概率$p_i$计算香农熵$S = -k_B \sum_i p_i \ln p_i$。或利用协方差矩阵计算Schlitter熵近似。方法依赖采样质量且对统计不足敏感，适合定性比较。对粗粒化自组装尤难准确定义，通常视为热力学态稳定性的补充说明。

#### 自由能估计方法

常见方法包括直接Boltzmann公式$F(\mathbf{x}) = -k_B T \ln P(\mathbf{x}$)，最大似然WHAM、umbrella采样、metadynamics（偏置势逼近$F$）及Replica Exchange（REST2）等。其中Metadynamics可生成多维CV下的FES（偏置收敛为$-F$），WHAM用于合并多窗口样本。各方法要求设计适当CV或窗口，计算量随维度增大急剧增长。对动态特征复杂、自组装路径冗长体系，增强采样尤为重要。

#### 指标对比表

| 指标 | 物理意义 | 敏感性/复杂度 | 软件工具支持 | 适用FES |
| --- | --- | --- | --- | --- |
| 最大团簇占比$f_{LCC}$ | 系统聚合程度（1=全聚集） | 低（简洁统计） | IMPULSE/GROMACS | 常用，易与其他CV联合，如$(f_{LCC}, R_g$) |
| 接触数$C_{AB}/C_{AA}$ | 组分间/同分子内相互作用强度 | 依距选敏感 | IMPULSE/MDAnalysis | 适合衡量混合态，可与$f_{LCC}$等联合绘制 |
| 混合度$\chi_{mix}$ | 异质混合程度（0分离，1混合） | 中等 | 需自定义 | 二维FES构建，如$(\chi_{mix}, f_{LCC}$) |
| 回转半径$R_g$ | 聚集物尺寸、紧凑性 | 中（形状敏感） | GROMACS/MDTraj | 常用CV，可与团簇大小/接触数联合绘制 |
| SASA/AP | 暴露表面积/结合面面积（聚集稳定性指标） | 高（计算量大） | GROMACS/MDTraj | 通常与$R_g$等结合，表征疏水/亲水效应 |
| 结构因子$S(q$) | 长程有序性（峰值反映周期结构，$q^*$首峰位） | 中等 | MDAnalysis/dynasor | 对局域团聚不敏感，常用于长程有序结构分析 |
| 径向分布$g(r$) | 粒子近邻分布（峰值与配位相关） | 高（需大量采样） | GROMACS | 常用结构表征，可基于第一个峰值定义配位数坐标 |
| 配位数分布 | 局部结构差异（粒子邻居数分布） | 中 | MDAnalysis | 辅助指示结构均匀度，可做直方图分析 |
| 网络指标（度、模块度等） | 体系连通性与社区结构（模块度高→分相） | 高（需计算图分） | NetworkX/IMPULSE | 可反映分相与混合，需配合其他指标综合评估 |
| 相互作用能分解 | 聚集驱动力（VDW/静电贡献） | 中（需二次计算） | GROMACS/PLUMED | 用于动力学分析，多与结构指标结合 |
| 构型熵估计 | 聚集态自由度大小 | 高（统计需求大） | MDAnalysis/PlaMO | 常作为自由能面的补充说明，不直接作CV |
| 自由能估计方法 | FES计算技术（-kTlnP、WHAM、MetaD、REST2等） | 高（计算密集） | PLUMED/GROMACS | 强调方法而非坐标，用于构建FES本身 |

每项指标需结合具体系统和需求评估：物理上是否能反映粒子稳定性或陷阱深度（例如$f_{LCC}$反映聚合程度，$R_g$/SASA反映紧凑度和暴露度）；对噪声/采样的敏感性（如$g(r$)和熵估计需大量采样，网络指标对小团簇波动敏感）；计算复杂度和可行性（简单几何量如$f_{LCC}$、$R_g$计算成本低，相互作用能需逐对累加）；是否可由现有软件直接输出（GROMACS/PLUMED自带RDF、$R_g$、能量分解；MDAnalysis可快速自定义计算）；以及能否用作联合CV绘制二维/三维自由能面（一般推荐2维组合，确保信号区分度较高且易于统计）。


---


### 二、近期文献进展：二维/多维FES的CV类型详解

本节聚焦2021-2025年明确构建聚合物（及聚合物-药物）纳米粒子自组装自由能景观的分子模拟研究，共分析10篇文献的CV设计、采样策略和FES构建方法。**这些文献体现了CV设计的核心思想：一个坐标反映“聚集到什么程度”，另一个坐标反映“构象是否紧密”**。

#### 1. Varner et al., 2025 – 二嵌段共聚物胶束链交换机制详解

**体系**：二嵌段共聚物胶束的链交换/链逃逸问题，强分凝条件下的单链逃逸过程。

**模拟方法**：结合粗粒化MD和增强采样，计算链逃逸过程的**二维自由能面**，并用forward-flux sampling研究稀有事件动力学。

**CV设计**：
- **distance-based CV**：推动链从胶束中逃逸
- **core block end-to-end distance**：确保链构象充分采样

**FES构建**：二维FES形式为$F(R, r) = -k_B T \ln P(R,r)$，其中$R$是链逃逸距离，$r$是core block end-to-end distance。作者从2D FES投影到1D自由能曲线：先由$F(R,r)$得到$P(R,r)=\exp[-\beta F(R,r)]$，再对构象变量积分得到$P(R)$，最后$\beta F(R)=-\ln P(R)$。

**关键发现**：二维FES揭示两条几乎简并的逃逸路径：
- 一条接近Halperin–Alexander的budding-like机制
- 另一条是链逐珠（bead-by-bead）伸展逃逸

计算不同core block长度下的自由能垒，发现其中一条路径的能垒满足$\beta\Delta F_{\rm barr}\sim N_{\rm core}^{2/3}$。

**方法学价值**：这篇最适合借鉴“一个进程坐标 + 一个构象坐标”的二维FES设计。如果只用一个聚集距离/团簇大小坐标，可能会把不同机制混在一起；最好再加一个能区分“紧密团聚、拉伸桥连、松散网络”的构象坐标，比如$R_g$、core compactness、异质接触数、端到端距离、链拉伸度或局部密度。

#### 2. Zhang & Meng, 2025 – 超分子二嵌段共聚物无序-有序转变详解

- **体系**：超分子二嵌段共聚物的disorder–order transition
- **模拟方法**：比较共价二嵌段共聚物和超分子二嵌段共聚物，使用smart Monte Carlo模拟动力学路径，再用**string method**构建**minimum free energy path**
- **CV设计**：FES可写成$F(S_{\rm order}, f_{\rm bond})$形式，其中：
  - $S_{\rm order}$：结构序参量
  - $f_{\rm bond}$：动态键比例
- **关键发现**：沿minimum free energy path讨论，得到transition state和free energy barrier，并将自由能分解为A–B interaction energy和association energy
- **方法学价值**：这篇说明自组装FES不一定要写成$F(r)$，也可以写成$F(S_{\rm order}, f_{\rm bond})$或沿minimum free energy path讨论。对应到二元纳米药物体系，可以借鉴成$F(S_{\rm assembly}, f_{\rm hetero\ contact})$，其中$S_{\rm assembly}$是整体有序/组装程度，$f_{\rm hetero\ contact}$是HA-OP、载体-药物或A-B异质接触比例。这比单纯距离更贴近“共组装是否可控”

#### 3. Gautham & Patra, 2022 – 聚合物接枝纳米粒子深度学习PMF详解

- **体系**：polymer-grafted nanoparticles
- **模拟方法**：从小规模polymer-grafted nanoparticle cluster的MD轨迹中学习pair interaction，然后用deep-learning PMF-based simulation预测大量接枝纳米粒子的3D自组装结构，包括percolating networks和bilayers
- **CV设计**：核心CV是颗粒间相对位置/距离，PMF形式为$W_{\rm eff}(\mathbf{R}_{ij}, \Omega_{ij}, \ldots)$
- **FES构建**：构建的是effective potential of mean force，而不是传统umbrella sampling得到的简单$F(q)$。用深度学习从小体系MD cluster轨迹中学习两颗polymer-grafted nanoparticles的有效相互作用，再把这个PMF放入更大规模的粒子模拟中
- **方法学价值**：不一定直接在全体系上构建高维FES，也可以先计算/学习“组装基元之间的PMF”，再用PMF预测大体系组装。这对大规模纳米药物自组装很现实，因为全体系$F(q_1,q_2,q_3)$采样困难；而基元-基元、载体-药物、HA-OP、OP-OP、HA-HA的pair/many-body PMF可以作为降维的热力学输入

#### 4. Wu, Pal & Keten, 2023 – 隐式链粒子模型详解

- **体系**：matrix-free polymer grafted nanoparticles，以PMMA的chemistry-specific coarse-grained MD为测试体系
- **模拟方法**：提出implicit chain particle model，核心是用strain-energy mapping framework和PMF计算建立粒子间有效相互作用
- **CV设计**：不是传统的单一两颗粒拉开距离PMF，而是把颗粒排列在close-packed lattice configuration中，通过bulk dilation/compression的strain-energy density匹配来推导有效相互作用。CV更接近于颗粒间距/晶格膨胀压缩程度
- **FES构建**：构建的是coarse-grained effective interaction/PMF，形式上类似$W_{\rm eff}(a) \leftrightarrow U_{\rm strain}^{\rm CG-MD}(a)$，其中$a$是晶格尺度或颗粒间距相关坐标
- **关键发现**：ICPM可将计算速度相对CG-MD提升约$10^5$–$10^6$倍
- **方法学价值**：适合借鉴“从显式链模型中抽取有效自由能相互作用”的思想。对二元聚合物体系，如果全体系太大，可以先做若干代表性小体系PMF：例如HA-OP、OP-OP、HA-HA、载体-药物之间的effective PMF，再把这些PMF作为coarse-grained self-assembly landscape的输入

#### 5. Munaò et al., 2018 – 原子级纳米颗粒PMF详解

- **体系**：atomistic silica/gold nanoparticles，包括bare gold nanoparticles和polyethylene-coated gold nanoparticles
- **模拟方法**：用atomistic MD计算纳米颗粒之间的PMF。先用silica nanoparticles对比Hamaker理论来验证过程，再计算bare与polyethylene-coated gold nanoparticles的有效相互作用
- **CV设计**：主要CV是两颗纳米颗粒之间的interparticle separation，即颗粒中心距离。对coated gold nanoparticles，还考察grafting density $\rho_g$对PMF的影响
- **FES构建**：构建一维PMF：$W(r) = -k_B T \ln P(r)+C$，或等价地由约束/平均力积分得到$W(r)$
- **关键发现**：
  - silica nanoparticles的PMF与粒径相关性不强，但较大颗粒出现明显surface interaction peak
  - bare gold nanoparticles作用较弱
  - polyethylene-coated情况下，有效相互作用随接枝密度增强。中等$\rho_g$下PMF类似Lennard-Jones型，而高$\rho_g$、小间距下逐渐变为更强排斥
- **方法学价值**：适合借鉴“表面聚合物层如何改变纳米颗粒PMF”的分析逻辑。对纳米药物载体来说，表面修饰密度、链长、亲疏水性、电荷状态都可以通过pair PMF表征其聚集倾向或抗聚集稳定性

#### 6. Egorov, 2011 – 立体稳定lock-and-key胶体详解

- **体系**：sterically stabilized lock-and-key colloids in polymer solution。key particle和lock cavity都设定为cylindrical shape，表面均匀接枝polymer chains，同时溶液中有free polymer chains
- **模拟方法**：使用self-consistent field theory，计算sterically stabilized lock-key particles在polymer solution中的PMF
- **CV设计**：由于假设key和lock都沿$z$轴同轴排列，PMF是单坐标函数$W(z)$，其中$z$是lock-key separation
- **FES构建**：先通过SCF理论得到不同$z$下的Helmholtz free energy $A(z)$，再定义$\beta W(z) = \beta A(z) - \beta A(\infty)$
- **关键发现**：lock-key interaction可通过几何匹配、接枝密度、自由链体积分数和焓相互作用调控。尺寸匹配时depletion attraction最强，聚合物steric stabilization可使binding-unbinding transition更尖锐
- **方法学价值**：这篇的重点不是动态轨迹采样，而是用SCF直接计算自由能面。它适合借鉴到“载体表面接枝层/聚合物刷/溶剂化层调控粒子间可逆结合”的场景。如果二元结合几何明确，比如HA与OP局部复合、载体表面基元与另一个颗粒/膜片段结合，可以把$z$或$r$作为PMF坐标，并扫描链长、接枝密度、溶剂质量、电荷状态

#### 7. Wang & Ferguson, 2017 – 环状聚合物拓扑约束

- **体系**：polyethylene ring polymers，包括trefoil knot、catenane、Borromean等拓扑状态
- **模拟方法**：用MD加nonlinear manifold learning，抽取低维自由能面
- **CV设计**：不是纳米颗粒自组装，但它是很好的“非预设CV的聚合物自由能面”参考。从多指标中学习低维坐标，再构建$F(\xi_1,\xi_2)=-k_B T \ln P(\xi_1,\xi_2)$，其中$\xi_1,\xi_2$是数据驱动的慢变量
- **FES构建**：这些FES揭示degree of polymerization和topological constraints如何影响可热访问构象、手性对称破缺、folding/collapse pathways
- **方法学价值**：如果不想手动限定CV为距离/$R_g$/contact number，可以用manifold learning、tICA、diffusion map、PCA/UMAP之类从多指标中学习低维坐标，再构建FES。这会更像“真实景观”，但解释性要靠事后把$\xi$与$R_g$、接触数、团簇大小、混合度相关联

#### 8. Sucerquia et al., 2022 – 银团簇ab initio metadynamics详解

- **体系**：$\ce{Ag5}$/$\ce{Ag6}$ clusters
- **模拟方法**：用ab initio metadynamics，通过PLUMED和ASE接口计算free-energy landscape
- **CV设计**：选用的CV是**radius of gyration**和**coordination number**，用它们比较planar/non-planar isomers的相对自由能
- **FES构建**：这对聚合物纳米粒子非常自然：$F(R_g, C_{\rm contact})$，低$R_g$、高contact number是紧密稳定颗粒；高$R_g$、低contact number是分散或松散网络；低$R_g$、低异质接触可能是单组分塌缩陷阱
- **方法学价值**：很直接展示了“非距离型二维CV”如何做纳米团簇FES。$(R_g)$和coordination/contact number对聚合物纳米粒子非常自然，物理含义很清楚

#### 9. Balestra & Semino, 2022 – ZIF-8自组装早期阶段详解

- **体系**：ZIF-8早期成核与热分解
- **模拟方法**：用all-atom well-tempered metadynamics，明确探索了一组physically relevant collective variables，选择合适子集
- **CV设计**：说明自组装FES的CV可以是**coordination/connectivity**、**cluster size**、**ring count**等，而不必是距离
- **关键发现**：结果包括Zn–N connectivity快速增加、小团簇蒸发并形成少数大团簇、$\ce{Zn(MIm)4^{2-}}$/$\ce{Zn(MIm)3^-}$复合物、4/5/6-membered rings等寿命差异
- **方法学价值**：虽然这是MOF，不是聚合物，但它说明自组装FES的CV可以是connectivity、cluster size、ring count。对二元聚合物体系，ring count不一定适用，但connectivity、largest cluster size、heterogeneous contact network是非常适用的

#### 10. Méndez & Semino, 2024 – ZIFs自组装热力学

- **体系**：ZIF-4自组装的early nucleation和late growth
- **模拟方法**：用reactive force field + well-tempered metadynamics
- **CV设计**：自由能分析聚焦金属离子配位变化、building block形成、ligand coordination saturation，以及不同晶面/多晶型增长的热力学差异
- **方法学价值**：对于多步自组装，CV可以按“化学连接/局部配位饱和度/生长单元加入程度”定义。对应到二元聚合物体系，可以类比为$F(n_{\rm AB\ contact}, n_{\rm core})$或$F(\text{hetero-coordination}, \text{cluster growth})$。这比单个距离更容易表达“成核—生长—稳定化”的过程

---

#### 文献CV类型总结

| CV类型 | 代表文献 | 典型形式 | 方法学价值 | CV设计思想 |
| --- | --- | --- | --- | --- |
| 距离 + 链构象 | Varner 2025 | $F(r, R_{ee})$ | 组装进程 + 链伸展/紧密度 | 过程坐标 + 构象坐标 |
| 结构序参量 + 动态键比例 | Zhang & Meng 2025 | $F(S_{\rm order}, f_{\rm supra})$或MFEP | 有序组装程度 + 异质复合比例 | 聚集程度 + 协同效应 |
| 颗粒间PMF | Gautham & Patra 2022; Munaò 2018 | $W(r)$或ML-learned PMF | A-A、B-B、A-B基元相互作用 | 简化为有效相互作用 |
| 压缩/膨胀自由能 | Wu 2023 | strain-energy mapped PMF | 纳米颗粒紧密堆积稳定性 | 体积变化 + 自由能响应 |
| 拓扑/数据驱动低维坐标 | Wang & Ferguson 2017 | $F(\xi_1, \xi_2)$ | 从多指标自动学习慢变量 | 无预设CV，数据驱动 |
| $R_g$ + coordination/contact number | Sucerquia 2022; ZIF metadynamics文献 | $F(R_g, CN)$, $F({\rm connectivity}, {\rm cluster\ size})$ | 最适合转译成聚合物纳米粒子稳定性景观 | 整体紧密度 + 局部连接度 |

**CV设计的核心原则**：成功的二维FES设计通常遵循“聚集程度+构象紧凑度”的逻辑——一个坐标描述“组装到什么程度”，另一个坐标描述“结构是否紧密”。这种设计能区分不同的组装机制（如紧密团聚vs松散网络，拉伸桥连vs塌缩成团）。

---

### 三、关键符号与公式的物理意义

为便于读者理解文献中的CV设计和自由能表达式，本节对常用符号和公式的物理意义进行解释。

#### 自由能面基本概念

文献中的$F$、PMF、free-energy surface/landscape大多不是“势能面”($U$)，而是沿某些集体变量($q$)投影后的有效自由能：

$$
F(q)=-k_B T \ln P(q)+C
$$

二维情况下：

$$
F(q_1,q_2)=-k_B T \ln P(q_1,q_2)+C
$$

其中$P$是体系在某个坐标区域出现的概率，$k_B T$是热能尺度，$C$是任意零点。更低的$F$代表该状态更常出现、更热力学稳定；能垒$\Delta F_{\rm barr}$代表从一个稳定态到另一个状态需要跨越的自由能代价。

#### 核心参数符号表

| 符号/表达式 | 出现场景 | 物理含义 | 稳定性解释 |
| --- | --- | --- | --- |
| $\beta=1/k_B T$ | 多数自由能文章 | 把自由能换算成热能单位 | $\beta\Delta F$越大，越难跨越 |
| $\Delta F_{\rm barr}$ | Varner、Zhang & Meng、Seeger等 | 初态到过渡态的自由能垒 | 能垒越高，动力学越慢、结构越kinetically stable |
| $N_{\rm core}$ | diblock micelle | 疏水核心block的聚合度 | core越长，链逃逸越难 |
| $\beta\Delta F_{\rm barr}\sim N_{\rm core}^{2/3}$ | Varner | budding-like过渡态的表面自由能尺度 | 胶束链交换能垒随core block长度亚线性增长 |
| $r$ | Seeger、Munaò、Egorov等 | 距离型反应坐标 | 描述结合/解离或链逃逸 |
| $R_{\rm ee}$ | Varner | core block端到端距离 | 区分塌缩逃逸和拉伸逃逸 |
| $S_{\rm order}$ | Zhang & Meng | 结构有序参数 | 越高越接近有序组装相 |
| $f_{\rm bond}$, $f_{\rm supra}$ | supramolecular copolymer | 动态键/超分子连接比例 | 表示可逆连接网络成熟程度 |
| MFEP | Zhang & Meng | 自由能面上的最低自由能路径 | 可识别transition state和pathway |
| $W(r)$ | PMF文献 | 距离$r$下的平均力势 | 低谷代表稳定结合，高峰代表排斥/能垒 |
| $W_{AA}$, $W_{BB}$, $W_{AB}$ | 对二元体系的类比 | 同质/异质组分PMF | 判断共组装还是自聚集 |
| $R_g$ | cluster/metadynamics文献 | 回转半径，整体紧密度 | 低$R_g$通常更紧密 |
| $CN$ | cluster/metadynamics文献 | coordination/contact number | 高CN表示局部连接更多 |
| $F(R_g,CN)$ | 纳米团簇FES | 紧密度+局部连接二维景观 | 可区分松散态、紧密态、亚稳态 |
| connectivity | MOF自组装文献 | 关键连接/配位数量 | 表示成核和网络形成程度 |
| cluster size | 自组装文献 | 团簇大小 | 表示成核、生长、并合 |
| $\xi_1$,$\xi_2$ | manifold learning FES | 数据驱动慢变量 | 可发现非人工预设的构象盆地 |

#### 典型公式的物理意义

**1. Varner的能垒标度律**

$$
\beta\Delta F_{\rm barr}\sim N_{\rm core}^{2/3}
$$

其中$\beta=1/k_B T$，$\Delta F_{\rm barr}$是链逃逸能垒，$N_{\rm core}$是疏水核心block的聚合度。$2/3$次方来自Halperin–Alexander budding-like机制的物理图像：逃逸链在过渡态中形成一个类似“球状芽”的globular transition state，其表面自由能随体积/链长的$2/3$次方增长。

**核心意义**：胶束稳定性不是简单随链长线性增加；在budding-like逃逸路径中，能垒近似随核心链段长度的表面积尺度增长。

**2. Zhang & Meng的二维FES**

$$
F(S_{\rm order}, f_{\rm bond})
$$

其中$S_{\rm order}$是结构有序参数，$f_{\rm bond}$或$f_{\rm supra}$表示形成supramolecular bonds/supramolecularly connected chains的比例。

**核心意义**：自组装路径不仅取决于结构是否有序，也取决于可逆连接/结合网络是否形成。动态键可能降低或改变能垒，但也可能引入中间态和路径复杂性。

**3. Egorov的SCF自由能定义**

$$
\beta W(z) = \beta A(z) - \beta A(\infty)
$$

其中$A(z)$是lock-key距离为$z$时的Helmholtz free energy，$A(\infty)$是两者相隔无限远时的自由能，$W(z)$是相对于无限远分离状态的PMF。

**核心意义**：聚合物刷、自由链耗竭作用和几何匹配共同决定颗粒识别/结合的自由能。

**4. 纳米团簇的二维FES**

$$
F(R_g, CN)
$$

其中$R_g$是回转半径，描述团簇整体尺寸/紧密度；$CN$是coordination number，描述原子之间的局部配位/接触程度。

**核心意义**：低$R_g$、高CN的盆地通常对应紧密稳定构型；高$R_g$、低CN对应松散构型；中间盆地可能对应亚稳态异构体。对聚合物纳米颗粒，$CN$可以替换成contact number：$F(R_g, C_{\rm contact})$。

**5. Minimum Free Energy Path (MFEP)**

MFEP是自由能面上从初态到终态最可能经过的低自由能路径。对于disorder–order transition，它大致表示：

$$
\text{disordered state} \rightarrow \text{transition state} \rightarrow \text{ordered state}
$$

沿MFEP可以定义：

$$
\Delta F_{\rm barr}=F_{\rm TS}-F_{\rm initial}
$$

其中$F_{\rm TS}$是过渡态自由能，$F_{\rm initial}$是初态自由能。这个能垒越高，转变越慢；中间有多个局部极小值，就说明路径上有metastable intermediates或kinetic traps。

#### CV选择的总原则

不同类型的CV适用于不同的自组装分析需求：

- **距离$r$**：适合描述“结合/解离”过程，如颗粒靠近/远离、链逃逸等
- **$R_g$、contact number、coordination、cluster size**：适合描述“颗粒是否紧密稳定”，反映聚集体整体紧密度和局部连接程度
- **$S_{\rm order}$、$f_{\rm bond}$、network/connectivity**：适合描述“是否形成有序、协同、可控的组装结构”，反映组装质量和协同效应
- **$\xi_1$,$\xi_2$**：适合在人工CV不确定时用数据驱动方式寻找自由能景观坐标，通过流形学习发现非预设的慢变量

这些符号和公式的核心是：通过合适的集体变量投影，将复杂的多维自组装过程降维到可理解、可计算的自由能景观，从而定量分析稳定性、动力学路径和可控性。

---

### 四、推荐坐标组合与计算流程

基于上述指标评估和文献分析，以下核心坐标组合可用于主图/备选图的自由能面构建。**CV设计的核心思想是：一个坐标反映“聚集程度”（如$f_{LCC}$、接触数），另一个坐标反映“构象紧凑度”（如$R_g$、链构象），这样可以区分不同的组装机制和路径**。

#### 组合1：$(f_{LCC}, R_g$)

$f_{LCC}$捕捉聚集程度，$R_g$表征整体尺寸，两者可区分紧密团簇与分散态。

**计算流程**：
- **数据提取**：逐帧用MDAnalysis或自定义脚本确定最大簇大小并归一化得$f_{LCC}(t)$；用GROMACS或MDTraj计算$R_g(t)$
- **直方化/KDE**：对结果进行二维直方化或核密度估计，计算每个格点概率$P(f, R_g)$
- **归一化**：归一化自由能$F = -k_B T \ln P$
- **误差估计**：建议使用细致的网格（bin宽视数据散布调整），通过多次Bootstrap估误差
- **增强采样**：缺采样区域可考虑温度扩展或Metadynamics增强（比如在$f_{LCC}$方向施加偏置）以填补低概率区间

#### 组合2：$(C_{AB}, C_{AA})$或$(C_{AB}, \chi_{mix})$

该组合直观衡量组分混合程度。

**计算流程**：
- **数据提取**：根据距离截断计算每帧$C_{AA}(t), C_{AB}(t)$；或计算混合度$\chi_{mix}(t)$
- **直方图构建**：构建二维直方图$P(C_{AA}, C_{AB})$（或$P(\chi_{mix}, C_{AB})$）
- **归一化**：归一化得$F$面
- **平滑处理**：由于接触统计可能波动大，需足够长轨迹并可适用滑动窗口平均平滑
- **增强采样**：增强采样建议对非混合态构造预偏置

#### 组合3：$(S(q^*), f_{LCC})$

适用于有序自组装体系，如纳米晶体。

**计算流程**：
- **结构因子计算**：对每帧计算小范围$q$的$S(q)$（用FFT或Dynasor库）
- **峰位提取**：提取主峰$q^*$或对应峰值$S(q^*)$
- **FES构建**：与$f_{LCC}$一起构建FES：$F(f_{LCC}, S(q^*))$

**方法学价值**：此组合将聚集程度与长程有序性结合，可显著区分混沌聚集与形成晶格结构的情况。

#### 组合4（备选）：$(\langle z \rangle, E_{int}$)

使用每帧平均配位数$\langle z \rangle$（通过积分$g(r$)首谷得到）与体系总相互作用能$E_{int}$，构建$F(\langle z \rangle, E_{int}$)面。此图可揭示形态相变中结构紧凑度与结合能的关系。

#### 组合5（备选）：$(Q, R_g$)

适合分析分相体系。先构造邻接图计算社团模块度$Q(t$)，然后与$R_g(t$)配对。$F(Q, R_g$)可显示不同分相（高$Q$小$R_g$）和混合状态（$Q \approx 0$大$R_g$）区域。

#### 流程示意

{% raw %}
```mermaid
graph TB
    Traj([轨迹文件]) --> Pre{{预处理}}
    Pre --> ComputeMetrics{计算指标}
    ComputeMetrics --> fLCC[f_{LCC}]
    ComputeMetrics --> RG[R_g]
    ComputeMetrics --> Contacts[C_{AB},C_{AA}]
    ComputeMetrics --> Struct[SASA, g(r)]
    ComputeMetrics --> Energy[E_{int}]
    fLCC --> Binbin[二维直方或KDE]
    RG --> Binbin
    Contacts --> Binbin
    Struct --> Binbin
    Energy --> Binbin
    Binbin --> FreeE["自由能计算$F=-k_BT\ln P$"]
    FreeE --> Plot[绘制二维等高或色图]
```
{% endraw %}

### 五、评分函数与实现建议

为量化自由能面“漏斗”特性，可定义评分函数综合考虑全局稳定性和结构多样性。例如建议的漏斗评分：

$$
S_{\rm funnel} = w_1 \Delta F_{\rm global} - w_2 N_{\rm trap} - w_3 {\rm Roughness} - w_4 \sigma_{\rm basin} + w_5 P_{\rm prod}
$$

其中：$\Delta F_{\rm global} = F_{\rm dispersed} - F_{\rm assembled}$为主井深度差（全局稳定性），$N_{\rm trap}$为能量陷阱（局部极小）数量，Roughness为势能面粗糙度度量（如主路径振荡总和），$\sigma_{\rm basin}$为主要基态宽度（椭圆拟合方差），$P_{\rm prod}$为从游离态演化到聚集态的“产物概率”。各权重$w_i$可根据系统需求调节（例如强调稳定性则增大$w_1$）。该函数可结合路径分析工具计算，参考蛋白折叠领域“foldability score”概念。

## 关键结论与批判性总结

- 综述了构建自组装体系自由能面所需的**完整指标体系**，包括物理意义、敏感性、软件支持和FES适用性的全面评估
- 推荐的**核心坐标组合**为$(f_{LCC}, R_g$)、$(C_{AB}, \chi_{mix}$)、$（S(q^*), f_{LCC}）$等二维设计，能够兼顾稳定性与形态区别
- 提出**漏斗评分函数**$S_{\rm funnel}$，综合考虑全局稳定性、陷阱深度与数量、产物生成概率
- 给出软件缺失功能的**补充方案**（PLUMED插件、WHAM工具、MDAnalysis脚本等）
- 梳理了近10-15年自组装FES相关文献的**6类CV设计**：距离+链构象、序参量+动态键比例、颗粒间PMF、压缩/膨胀自由能、拓扑/数据驱动坐标、$R_g$+配位数

### 局限性

- 现有文献中**明确构建自组装二维FES的工作仍不多**，许多自组装模拟只报告micelle、vesicle、lamella、cluster size、$R_g$、$S(q$)、morphology diagram，并未真正构建$F(q) = -k_B T \ln P(q$)或PMF
- 漏斗评分函数的**权重选择**目前缺乏系统指导，需要根据具体体系调节
- **增强采样方法的选择**在自组装体系中尚无统一标准，不同方法各有优劣
- 自组装FES的**实验验证**仍然困难，特别是动力学陷阱和亚稳态结构分布的对应关系
- 粗粒化模型的选择会显著影响FES结果，不同分辨率的模型可能给出不同的漏斗形貌

### 后续工作优先级

1. 对推荐坐标组合做实际模拟验证，检查是否合理区分聚集态
2. 实现并测试上述评分函数，评估能否量化漏斗质量
3. 针对软件功能空缺，开发补充脚本或PLUMED模块
4. 使用增强采样（如并行Metadynamics、REST2等）提高FES可靠度
5. 基于生成的FES提出可实验验证的预测（如体系在不同参数下的相行为）

以上结论均基于现有文献与官方文档所述原理。未来可持续关注相关软件更新和新的案例研究。对于自组装体系自由能面构建的进一步工作，建议优先关注**二维/多维FES设计**和**数据驱动CV抽取**两个方向，它们将是未来改善现有方法学局限的关键。
