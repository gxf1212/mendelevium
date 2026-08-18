---
title: "MD怎么说明激酶活性变化？多副本轨迹看活性标记的概率分布偏移"
date: "2026-08-17"
last_modified_at: 2026-08-17
tags: [kinase, molecular-dynamics, conformational-ensemble, allostery, apo-holo, population-shift, enhanced-sampling, structure-dynamics-function]
description: "用MD模拟论证激酶活性调控的方法论总结：多副本MD比较apo和holo体系的构象集合偏移，从αC螺旋、DFG基序、R-spine、ATP口袋等结构标记的概率分布变化看活性状态转换"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-e7z65r.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-e7z65r.jpg"
author: Xufan Gao
lang: zh-CN
---

# MD怎么说明激酶活性变化？多副本轨迹看活性标记的概率分布偏移

> AI总结，请自行甄别，仅供参考。

## 背景

激酶活性的结构基础已经研究了几十年：αC螺旋内旋（αC-in）、DFG-in、A-loop有序化、R-spine堆叠连续，这些结构特征对应活性态，它们的反向对应非活性态。但晶体结构有个根本问题：低温和单一构象，往往捕捉的是最稳定的一个状态，而激酶在溶液中实际上是在多个亚态之间采样的。

所以，当你想说明“某个配体或条件让激酶变活跃了（或失活了）”时，单靠晶体结构是说不清楚的。你需要回答的是：这个配体是否改变了构象集合的概率分布？

MD模拟可以回答这个问题，但前提是使用合适的方法。很多研究用MD说激酶活性，最后只给了RMSD、RMSF和“配体稳稳地结合在口袋里”——这些只能证明模拟是合理的，不能证明机制。

> 这里可能需要叠一下甲：本文里面的MD模拟方法都是普通MD适用于没有资源和时间去做复杂的增强采样、自由能计算等更严格的验证方法

## 核心方法：多副本apo vs holo的构象集合比较

### 模拟设置

最基本的设置是对比两组模拟：

- **apo组**：没有配体的激酶，至少3–5条独立轨迹，每条100–500 ns
- **holo组**：结合配体后的激酶，同样的轨迹数和长度

关键点是**轨迹必须相互独立**：不同的初始速度、从不同构象出发，这样才能估计概率分布。如果只跑一条轨迹然后说“αC从out变成了in”，无法判断这是系统性的偏移还是随机涨落。更稳妥的做法是**从多个可信的实验构象作为补充起点**，检查结论是否依赖初始结构。PAK1那篇用了3个独立400 ns轨迹，是很好的reviewer-friendly先例。

多replica的价值在于**评估结果的重复性和不确定性**，但不能替代足够长的单轨迹采样。如果关键的αC/DFG转换需要超过单条replica的长度，所有轨迹可能都困在原来的basin。

> 轨迹长度应覆盖目标构象转换的时间尺度。EGFR那篇用的就是 3 × 1 μs，而不是用很多短轨迹代替长采样。

### 常用的一组结构标记

每条轨迹中，持续追踪以下结构量，最后画出**分布直方图**或**自由能面**进行比较。

#### 1. αC螺旋的位置

αC-helix有两种基本状态：**αC-in**（朝向催化位点内旋）和 **αC-out**（外旋）。量化的关键是**根据具体激酶的active/inactive参考结构定义合适的CV**，如αC-helix的RMSD到active reference，或αC orientation的几何参数。

> 如果holo组的分布明显向短距离方向偏移，说明配体促进了αC-in状态，这是活性构象的一个关键特征。

#### 2. DFG基序的取向

DFG（Asp-Phe-Gly）基序决定催化活性中心的组装。**DFG-in**（Asp朝向ATP结合位点）对应活性态，**DFG-out** 对应非活性态，中间态也存在。

DFG-in/out的判断基于**DFG-Phe是否进入R-spine** 和 **DFG-Asp是否定位到催化位点**，这是两个独立的几何事件。具体方法：用Phe侧链的χ二面角或centroid位置、DFG-Asp与Mg<sup>2+</sup>/ATP的距离，或者DFG整体相对于参考结构的构象CV来量化。相邻残基的Cα距离变化不大，不能表征DFG flip。

> 注意 **DFG-in不等于一定active**。某些Type-Iα抑制剂结合的inactive状态也是DFG-in，必须联合αC、R-spine、A-loop一起判断。

#### 3. A-loop的开放程度

催化环（activation loop，A-loop）的构象与活性密切相关。在活性态中，A-loop通常呈有序构象，覆盖在底物上方；不同激酶对A-loop状态的区分方式不同。p38α等研究用**active/inactive contact-map CV**来量化，而不是简单的“开/关”二值判断。用A-loop RMSD（关于active reference）或关键contact/distance来量化。

#### 4. 盐桥的形成（β3-Lys–αC-Glu）

αC-in构象的一个关键稳定因素是**β3 链上的Lys残基和αC-helix上的Glu残基**之间形成盐桥。用这两个残基的Cα距离或侧链重原子距离来量化。

> `distance < 4 Å` 通常认为是盐桥形成。可以比较apo和holo中**盐桥形成的时间占比**。

#### 5. R-spine的连续性

R-spine（regulatory spine）是由 **β4 疏水残基–αC疏水残基–DFG-Phe–HRD-His/Tyr** 组成的斜向疏水堆叠，并由 **αF-helix附近的保守Asp** 锚定，贯穿激酶的核心。**R-spine连续对应活性态**（各残基侧链紧密堆积），断裂（某个残基移出堆叠）对应非活性态。

R-spine描述的是 **side-chain hydrophobic packing**，因此用侧链接触、堆积几何或相对active reference的R-spine RMSD来量化，而不是简单规定Cα距离阈值。

#### 6. ATP口袋的预组织程度

这个指标比较特殊，因为它衡量的是**配体结合前的口袋状态**。如果配体结合后ATP口袋已经预先组织成结合构象（preorganized），说明结合过程不需要付出大的构象熵代价。

具体的量化方式：定义ATP口袋的若干关键残基（如P-loop上的残基、αC-Glu、 hinge区域残基），计算它们到ATP/配体上对应接触点的距离，然后和参考的晶体结构比较RMSD。**RMSD越小说明口袋越预组织**。

#### 7. 催化前几何

这个指标最直接，因为它衡量的是**磷酸转移反应能否发生**。如果模拟体系中有Mg<sup>2+</sup>-ATP和底物，可以用底物上待磷酸化的Ser/Thr/Tyr的氧原子到ATPγ-磷酸的距离来量化。

### 看分布，不是看平均值

这是容易出问题的地方：不能只报告平均距离。如果apo组的分布是双峰的（αC-in和αC-out两个峰），而holo组变成了单峰（主要是αC-in），那么平均距离当然会变，但更有意义的信息是峰的相对高度变化，也就是各状态的概率变化。

> 比如，apo组是双峰（两个态共存），holo组变成单峰（一个态占主导）。平均值的移动可能掩盖了分布形状的变化。正确的做法是画出两组的直方图或自由能面（$-k_B T \ln P(x)$），直接展示分布形状的差异。

### 变构路径追踪：从结合位点到活性位点

如果配体结合位点不在ATP口袋（比如是一个远端的变构位点），那么还需要追踪信号如何从结合位点传导到活性位点。

> 完整的传导路径是：**配体结合位点 → 局部结构变化 → αC/DFG/A-loop变化 → 催化几何变化**

每一步都需要结构证据支撑，形成一个**因果链**。文献中常用的追踪手段包括动态交叉相关矩阵（DCCM）来观察运动模式的变化，以及接触图（contact map）来观察哪些残基对之间的相互作用出现了系统性变化。

但这只是辅助手段。DCCM显示相关运动改变了，不等于功能改变了，**最终还是要回到结构标记的分布偏移上**。

## 文献中的具体案例

先说明一下后文反复出现的“MD追踪”到底在追踪什么。这里不是泛指“看轨迹”，而是持续统计一组可复现的观测量：**几何量**（关键距离、角度、二面角、RMSD）、**状态量**（如αC-in/αC-out、DFG-in/DFG-out、open/closed的帧占比）、**相互作用量**（盐桥/氢键/疏水接触占有率）和**动力学耦合量**（DCCM、接触网络重排）。判断是否“有机制”看的是这些指标在apo与holo之间是否出现方向一致、可重复的分布偏移，而不是某一条轨迹里出现一次构象切换。

### PAK1：自抑制片段释放与DFG-in状态（Cell 2026）

He等人（Cell, 2026, [DOI: 10.1016/j.cell.2026.03.008](https://doi.org/10.1016/j.cell.2026.03.008)）发现PAK1-A1小分子直接结合到自抑制区和kinase domain交界面的“autoinhibition-release site”，提高PAK1 kinase活性。PAK1这篇解读可以参考：[肽段引导策略发现PAK1变构激活剂（上）](https://mp.weixin.qq.com/s/qe0ko5E0M2i1YsDtTn6f2A)

**3条独立400 ns MD**（OpenMM，Amber ff14SB + GAFF）追踪三个量，三类量独立一致地指向同一个结论：**PAK1-A1稳定了DFG-in/αC-in的活性构象集合**。

- **DFG motif扭转角分布**：PAK1-A1保持**DFG-in**，NVS-PAK1-1（抑制剂对照）锁定**DFG-out**；
- **αC螺旋Glu315到activation loop Phe408的距离**：PAK1-A1相比NVS-PAK1-1有差异化的距离分布；
- **蛋白内接触图**：PAK1-A1占据DEK位点后KIS的接触模式改变，远离催化结构域，自抑制界面被破坏。

实验验证分三个层面：**SPR**直接测量结合亲和力，结合位点突变（K141D、V318D、N383L、V385D、D407K）显著降低SPR信号和FRET构象响应，验证了MD预测的配体接触；**FRET和HDX-MS**从动力学角度独立验证了自抑制片段释放和activation loop构象变化。同时，**ADP-Glo和RapidFire-MS**酶活实验确认这些小分子确实提高了PAK1的催化活性。

### MEK1：催化前几何直接解释磷酸化抑制

Mudedla等人（ACS Omega, 2024, [DOI: 10.1021/acsomega.4c03615](https://doi.org/10.1021/acsomega.4c03615)）研究MEK1/2结合selumetinib、trametinib、cobimetinib等变构抑制剂后的构象变化，用常规MD模拟，核心看一个量：**activation-loop上Ser222的氧原子到ATPγ-磷酸的距离分布**。

> 无抑制剂时，Ser222可以进入约5 Å的“可磷酸化区域”；加入变构抑制剂后，activation loop的**柔性被配体相互作用限制**，Ser222接近γ-磷酸的构象占比消失，直接解释了为什么抑制剂能阻止磷酸化。

分析方法为GROMACS常规MD + PCA（AMBER99SB-ILDN力场）。一个和催化反应直接挂钩的**distance observable**就建立了机制解释。

### EGFR L747P：一组突变的构象标记解释药物敏感性

Yoshizawa等人（npj Precision Oncology, 2021, [DOI: 10.1038/s41698-021-00170-7](https://doi.org/10.1038/s41698-021-00170-7)）用**3×1 μs MD加MP-CAFEE自由能计算**研究EGFR L747P突变对多代TKI的敏感性差异。L747P位于P-loop，MD追踪了三个活性构象标记的分布变化：**P-loop RMSF**反映loop的柔性改变，**αC螺旋取向**反映活性构象的偏好偏移，**K745–E762盐桥占有率**反映αC-in的稳定性。同时，MP-CAFEE计算的结合自由能解释了配体亲和力差异。构象标记和热力学两条证据线一致地解释了L747P对不同代TKI敏感性差异的结构基础。

> EGFR L747P案例的设计是**多副本MD → 活性构象标记分布 → 自由能计算 → 实验验证**四条线，构象和热力学互相独立地验证了同一个结论。

### LRRK2：抑制剂的open/closed状态调控

Schmidt等人（PNAS, 2021, [DOI: 10.1073/pnas.2100844118](https://doi.org/10.1073/pnas.2100844118)）和Weng等人（ACS Chemical Biology, 2023, [DOI: 10.1021/acschembio.2c00868](https://doi.org/10.1021/acschembio.2c00868)）的工作用MD + HDX-MS研究LRRK2结合抑制剂后的动力学变化。两个抑制剂作为**perturbation controls**：type-I抑制剂MLi-2和type-II抑制剂rebastinib稳定不同的kinase conformational state，MD追踪了kinase domain的“breathing”运动、**open/closed构象占比**（按domain间距离或开合角将每帧分为open或closed后统计比例）、**αC/activation-segment架构变化**（看αC-in/αC-out取向、β3-Lys–αC-Glu盐桥占有率和A-loop接触模式）和interaction networks。不同分子给出不同的**dynamical fingerprint**。

> “active-like conformation”不能直接等价于“活性升高”。功能的定义必须由实验数据决定。一个构象看起来像active不等于它真的催化效率高。

### p38α：用contact map代替RMSD描述活性状态

Kuzmanic等人（eLife, 2017, [DOI: 10.7554/eLife.22175](https://doi.org/10.7554/eLife.22175)）的研究（用ATP–Mg<sup>2+</sup>体系）的核心创新是用**active/inactive activation-loop的contact-map作为CV**来构建自由能面，而不是只报一个RMSD。具体来说，分析了**D112–ATP的催化几何**、**R49–D112和K53–E71盐桥**等active-state markers的概率分布，把A-loop的构象变化和DFG状态联系起来。采样上用Metadynamics加速采样。

### PKCα：从G-loop到远端调控位点的传播

Lippert等人（Journal of Biological Chemistry, 2021, [DOI: 10.1074/jbc.RA120.016470](https://doi.org/10.1074/jbc.RA120.016470)）研究ATP位点的抑制剂BimI和sotrastaurin如何破坏PKCα的远端调控。MD追踪**G-loop位移**、**domain间相互作用**和**局部/远端构象变化**，对应的可计算指标可以写成：G-loop相对口袋的质心距离或关键残基距离分布、domain界面的接触占有率/最小重原子距离、以及远端残基的RMSF与DCCM变化，并结合FRET传感器实验验证。ATP位点的抑制剂推开通常覆盖ATP口袋的Gly-rich loop，由此沿结构相互作用网络传播到远端调控界面，传导路径是**配体 → G-loop → 远端调控界面**。信号传播路径不一定沿着最短距离，而是沿着已有的结构相互作用网络。

### PDK1：同一个口袋，不同配体推向不同调控状态

Schulze等人（Cell Chemical Biology, 2016, [DOI: 10.1016/j.chembiol.2016.06.017](https://doi.org/10.1016/j.chembiol.2016.06.017)）研究PDK1的ATP位点和PIF pocket之间的变构关系。用了**unbiased MD + PT-WTE增强采样**，追踪**αB/αC动力学**、**hinge运动**、**Lys144盐桥切换**和**PIF pocket的柔性差异**。不同ATP位点配体可以增强或削弱PIF pocket的协同性，同一个口袋被不同分子推向相反方向的调控状态。**配体并非只是待在ATP口袋，而是在重新分配整个构象集合**。

## 突变验证：一种强有力的因果性检验

模拟结果显示holo组的构象偏移了，但这只是相关性。如果能进一步做突变模拟，因果链会更完整。

> 构象偏移 + 突变后偏移消失 + 实验活性变化一致，三者构成完整的因果链。突变MD本身仍然是计算证据，最理想的是同时有实验活性数据。

方法：将配体的关键接触残基突变（比如把形成氢键的Asp突变成Ala），重新跑holo模拟。如果突变后构象分布回到apo的水平，说明原来的偏移确实是由这些相互作用驱动的。但突变MD本身仍然是计算证据，最理想的验证是同时有实验活性数据，PAK1那篇正是这样做的。

这个策略和文献中的药物抗性突变分析是一致的：如果某个突变导致抑制剂失效，模拟应该显示突变后的构象分布回到apo状态。

## 证据强度排名

从最强到最弱：

1. **结构-动态-功能三联证据**：晶体结构确认结合模式 + MD显示构象偏移 + 体外活性数据验证，三者一致
2. **构象集合偏移 + 突变验证**：MD显示偏移，且突变后偏移消失，形成因果链
3. **构象集合偏移（多标记一致）**：多个活性标记同时向活性方向偏移
4. **单标记偏移**：只有一个标记变化，其他没变，说服力有限
5. **RMSD稳定、配体没解离**：只说明模拟合理，不构成机制证据

## 关键结论与批判性总结

### 主要影响

- **用构象分布而非单点结构来说活性**，这个思路已经在这个领域被充分验证，从PAK1到MEK1到EGFR都有对应的结构-动力学-功能证据链
- 变构效应是远端构象变化，晶体结构看不到，MD可以追踪
- **证据链的完整性**比单个图更重要：结合模式 + 构象偏移 + 突变验证 + 实验验证，四条线都要有。

### 存在的局限

- **采样不足**是MD最大的问题。100–500 ns的模拟对于αC和DFG的大尺度转变可能不够，这些转变的时间尺度可能是μs甚至更长。多副本可以部分缓解，但不能替代增强采样
- 分布偏移的方向和活性之间的**定量关系**不是通用的。不同激酶可能有更复杂的构象-活性映射，不能简单套用
- 这些结构标记的**阈值是经验性的，不同研究可能用不同的标准**，比较时要注意口径一致

### 实用建议

- 先跑apo的基线模拟，了解这个激酶本身的构象分布，再跑holo做比较，不能跳过apo直接看holo
- 每个体系至少 3 条独立轨迹，不要只跑一条
- 结果展示用分布图（直方图或自由能面），不要用单个数值
- 报告偏移量时，同时给出apo的分布宽度作为参考。如果偏移量小于apo的分布宽度，说明这个偏移可能是噪声
- **每个replica单独算分布，再报告均值 ± 置信区间**。MD帧之间高度时间相关，不能把10万个帧当成10万个独立测量。用block-averaging估计误差
- PCA和tICA可以用于可视化整体构象空间的变化，但**apo和holo的轨迹要一起分析（pooled）**，不要分开做PCA然后说“第一主成分变了”。分开做的话主成分方向本身就没法比较
