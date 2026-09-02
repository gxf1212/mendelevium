---
title: "Basit 2022：ECCR 电荷缩放 + MM-PBSA 计算钙调蛋白 Ca²⁺ 结合亲和力（技术细节）"
date: "2026-09-02"
last_modified_at: 2026-09-02
tags: [calmodulin, calcium, ECCR, charge-scaling, MM-PBSA, binding-free-energy, force-field, EF-hand, method-details]
description: "Basit 等（2022）用平均场电荷缩放（ECCR）引入极化、配合 MM-PBSA 与多重短轨迹预测钙调蛋白 17 个突变体的 Ca²⁺ 结合亲和力。本文逐节拆解其方法细节、电荷补偿规则、MM-PBSA 协议与结果，并对照 SI 核实关键数字。"
author: Xufan Gao
lang: zh-CN
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/Wallpaper_compressed/cyber-3400789_1920.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/Wallpaper_compressed/cyber-3400789_1920.jpg"
---

# Basit 2022：ECCR 电荷缩放 + MM-PBSA 计算钙调蛋白 Ca²⁺ 结合亲和力（技术细节）

## 方法概括

Basit等2022年的目标不是用严格的ABFE直接计算Ca²⁺结合自由能，而是建立一个**低成本的CaM突变体亲和力预测模型**。体系为human calmodulin，初始结构采用Ca²⁺饱和的晶体结构PDB 1CLL；研究重点是C-lobe的两个EF-hand，即loop-3和loop-4。所有突变体直接由1CLL对应残基突变得到，质子化状态由PDB2PQR在中性pH下确定。

### 1. 基础MD体系

- AMBER 16
- 蛋白：ff14SB
- 水：TIP3P
- 水盒边界距离：10 Å
- 盐：0.15 M KCl
- K⁺/Cl⁻基础参数：Joung–Cheatham
- Ca²⁺基础LJ参数：Li–Merz 2013
- 温度：298 K
- Langevin thermostat，collision frequency = 2 ps⁻¹
- Berendsen barostat，1 bar
- PME处理长程静电
- SHAKE约束含H键，时间步长2 fs

他们比较了2、20和200 ns三种采样长度。最初WT、N97S和F141L跑200 ns，并分析后100 ns；随后发现10条20 ns与10条2 ns轨迹给出的平均MM-PBSA结果接近，因此最终17个突变体主要采用**10条独立2 ns轨迹**。

这个设计本质是“multiple short replicas”，目的是降低平均值的standard error，而不是充分探索CaM的慢构象变化。作者自己也承认2 ns不能充分描述CaM构象涨落。

---

## 2. ECCR到底改了什么？

这篇并没有重新拟合新的Ca²⁺ LJ参数。

它以Li–Merz Ca²⁺参数为基础，只修改静电部分：

| 离子 | 标准电荷 | ECCR |
|---|---:|---:|
| Ca²⁺ | +2.00 | +1.50 |
| K⁺ | +1.00 | +0.75 |
| Cl⁻ | −1.00 | −0.75 |

即统一乘0.75，用平均场方式近似电子极化。

### 蛋白如何补偿少掉的+0.5e？

这里SI给得很具体。

一个Ca²⁺由+2.0变成+1.5，相当于体系少了+0.5e。因此作者把每个EF-hand中Ca附近的一组O原子**统一调得稍微不那么负**，总共补回+0.5e。

例如loop-3选了：

- Asp93、Asp95：两个羧酸O均修改
- Asn97：侧链O
- Tyr99：O
- Glu104：两个羧酸O均修改

一共8个O，因此每个O增加约：

+0.5 / 8 = +0.0625e

所以：

- Asp O：−0.8014 → −0.7389
- Asn O：−0.5931 → −0.5306
- Tyr O：−0.5679 → −0.5054
- Glu O：−0.8188 → −0.7563

Loop-2也是8个O，每个+0.0625e。

Loop-1和loop-4各选择9个O，因此每个增加约：

+0.5 / 9 ≈ +0.0556e

例如loop-4：

- Asp O：−0.8014 → −0.7458
- Gln O：−0.5679 → −0.5123
- Glu O：−0.8188 → −0.7632

所以这篇的“charge compensation”规则其实非常简单：

**Ca²⁺降低0.5e → 把+0.5e平均分摊给该EF-hand预先指定的Ca附近O原子**。

而且Asp/Glu的**两个羧酸O都修改**，并不只修改晶体结构中直接朝向Ca²⁺的那个O。完整逐原子数值就在SI Table S1。

![图5 loop-3与loop-4结合位点的残基示意](./Basit-2022-ECCR-MM-PBSA_figs/main_p7_i1_1750x1016.png)

图5（论文 Figure 5）：C-lobe 两个结合位点的残基排布（单字母缩写 + 残基编号，W 为水分子）；本文做电荷补偿的氧原子，即图中直接配位 Ca²⁺ 的那些侧链氧（Asp/Glu 的羧酸氧、Asn/Gln 的酰胺氧、Tyr/Thr 的酚/醇氧）。

这点对你的体系很重要：它不是一个通用“看到Ca附近多少个O就自动重新拟合”的算法，而是人为定义一组EF-hand邻近O，然后做总电荷补偿。

---

## 3. MM-PBSA怎么计算？

使用AMBER的MMPBSA.py，采用single-trajectory protocol：receptor和Ca²⁺轨迹都从同一条CaM–Ca²⁺复合物轨迹中提取。

基本形式是：

ΔG_bind ≈ ΔE_MM + ΔG_solv − TΔS

但这篇**直接忽略TΔS**。

能量包括：

ΔE_MM = ΔE_vdW + ΔE_el

溶剂化部分包括：

ΔG_solv = ΔG_PB + ΔG_cavitation + ΔG_dispersion

PB计算：

- ionic strength = 0.15 M
- solvent dielectric = 78.35
- **solute/internal dielectric = 8**

其中ε_in = 8不是无先验选择。作者明确说，他们在之前的CaM工作中比较后发现8效果最好，所以本论文沿用8。

因此这套protocol已经包含了一项针对CaM体系的经验校准。

---

## 4. 两个EF-hand的协同性怎么处理？

实验给的是CaM C-lobe两个位点的macroscopic binding affinity。

作者分别计算四种site-specific过程：

- ΔG₃：loop-4为空时，Ca结合loop-3
- ΔG₄,₃：loop-3已有Ca时，Ca结合loop-4
- ΔG₄：loop-3为空时，Ca结合loop-4
- ΔG₃,₄：loop-4已有Ca时，Ca结合loop-3

最终按：

ΔG_bind = (ΔG₃ + ΔG₄,₃ + ΔG₄ + ΔG₃,₄) / 4

得到每个位点平均的宏观结合自由能（即论文里的 ΔG_bind/2，与实验报告的每位点平均值 −7.64 kcal/mol 直接可比；论文的 ΔG_bind 全量是它的两倍）。SI Figure S1 给出了完整的计算流程。

![图2 钙调蛋白C-lobe两个EF-hand（loop-3、loop-4）的两步Ca²⁺结合自由能图](./Basit-2022-ECCR-MM-PBSA_figs/main_p4_i1_996x1018.png)

图2（论文 Figure 2）：白色圆圈为空位点（左为 loop-3），红色圆圈为已结合 Ca²⁺ 的位点；ΔGx,y 表示当位点 y 已被占据时 Ca²⁺ 结合位点 x 的自由能变化。

![图S1 钙调蛋白Ca²⁺结合自由能计算流程图](./Basit-2022-ECCR-MM-PBSA_figs/si_p3_i1_1338x947.jpeg)

图S1（论文 Supporting Information Figure S1）：从复合物轨迹提取受体与配体、四种 site-specific 过程到宏观结合自由能的完整流程。

---

## 结果到底有多好？

这里必须把**原始MM-PBSA结果**和**后面的回归拟合结果**分开。

## 1. ECCR确实明显改善了普通ff14SB

以10 × 2 ns结果为例：

| 系统 | STD | ECCR | 实验 |
|---|---:|---:|---:|
| WT | −11.04 | **−8.10** | −7.64 |
| N97S | −9.90 | **−8.46** | −6.81 |
| F141L | −11.61 | **−8.19** | −6.62 |

单位均为kcal/mol。

所以对这几个体系，标准+2e Ca明显overbinding；简单charge scaling确实把结果往实验方向拉了很多。

SI Table S7给出了全部体系的ECCR结果。例如：

| 系统 | ECCR ΔG | 实验 |
|---|---:|---:|
| WT | −8.10 | −7.64 |
| D95V | −6.14 | −6.09 |
| G113R | −7.28 | −6.95 |
| D131E | −8.77 | −5.89 |
| Q135P | −9.25 | −6.58 |
| E140K | −4.61 | −5.62 |



我直接根据Table S7重新统计全部WT+17 mutants，**在任何后续回归拟合之前**：

- absolute ΔG MAE ≈ **1.20 kcal/mol**
- absolute ΔG RMSE ≈ **1.42 kcal/mol**
- 平均bias ≈ **−1.09 kcal/mol**

也就是说，ECCR+MM-PBSA本身其实已经不算差，但整体仍有约1 kcal/mol的overbinding，而且D131E、Q135P、D131V等有2–3 kcal/mol的大误差。

---

## 2. 真正比较突变效应时，原始结果没有图上那么漂亮

论文真正关心的是：

ΔΔG = ΔG_mutant − ΔG_WT

因为实验主要研究突变导致的亲和力改变。

**拟合之前**，Figure 9的原始ECCR+MM-PBSA结果是：

- RMSE = **1.15 kcal/mol**
- Pearson r = **0.51**
- Spearman r = **0.55**
- 17个突变中有**4个连ΔΔG正负号都预测错了**

作者自己也明确称此时相关性只是moderate。

所以如果你看到后面的：

- RMSE ≈ 0.3 kcal/mol
- r ≈ 0.8

那已经不是原始MM-PBSA结果。

---

## 后面确实拟合了，而且拟合得相当明显

作者把MM-PBSA分解得到的：

- ΔΔE_vdW
- ΔΔE_el
- ΔΔG_PB
- ΔΔG_cavitation
- ΔΔG_dispersion

拿去对实验ΔΔG做**多元线性回归**。

五描述符模型的相关性一下提高到：

- Spearman r = 0.82
- Pearson r = 0.84

但作者自己检查后发现，除了截距外，五个descriptor的系数都不显著，因此明确说这个模型**不能作为predictive model**。

之后逐步删descriptor。

最后得到的“最佳模型”甚至只剩一个量：

**ΔΔG_model = 0.9863 − 0.3083 × ΔΔE_vdW**

也就是说，最后那个很漂亮的模型，本质上已经不是MM-PBSA自由能，而是：

> 用MD得到的vdW能量变化作为descriptor，再对实验亲和力重新拟合一个线性模型。

它得到：

- RMSE ≈ **0.34 kcal/mol**
- Spearman r ≈ **0.80**
- Pearson r ≈ **0.80**

SI中的多描述符模型甚至有RMSE≈0.31 kcal/mol。

作者也做了leave-3/4/5-out cross-validation：

| CV | test RMSE | test Spearman | test Pearson |
|---|---:|---:|---:|
| leave-5-out | 0.33 | 0.75 | 0.75 |
| leave-4-out | 0.32 | 0.73 | 0.72 |
| leave-3-out | 0.32 | 0.70 | 0.69 |

所以不是完全裸拟合，但仍然是在**同一批CaM点突变数据**内部做cross-validation，没有真正独立的外部蛋白或另一批CaM突变测试集。

---

## 我怎么评价这篇的严谨性？

我的评价是：**作为“CaM突变体快速预测模型”，做得比较认真；作为“证明ECCR+MM-PBSA能准确计算Ca²⁺物理结合自由能”，证据远没有摘要看起来那么强**。

### 做得好的地方

1. **没有只跑一条轨迹**。  
   最终采用10条独立trajectory，并专门比较2、20和200 ns方案。

2. **STD和ECCR正面对照**。  
   可以清楚看到+2e固定电荷模型确实严重overbind，而ECCR明显改善。

3. **charge modification完全公开**。  
   SI甚至给到逐原子partial charge，复现性很好。

4. **考虑了两个EF-hand的occupancy/cooperativity**。  
   不是简单把一个Ca删除然后算一次PBSA。

5. **作者没有掩盖原始结果**。  
   Figure 9明确告诉你原始ECCR只有r≈0.5、RMSE≈1.15，而且有4个sign error。

6. **回归部分至少做了显著性检验和cross-validation**。  
   他们甚至主动否掉了看起来相关性更高、但统计上不显著的五参数模型。

这些都算比较规范。

---

## 但有几个明显限制

### 1. “RMSE 0.3 kcal/mol”不能拿来证明MM-PBSA准

这是最重要的一点。

真正未经实验拟合的ECCR+MM-PBSA是：

**RMSE 1.15 kcal/mol，r≈0.5**。

0.3 kcal/mol那个结果来自**用实验ΔΔG训练后的回归模型**。

因此不能写：

> ECCR-MM/PBSA predicts Ca binding affinity with an RMSE of 0.3 kcal/mol.

更准确应该写：

> Raw ECCR/MM-PBSA achieved an RMSE of about 1.15 kcal/mol for relative binding free energies; subsequent empirical regression against experimental data reduced the cross-validated error to approximately 0.3 kcal/mol.

两件事不是一个概念。

### 2. ε_in = 8本身已经经过CaM体系校准

他们没有blindly采用一个标准PB参数。

论文明确说ε_in = 8是之前在类似CaM体系比较后得到的最佳值。

因此“原始MM-PBSA”其实也不是完全无经验调整的ab initio prediction。

### 3. 最终模型只剩vdW descriptor，有点反常

最后最成功的预测式是：

ΔΔG_model = 0.9863 − 0.3083 × ΔΔE_vdW

Ca²⁺结合明明是高度electrostatic的问题，最后统计模型却只保留了vdW变化。

这不代表“Ca²⁺结合由vdW主导”，只是说明在这一小组高度相似的CaM mutants中，ΔE_vdW恰好是最稳定的经验descriptor。

所以它更像**QSAR式校准模型**，不能直接赋予强物理解释。

### 4. 2 ns轨迹很短

10 × 2 ns确实能让standard error变小，但：

**很多短trajectory ≠ 已经采样到慢构象转换**。

作者自己也承认2 ns不能充分描述CaM结构涨落。

对这些从同一个holo晶体结构出发的点突变，可能还能工作；换成一个磷酸化可能改变EF-hand构象ensemble的体系，风险明显更大。

### 5. single-trajectory MM-PBSA没有真正模拟结合过程

蛋白、Ca和complex都从同一个bound trajectory提取，因此大量构象重组代价被假定相消。

它实际上问的是：

> 在已有holo-like构象附近，这个Ca的endpoint interaction energetics怎样？

不是严格意义上的：

> Ca从bulk water进入EF-hand的标准结合自由能是多少？

所以不能和umbrella sampling、double-decoupling ABFE等严格自由能方法等价。

### 6. 没算构象熵

作者认为WT和mutants相近，因此希望熵项相消。对于普通point mutants可能有一定合理性。

但如果你的磷酸化真正改变EF-hand opening、两臂运动或apo/holo构象平衡，这个假设就更危险。

---

## 对这篇论文最合适的定位

我觉得它真正证明的是：

> **把Ca²⁺从+2e缩到+1.5e，同时对EF-hand附近O做局部charge compensation，确实能显著缓解ff14SB/TIP3P中Ca²⁺的过强静电结合，并改善CaM突变体之间的相对亲和力预测**。

它没有证明：

> **ECCR + MM-PBSA已经可以无参数、定量预测任意EF-hand的绝对Ca²⁺结合自由能**。

未经实验回归时，它的表现其实是一个很现实的水平：absolute ΔG大约**1.2 kcal/mol MAE**，relative ΔΔG **1.15 kcal/mol RMSE、r≈0.5**。这已经值得作为你的一个benchmark方法，但我不会把最终0.3 kcal/mol那个数字拿来给ECCR力场本身背书。

对你的体系反而有个好消息：如果pSer/pThr确实离Ca第一配位层较远，你不需要照这篇去改磷酸基。真正需要决定的是，是否像他们一样把**每个EF-hand第一配位环境附近的O统一分摊+0.5e**。这个细节比“Ca设成+1.5”本身更值得谨慎复现。