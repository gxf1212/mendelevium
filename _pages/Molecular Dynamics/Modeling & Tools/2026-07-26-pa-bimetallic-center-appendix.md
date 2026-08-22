---
title: "PAN双金属中心参数化：方法细节"
date: "2026-07-26"
last_modified_at: "2026-07-26"
tags: [molecular-dynamics, force-field, parameterization, metal, bimetallic, amber, pa-endoribonuclease, baloxavir]
description: "PAN双金属中心参数化方法细节，包括四类金属参数化方案、MD模拟协议与分析方法"
image: "/assets/img/thumbnail/bricks.webp"
thumbnail: "/assets/img/thumbnail/bricks.webp"
author: Xufan Gao
lang: zh-CN
---

# PAN双金属中心参数化：方法细节

## 方法设计：金属参数化策略与模拟协议

### 整体工作流

本文采用"基准测试→应用"的二阶段策略：

1. **WT PAN基准测试**：在apo和BXA结合态下，分别测试`12-6`非键、`12-6-4`非键、键合、杂化四类参数化方案，评估500 ns模拟中的结构稳定性、配位几何保持、配体动力学等
2. **突变体应用**：基于基准测试结果，对I38T/F/M、A36V、E23K五种突变体分别在apo和BXA结合态进行3副本500 ns MD模拟（共30条轨迹），通过构象采样和相互作用分析揭示耐药机制

### 系统构建

所有系统基于A/California/04/2009（H1N1）毒株的PAN序列，使用Modeller 10.6重建缺失的51-72号柔性loop（基于PDB 4LN7、4M4Q、9DOJ模板）。缺失loop结构对参数化基准测试影响很大，因为它是apo态骨架RMSD的主要波动源——文章通过"排除loop"和"包含loop"两种RMSD计算方式分别报告。

### 四类金属参数化方案

#### `12-6`非键模型

将金属-配体相互作用完全交给经典非键项：`12-6` LJ描述范德华，库仑描述静电。$\mathrm{Mn^{2+}}$参数采用针对Amber TIP3P水模型优化的`12-6` normal usage（compromised）参数集。这是**最简单**的方案，物理上假设金属-配体相互作用不涉及显著的电荷转移或诱导偶极。

#### `12-6-4`非键模型

在`12-6` LJ基础上引入 $C_4/r^4$ 项，描述离子-诱导偶极相互作用。$\mathrm{Mn^{2+}}$与水、咪唑氮、羧基氧的 $C_4$ 项来自Li/Merz组和Jafari组等专门针对TIP3P开发的工作。$\mathrm{Mn^{2+}}$与$\mathrm{OH^-}$的 $C_4$ 项基于 $\log K = 3.452$ 的伞形采样拟合。这一方案在物理上比`12-6`更准确，特别是对二价金属的强极化效应。

#### 键合模型

对金属-配位原子施加完整的谐振势（键、角、二面角），通过MCPB.py工作流构建：Seminario方法从QM频率计算获得键和角力常数，RESP拟合从大模型获得部分电荷。为避免两个配位水对同一$\mathrm{Mn^{2+}}$产生的1-4相互作用问题，**保留谐振约束而非固定键合**，以维持配位距离在平衡位置附近的同时允许几何柔性。

#### 杂化模型

仅对金属-配位氨基酸施加键合约束，**配位水仍保持完全非键**且可与溶剂水交换。这一策略源于一个关键观察：BXA结合态需要保持金属-氨基酸的几何，配位水则可灵活交换。文章对$\mathrm{Mn^{2+}}$和配位氨基酸采用RESP拟合（基于大模型ESP），而BXA和水的部分电荷固定为非键模型值，并移除金属-配位水的谐振约束。

### MD模拟协议

所有系统在Amber ff14SB力场下进行：

- **等距TIP3P水盒溶剂化**：溶质到盒子边缘最小距离13 Å
- 添加$\mathrm{Na^+}$中和并补充0.15 M NaCl，模拟生理离子强度
- **50,000步能量最小化**：约束从50降到 $0\,\mathrm{kcal\cdot mol^{-1}\cdot \mathring{A}^{-2}}$
- **350 ps NVT 升温至310 K**：约束从20降到 $0\,\mathrm{kcal\cdot mol^{-1}\cdot \mathring{A}^{-2}}$
- **1 ns NPT 平衡**
- **500 ns NPT 生产模拟**：每系统3独立副本
- **长程静电**：PME，截断10 Å
- **温度耦合**：Langevin thermostat，碰撞频率 $1.0\,\mathrm{ps^{-1}}$
- **键约束**：SHAKE

### 分析方法

通过cpptraj模块分析：

- **配位几何**：径向分布函数（RDF，0.01 Å分辨率）积分第一峰获得配位数，以第一极小值作为配位边界
- **结构波动**：backbone RMSD、residue-level RMSF
- **构象采样**：主成分分析（PCA）+ KDE等高线，mean-shift聚类识别代表构象
- **相互作用**：通过ProLIF编码为指纹图谱，统计hydrophobic、π-stacking、HBA/HBD、π-cation、VdWContact等相互作用occupancy
- **口袋体积**：POVME 3.0
