---
title: "用 PyMOL 和 UCSF ChimeraX 能补建蛋白二硫键、优化侧链吗？"
date: "2026-07-30"
last_modified_at: 2026-07-30
tags: [pymol, chimerax, protein-structure, disulfide-bond, rotamer, structural-editing]
description: "PyMOL 和 UCSF ChimeraX 在蛋白结构编辑中的能力边界：哪些步骤两个软件都能胜任，哪些必须交给 Rosetta、ISOLDE 或分子动力学"
image: "/assets/img/thumbnail/bricks.webp"
thumbnail: "/assets/img/thumbnail/bricks.webp"
author: Xufan Gao
lang: zh-CN
---

# 用 PyMOL 和 UCSF ChimeraX 能补建蛋白二硫键、优化侧链吗？

## 1. 链间天然二硫键的补建

### PyMOL

PyMOL 可以完成指定残基之间的共价键构建，因此能够用于补建同源二聚体蛋白特异性的链间二硫键。以两条链对应的 Cys89 和 Cys98 为例，通过连接两个半胱氨酸侧链的硫原子（SG）建立 S–S 键：

```
bond chain A and resi 89 and name SG, chain B and resi 89 and name SG
bond chain A and resi 98 and name SG, chain B and resi 98 and name SG
```

随后需要将参与二硫键形成的半胱氨酸残基**由还原态 CYS 修改为氧化态 CYX**，以匹配多数分子模拟力场中的二硫键残基定义。

> PyMOL 的 `bond` 命令**只建立拓扑连接，不优化几何**。补建二硫键后必须自行检查 S–S 键长（约 2.0–2.1 Å）、Cβ–Sγ–Sγ–Cβ 二面角，以及两条链之间是否存在新的原子碰撞。

### UCSF ChimeraX

UCSF ChimeraX 也支持直接编辑蛋白共价连接，`bond` 命令连接两个半胱氨酸 SG 原子：

```
bond #1/A:89@SG #1/B:89@SG
bond #1/A:98@SG #1/B:98@SG
```

也可以通过 Structure Editing 工具在图形界面中选择两个硫原子并建立新键。

UCSF ChimeraX 在结构编辑上提供了比 PyMOL 更完善的交互支持，**可以显示新增共价键、调整键长、检查局部空间关系**。如果使用 ISOLDE 扩展，还可以基于分子力学约束实时调整结构，使二硫键几何参数更加合理。

## 2. 链间结构拼接后的侧链优化

### PyMOL

PyMOL 可以通过 Mutagenesis Wizard 进行单个残基的侧链重建和 rotamer 替换。

该工具根据 Dunbrack rotamer library 提供不同侧链构象，并计算侧链与周围原子的空间冲突，帮助选择 **steric clash 较小的构象**。因此，如果拼接区域仅存在少数残基的侧链冲突，可以使用 PyMOL 手动选择更合理的 rotamer。

> 但 PyMOL 不具备自动化的蛋白界面侧链重新堆积（side-chain repacking）能力，也不进行全局能量最小化——**无法自动寻找整个接口区域的最低能构象**，无法系统解决 backbone 偏移导致的局部结构缺陷，也无法替代 Rosetta Relax、ISOLDE 或分子动力学优化。

### Rotamers 工具

UCSF ChimeraX 内置 Rotamers 工具，可用于蛋白侧链构象优化，**在链间拼接场景下比 PyMOL 更适合作为侧链修正方案**。工作流程为：选择存在空间冲突的残基 → 调用 rotamer library 生成候选构象 → 根据 steric clash、氢键作用和局部几何关系筛选 → 替换原始侧链。

对于拼接区域中二硫键附近的侧链碰撞、疏水残基堆积或极性残基氢键网络异常，可以使用 Rotamers 进行局部修正。**这些不涉及 backbone 调整**。较大的构象缺陷仍需交给 Rosetta FastRelax / InterfaceRelax、ISOLDE 实时分子动力学优化，或 OpenMM / GROMACS 局部能量最小化。

## 综合评价

| 功能 | PyMOL | UCSF ChimeraX |
| ----- | ----- | -------------- |
| 指定 Cys–Cys 二硫键建立 | ✓ | ✓ |
| 修改 CYS/CYX 状态 | ✓ | ✓ |
| 检查 S–S 几何合理性 | 手动 | 更方便 |
| 单残基 rotamer 优化 | ✓ | ✓（更推荐） |
| 多残基界面侧链重打包 | ✗ | ✗ |
| 自动消除空间冲突 | 有限 | 有限 |
| 局部结构 relaxation | ✗ | 需 ISOLDE 扩展 |
| 全局能量优化 | ✗ | ✗ |
