---
title: "【笔记整理|2024-07】力场与分子建模：从Martini 3粗粒化到OPLS-AA全原子力场"
date: "2024-07-01"
tags: [force-fields, martini, opls-aa, molecular-modeling, coarse-graining, parameterization, technical-notes]
description: "全面解析Martini 3粗粒化力场和OPLS-AA全原子力场的参数化策略、珠子类型系统、介电常数设置和实际应用，为分子建模提供力场选择指南"
thumbnail: "/assets/img/thumbnail_mine/wh-dp5x3l.jpg"
image: "/assets/img/thumbnail_mine/wh-dp5x3l.jpg"
---

# 【笔记整理|2024-07】力场与分子建模：从Martini 3粗粒化到OPLS-AA全原子力场

## 引言

分子力场是分子动力学模拟的基石，不同的力场适用于不同的研究目的和应用场景。本文整理了从技术讨论中提取的关于Martini 3粗粒化力场、OPLS-AA全原子力场以及其他重要力场系统的关键知识和实用技巧，涵盖力场选择、参数化策略和应用实践。

## Martini 3粗粒化力场

### Martini 3设计理念

Martini 3是目前最先进的粗粒化力场之一，其设计理念基于系统性的参数化策略：

> The OPLS-AA force field has followed a consistent philosophy throughout the course of its development. Nonbonded parameters are optimized to reproduce experimental liquid phase properties, and torsional parameters are fit to available experimental or quantum chemical data.

> The Martini FF has been parametrized towards dielectric screening constant of 15, part of the electrostatic interactions have been included in the LJ parameters, therefore changing the screening constant would mean that you would also have to parametrize the LJ interactions. In short we would not advice fiddling with the screening.

### Martini 3珠子类型系统

Martini 3采用了系统的珠子类型命名和分类方法：

> 能否用图片和表格详细总结martini3的beads类型、命名方式的传统，主要参考Supporting information for: Martini 3: A General Purpose Force Field for Coarse-Grained Molecular Dynamics，其中的Supplementary Notes A) Description of beads and Martini 3 model和Supplementary Notes B) Parametrization strategy of the beads，还有其他的可能的资料，如：Table 1: Small molecule building blocks. CG particle type, corresponding chemical building block, and examples of molecules in which such a block appears.

**珠子类型示例：**
> The magnesium ion is represented by one TQ3p bead with a charge of +1

### Martini 3参数化资源

Martini 3提供了丰富的参数化资源和数据库：

> https://github.com/Martini-Force-Field-Initiative/M3-Sterol-Parameters/blob/main/martini_v3.0_sterols_v1.0.itp

> https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters

> https://github.com/ricalessandri/Martini3-small-molecules/tree/main

### Martini 3蛋白质-配体结合模拟

Martini 3在蛋白质-配体结合模拟方面具有独特优势：

> CHAPTER 1 A PRACTICAL INTRODUCTION TO MARTINI 3 AND ITS APPLICATION TO PROTEIN-LIGAND BINDING SIMULATIONS

### Martini 3介电常数

Martini 3的介电常数设置是其重要特征：

> There is actually an option in the mdp file to change the dielectric

## OPLS-AA全原子力场

### OPLS-AA设计哲学

OPLS-AA力场具有明确的参数化哲学和一致性原则：

> The OPLS-AA force field has followed a consistent philosophy throughout the course of its development. Nonbonded parameters are optimized to reproduce experimental liquid phase properties, and torsional parameters are fit to available experimental or quantum chemical data.

### OPLS-AA参数转换

OPLS-AA力场的参数在转换为GROMACS格式时需要注意一些细节：

> https://github.com/leelasd/OPLS-AAM_for_Gromacs/tree/master

> parmed CharmmParameterSet, all bonds,angles,dihedrals have two copied, where atom names are reversed, so we don't need to sort?

### PolyParGen聚合物参数化

PolyParGen为聚合物和大分子提供OPLS-AA和Amber力场参数：

> PolyParGen provides OPLS-AA and Amber force field parameters for polymers or large molecules. In the case that PolyParGen generates OPLS-AA parameters...

## 分子力场参数化

### 参数化策略

不同力场采用不同的参数化策略，需要根据研究需求选择：

> We can use mols2grid to display and scroll through the cluster samples

### 力场参数文件格式

力场参数文件的格式和结构对于正确使用力场至关重要：

> vmd modeling top_opls_aam.inp problematic IC: VAL, ILE, MET, CYS, PRO....

> vdwGeometricSigma yes

### 排除约束设置

合理的排除约束设置是力场配置的重要部分：

> For the [ exclusions ] section:

> For the [ constraints ] section:

> Extra exclusions within a molecule can be added manually in a [ exclusions ] section. Each line should start with one atom index, followed by one or more atom indices. All non-bonded interactions between the first atom and the other atoms will be excluded.

## 特殊相互作用与拓扑处理

### 质子海绵效应

质子海绵效应在分子模拟中是一个特殊的现象：

> proton sponge effect

### 受限弯曲势能

受限弯曲势能用于模拟特殊的分子结构：

> https://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html#restricted-bending-potential

### 虚拟位点

虚拟位点是分子力场中用于简化计算的重要技术：

> https://manual.gromacs.org/current/reference-manual/functions/interaction-methods.html#virtualsites

## 力场兼容性与转换

### 不同力场的兼容性

不同力场之间的兼容性是混合模拟中的关键问题：

> WARNING 3 [file ../../mdps_cg_78.4_mem/em.mdp]:

> ERROR 1 [file ../../mdps_cg_78.4_mem/nvt_neutral.mdp]:

### 力场参数验证

力场参数的验证确保模拟的可靠性：

> WARNING 4 [file system.top, line 13]:

### 力场组合使用

在某些情况下，需要组合使用不同的力场：

> 36 1 makes vmd output "psfgen) Created by CHARMM version 36 1"

> not useful in FEbuilder

## 分子建模工具与技术

### SMARTS模式匹配

SMARTS模式匹配是分子结构识别的重要工具：

> SMARTS matching

> emm, cannot ensure won't cause the same problem as rdkit

### 分子体积计算

分子体积计算是分子表征的重要参数：

> https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html#rdkit.Chem.AllChem.ComputeMolVolume

> from rdkit.Chem import rdMolDescriptors

### 分子表示与立体化学

立体化学的正确表示对分子模拟至关重要：

> Stereogenic centers belonging to an AND n group (e.g. AND1) represents a mixture of two enantiomers: the structure as drawn AND the epimer in which the stereogenic centers have the opposite configuration. (Note, that it is not a racemic mixture, but a mixture of the enantiomers of any ratio. Of course, a 1:1 mixture (racemic mixture) is included in this sense.)

## 特殊分子系统

### 膜蛋白与去垢剂

膜蛋白的模拟需要特殊的去垢剂处理：

> In addition, many proteins (especially membrane proteins) would aggregate if the SDS were simply washed out, this could lead to loss of activity. Non-ionic detergents like Triton solubilise proteins gently, often maintaining its activity.

### 荧光染料特性

荧光染料在生物物理研究中具有广泛应用：

> FITC reacts with a primary amine on the protein to form a covalent amide bond.

> Hoechst dyes are cell membrane-permeant, minor groove-binding blue fluorescent DNA stains. These dyes are widely used in cell cycle and apoptosis studies as nuclear counterstains.

### 圆二色谱计算

圆二色谱（CD）是研究蛋白质二级结构的重要技术：

> The DichroCalc web server [38] was used to calculate CD spectra from molecular

## 自由能计算与力场应用

### 软核相互作用

自由能计算中的软核相互作用避免奇点问题：

> https://manual.gromacs.org/current/reference-manual/functions/free-energy-interactions.html#soft-core-interactions-beutler-et-al

### 自由能计算工具

专业的自由能计算工具提高了模拟效率：

> https://github.com/delphi001/DelphiPka

> https://rowansci.com/tools/pka

> https://github.com/mms-fcul/PypKa

> https://valdes-tresanco-ms.github.io/gmx_MMPBSA/v1.5.5/command-line/

### 自由能计算标准流程

标准化的自由能计算流程确保结果的可比性：

> https://alchemistry.org/wiki/Exponential_Averaging

## 力场发展与前沿趋势

### 新兴力场系统

力场技术不断发展，出现了许多新兴的力场系统：

> https://www.bohrium.com/notebooks/38543442597

### 开源力场项目

开源力场项目促进了力场技术的普及和发展：

> https://github.com/OpenFreeEnergy/openfe-benchmarks

> https://github.com/drazen-petrov/SMArt

> https://github.com/OpenFreeEnergy/konnektor

### 商业力场软件

商业力场软件提供了专业的技术支持和服务：

> NVIDIA NIM for Boltz-2

> https://qsimulate.com/documentation/fep_tutorial/fep_tutorial.html

## 力场验证与质量控制

### 力场验证标准

力场验证是确保模拟结果可靠性的关键步骤：

> math font still use normal

### 力场参数数据库

力场参数数据库为研究人员提供了丰富的资源：

> https://www.wiredchemist.com/chemistry/data/metallic-radii

### 力场性能评估

力场性能评估帮助选择最适合的力场：

> https://www.r-ccs.riken.jp/labs/cbrt/tutorial/remd-tutorials/tutorial-2-1/

> https://manual.gromacs.org/current/reference-manual/analysis/correlation-function.html

## 总结与最佳实践

1. **力场选择**：根据研究目的选择合适的力场系统，Martini 3适合大系统长时间尺度，OPLS-AA适合高精度全原子模拟
2. **参数化策略**：理解不同力场的参数化哲学，确保参数的一致性和可靠性
3. **兼容性考虑**：在混合力场模拟中，充分考虑不同力场之间的兼容性问题
4. **验证流程**：建立完善的力场验证流程，确保模拟结果的可靠性
5. **工具使用**：熟练使用各种力场建模和分析工具，提高研究效率
6. **前沿跟踪**：关注力场技术的最新发展，及时更新知识体系
7. **质量控制**：建立严格的质量控制标准，确保研究成果的可重复性
8. **社区参与**：积极参与开源力场项目，促进力场技术的发展

通过这些力场知识和建模技巧的掌握，可以显著提高分子动力学模拟的质量和效率。

## 参考资源

- [Martini 3固醇参数](https://github.com/Martini-Force-Field-Initiative/M3-Sterol-Parameters/blob/main/martini_v3.0_sterols_v1.0.itp)
- [Martini 3脂质参数](https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters)
- [Martini 3小分子参数](https://github.com/ricalessandri/Martini3-small-molecules/tree/main)
- [OPLS-AA for GROMACS](https://github.com/leelasd/OPLS-AAM_for_Gromacs/tree/master)
- [gmx_MMPBSA手册](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/v1.5.5/command-line/)
- [DelphiPKa](https://github.com/delphi001/DelphiPka)
- [PypKa](https://github.com/mms-fcul/PypKa)
- [RowanSci pKa工具](https://rowansci.com/tools/pka)
- [自由能计算指数平均方法](https://alchemistry.org/wiki/Exponential_Averaging)
- [限制性弯曲势能文档](https://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html#restricted-bending-potential)