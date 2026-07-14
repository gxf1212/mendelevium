---
title: "PSP与RadonPy对比——两种聚合物高通量建模工具的异同"
date: "2026-07-14"
last_modified_at: "2026-07-14"
tags: [PSP, RadonPy, polymer-modeling, high-throughput, force-field, GAFF2, OPLS-AA, automated-simulation, polymer-informatics]
description: "对比分析Polymer Structure Predictor (PSP)和RadonPy两种开源聚合物高通量建模工具的功能边界、力场支持、链长控制和实际使用体验"
thumbnail: "/assets/img/4K_1080P_compressed/10201227SUn.jpg"
image: "/assets/img/4K_1080P_compressed/10201227SUn.jpg"
author: Xufan Gao
lang: zh-CN
---

# PSP与RadonPy对比——两种聚合物高通量建模工具的异同

## 本文信息

### PSP（Polymer Structure Predictor）

- **标题**：Polymer Structure Predictor (PSP): A Python Toolkit for Predicting Atomic-Level Structural Models for a Range of Polymer Geometries
- **作者**：Harikrishna Sahu, Kuan-Hsuan Shen, Joseph H. Montoya, Huan Tran, Rampi Ramprasad
- **发表期刊**：Journal of Chemical Theory and Computation
- **发表时间**：2022年3月4日
- **单位**：Georgia Institute of Technology（美国）
- **引用格式**：Sahu, H., Shen, K.-H., Montoya, J. H., Tran, H., & Ramprasad, R. (2022). Polymer Structure Predictor (PSP): A Python Toolkit for Predicting Atomic-Level Structural Models for a Range of Polymer Geometries. *J. Chem. Theory Comput.*, 18(4), 2737–2748. https://doi.org/10.1021/acs.jctc.2c00022
- **代码与数据**：https://github.com/RamprasadGroup/PSP（MIT许可证）

### RadonPy

- **标题**：RadonPy: automated physical property calculation using all-atom classical molecular dynamics simulations for polymer informatics
- **作者**：Yoshihiro Hayashi, Junichiro Shiomi, Junko Morikawa, Ryo Yoshida
- **发表期刊**：npj Computational Materials
- **发表时间**：2022年10月19日
- **单位**：The Institute of Statistical Mathematics（ISM）、The University of Tokyo（日本）
- **引用格式**：Hayashi, Y., Shiomi, J., Morikawa, J., & Yoshida, R. (2022). RadonPy: automated physical property calculation using all-atom classical molecular dynamics simulations for polymer informatics. *npj Computational Materials*, 8, 222. https://doi.org/10.1038/s41524-022-00906-4
- **代码与数据**：https://github.com/RadonPy/RadonPy（BSD-3许可证）

## 摘要

### PSP摘要

> PSP是一个Python工具包，能够从聚合物重复单元的SMILES字符串出发，生成**层级化的聚合物结构模型**，包括寡聚物、无限链、晶体和无定形结构。PSP输出**GAFF2/OPLS-AA力场文件**和多种格式（XYZ/PDB/POSCAR/LAMMPS），可直接用于下游DFT和MD模拟。PSP在1000个不同聚合物SMILES的测试中，二聚体模型生成成功率达到**97.2%（单核）**和**98.5%（多核）**，其中93.5%的模型在15秒内完成构建。PSP的核心优势在于**从化学结构一步到位的建模能力**和**对晶体结构预测的独特支持**。

### RadonPy摘要

> RadonPy是一个开源Python库，能够**全自动计算聚合物的一系列物理性质**，包括密度、热容、热导率、折射率、热膨胀系数等**62种性质**（当前版本）。用户提供重复单元的SMILES、聚合度、链数和温度，RadonPy自动完成构象搜索、DFT电荷计算、链生成、填充分子、平衡MD和性质计算的全流程。RadonPy在**1070种无定形聚合物**上验证了自动化MD流程，密度预测R²=0.890、热导率R²=0.490、折射率R²=0.809。RadonPy是**首个用于聚合物信息学的全流程自动化MD工具**。

### 核心结论

- **PSP和RadonPy定位差异显著**：PSP是**结构模型生成器**（输出模型+力场文件），RadonPy是**全自动性质计算流水线**（输出性质数值）
- **PSP更开放灵活**：用户可调参数多，输出格式丰富（XYZ/PDB/POSCAR/LAMMPS），适合定制化研究
- **RadonPy更黑箱高效**：内置21步压缩/解压平衡协议和自动平衡判定，适合高通量数据库构建
- **两者都不直接输出GROMACS的itp格式**（vs PEMD原生支持），力场选择也有所不同
- **PSP不支持共聚物**（仅均聚物），**RadonPy完整支持均聚物+交替/无规/嵌段共聚物+支化聚合物**

## 背景

全原子MD模拟在聚合物科学中的应用越来越广，但**从化学结构到可运行的MD模型**这一步始终是瓶颈。研究者需要：将化学式转为3D结构、分配力场参数、构建合理的无定形态模拟盒子、平衡体系、跑生产模拟、分析性质——每一步都可能卡住非专家用户。

现有的自动化思路大致有两个方向：

**结构生成方向**：工具聚焦于把化学式变成3D结构和对的力场文件，后续模拟交给用户自己处理。典型代表是PSP、mBuild。它们的优势是灵活，用户能控制每一步的参数；代价是自己需要跑MD和分析。

**全流水线方向**：工具把建模到模拟到分析全部封装成一条自动化链路，用户只需提供化学式，拿到最终性质。典型代表是RadonPy。优势是省事、可重复，适合大数据量筛选；代价是灵活度降低。

之前介绍的PEMD（2026年）介于两者之间：[PEMD——固体聚合物电解质高通量模拟与分析框架](/molecular%20dynamics/modeling%20&%20tools/polymer/2026-07-12-pemd-spe.html)。它是一条端到端流水线，但针对的是**电解质输运和电化学稳定性**这一特定科学问题，力场只用OPLS-AA、MD引擎只用GROMACS。PSP和RadonPy相比之下是**更通用的聚合物工具**。

### 关键科学问题

- **对聚合物研究新手来说，最大的障碍往往不是跑MD本身，而是建不出一个化学结构正确、力场参数合理的初始模型**。PSP和RadonPy都尝试降低这一门槛，但选择了不同的切入点：PSP让用户能看到和调整中间结构，RadonPy则把结构生成封装为流水线的一个环节。
- **高通量筛选需要标准化流程**。RadonPy成功在1000多种聚合物上跑了自动化计算，证明了全自动流水线的可行性；但这种一统配置的方法，对性质预测的精度有一定折扣（如膨胀系数的R²仅为0.178-0.217）。
- **GROMACS用户群体的需求未被充分满足**。PSP和RadonPy都默认输出LAMMPS格式。对于国内GROMACS占主流的MD用户，这意味着需要额外转换步骤。PEMD的出现填补了这一缺口。

---

## 研究内容

### PSP：层级结构模型生成器

PSP的设计哲学是**从SMILES出发，按层级构建越来越复杂的聚合物模型**。代码包由四个模块组成：MoleculeBuilder、ChainBuilder、CrystalBuilder、AmorphousBuilder。

![fig2](psp_figs/fig2.png)

**图：PSP工作流概览**。从聚合物SMILES出发，MoleculeBuilder和ChainBuilder构建寡聚物和聚合物链，CrystalBuilder和AmorphousBuilder分别构建晶体和无定形模型。PSP提供从单体到晶体的层级化建模能力。

#### MoleculeBuilder：有限链（寡聚物）

从SMILES字符串（用`[*]`标注两个连接位点）出发，首先生成目标聚合度的线性SMILES，然后利用RDKit做构象搜索和UFF几何优化。如果RDKit生成的结构缺少链内非键相互作用，PSP可以可选用GAFF2做一轮15ps的NVT短模拟加优化，得到多个构象的3D模型。

输入通过CSV文件传入：重复单元SMILES、左右封端SMILES和聚合度列表。PSP支持线性链和环形链两种模式。

#### 性能测试

PSP在1000个不同聚合物SMILES上测试了二聚体模型生成能力。结果显示，单核处理器成功率为**97.2%**，多核处理器成功率为**98.5%**。多核处理器的计算中，93.5%的模型在**15秒内完成**构建，包括大于150个原子的寡聚物。

![fig3](psp_figs/fig3.png)

**图：PSP生成的分子结构及性能测试**。（a）三种代表性分子的几何结构，碳、氢、氮、氧分别以棕色、浅棕色、青色和红色显示。（b，c）使用单核和多核处理器从1000个SMILES字符串生成二聚体模型的性能测试。PSP展示了高成功率和快速的计算能力。

#### ChainBuilder：无限链

PSP的核心卖点之一——构建**周期性无限链**模型，用于后续的晶体结构预测或第一性原理计算。核心算法是用模拟退火优化一个目标函数：

$$
E = E_{\text{UFF}} + E_{\text{UFF}}(1 - \alpha/180) + E_{\text{connectivity}}
$$

三个项分别对应UFF能量、端基反平行度（$\alpha$接近180°为最优），以及连接正确性罚项。如果单体的模拟退火失败，PSP会自动尝试用二聚体替代单体，增加骨架柔性。

> **关键限制**：对于骨架中有双键的聚合物，ChainBuilder无法生成寡聚物（因为氢原子不能饱和连接位点的价态）。用户需要改写SMILES，让连接位点用单键连接。

#### CrystalBuilder与AmorphousBuilder

晶体结构的预测采用**刚体采样**方式：将一个链拷贝沿z轴放置、旋转和平移，确保原子间距不低于2Å，产生多个候选晶体结构。无定形结构通过PACKMOL将指定数量的寡聚链填充到模拟盒子中，然后由用户选择OPLS-AA（经LigParGen）或GAFF2（经PySIMM/Antechamber）参数化。

OPLS-AA参数化依赖LigParGen Web Server，**单个分子上限200个原子**（LigParGen的限制）。对于分子量较大的聚合物链，需要通过短链生成参数再匹配的思路扩展。

### RadonPy：全自动物理性质计算流水线

RadonPy的设计哲学是**从SMILES一步到性质**。用户提供重复单元的SMILES、聚合度、链数和温度，剩下的全部自动化。

![fig1](radonpy_figs/fig1.jpeg)

**图：RadonPy自动化MD计算工作流**。RadonPy能够自动化每个过程以执行全原子经典分子动力学模拟。多个聚合物在超级计算机的许多计算节点上独立并行计算。RadonPy从SMILES出发，经过构象搜索、DFT电荷计算、链生成、填充分子、平衡MD，最终输出62种物理性质。

#### 自动化工作流总览

整个流水线包括：

1. **构象搜索**：RDKit ETKDG v2生成1000个构象 → MMFF预筛 → 聚类选最优4个 → **DFT优化**（ωB97M-D3BJ/6-31G(d,p)）确定最稳定构象
2. **电子性质计算**：RESP/ESP/Mulliken/Löwdin/Gasteiger五种电荷方法可选，以及HOMO/LUMO/偶极矩/极化率计算
3. **聚合物链生成**：**自避随机游走算法**构建链，控制约1000个原子/链（约10条链/盒子，总计约10000原子），支持**等规/间规/无规立构**控制
4. **填充分子**：10条链随机放置于盒子 → 0.05 g/cm³极低初始密度 → 填充分子模拟逐步压缩
5. **21步压缩/解压平衡协议**：温度300K到600K来回切换，压力最高达50000 atm再回常压——这个协议来自Larsen等的工作，目的是快速获得松弛的无定形态
6. **自动平衡判定**：每5ns检查能量、密度和回旋半径的RSD是否低于阈值，最长等50ns
7. **NEMD计算热导率**：Müller-Plathe反向NEMD方法
8. **性质提取与存储**：62种性质存CSV，轨迹存LAMMPS dump，最终状态存pickle

#### 力场与参数化

RadonPy当前版本支持多种力场（代码包`radonpy/ff/`目录下可验证）：

- **GAFF/GAFF2/GAFF2_mod**（标准有机聚合物）
- **Dreiding/Dreiding_UT**（通用力场）
- **Amber（ff19SB）/GLYCAM_06j**（生物大分子，如肽、多糖、水模型）

对含氟聚合物使用**Träg和Zahn修正参数**（GAFF2_mod），因为标准GAFF2在高密度区表现不佳。缺少的键角参数按GAFF2类似规则**经验估计**（不是留空报错）。原子电荷可以来自**DFT-level RESP拟合**（Psi4），也可以用更快的ESP、Mulliken、Löwdin或Gasteiger方法。这说明RadonPy在精度和自动化之间提供了**多种精度档位的选择**，代价是DFT步骤会显著增加单链的计算时间（论文中约6%聚合物因DFT不收敛而失败）。

#### 共聚物与支化支持

RadonPy是三个工具中对共聚物支持最完整的。代码包（`radonpy/core/poly.py`）原生实现了交替共聚物（`copolymerize_mols`）、无规共聚物（`random_copolymerize_mols`）、嵌段共聚物（`block_copolymerize_mols`）和支化聚合物的自避随机游走构建。这些功能在2022年论文发表时已在后续版本中实现，README已明确列出。

### 功能对比总表

| 维度 | PSP | RadonPy | PEMD（参考） |
| --- | --- | --- | --- |
| 核心定位 | 层级结构模型生成 | 全自动性质计算流水线 | 电解质专用高通量框架 |
| 输入格式 | 重复单元SMILES（`[*]`标注两个连接位点） | 重复单元SMILES（`*`标注连接位点） | 单体SMILES（`*`连接位点）+JSON配置 |
| 结构输出格式 | XYZ/PDB/POSCAR/LAMMPS data | LAMMPS dump文件、pickle | GROMACS拓扑（itp/gro） |
| 支持的力场 | **GAFF2 + OPLS-AA**（二选一） | **GAFF、GAFF2、GAFF2_mod、Dreiding、Dreiding_UT、Amber（ff19SB）、GLYCAM_06j** + 水模型（TIP3P/TIP4P/TIP5P） | OPLS-AA（RESP/RESP2可选） |
| 电荷来源 | CM1A（OPLS-AA）/ Antechamber（GAFF2） | **DFT RESP、ESP、Mulliken、Löwdin、Gasteiger** | CM1A-LBCC / 多构象RESP/RESP2 |
| 线性交替共聚物 | **不支持**（可写为超重复单元，但序列不可控） | **支持**（原生模块，随机游走算法实现） | 支持（均聚/交替/无规/嵌段） |
| 生成任意链长 | **支持**（聚合度参数，寡聚物长度任意） | **支持**（以约1000原子/链为目标，可配置聚合度） | 支持（无长度限制） |
| 带电单体（+1/-1） | **未明确说明**（官方无专门说明） | **支持**（RESP/ESP自动处理质子化带电状态） | 极强（内置Psi4/Multiwfn一键RESP/RESP2） |
| 原生GROMACS itp | **不支持**（LAMMPS格式为主） | **不支持**（LAMMPS格式） | **原生支持** |
| 自带性质计算 | **无**（只生成结构，留给用户计算） | **62种性质**（密度/热导/折射率/热容/Tg/溶解度参数等） | 离子电导率/输运数/ESW/溶剂化结构 |
| MD引擎 | LAMMPS（PySIMM包装） | LAMMPS（直接调用） | GROMACS |
| DFT引擎 | 无（结构优化用UFF/RDKit） | Psi4（DFT构象优化+RESP电荷） | Gaussian 16 |
| 运行模式 | 本地/Colab | **超算/集群**（多节点并行） | 本地/集群 |
| Python版本 | 2.7/3.x（老旧） | 3.9-3.13 | 3.x |
| 许可证 | MIT | BSD-3 | MIT |

### 关键选型差异

#### 什么时候用PSP？

- **需要定制化建模范式**：不想要默认的无定形填料方式，想用自己的方法生成初始结构
- **需要第一性原理输入**：PSP输出VASP POSCAR格式，适合直接接DFT计算
- **需要晶体结构预测**：这是PSP独有的能力（RadonPy和PEMD都不做）
- **能接受LAMMPS工作流**：PSP的力场输出是LAMMPS data格式

#### 什么时候用RadonPy？

- **高通量性质筛选**：RadonPy提供了标准化流程，1000+聚合物的自动运行已验证可行
- **对热物理性质感兴趣**：密度、热导率、折射率、热容、Tg、膨胀系数等62种性质直接输出
- **有超算/集群资源**：RadonPy设计为多节点并行（论文用Fugaku超算）
- **需要共聚物支持**：原生支持交替/无规/嵌段共聚物和支化聚合物

#### 什么时候用PEMD？

- **研究聚合物电解质**：PEMD是唯一原生支持离子电导率、输运数、溶剂化结构和ESW分析的框架
- **需要GROMACS itp输出**：PEMD的力场参数化直接输出GROMACS拓扑文件，不需要格式转换
- **需要强带电单体支持**：PEMD内置Psi4/Multiwfn工作流，支持一键式高精度RESP/RESP2电荷拟合

> 选型总结：如果需要**晶体结构**或**VASP输入**，用PSP。如果只需要**标准无定形态性质**且**有超算资源**，用RadonPy。如果研究**电解质**且用**GROMACS**，用PEMD。如果三者都不完全满足——说明你可能需要组合使用多种工具。

### 实操感受：安装与上手难度

**PSP的安装依赖链最长**。根据README，需要安装PySIMM、LAMMPS、AmberTools（含ANTECHAMBER）、PACKMOL、LigParGen（含BOSS可执行文件），还有RDKit和Open Babel。RDKit和AmberTools可以通过conda安装，但PySIMM和LigParGen的配置比较繁琐，LigParGen还需要单独下载BOSS可执行文件并配置`$BOSSdir`环境变量。论文提供了一个Colab Notebook作为快速体验方式，这对新手比较友好。

**RadonPy的安装相对简单**。根据README，主要依赖Psi4（DFT引擎）和LAMMPS（MD引擎），两者都有conda包。安装命令清晰：先创建conda环境，然后依次`conda install -c conda-forge rdkit psi4 dftd3-python resp mdtraj psutil scipy pandas matplotlib`，再`conda install -c conda-forge lammps`。Python版本要求3.9到3.13，覆盖较新环境。RadonPy已提供PyPI包（`pip install radonpy-pypi`），最小化安装只需`pip install radonpy-pypi`。但RadonPy主要为超算环境设计，单机跑1000+聚合物不太现实。

**PEMD的安装适中**。核心依赖GROMACS和Gaussian 16。PDF提到支持conda安装，依赖链比PSP短。

---

## 关键结论与批判性总结

- **PSP和RadonPy代表了聚合物模拟自动化的两种互补路径**。PSP是自底向上的结构生成器，用户在每个阶段都有控制权；RadonPy是自顶向下的流水线，优先保证标准化和吞吐量。选择取决于研究需求是探索性建模还是批量筛选。
- **两者的共同短板是GROMACS生态支持不足**。在LAMMPS占主流的美国课题组和GROMACS占主流的中国/欧洲课题组之间，存在工具链偏好差异。PSP和RadonPy默认输出LAMMPS格式，而PEMD选择支持GROMACS。对于国内大多数MD用户来说，转换一步是绕不开的。
- **共聚物支持的差距在拉大**。RadonPy当前版本已原生支持交替/无规/嵌段共聚物和支化聚合物（代码包可验证），而PSP仍局限于均聚物，论文2022年版本明确说only applicable to homopolymers。对于需要交替共聚物（如AB交替的离子聚合物）的研究场景，RadonPy目前更合适。
- **带电单体处理需要谨慎验证**。RadonPy的RESP/ESP方法能自动处理质子化带电状态，但论文数据和示例主要针对中性聚合物。PSP没有专门说明永久带电聚合物的支持情况。对于带电聚合物电解质研究，PEMD的强电荷拟合支持（RESP/RESP2+缩放）更可靠。
- **全自动流水线的最后一公里问题**：RadonPy能在1000+聚合物上自动化跑完MD，但论文验证显示膨胀系数等性质的实验-计算相关性很低（R²约0.2）。高通量生成的数据量不等于数据质量，需要用实验数据或迁移学习来校准。