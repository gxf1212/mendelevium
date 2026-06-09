---
title: "PUCHIK：非球形纳米粒子界面分析的Python工具包"
date: "2026-06-08"
tags: [puchik, nanoparticle, aspherical, interface-analysis, molecular-dynamics, alpha-shape, convex-hull, software-tools]
description: "PUCHIK为非球形纳米粒子MD模拟提供界面表征、密度计算和体积分析，支持alpha shape和convex hull两种方法，适用于胶囊状、棒状等复杂形貌"
thumbnail: "/assets/img/thumbnail_mine/empty.jpg"
image: "/assets/img/thumbnail_mine/empty.jpg"
author: Xufan Gao
lang: zh-CN
---

# PUCHIK工具包：非球形纳米粒子界面、密度与体积的自动化分析

## 本文信息

- **标题**：PUCHIK：用于分析非球形纳米粒子分子动力学模拟的Python工具包
- **作者**：Hrachya Ishkhanyan，Alejandro Santana-Bonilla，Christian D. Lorenz
- **发表期刊**：*Journal of Chemical Information and Modeling*
- **发表时间**：2025年（Volume 65，Pages 1694-1701）
- **DOI**：https://doi.org/10.1021/acs.jcim.4c02128
- **单位**：英国伦敦国王学院（King's College London）物理系与工程系；亚美尼亚国家科学院信息学与自动化学研究所
- **引用格式**：Ishkhanyan, H.; Santana-Bonilla, A.; Lorenz, C. D. (2025). PUCHIK: A Python Package To Analyze Molecular Dynamics Simulations of Aspherical Nanoparticles. *J. Chem. Inf. Model.*, 65, 1694-1701. https://doi.org/10.1021/acs.jcim.4c02128
- **代码与数据**：PUCHIK代码以开源形式提供，遵循标准的Python包结构与文档规范。

## 摘要

> 准确描述纳米粒子的界面对于理解其内部结构、界面性质乃至最终功能至关重要。虽然当前计算方法对球形和准球形纳米粒子提供了合理的描述，但针对胶囊状和棒状体系等**非球形结构**的有效模型仍然存在需求。本工作引入了**Python Utility for Characterizing Heterogeneous Interfaces and Kinetics（PUCHIK）**，这是一种为描述**球形和非球形纳米粒子**而开发的新算法。通过准确描述纳米粒子界面的位置，该算法允许计算各种重要物理量（例如不同原子/分子类型相对于界面的密度、纳米粒子体积等）。PUCHIK基于**SciPy、MDAnalysis和Cython**构建，提供了经过优化的Python实现，执行时间与粒子数呈线性关系。PUCHIK能够可靠地表征纳米粒子界面，为纳米科学和纳米技术中的**in silico材料设计**提供了强大工具。

![abs](puchik_figs/abs.png)

**摘要图：PUCHIK的核心工作流程**——从MD结构到原子点集、再到Convex hull和Alpha shape两种界面建模方法的完整流程。Convex hull形成凸形包络，Alpha shape则生成贴合粒子实际形貌的凹形界面。

### 核心结论

- PUCHIK首次实现了对**非球形纳米粒子**（胶囊状、棒状等）界面的准确表征，弥补了传统工具仅适用于球形粒子的局限
- 采用**alpha shape**和**convex hull**两种方法定义界面，通过Cython优化后实现与粒子数呈**线性关系**的计算复杂度
- 在TX100胶束和吲哚美辛共溶剂体系的对比测试中，PUCHIK成功避免了nanoCISC算法的水密度虚高问题
- 开源、易用（3行代码完成密度计算），定位是纳米粒子界面分析的**标准化工具**

## 背景

纳米粒子的**界面表征**是理解其结构-性质关系的核心。传统的密度分析方法（如以质心为基准的径向密度分布）对球形粒子效果良好，但对**非球形粒子**（如胶囊状、棒状、不对称胶束）会产生严重误判。现有工具如**nanoCISC**虽能处理部分复杂形貌，但在计算密度时可能出现**水密度虚高**、**组分密度分布不合理**等问题。PUCHIK通过**计算几何方法**（alpha shape和convex hull）精确定义纳米粒子的核心-壳界面，进而计算相对于界面的密度分布和体积。

### 配套资源

- **算法依赖**：SciPy（ConvexHull）、MDAnalysis（轨迹/拓扑管理）、Cython（性能优化）
- **计算复杂度**：$O(mN)$，其中m为凸包顶点数，N为粒子数，实测执行时间与N呈线性关系
- **优化策略**：支持Python单进程（SP）、多进程（MP）以及Cython加速，MP模式可将单帧计算时间从0.40秒降至0.13秒
- **适用体系**：固体、空心、介孔材料，以及表面活性剂胶束、药物纳米载体等软物质体系

对于涉及**非球形纳米粒子、表面活性剂自组装、药物纳米载体**等体系的MD研究者，PUCHIK提供了一个**专业但易用**的界面分析工具，建议作为标准工作流的一部分。

## 创新点

- **alpha shape界面定义**：首次将alpha shape方法引入纳米粒子界面分析，能够精确表征凹形界面，避免convex hull的过度包裹问题
- **线性时间复杂度**：通过Cython优化和多进程并行，实现与粒子数呈线性关系的执行时间，显著优于传统方法
- **非球形体系适用性**：专门针对胶囊状、棒状等非球形纳米粒子设计，突破了球形假设的局限
- **开源易用**：提供简洁的API（3行代码完成密度计算），配套完整文档和示例

### 工具能力速览

| 功能类 | 代表方法 | 核心功能 | 适用场景 |
| --- | --- | --- | --- |
| 界面定义 | Convex Hull、Alpha Shape | 计算几何方法精确界定界面 | 球形/非球形纳米粒子 |
| 密度计算 | `Interface.density()` | 相对界面的径向密度分布 | 核-壳结构、胶束、纳米载体 |
| 体积计算 | `Interface.volume()` | 界面包围的体积 | 多孔材料、空心粒子 |
| 并行优化 | Python MP、Cython | 多进程加速和编译优化 | 大体系、长轨迹 |

> **对化学无关设计的强调**：虽然论文主要展示表面活性剂体系，PUCHIK的设计目标是适应**任何纳米粒子体系**——这种"通用化分析"思路是软物质模拟分析工具的发展方向之一。

---

## 研究内容

### 一、方法学设计

PUCHIK的命名来自亚美尼亚语的"气球"，寓意其能适应各种形状的纳米粒子。整个包主要建立在以下组件之上：

- **SciPy**：通过`ConvexHull`类构建凸包界面（Qhull库的Python封装）
- **CGAL（Computational Geometry Algorithms Library）**：alpha shapes在C++层面通过CGAL实现
- **MDAnalysis**：负责读取轨迹和拓扑，提供类似CHARMM和VMD的选择语法
- **Cython**：优化计算密集型部分，实现显著性能提升

PUCHIK的包结构分为两个子包：`core`（主要类`Interface`）和`utilities`（辅助工具）。用户只需3行代码即可完成密度计算：

```python
from puchik.core import Interface
interface = Interface(universe, selection='name CA')
density = interface.density()
```

`Interface`类的设计遵循**关注点分离原则**：拓扑和轨迹读取由MDAnalysis负责，凸包构建通过SciPy的`ConvexHull`实现，alpha shape通过CGAL（C++库）实现，密集型计算通过Cython优化。这种**模块化设计**使得每个组件都可以独立替换或升级，而不影响整体功能。

整套工具采用**化学无关设计**。虽然论文主要展示表面活性剂体系，算法可应用于**任何纳米粒子体系**，包括**金属纳米粒子、聚合物胶束、药物纳米载体**等。

### 二、界面定义：Convex Hull vs Alpha Shape

PUCHIK提供两种界面定义方法：**convex hull**（凸包）和**alpha shape**（α形状）。Convex hull是包含所有点的最小凸集，而alpha shape通过调节α参数可以生成凹形界面，更精确地贴合非球形粒子的实际形貌。

![fig1](puchik_figs/fig1.png)

**图1：标准几何体测试**——用圆柱和球形验证PUCHIK的密度计算准确性。

- **图1a**：测试结构——左为圆柱（半径和半高均为2.9 nm），右为球形（半径2.9 nm）
- **图1b**：标准方法（左，以质心为基准）与PUCHIK算法（右，以convex hull界面为基准）的密度对比

PUCHIK计算的密度与理论值（0.0375 Å⁻³）吻合良好。

### 三、非球形胶束案例分析：TX100体系

论文以**Triton X-100（TX100）表面活性剂胶束**为例，对比PUCHIK与现有工具nanoCISC在非球形体系中的表现。该胶束因包裹吲哚美辛药物分子而呈现**不对称形状**（由6750个重原子组成，尺寸约110 Å × 84 Å × 74 Å）。

![fig2](puchik_figs/fig2.png)

**图2：TX100胶束的密度计算对比**——展示PUCHIK在真实非球形体系中的优势。

- **图2a**：拉长的TX100胶束的快照
- **图2b**：nanoCISC算法计算的水密度（蓝色）、疏水尾（粉色）和亲水头（青色）密度分布——水密度虚高（>0.033 Å⁻³），且头基密度在核内高于尾基，不符合稳定核-壳模型
- **图2c**：PUCHIK算法计算的密度分布——水密度在核内接近0，PEO密度在r=0处达峰值后降为0，符合物理预期

nanoCISC的主要问题在于两点：**水密度虚高**（计算得到的水密度远高于体相水密度~0.033 Å⁻³）和**结构不合理**（头基密度在核内高于尾基密度，不符合典型核-壳胶束的分布）。相比之下，PUCHIK通过准确界定界面，得到了**物理上合理**的密度分布。

### 四、Alpha Shape的优势：处理凹形界面

对于具有凹陷或复杂形貌的纳米粒子，convex hull会过度包裹，导致密度计算出现偏差。Alpha shape方法通过调节α参数，能够生成更贴合实际形貌的凹形界面。

![fig3](puchik_figs/fig3.png)

**图3：Convex Hull vs Alpha Shape对比**——同一表面活性剂纳米粒子的两种界面建模方法。

- **图3a**：Convex hull建模——红色区域虽属于凸包，但几乎不含粒子原子，被水分子填充
- **图3b**：Alpha shape建模——形成凹形界面，更贴合纳米粒子的整体形状
- **图3c**：使用convex hull计算的密度——水密度在核内显著偏高
- **图3d**：使用alpha shape计算的密度——水密度明显降低，更符合物理现实

Alpha shape不仅提高了界面定义的准确性，还通常**包裹更小的体积**，从而避免了convex hull的过度包裹问题。这意味着基于alpha shape计算得到的密度分布更贴近真实物理情况，**对于研究界面附近水分子分布、活性位点可及性等问题尤为关键**。

### 五、计算性能：线性时间复杂度

PUCHIK通过Cython优化和多进程并行，实现了与粒子数呈线性关系的执行时间。论文测试了168,989个原子的体系（其中含约51,000个水分子和约1,100个界面原子），结果显示：

![fig4](puchik_figs/fig4.png)

**图4：执行时间与粒子数的线性关系**——展示PUCHIK的可扩展性。

**表1：不同优化技术的单帧执行时间对比**

| 优化技术 | 执行时间（秒/帧） | 加速比（基于单进程Python） |
| --- | --- | --- |
| Python SP（单进程） | 0.40 | 1.0× |
| Python + Cython SP | 0.37 | 1.1× |
| Python MP（多进程） | 0.13 | 3.1× |
| Python + Cython MP | 0.12 | 3.3× |

注：加速比基于原文表1的执行时间计算（0.40/0.40=1.0、0.40/0.37≈1.1、0.40/0.13≈3.1、0.40/0.12≈3.3）。原文表述为"**The rate of change in computational time is ∼200% less for parallelized code than for the single process**"，即并行代码比单进程快约**3倍**。

多进程模式带来**约3倍加速**，使PUCHIK能够高效处理大规模体系。线性时间复杂度保证了算法在**大体系、长轨迹**分析中的可扩展性。对于**含数十万原子**的药物-蛋白复合物或胶束-膜相互作用体系，PUCHIK能在合理时间内完成全轨迹的界面分析。

---

## 关键结论

- PUCHIK为非球形纳米粒子的界面表征提供了**准确且高效**的解决方案。通过alpha shape和convex hull两种方法，PUCHIK能够精确界定界面，进而计算相对界面的密度分布和体积。在TX100胶束测试中，PUCHIK避免了nanoCISC的水密度虚高问题；在alpha shape对比中，成功降低了convex hull带来的过度包裹误差。
- PUCHIK的核心优势在于**线性时间复杂度**和**物理上合理的结果**。多进程模式带来**约3倍加速**，使其能够高效处理大规模体系，**大体系、长轨迹**分析的可扩展性得以保证。
- 论文明确指出PUCHIK的设计目标是"**in silico材料设计的强大工具**"，为纳米科学和纳米技术提供支持。

### 局限性

- Alpha shape的α参数需要用户根据体系特点调整，缺乏自动选择策略
- 论文主要在表面活性剂体系验证，对**金属纳米粒子、无机材料**等硬物质的迁移性有待进一步测试
- PUCHIK目前不支持命令行执行，必须在Python解释器中运行，对不熟悉Python的用户有一定门槛
- 密度计算的性能优化主要针对**单帧分析**，对长轨迹的批处理效率未详细讨论
- 与PySoftK等软物质分析工具的**集成性**未在论文中展示，用户需要自行编写接口代码

对于涉及**非球形纳米粒子、表面活性剂自组装、药物纳米载体**等体系的MD研究者，PUCHIK提供了一个**不可或缺**的界面分析层，建议作为标准工作流的重要组成部分。结合PySoftK（同一团队开发，聚焦软物质聚集体分析），二者可覆盖从**简单球形胶束**到**复杂胶囊状纳米载体**的界面分析需求。这一系列工作反映了软物质/纳米粒子MD分析领域正在**从零散脚本向标准化工具集**演进的发展趋势。
