---
title: "机器学习赋能药物发现：四款实用工具与方法全解析"
date: "2025-10-21"
tags: [machine-learning, drug-discovery, molecular-dynamics, data-management, force-field, interpretability, SHAP, dpdata, gmx-ffconv]
description: "介绍四项最新机器学习与计算化学工具：dpdata数据管理、gmx_ffconv力场转换、SHAP模型解释优化、图论机器学习预测抗病毒药物性质"
image: "/assets/img/thumbnail_mine/wh-e7z65r.jpg"
thumbnail: "/assets/img/La-Mancha.jpg"
author: Xufan Gao
lang: zh-CN
---

# 机器学习赋能药物发现：四款实用工具与方法全解析

## 引言

机器学习势能（MLP）和人工智能正在深刻改变药物发现和材料科学领域。从分子动力学模拟到虚拟筛选，从力场开发到模型可解释性分析，研究者们不断开发新工具来**提升计算效率、增强预测准确性、改善模型透明度**。本文将介绍四项近期发表的重要工作，涵盖数据管理、力场转换、模型优化和性质预测等多个关键环节。

---

## 一、dpdata：可扩展的原子机器学习数据集工具包

### 本文信息
- **标题**: dpdata: A Scalable Python Toolkit for Atomistic Machine Learning Data Sets
- **作者**: Jinzhe Zeng, Xingliang Peng等（中国科学技术大学、北京大学等）
- 发表时间: 2025年
- **单位**: 中国科学技术大学、北京大学、AI for Science Institute等
- **引用格式**: Zeng, J., Peng, X., Zhuang, Y.-B., et al. (2025). dpdata: A Scalable Python Toolkit for Atomistic Machine Learning Data Sets. *J. Chem. Inf. Model.* https://doi.org/10.1021/acs.jcim.5c01767
- **源代码**: https://github.com/deepmodeling/dpdata

### 核心问题

机器学习势能的成功高度依赖于**大规模、高质量的参考数据集**。然而，不同软件包采用异构的文件格式和数据模式，严重阻碍了互操作性：

- 电子结构和分子动力学软件使用各自的输入/输出格式
- MLP训练数据通常采用extended XYZ、NumPy数组、pickle、ASE数据库、HDF5等格式
- 即使格式相同，各软件包的数据模式和单位约定也常常不同

### dpdata的解决方案

**核心设计理念**

dpdata是一个开源Python库，采用**灵活的插件式架构**，支持在广泛的文件格式之间进行读取、写入和转换。与ASE等工具不同，dpdata设计为在**系统级别**而非逐个构型操作数据，显著提升了内存效率和推理速度。

**关键功能**

1. **格式支持广泛**：
   - MLP包：DeePMD-kit、QUIP GAP、MACE、NequIP、GPUMD、n2p2
   - MD软件：LAMMPS、AMBER、GROMACS
   - 量子化学：ABACUS、OpenMX、Gaussian、FHI-aims、VASP、Quantum ESPRESSO等
   - 通用格式：XYZ、MOL、SDF、ASE、Pymatgen

2. **数据处理工具**：
   - 自动train-test分割
   - 坐标扰动（用于主动学习）
   - 异常能量去除
   - Δ-learning数据集生成
   - 误差指标计算
   - 单位转换

3. **插件扩展性**：
   - 用户可定义自定义数据类型、格式、驱动和最小化器
   - 示例：dpdata_abinit、cp2kdata、dpdata_ani

**性能优势**

**内存效率对比**：加载QDπ数据集（1,460,161个构型，1.85 GB）
- dpdata: 1.93 GB
- ASE: 7.47 GB（约4倍差距）

**推理加速对比**（dpdata driver vs ASE calculator）
- Water数据集: 4-8倍加速
- Copper数据集: ~6倍加速
- HEA数据集: ~4倍加速

dpdata的系统级设计允许并行处理多个构型，而ASE按顺序逐个处理。

### 实际应用

dpdata已被多项研究用于：
- **格式转换**：将DFT/AIMD输出转换为MLP所需格式
- **数据存储**：以dpdata兼容格式共享数据
- **坐标扰动**：丰富训练集多样性
- **项目集成**：DP-GEN、ChecMatE、PFD-kit、CatFlow、APEX、PyHEA等

---

## 二、gmx_ffconv：GROMACS全原子力场快速转换工具

### 本文信息
- **标题**: gmx_ffconv: A Fast, User-Friendly Semi-Automated All-Atom Force Field Converter for GROMACS
- **作者**: Jasmine E. Aaltonen（Lancaster大学）
- 发表时间: 2025年
- **单位**: Lancaster大学化学系（英国）
- **引用格式**: Aaltonen, J. E. (2025). gmx_ffconv: A Fast, User-Friendly Semi-Automated All-Atom Force Field Converter for GROMACS. *J. Chem. Inf. Model.*, 65, 9850-9855. https://doi.org/10.1021/acs.jcim.5c02200
- **源代码**: https://github.com/Jassu1998/gmx_ffconv

### 核心问题

GROMACS力场转换通常是**耗时且易错**的过程：

- 不同力场采用各自的命名约定和原子排序
- GROMACS要求坐标文件中的原子顺序必须与拓扑文件严格匹配
- 即使像DPPC这样的标准脂质，也无法直接通过pdb2gmx从AMBER Lipid21转换到CHARMM36

现有工具的局限：
- **CHARMM-GUI Force Field Converter**：需要CHARMM输入文件，仅支持AMBER和CHARMM
- **pdb2gmx**：需手动修改残基拓扑文件（.rtp），确保坐标文件语法匹配

### gmx_ffconv的解决方案

**工作原理**

gmx_ffconv通过**分子图匹配**解决原子排序和命名不匹配问题，包含两个核心工具：

1. **ffmap**：通过图同构找到两个力场间的映射
   - 从ITP文件读取原子和键信息
   - 根据原子质量识别化学元素（误差容忍度±0.3 amu）
   - 构建标记图（原子=节点，键=边）
   - 使用NetworkX的VF2算法进行图同构匹配

2. **groconv**：根据映射重新排列坐标文件
   - 读取原始GRO文件
   - 按用户指定的分子类型和数量重组
   - 自动重命名残基和重新编号以匹配新力场
   - 输出重排的GRO文件

**验证系统**

| 系统 | 分子类型 | 分子数 | 总原子数 |
|------|---------|--------|---------|
| 苯乙酸 | BZAA | 1 | 18 |
| 病毒膜 | CHL, DPPC等 | 675,234 | 2,270,122 |
| 人血清白蛋白(HSA) | PROA, PROB | 2 | 18,246 |
| 糖基化SARS-CoV-2刺突蛋白 | PROA-C | 3 | 72,990 |

**性能表现**

时间成本（秒）：

| 分子 | CHARMM → AMBER | AMBER → CHARMM |
|------|---------------|---------------|
| BZAA | 0.10 | 0.10 |
| CHL | 0.10 | 0.10 |
| DPPC | 65.48 | **0.11** |
| DOPE | 60.02 | 0.33 |

注意：某些方向的转换可能快数百倍（如DPPC），这取决于节点排序如何影响VF2算法的搜索过程。

**病毒膜系统转换**：
- ffmap总时间（顺序）: 207.92秒
- ffmap总时间（并行）: 71.31秒
- groconv时间: 4.47秒

### 使用场景

- **力场验证**：使用相同起始坐标比较不同参数化或力场
- **系统转换**：轻松转换文献中的预平衡系统到偏好力场
- **一致性名称**（v1.0.3+）：通过CSV文件确保原子名称在力场间一致

**局限性**
- 不支持水模型转换（3点 ↔ 4点模型）
- 质子化状态必须一致（不支持互变异构体）
- 双硫键等特征仅在两个拓扑都存在时支持

---

## 三、通过SHAP和特征分析改进机器学习分类预测

### 本文信息
- **标题**: Improving Machine Learning Classification Predictions through SHAP and Features Analysis Interpretation
- **作者**: Leonardo Bernal, Giulio Rastelli, Luca Pinzi（Modena and Reggio Emilia大学）
- 发表时间: 2025年
- **单位**: 意大利Modena and Reggio Emilia大学生命科学系
- **引用格式**: Bernal, L., Rastelli, G., Pinzi, L. (2025). Improving Machine Learning Classification Predictions through SHAP and Features Analysis Interpretation. *J. Chem. Inf. Model.* https://doi.org/10.1021/acs.jcim.5c02015

### 核心问题

树基机器学习算法（ET、RF、GBM、XGBoost）在早期药物发现中广泛应用，但常面临：

1. **误分类问题**：假阳性/假阴性影响虚拟筛选效率
2. **可解释性不足**：难以理解预测背后的化学机制
3. **传统置信度过滤的局限**：
   - predict_proba阈值过滤会丢弃大量化合物
   - 无法检测到具有高置信度但实际错误的"局部误分类"

### 创新方法：SHAP与特征值联合分析

**研究设计**

在三个前列腺癌细胞系（PC3、DU-145、LNCaP）的ChEMBL抗增殖数据上开发分类器：
- 算法：ET、RF、GBM、XGBoost
- 特征：RDKit描述符、MACCS keys、ECFP4指纹、custom-fragments

**最佳模型性能**

| 数据集 | 最佳模型 | MCC | F1-score |
|--------|---------|-----|----------|
| DU-145 | ET/GBM-RDKit | 0.60 | 0.83 |
| PC3 | XGB-ECFP4 | 0.64 | 0.86 |
| LNCaP | GBM/XGB-RDKit | 0.62 | 0.88 |

**误分类检测框架**

研究发现：**误分类化合物的特征值（"RAW"）和SHAP值常落在相反类别的范围内**。

基于此，开发了四种标记规则：

1. **"RAW"规则**：化合物的RAW特征值落在相反类别范围内的数量超过阈值
2. **"SHAP"规则**：SHAP值落在相反类别范围内的数量超过阈值
3. **"RAW OR SHAP"**：满足任一条件即标记（高灵敏度）
4. **"RAW AND SHAP"**：同时满足两个条件才标记（高精度）

**阈值定义**：采用分层分位数方法

$$
T_{\text{glob}}(M) = \text{quantile}_p(M_{\text{correct}})
$$

$$
T_C(M) = \text{quantile}_p(M_{\text{correct in C}}), \quad \text{if } |C| \geq 3
$$

其中 $M$ 是"相反类别范围内的特征数量"，$p$ 通常选择80-th或85-th分位数。

**检测效果**

在50%预测置信度下检测到的误分类化合物百分比：

| 数据集 | RAW | SHAP | RAW OR SHAP | RAW AND SHAP |
|--------|-----|------|------------|--------------|
| LNCaP | 48.6% | 46.2% | **63.6%** | 31.2% |
| PC3 | 19.0% | 7.5% | **20.7%** | 5.8% |
| DU-145 | 21.5% | 21.7% | **24.9%** | 18.3% |

**与置信度阈值协同**

随着predict_proba阈值从50%提升到90%，标记规则的效果进一步增强：

- PC3（RAW OR SHAP）：移除误分类从21% → 29%
- DU-145（RAW OR SHAP）：24.9% → 41.9%
- LNCaP（RAW OR SHAP）：63.6% → 70.4%

### 实际意义

- **虚拟筛选优化**：在大型化合物库筛选中，最大化灵敏度以识别边界化合物
- **二次筛选精炼**：在聚焦筛选中，使用高精度规则保留真阳性
- **特征可解释性**：误分类化合物显示的关键描述符（如"EState_VSA1"、"SMR_VSA6"）为结构优化提供洞察

---

## 四、图论+机器学习：用拓扑指数预测抗病毒药物性质

### 本文信息
- **标题**: A Graph-Based Machine Learning Framework for Predicting Physicochemical Properties of Antiviral Drugs via Topological Indices
- **作者**: Irfan Haider, Muhammad Ahsan等（巴基斯坦COMSATS大学等）
- 发表时间: 2025年
- **单位**: COMSATS大学（巴基斯坦）、印度中央大学、中东技术大学（塞浦路斯）等
- **引用格式**: Haider, I., Ahsan, M., Siddiqui, M. K., et al. (2025). A Graph-Based Machine Learning Framework for Predicting Physicochemical Properties of Antiviral Drugs via Topological Indices. *J. Chem. Inf. Model.* https://doi.org/10.1021/acs.jcim.5c00117
- **源代码**: https://github.com/IrfanHaider/graph_based_antiviral_drugs.git

### 创新框架：两阶段机器学习

传统QSPR方法直接从分子结构预测性质，本研究引入**拓扑指数作为中间桥梁**：

**阶段一：SMILES → 拓扑指数**
- 输入：SMILES字符串
- 输出：六种拓扑指数（M1、M2、ABC、Randić、Harmonic、Forgotten）
- 方法：RDKit解析分子图，ML模型预测指数

**阶段二：拓扑指数 → 理化性质**
- 输入：预测的拓扑指数
- 输出：六种性质（摩尔折射率、极性表面积、极化率、摩尔体积、分子量、复杂度）
- 模型：四种ML算法比较

**拓扑指数定义**

1. **First Zagreb (M1)**：

$$
M_1(G) = \sum_{v \in V(G)} d_v^2
$$

反映分子的**整体连接性和分支度**。

2. **Second Zagreb (M2)**：

$$
M_2(G) = \sum_{uv \in E(G)} d_u d_v
$$

捕捉**相邻原子的连接特征**。

3. **ABC指数**：

$$
\mathrm{ABC}(G) = \sum_{uv \in E(G)} \sqrt{\frac{d_u + d_v - 2}{d_u d_v}}
$$

与**分子稳定性和应变能**相关。

4. **Randić指数**：

$$
R(G) = \sum_{uv \in E(G)} \frac{1}{\sqrt{d_u d_v}}
$$

反映分子的**分支程度**。

5. **Harmonic指数**：

$$
H(G) = \sum_{uv \in E(G)} \frac{2}{d_u + d_v}
$$

与分子的**电子性质**相关。

6. **Forgotten指数**：

$$
F(G) = \sum_{v \in V(G)} d_v^3
$$

对高度顶点赋予更大权重，适用于**复杂结构分子**。

### 预测性能

**阶段二：理化性质预测**

| 性质 | 最佳模型 | $R^2$ |
|------|---------|-------|
| 分子量（MW） | XGBoost | **0.9950** |
| 极化率（P） | 神经网络 | **0.9891** |
| 摩尔折射率（MR） | 线性回归 | 0.9863 |
| 摩尔体积（MV） | 随机森林 | 0.9732 |

**关键发现**
- **M1和Forgotten与MW、P、MR的相关系数超过0.95**
- XGBoost和随机森林显著优于线性回归
- 极性表面积（PSA）预测较难（$R^2$=0.4242）

### 优势与局限

**优势**
- **降低复杂度**：每阶段输入输出维度低
- **提高可解释性**：拓扑指数有明确化学意义
- **模块化设计**：两阶段可独立优化
- **计算效率**：相比量子化学计算极低成本

**局限性**
- **数据集规模小**：59个样本限制泛化能力
- **缺乏3D信息**：忽略立体化学和构象效应
- **PSA预测不佳**：度基指数对极性特征表征能力有限

---

## 总结与展望

本文介绍的四项工作展示了机器学习和计算化学工具链的不同环节：

### 工具定位

| 工具 | 功能 | 适用场景 |
|------|------|---------|
| **dpdata** | 数据管理与转换 | MLP开发、大规模数据处理 |
| **gmx_ffconv** | 力场快速转换 | 比较模拟、系统迁移 |
| **SHAP+特征分析** | 模型优化与误分类检测 | 虚拟筛选、模型可解释性 |
| **图论ML框架** | 性质预测 | 抗病毒药物设计、QSPR建模 |

### 共同趋势

1. **效率优先**：dpdata实现4倍内存节省，gmx_ffconv秒级转换复杂系统
2. **可解释性**：SHAP分析不仅解释模型，还能主动改进预测
3. **插件化设计**：dpdata和gmx_ffconv均支持用户扩展
4. **实用导向**：所有工具均开源，提供详细文档和示例

### 未来方向

- **工具整合**：将dpdata用于MLP数据管理，gmx_ffconv用于多力场验证，SHAP用于模型诊断
- **深度学习融合**：图神经网络替代ECFP4以减少比特碰撞，提升拓扑指数预测
- **主动学习**：结合SHAP标记和dpdata坐标扰动，优化训练集采样
- **跨尺度建模**：从拓扑指数到全原子MD，再到粗粒化模拟的无缝衔接

**参考资源**
- dpdata文档：https://docs.deepmodeling.com/projects/dpdata
- gmx_ffconv教程：https://github.com/Jassu1998/gmx_ffconv
- SHAP官方文档：https://shap.readthedocs.io

这些工具的出现标志着计算化学和药物发现正在向**自动化、智能化、可解释化**方向发展，为研究者提供了更高效的武器库。

