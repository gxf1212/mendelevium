---
title: "SARS-CoV-2复制酶多蛋白1ab (CHEMBL4523582) 数据集：机器学习药物设计的真实用例"
date: "2026-01-14"
tags: [sars-cov-2, replicase, pp1ab, rdRp, chembl, dataset, qsar, machine-learning, drug-design]
description: "全面整理CHEMBL4523582（SARS-CoV-2复制酶多蛋白1ab）活性数据集，聚焦RdRp定量数据来源、处理策略与机器学习应用案例，并单列3CLPro研究。"
image: "img/thumbnail_mine/wh-g83qpe.jpg"
thumbnail: "img/thumbnail_mine/wh-g83qpe.jpg"
author: Xufan Gao
lang: zh-CN
---

# 从101条IC50到千万级虚拟筛选：SARS-CoV-2 RdRp数据集与机器学习实战指南

## 概要

本文系统整理SARS-CoV-2 RNA依赖的RNA聚合酶（RdRp/nsp12）的IC50/EC50数据集，并提供数据稀缺场景下的机器学习药物设计实施指南。

- **数据资源与规模**：CHEMBL4523582（SARS-CoV-2 pp1ab）包含约1,455条IC50/EC50记录，但真正针对nsp12/RdRp的连续活性数据仅**101条**（72条IC50 + 9条EC50 + 19条Inhibition），其余多为3CLPro或其他非结构蛋白测定，需从聚合靶点中通过关键词筛选分离出RdRp专属记录。
- CAS汇编的**1,209条**RdRp抑制剂（Active 464/Inactive 745，阈值IC50≤10 µM）是目前最完整的二分类数据集。SARS-CoV-1（CHEMBL5118）提供**27条**nsp12记录及完整数据集约**2,410条**记录，可用于迁移学习。其他数据源包括BindingDB（36条）、PubChem（78条）、**NCATS OpenData RdRp Enzymatic Activity**（2,678个化合物的4点剂量反应；375条有AC50拟合值）、VDDB（约300万条多靶点）、TDC和Guide to PHARMACOLOGY。
- **模型方法**：总结**5种数据分割策略**（随机、互异、骨架、时间、分层）、**连续值处理方法**（pIC50转换、阈值二分类IC50≤10 µM）、**重复测定处理**（聚合、去重、保守值）及**应对数据稀缺**的方案（迁移学习、多任务学习、半监督学习、数据增强）。小样本（<200条）优先二分类（RF、SVM、CORAL），中等规模（200-500条）可尝试pIC50回归，大规模（>500条）推荐深度学习（DNN、GNN）。
- **应用案例**：涵盖**10个ML/药设案例**，包括Ivanov等的AutoML QSAR（1,209条，RBF-SVM最优）、Bazzi-Allahri等的CORAL SMILES模型（10组互异划分，验证集MCC最高0.89）、Tsuji的量子化学联动虚拟筛选（1,838,257个分子）、Sarma的药效团+MD筛选（100 ns验证）和生成式分子设计（75个SMILES候选）。QSAR+药效团+对接三段式虚拟筛选可在6,020万分子中缩减至数十个候选。
- **结构与文献**：提供Cryo-EM结构（EMD-11692/7AAP含favipiravir-RTP、EMD-30794/7DOI含penciclovir）、对接模板（PDB 6Y2F/6NUR）及**28篇关键文献**的DOI和详细说明，涵盖数据集构建、QSAR建模、虚拟筛选、3CLPro研究和核苷类似物机制。

## 数据现状与测定来源

### RdRp数据源对照表

#### 能用的（含IC50/EC50/AC50连续值数据集）

| 数据源 | 规模/时间 | 连续label数据量 | 采集说明 | 活性类型 |
| --- | --- | --- | --- | --- |
| ChEMBL (CHEMBL4523582) | SARS-CoV-2 pp1ab，持续增长 | **101条**（站内脚本筛选所得） | ChEMBL数据库SARS-CoV-2 pp1ab靶点（17,929行原始数据），需结合assay描述和关键词从聚合靶点中分离出RdRp专属记录。 | 酶学IC50、细胞EC50 |
| ChEMBL (CHEMBL5118) | SARS-CoV-1 pp1ab，约2,410条总记录 | **27条**（站内脚本筛选所得） | ChEMBL数据库SARS-CoV-1 pp1ab靶点（完整数据集约532条），可同时筛选两个靶点并标注Virus_Type用于迁移学习。 | IC50、Ki |
| PubChem BioAssay | 零散更新；2026-01-09抓取1,177个候选AID | **78条**（5个AID，60条IC50 + 18条EC50，均为nM） | PubChem生物测定数据库，通过protein name接口抓取RdRp相关assay，经关键词过滤与表格解析抽取连续值。 | IC50、EC50 |
| NCATS OpenData Portal | SARS-CoV-2 Assays（RdRp Enzymatic Activity） | **2,678条**（4点剂量反应）；其中**375条**有AC50拟合值（µM） | NCATS高通量筛选数据，生化酶学RdRp测定（Biochemical + Fluorescence），4点剂量反应曲线拟合AC50值。 | AC50（生化） |
| BindingDB COVID-19集合 | 2020年打包发布 | **36条** | COVID-19 target set包含36条RdRp IC50（0.001–10 µM），含部分与ChEMBL重叠的测定，可做交叉验证。 | IC50 |

**预计去重后**：约**1,800-2,000个独立化合物**（含SARS-CoV-1和重复），其中约**500-600条有连续IC50/EC50/AC50值**

#### 不适用的（二分类/对接/多靶点难拆）

| 数据源 | 规模/时间 | 连续label数据量 | 采集说明 | 活性类型 | 不适用原因 |
| --- | --- | --- | --- | --- | --- |
| CAS汇编（Ivanov等，2020） | 1,209个化合物 | **1,209条** | Active 464 / Inactive 745（阈值IC50≤10 µM），覆盖SARS-CoV-2、SARS-CoV-1、HCV等RdRp抑制剂，SMILES与标签在ACS Omega论文附录开放。 | 二分类IC50标签 | **二分类**（Active/Inactive），适合分类基线，不是连续值补充 |
| DockCoV2 | 2021年发布（NAR） | **7个靶点**（含RdRp）的对接数据 | 针对SARS-CoV-2的药物数据库，提供RdRp等7个蛋白的分子对接分析结果，包含FDA批准药物、临床试验化合物和天然产物。 | 对接分数、结合能 | **对接结果**，不是实验连续活性 |
| VDDB（Tao 2022/2023） | 39种病毒、>71万小分子、约300万条活性 | **约300万条**（117靶点合计，未拆分到单靶点） | 自带117个靶点级机器学习模型（含RdRp），支持下载与在线预测，可用于外部验证或迁移学习。 | 靶点预测、QSAR输出 | **多靶点混合**，未拆分到单靶点 |
| Therapeutic Data Commons (TDC) | 核心“TDC-1”提供66个标准化数据集（22个任务） | **多模态**（小分子、蛋白、单细胞等；可筛选COVID/RdRp相关任务） | 2024-06-22发布“TDC-2”，增加>10种模态、>1,000个多模态数据集与Model Hub，Python API交付ML-ready拆分与评测脚本。 | 多种活性类型（QSAR、生成、单细胞） | **当前常用COVID任务不对齐RdRp连续活性**：如Touret为细胞感染二分类（不区分靶点），Diamond为3CLPro二分类 |
| Stanford CoV-RDB | 2020年发布，持续更新 | - | Coronavirus Antiviral Research Database：涵盖SARS-CoV-2、SARS-CoV、MERS-CoV的体外/动物/临床数据。 | IC50、EC50 | **抗性数据库**：聚焦于单抗、恢复期血浆和疫苗接种后的敏感性，以及mAbs、3CLpro抑制剂、RdRP抑制剂的抗性突变，不提供小分子抑制剂IC50/EC50数据集 |
| dbSCI | 2022年发布 | **5条** Nsp12数据（164条细胞实验中筛选） | 手工注释的SARS-CoV-2抑制剂数据库，下载页提供`dbSCI_Cell_experiments.xlsx`（164条细胞实验）、`dbSCI_Virtual.xlsx`（64条临床试验）。 | IC50、EC50 | **数据量太少且质量差**：Nsp12数据仅5条（Ribavirin、Ivermectin、Enisamium Iodide、Remdesivir triphosphate、VR17-04），IC50值从2.5μM到40.7mM单位混杂，1条缺失SMILES；已提供`extract_dbsci_rdrp.py`脚本提取数据，但不适用于机器学习建模 |

不适用的数据库：

- https://tdcommons.ai/single_pred_tasks/hts/#sars-cov-2-in-vitro-touret-et-al
- Stanford CoV-RDB: https://covdb.stanford.edu/（抗性数据库）
- dbSCI: https://bioinfo.henu.edu.cn/COVID/COVIDIH.html（下载页：https://bioinfo.henu.edu.cn/COVID/download.html，仅5条Nsp12）
- AntiviralDB: https://www.antiviraldb.com/（无法筛选SARS-CoV-2 RdRp数据）
- COVID-19 MolSSI Hub: https://covid.molssi.org/（治疗剂数据很少）

不适用的文献：

- **ACS Infect Dis 2025**（Dual-Site Inhibition）：仅4个IC50/EC50有完整数据，不足以作为数据集
- **核苷类似物综述（2023）**：综述为主，非结构化数据集
- **[SARS-CoV-2 cytopathicity dataset](https://doi.org/10.1038/s41597-021-00848-4)**（Scientific Data, 2021，PMC7910569）
  - **规模**：筛选5632个化合物，发现258个hits（>75%抑制），对67个最活跃分子进行8点剂量反应测试
  - **问题**：细胞水平表型筛选（Caco-2细胞病变效应），**靶点不明确**，可能是任何抗病毒机制（膜融合抑制、蛋白酶抑制、RdRp抑制等），不能确定为RdRp专属数据
- **ChemRxiv 2023**（de novo molecule generation）：只有预测活性，非实验IC50
- **Digging for Discovery 2021**：虚拟筛选，无实验IC50
- **DrugDevCovid19 Atlas**：非RdRp专属，含多靶点混合数据
- **Kaggle COVID-19 Drug Discovery Data**：非RdRp专属通用数据集
- **COVID-19 Molecular Structure and Therapeutics Hub**（https://covid.molssi.org/）：社区驱动数据仓库，**治疗剂数据很少**，不足以作为RdRp数据集
- **AntiviralDB**（https://www.antiviraldb.com/）：2025年发布的专家库集抗病毒药物数据库，声称包含冠状病毒抗病毒化合物数据，但**无法通过Virus: SARS-CoV-2 + Target: RdRp筛选出任何数据**，靶点标注不完整或缺失
- **PNAS 2024**（Structural basis of polymerase inhibition，PMC11912441）：只有1个化合物HeE1-2Tyr（IC50=5μM）

#### 存疑的（需手动确认是否可导出RdRp连续值）

请手动访问以下链接，确认是否可按靶点（Polymerase/RdRp）导出IC50/EC50连续值数据集：

| 数据源 | 规模/时间 | 连续label数据量 | 采集说明 | 活性类型 |
| --- | --- | --- | --- | --- |
| COVID19 Drug Repository | 2021年发布（NAR） | **约460项**（184个批准药 + 384个研究药物） | 通过文献挖掘构建的COVID-19治疗剂数据库，包含小分子抑制剂数据，可补充RdRp抑制剂的化学空间。 | 多种活性类型 |

**注意**：上表"存疑"表示需要你手动确认是否能导出RdRp连续活性表（IC50/EC50/AC50等）；确认后我再把它们纳入可用数据源并补抽取脚本/统计。

#### 文献SI（需手动验证）

1. **[High-throughput screening identifies non-nucleoside inhibitors](https://www.slas-discovery.org/article/S2472-5552(25)00082-6/fulltext)**（SLAS Discovery, 2025）
   - 发现非核苷类RdRp抑制剂，使用全长nsp12（1–932位）
   - **需确认**：SI是否提供完整IC50表格


2. **[Current understanding of nucleoside analogs](https://doi.org/10.1016/j.csbj.2023.09.001)**（Comput Struct Biotechnol J, 2023，PMC10498173）
   - 核苷类似物抑制RdRp的综合综述
   - **需确认**：是否包含结构化IC50数据表

### 浅总结

#### 为什么找不到“5,000条RdRp活性数据”？

1. **聚合式靶点设计**：CHEMBL将pp1ab映射到单一目标（CHEMBL4523582），任何nsp5、nsp12、nsp13、nsp15的测定结果（无论是酶学还是细胞水平）都会归到同一target，导致靶点标签混乱。**关键问题不是测定类型（酶学vs细胞），而是靶点归属不清**。
2. **Guide-to-Pharmacology互链误导**：在Ligand Activity Chart中，JAK3 inhibitor IV、GC-376、SU-3327、PBIT等被标注为“RNA-dependent RNA polymerase”抑制剂，但原始测定是3CLPro荧光FRET，靶点完全错误。
3. **数量级差异的客观存在**：依据Ashraf 等（2023, PLOS ONE），CHEMBL在2020–2021年仅能提供161条pp1ab IC50，其中<1 µM活性仅5条，失活>10 µM有156条；到2024年增至约1,455条，但真正属于RdRp的仍属少数。
4. **COVID Moonshot项目聚焦3CLPro**：COVID Moonshot是开放科学药物发现项目，但其主要靶点是SARS-CoV-2主蛋白酶（Mpro/3CLpro），而非RdRp，因此该项目的大规模数据集对RdRp研究的直接贡献有限。

#### 数据增长时间线

- **2020年初**：CHEMBL几乎无SARS-CoV-2 RdRp记录，研究者多使用SARS-CoV（CHEMBL5118）或HCV RdRp数据做迁移。
- **2020年10月**：Ivanov 等整合了1,209条RdRp抑制剂，首次提供系统的Active/Inactive标签。
- **2021年**：Eur J. Med. Chem. 报道的HEK293T转染实验提供了十余条nsp12 (±nsp14/nsp10) 细胞系IC50。
- **2022–2023年**：BindingDB、PubChem以及多项分子设计研究陆续上传新的酶学数据，ChEMBL条目累计到千级，但3CLPro/Hydrolase等测定增长更快。
- **2024年**：Bazzi-Allahri 等使用2,377条分子（3CLPro+RdRp）在CORAL中构建SMILES-QSAR模型，并开展6,020万规模的虚拟筛选。

#### 结构与表征资源

- **Cryo-EM与PDB**：2020年以来已发布多份含nsp12-nsp8-nsp7-RNA复合物的高分辨率结构，为分子动力学、量子化学与虚拟筛选提供了可靠模板：
  - **EMD-11692 / 7AAP**：掺入favipiravir-RTP的RdRp复合物结构
  - **EMD-30794 / 7DOI**：掺入penciclovir的RdRp复合物结构
  - **PDB 6Y2F / 6NUR**：常用于分子对接的RdRp结构模板
- **官方知识库**：Guide to PHARMACOLOGY将RdRp列为“SARS-CoV-2 polymerase complexes”并提供CSV下载链接，适合作为ChEMBL以外的交叉验证。

#### 连续IC50/EC50数据的特殊挑战

与3CLPro相比，SARS-CoV-2 RdRp的连续活性数据面临以下独特挑战：

1. **测定复杂度高**：RdRp是多亚基复合物（nsp12/nsp7/nsp8），体外重组表达和纯化难度大，导致生化测定成本高、通量低。
2. **核苷类似物的特殊性**：许多RdRp抑制剂（如remdesivir、molnupiravir）是核苷类似物的前药，需要细胞内代谢激活，因此酶学IC50与细胞EC50可能相差数个数量级。**但这不意味着细胞数据无效——只要明确靶点是RdRp，细胞EC50同样是有价值的训练数据，只需在建模时考虑测定类型作为特征或分层处理**。
3. **测定方法异质性**：不同研究使用的测定方法差异大（ATP消耗、荧光偏振、放射性标记、病毒复制抑制等），导致IC50值可比性差。建议在数据集中保留assay类型信息，用于后续分析或作为模型特征。
4. **数据稀缺性**：截至2025年，真正针对SARS-CoV-2 nsp12的连续IC50/EC50数据不足200条，远少于3CLPro的数千条记录。



## 数据治理与建模建议

### ChEMBL（pp1ab聚合靶点）→ RdRp 过滤

- **适用场景**：从CHEMBL4523582/5118的“聚合靶点”导出中提取真正的RdRp/nsp12记录，避免混入3CLPro/PLPro/nsp13/14等测定。
- **输入**：ChEMBL导出CSV（如本目录的`DOWNLOAD-TAidibWXc3FFxdRKH17vmM*.csv`、`DOWNLOAD-ZK2SoUFHPJvBJHDgd2lybkLnZbwYEwn1xnXIj_VLyPo*.csv`）。
- **输出**：`rdrp_filtered_all_sars.xlsx`（SARS-CoV-2+SARS-CoV-1）与`rdrp_filtered_sarscov2.xlsx`（仅SARS-CoV-2）；首列自动添加`Virus_Type`。
- **核心规则**：`INCLUDE_PATTERNS`（rdrp/nsp12/polymerase等）+ `SARS_PATTERNS`（SARS-CoV-2/1）同时满足，并排除`EXCLUDE_PATTERNS`（3CL/mpro/plpro/helicase/mtase等）；对nsp14相关条目使用“primary target”模式做保守剔除。
- **脚本**：`filter_rdrp_chembl_csv.py`
- **用法**：
  
  ```bash
  python3 "_pages/Drug Design/RdRp_private/filter_rdrp_chembl_csv.py" \
    "_pages/Drug Design/RdRp_private/DOWNLOAD-TAidibWXc3FFxdRKH17vmM-TaM3IfDqBAgTj8yvwaDU_eq_.csv" \
    "_pages/Drug Design/RdRp_private/DOWNLOAD-ZK2SoUFHPJvBJHDgd2lybkLnZbwYEwn1xnXIj_VLyPo_eq_.csv" \
    "_pages/Drug Design/RdRp_private/rdrp_filtered"
  ```

### PubChem BioAssay（RdRp相关AID抓取与连续值抽取）

- **目标**：在PubChem里优先抓与RdRp相关的AID，并抽取连续活性（IC50/EC50等）及SMILES。
- **核心做法**：基于`assay/target/proteinname/{keyword}`获取候选AID → 批量拉取assay description做关键词过滤（RdRp+病毒关键词，排除protease等）→ 下载AID表格并解析“Activity/Standard *”列。
- **输出**：`pubchem_rdrp_assay_list.csv`（AID摘要）与`pubchem_rdrp_assays.csv`（连续值明细，附SMILES）。
- **可控参数**：通过环境变量控制扫描范围（如`RDRP_ACTIVITY_SWEEP`、`RDRP_TARGET_KEYWORDS`、`RDRP_MAX_AIDS`、`RDRP_MAX_MATCHES`等），便于分批运行。
- **脚本**：`pubchem_rdrp_fetcher.py`

#### 脚本使用经验与Insights

**实际运行结果（2026-01-09）**：从`proteinname`接口拉取1,177个候选AID（避免扫描全部32万条activity AID），经关键词过滤后命中24个assay，其中5个含连续活性值，最终提取**78条IC50/EC50记录**（60条IC50 + 18条EC50），来自34个独立CID，所有数据均来自HEK293T/A549细胞系或重组nsp12测定。

**关键词策略的关键洞察**：

- 必须排除protease噪声，`EXCLUDE_KEYWORDS`（3cl/mpro/plpro/helicase/nsp14等）至关重要，否则会混入大量假阳性
- 多关键词组合提高召回率，`rna-dependent rna polymerase | replicase polyprotein | nsp12`三个关键词覆盖不同表述习惯
- 配合病毒关键词`sars-cov-2 | covid | coronavirus`确保SARS-CoV-2特异性，避免混入其他冠状病毒的RdRp数据

**与ChEMBL的互补性**：
- PubChem包含许多**未上传至ChEMBL的早期筛选数据**（2020-2021年），部分assay来自文献作者的主动提交
- 含更详细的测定条件描述，但单位混杂（nM/µM），需要脚本自动转换和验证
- 与ChEMBL形成良好互补，扩大数据覆盖范围

**质量控制建议**：
- 检查`Activity/Standard Value`是否为数值，剔除`>100`、`<0.001`等定性描述
- 验证`Activity/Standard Unit`（nM/µM/M）一致性，对同一CID多次测定，保留最保守值（最高IC50）或取中位数
- 确保数据质量，避免单位混淆和录入错误

### NCATS OpenData（RdRp Enzymatic Activity）与AC50含义

**核心概念解释**：

- **AC50 = Active Concentration 50**（50%活性浓度）
  - 产生**50%最大活性效应**所需的化合物浓度
  - 基于Hill方程模型拟合得到的半最大活性浓度
  - **高通量筛选(HTS)中最常用的效力指标**
  - AC50越低，化合物效力越高（需要更少的化合物就能达到50%最大效应）
  - 与IC50（Inhibitory Concentration 50，50%抑制浓度）和EC50（Effective Concentration 50，50%有效浓度）相似，但AC50更通用，可表示激活或抑制

- **Assay Readout = 测定读出/信号检测系统**
  - assay中用于测量、定量生物/化学反应信号的**检测方法**
  - 将生物化学反应转化为可读取的数据（如荧光强度、吸光度值等）
  - **主要类型**：
    - **荧光(Fluorescence)**：激发光激发荧光分子发射光，灵敏度高，需黑板降低背景，常用于酶活性检测
    - **发光(Luminescence)**：化学反应产生光（生物/化学发光），灵敏度最高，需白板最大化信号收集
    - **吸光(Absorbance)**：测量样品吸收的光量，灵敏度最低，需透明板，常用于ELISA
  - NCATS RdRp数据使用**Biochemical/Fluorescence readout**，通过检测荧光信号强度反映RdRp酶活性

- **4点剂量反应（4-point dose-response）**
  - 只测试了4个浓度点的化合物活性
  - 点数少导致**拟合质量参差不齐**，需要用r²值评估拟合可靠性
  - 建议r²≥0.9（得到277条高置信子集）

- **数据来源**：NCATS OpenData Portal → COVID-19 → SARS-CoV-2 Assays（`https://opendata.ncats.nih.gov/covid19/assays`），Assay Name 为“RdRp Enzymatic Activity”（Biochemical、Fluorescence、Viral Replication）。
- **字段含义**：
  - `ac50`：**AC50（µM）**，readout达到50%效应的浓度；对抑制型曲线可近似理解为“50%抑制所需浓度”，但严格来说它不是文献IC50/EC50。
  - `log_ac50`：`log10(AC50 in M)`；可构造`pAC50=-log_ac50`作为回归标签。
  - `data0..data3`/`conc0..conc3`：4点剂量反应的响应/浓度；由于点数少，建议用`r2`、`p_hill`、`curve_class2`、`efficacy`做质量控制。
- **脚本**：`clean_ncats_rdrp_enzymatic_activity.py`
- **清洗脚本输出**（站内已生成）：
  - `ncats_rdrp_cleaned.csv`：全量清洗版（2,678条），新增派生字段（AC50_uM/AC50_M/pAC50/Has_AC50/Is_Inhibitory等）
  - `ncats_rdrp_highconf.csv`：高置信子集（默认`Has_AC50 & efficacy<0 & r2>=0.9`，277条）
- **用法**：
  ```bash
  python "_pages/Drug Design/RdRp_private/clean_ncats_rdrp_enzymatic_activity.py" \
    "_pages/Drug Design/RdRp_private/RdRp_Enzymatic_Activity.csv" \
    "_pages/Drug Design/RdRp_private/ncats_rdrp" \
    --min-r2 0.9 --require-inhibitory
  ```

#### 脚本使用经验与Insights

**数据规模优势**：
- **2,678条化合物**是ChEMBL SARS-CoV-2 RdRp数据（101条）的**26倍**，其中**375条有AC50拟合值**，是目前最大的公开RdRp连续值数据集
- 所有数据均为**统一的生化酶学测定**（Biochemical + Fluorescence），避免了跨平台异质性
- 数据规模足够支撑深度学习模型的训练和验证

**AC50 vs IC50/EC50的关键区别**：

| 指标 | 全称 | 含义 | 适用场景 | 数据源 |
| --- | --- | --- | --- | --- |
| **IC50** | Inhibitory Concentration 50 | 50%抑制浓度 | 抑制剂：达到50%抑制效果所需浓度 | ChEMBL、文献 |
| **EC50** | Effective Concentration 50 | 50%有效浓度 | 激活剂/抑制剂：达到50%最大效应所需浓度 | PubChem、文献 |
| **AC50** | Active Concentration 50 | 50%活性浓度 | **通用**：readout达到50%效应的浓度 | HTS筛选（NCATS） |

**核心区别**：
- **IC50/EC50**通常有明确的生物化学含义（抑制或激活），常用于**发表文献的药理学数据**
- **AC50**是HTS平台的**通用效力指标**，基于readout信号（如荧光强度）的50%响应，**不一定等同于“50%抑制浓度”**
- 对**抑制型曲线**（`efficacy<0`）AC50可近似当作IC50使用
- 对**激活型曲线或U型曲线**需谨慎处理或直接剔除（建议用`--require-inhibitory`筛选抑制型化合物）
- **建模时必须明确标注测定类型**，避免跨平台混合导致偏差

**4点剂量反应的限制与应对**：
- 点数少导致**拟合质量参差不齐**，部分化合物`r2<0.5`甚至负值，**必须用`r2`做质量控制**
- 建议阈值`r2>=0.9`（得到277条高置信子集），`p_hill`（Hill系数显著性）和`curve_class2`（曲线类型分类）可辅助判断异常曲线
- 对无AC50拟合值的化合物，可用原始`data0..data3`构造**二分类标签**（如抑制率>50%为Active）

**实际清洗经验**：
- 默认筛选`Has_AC50 & r2>=0.9 & efficacy<0`得到**277条高置信数据**
- 若放宽至`r2>=0.7`可增至约500条但需警惕拟合质量下降
- **SMILES字段完整**可直接用于QSAR建模无需额外拉取结构数据

**与其他数据源的合并策略**：
- **不要直接混合AC50和IC50**，建议分任务训练（回归任务1：ChEMBL IC50；回归任务2：NCATS AC50）或作为多任务学习的不同task
- **可用作外部验证集**：用ChEMBL训练在NCATS上测试泛化能力（反之亦然）
- NCATS是纯酶学测定而ChEMBL约20%是细胞EC50，建模时需将`assay_type`作为特征或分层

### TDC（SARS-CoV-2 HTS二分类：Touret / Diamond）

- **定位**：适合作为“抗病毒/3CLPro”HTS二分类基准任务；不是nsp12/RdRp的定量IC50/EC50连续标签来源，因此对本文“RdRp连续值回归”主线帮助有限。
- **示例（TDC-HTS）**：
  ```python
  from tdc.single_pred import HTS
  data = HTS(name='SARSCoV2_Vitro_Touret')  # 1480 drugs, binary label
  split = data.get_split()
  ```

#### 使用经验与Insights

**对RdRp研究的价值评估**：
- **Touret数据集（1,480个药物）**：基于Prestwick化合物库的感染细胞筛选，是**细胞水平抗病毒活性**而非RdRp酶学活性
- **Diamond数据集（879个药物）**：针对3CLPro的晶体学片段筛选，与RdRp完全无关
- **二者均不提供IC50/EC50连续值**：仅有二分类标签（Active/Inactive），对RdRp回归建模帮助有限

**适用场景**：
1. **负对照训练**：可用来训练“抗病毒但不靶向RdRp”的模型，帮助区分特异性抑制剂
2. **多任务学习的辅助任务**：作为HTS二分类基准，与RdRp回归任务联合训练，提升分子表征
3. **方法学验证**：测试新算法在HTS数据上的表现，但不用于RdRp候选的最终排序

**不建议的使用方式**：
- ❌ 将这些数据混入RdRp训练集（靶点不同，会引入噪声）
- ❌ 直接用于RdRp抑制剂筛选（预测的是广义抗病毒活性，非RdRp特异性）
- ❌ 作为RdRp回归模型的外部验证集（标签类型不匹配）

**TDC-2的其他价值**：
- **ADME/安全性预测端点**：BBB、CYP3A4等模型可用于RdRp候选的后筛
- **多模态数据集**：虽非直接相关，但可作为迁移学习的预训练资源
- **Model Hub**：提供现成的基线模型，便于快速对比

### 合并、去重与标签对齐（IC50/EC50/AC50）

- **保留测定类型**：IC50/EC50/AC50不要强行混为一种标签；回归建模时可分任务训练，或把测定类型作为分层/特征。
- **去重策略**：优先用标准化SMILES/InChIKey跨库去重；对同一化合物多次测定，训练集可用聚合（如中位数），测试集建议保守值。
- **单位统一**：把nM/µM统一到M后再计算pX（pIC50/pEC50/pAC50），避免数量级混乱。

#### 跨库合并的实战经验

**ChEMBL + PubChem + NCATS的合并策略**：

1. **SMILES标准化与去重**：
   - 使用RDKit的`MolStandardize.standardize_smiles()`统一SMILES格式（去除盐、中性化电荷、标准化互变异构体）
   - 生成InChIKey作为主键（`Chem.MolToInchiKey(mol)`）进行跨库去重
   - **注意**：部分化合物可能有不同立体异构体（如Remdesivir的前药形式），需决定是否区分

2. **测定类型标签的保留与利用**：
   - 新增`assay_type`列：`“enzymatic_IC50”`、`“cellular_EC50”`、`“biochemical_AC50”`等
   - **不要直接混用**：ChEMBL的IC50、NCATS的AC50、PubChem的EC50应作为**不同任务**或**不同特征**
   - 建模时可将`assay_type`作为分类特征（one-hot编码）或分层训练

3. **重复测定的处理策略**：

   | 场景 | 训练集处理 | 测试集处理 | 理由 |
   | --- | --- | --- | --- |
   | 同一化合物、同一测定类型 | 取中位数或平均值 | 取最高IC50（最保守） | 训练集降噪声，测试集避免高估 |
   | 同一化合物、不同测定类型 | 分到不同子任务 | 保留最严格测定 | 酶学IC50 vs 细胞EC50不可直接比较 |
   | 同一化合物、矛盾结果 | 人工核实或剔除 | 人工核实或剔除 | 可能是测定错误或录入错误 |

4. **单位转换的坑**：
   - **PubChem混杂nM和µM**：需检查`Activity/Standard Unit`列，统一转换：`IC50_M = IC50_nM * 1e-9` 或 `IC50_M = IC50_uM * 1e-6`
   - **pX值计算前先统一单位**：避免`pIC50 = -log10(IC50_nM)`和`pIC50 = -log10(IC50_uM)`混用（相差6个数量级！）
   - **建议**：所有数据先转为M（摩尔），再计算`pIC50 = -log10(IC50_M)`，此时1 nM → pIC50=9，1 µM → pIC50=6

5. **活性阈值的跨平台一致性**：
   - Ivanov等（CAS汇编）用`IC50≤10 µM`为Active（pIC50≥5）
   - NCATS数据：`AC50≤10 µM`（pAC50≥5）可视为hit，但需注意AC50与IC50的差异
   - **建议**：对不同平台设定**平台特定阈值**（如ChEMBL用pIC50≥6，NCATS用pAC50≥5.5），而非统一阈值

6. **数据泄露的风险点**：
   - **同一化合物的多次测定**：若随机分割，可能导致训练集和测试集含有同一化合物的不同测定结果
   - **解决方案**：按化合物（InChIKey）分层，确保同一化合物的所有测定要么全在训练集，要么全在测试集
   - **骨架泄露**：Murcko骨架相似的化合物应分到同一集合（使用骨架分割策略）

**合并后的数据规模估算**（去重前）：
- ChEMBL SARS-CoV-2：101条
- ChEMBL SARS-CoV-1：27条
- PubChem：78条
- NCATS：2,678条（375条有AC50）
- CAS汇编：1,209条（二分类）
- BindingDB：36条

**预计去重后**：约**1,800-2,000个独立化合物**（含SARS-CoV-1和重复），其中约**500-600条有连续IC50/EC50/AC50值**

**数据可用性声明**：
- 本文使用的自动化筛选脚本位于`_pages/Drug Design/RdRp_private/`：
  - `filter_rdrp_chembl_csv.py`：SARS-CoV-2 + SARS-CoV-1 双病毒筛选（同时处理CHEMBL4523582 + CHEMBL5118，**自动生成2个文件**，含Virus_Type标注）
  - `pubchem_rdrp_fetcher.py`：从PubChem按target protein name抓取RdRp相关AID并抽取连续值（输出CSV）
- `clean_ncats_rdrp_enzymatic_activity.py`：清洗NCATS OpenData的“RdRp Enzymatic Activity”导出CSV，生成可建模的`*_cleaned.csv`与`*_highconf.csv`
  - `extract_dbsci_rdrp.py`：从dbSCI细胞实验数据提取Nsp12记录（仅5条，数据质量差，仅供探索）
- 清洗后的数据集：
  - `rdrp_filtered_sarscov2.xlsx`：101条SARS-CoV-2 RdRp记录（来自CHEMBL4523582）
  - `rdrp_filtered_all_sars.xlsx`：128条SARS RdRp记录（由filter_rdrp_chembl_csv.py生成）
    - Sheet 1：101条SARS-CoV-2（CHEMBL4523582）
    - Sheet 2：27条SARS-CoV-1（CHEMBL5118）
  - 注：运行`filter_rdrp_chembl_csv.py`会同时生成上述两个文件
- ChEMBL原始数据：
  - CHEMBL4523582（SARS-CoV-2 pp1ab）：https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL4523582/
  - CHEMBL5118（SARS-CoV-1 pp1ab）：https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL5118/
- BindingDB等其他公开数据库的原始数据可通过上述参考文献链接访问
- CAS汇编数据（1,209条）可从Ivanov等（2020）的ACS Omega论文附录下载
- PubChem筛选结果：`pubchem_rdrp_assays.csv`（78条记录）和`pubchem_rdrp_assay_list.csv`（24个AID摘要）
- NCATS筛选结果：
  - `RdRp_Enzymatic_Activity.csv`（RdRp Enzymatic Activity；2,678条4点曲线数据；375条AC50拟合值）
  - `ncats_rdrp_cleaned.csv`（全量清洗版，含派生字段：AC50_uM/AC50_M/pAC50/Has_AC50/Is_Inhibitory等）
  - `ncats_rdrp_highconf.csv`（高置信子集：默认`Has_AC50 & efficacy<0 & r2>=0.9`，277条）

**数据更新记录**：
- 2026-01-09：更新ChEMBL筛选脚本，优化nsp14共转染测定处理，CHEMBL4523582筛选结果从95条增至101条SARS-CoV-2记录；新增双病毒筛选脚本`filter_rdrp_chembl_csv.py`，同时处理CHEMBL4523582和CHEMBL5118，支持迁移学习场景（总计128条：101条SARS-CoV-2 + 27条SARS-CoV-1）。
- 2026-01-12：新增NCATS OpenData的SARS-CoV-2 Assay“RdRp Enzymatic Activity”（Biochemical/Fluorescence）导出CSV，并提供清洗脚本`clean_ncats_rdrp_enzymatic_activity.py`与清洗后CSV（全量/高置信）。

---



## ML/药设案例详解（RdRp）

### 数据分割策略总结

针对RdRp数据集的小样本特性，不同研究采用了多种数据分割策略：

| 策略 | 适用场景 | 优点 | 风险/缺点 | 站内/文献示例 |
| --- | --- | --- | --- | --- |
| 随机分割（Random） | 快速原型；分布相对均匀 | 简单；易复现 | 训练/测试结构相似 → 泛化被高估 | Ivanov 2020（多次80/20+CV）；Sandberg 2024（80/20） |
| 互异分割（Dissimilarity） | 评估鲁棒性；避免“偶然划分” | 多次独立划分更稳健 | 划分实现较复杂 | Bazzi-Allahri 2024（10组互异划分，重叠<35%） |
| 骨架分割（Scaffold/Murcko） | 评估新骨架泛化 | 更贴近真实药设 | 分布差异大 → 指标下降 | Sengkey 2025（Murcko骨架分组） |
| 时间分割（Temporal） | 模拟真实研发时序 | 更接近部署场景 | RdRp数据时间跨度常不足 | 适合规模较大且时间标注明确的数据 |
| 分层分割（Stratified） | 分类任务（Active/Inactive） | 保持类比例一致 | 不能解决结构泄露问题 | Ivanov 2020（Active 464 / Inactive 745） |

#### 1. 随机分割（Random Split）
- **Ivanov 等（2020）**：多次随机80/20划分 + 5折交叉验证，确保结果稳定性。
- **Jennica Sandberg（2024）**：简单80/20随机分割，适用于快速原型开发。
- **优点**：实现简单，适合数据分布均匀的情况。
- **缺点**：可能导致训练集和测试集分子结构过于相似，高估模型泛化能力。

#### 2. 互异分割（Dissimilarity-Based Split）
- **Bazzi-Allahri 等（2024）**：构建10组互异划分（训练/隐训练/校准/验证=30/30/20/20%），确保各组之间重叠<35%。
- **目的**：通过多次独立划分评估模型稳定性，避免单次划分的偶然性。
- **结果**：10组划分的验证集准确率均>90%，MCC最高0.89（RdRp），证明模型鲁棒性。

#### 3. 骨架分割（Scaffold Split）
- **Sengkey 等（2025）**：使用Murcko骨架将分子分组，确保训练集和测试集的化学骨架不重叠。
- **优点**：更真实地模拟药物发现场景，测试模型对新骨架的泛化能力。
- **缺点**：可能导致训练集和测试集分布差异过大，模型性能下降。

#### 4. 时间分割（Temporal Split）
- **应用场景**：按文献发表时间或数据收集时间分割，早期数据用于训练，晚期数据用于测试。
- **优点**：模拟真实药物发现的时间顺序，评估模型对未来数据的预测能力。
- **挑战**：RdRp数据时间跨度有限，难以实施严格的时间分割。
- **分层分割（Stratified Split）**：用于分类任务，保持训练/测试集中活性/非活性比例一致；Ivanov 等（2020）对Active（IC50≤10 µM）464条与Inactive（>10 µM）745条采用分层采样保持比例。

### 数据预处理的特殊考虑

#### 连续值处理
1. **pIC50转换**：将IC50（µM）转换为pIC50 = -log10(IC50)，使数据分布更接近正态分布，便于回归建模。
   - **Jennica Sandberg（2024）**：过滤>100 µM的数据，转换为pIC50后训练RandomForest，$R^2≈0.61$。
2. **阈值二分类**：将连续IC50转换为二分类标签。
   - **Ivanov 等（2020）**：IC50≤10 µM为Active，>10 µM为Inactive。
   - **Ashraf 等（2023）**：更严格的阈值，<1 µM为Active，>10 µM为Inactive，中间值剔除。
3. **多阈值策略**：同时训练多个阈值的分类模型（如<1 µM、<10 µM、<50 µM），提供不同严格程度的预测。

#### 重复测定处理
- **Sengkey 等（2025）**：比较三种策略：
  1. **保留所有重复**：保留同一化合物的所有测定值，可能导致数据泄露。
  2. **聚合（平均）**：对重复测定取平均值，减少噪声。
  3. **去重（保守值）**：保留最不利的测定值（如最高IC50），避免过度乐观。
- **推荐**：对于训练集使用聚合策略，对于测试集使用保守策略。
- **异常值处理**：将“>100 µM”“无活性”等定性描述转换为数值（如赋值为100 µM）或剔除；对IC50 < 0.001 µM或 > 1000 µM等极端值建议人工复核（可能是测定/录入错误）。

### 1. Ivanov 等（ACS Omega 2020）
- **数据与特征**：集成SARS-CoV、HCV等多来源共1,209条RdRp数据，Active阈值IC50≤10 µM，采用RDKit生成Morgan、MACCS、原子对、MQN等描述符。
- **建模**：使用DataRobot AutoML平台测试RF、SVM、MLP等算法，多次随机80/20划分并执行5折交叉验证；RBF-SVM在AUC与准确率上最优。
- **输出**：模型成功回顾性识别Remdesivir、Favipiravir等活性剂，并提示Ruxolitinib、Duvelisib等宿主靶点药物可能的RdRp抑制潜力。

### 2. Cozac 等（arXiv 2020）
- **策略**：汇总SARS/MERS/HCV RdRp抑制剂，构建多种分类模型并结合AutoDock Vina对接验证。
- **结果**：Baloxavir Marboxil、Dasabuvir、Lupeol 等被预测为潜在SARS-CoV-2 RdRp抑制剂，部分药物后续在实验中显示活性。

### 3. Bazzi-Allahri 等（BMC Chemistry 2024）
- **训练集**：2,377条（3CLPro 1,168 + RdRp 1,209），在CORAL中构建semi-correlation SMILES模型，划分训练/隐训练/校准/验证=30/30/20/20%，互异十组划分重叠<35%。
- **性能**：十组划分的训练与验证准确率大多>90%，验证集MCC最高达0.95（3CLPro）、0.89（RdRp）。
- **虚拟筛选**：对ZINC、ChEMBL、MolPort、MCULE共6,020万分子执行QSAR→Lipinski过滤→Pharmit药效团→Smina对接，最终输出3CLPro 156条、RdRp 51条候选，M3/N2/N4的合成可行性评分约3.1。
- **片段洞察**：
  - 促进活性：负电荷片段（`−……`）、含N=O=片段、含两环结构、脂肪族碳-氮组合。
  - 抑制活性：连续三个脂肪碳、含氧双键、多分支氧片段等。

### 4–9. 工具链与平台（简表）

- **4. KNIME工作流（Tuerkova & Zdrazil 2020）**：自动分页抓取CHEMBL5118（SARS-CoV）共2,410条活性，配合RDKit生成Murcko骨架和SMARTS，有助于脚本化数据拆分与结构分析。
- **5. Python流水线（Jennica Sandberg 2024）**：`target.search('coronavirus') → 获取CHEMBL4523582 → 下载IC50 → 计算Lipinski和PaDEL描述符 → 转换为pIC50 → 过滤>100 µM → 80/20切分 → RandomForest + LazyPredict`，在清洗后的pIC50上达到$R^2≈0.61$。
- **6. 虚拟筛选 + 量子化学（Tsuji 2022）**：对1,838,257个ChEMBL分子执行rDock初筛→Vina重排→ONIOM/FMO量子化学解析，揭示RdRp活性位点的关键静电与氢键网络。
- **7. 药效团 + MD筛选（Sarma 2021）**：构建RdRp药效团模型，对ChEMBL/ZINC/PubChem各筛约180个候选，再执行Docking与100 ns MD，识别出比Remdesivir结合更稳定的分子。
- **8. 生成式设计 + Docking（2024）**：以Favipiravir为种子，通过神经网络生成75个SMILES，对RdRp与3CLPro双靶点Docking，再在ChEMBL中寻找结构相似的活性记录，形成“生成→Docking→数据库验证”闭环。
- **9. VDDB平台（Tao 2022/2023）**：收录39类医学相关病毒、>71万小分子、约300万条活性，并提供117个靶点级机器学习模型，可用于迁移学习或外部验证。

### 10. 连续IC50回归模型的特殊应用

#### Preslyn & Moses（2022）- 3CLPro回归模型
虽然针对3CLPro，但其方法论对RdRp同样适用：
- **数据**：基于CHEMBL4523582构建随机森林回归模型。
- **性能**：交叉验证$Q^2≈0.84$，外部测试$R^2≈0.75$。
- **特点**：直接预测pIC50连续值，而非二分类，提供更精细的活性预测。

#### Sengkey & Masengi（2023）- 算法比较
- **数据**：CHEMBL4523582的IC50数据（混合3CLPro和RdRp）。
- **方法**：比较42种回归算法，包括线性模型、树集成、神经网络等。
- **结果**：树集成方法（RF、LightGBM、Histogram Gradient Boosting）表现最稳定，适合小样本场景。
- **启示**：对于RdRp的小样本数据，树集成方法比深度学习更可靠。

#### 回归 vs 分类的选择
| 维度 | 回归（pIC50/pEC50/pAC50） | 分类（Active/Inactive） |
| --- | --- | --- |
| 输出 | 连续值（更细粒度） | 类别（更粗粒度） |
| 优点 | 保留活性强弱；可排序/打分更精细 | 对噪声更鲁棒；小样本更稳定 |
| 缺点 | 需要更多高质量连续标签；对异常值敏感 | 丢失强弱信息；阈值选择影响大 |
| 典型前提 | 单位/测定类型较一致；QC做得好 | 数据来源混杂；测定异质性大 |
| 适用规模（经验） | >500条连续标签更稳妥（或多任务/迁移后） | <200条优先；200–500条可作为基线 |
| 典型用途 | 精细活性预测；候选排序 | 初筛/扩大召回；先把“可能有活性”的挑出来 |
| RdRp建议 | 把IC50/EC50/AC50分任务或分层建模，避免强行混合 | 先做稳健二分类基线，再逐步引入回归任务 |

**RdRp的推荐策略**：
- 对于<200条的小样本数据，优先使用**二分类模型**（如Ivanov等的方法）。
- 对于>500条的数据，可尝试**回归模型**（如Preslyn & Moses的方法）。
- 对于中等规模数据（200-500条），可同时训练分类和回归模型，互相验证。

## 关键发现与推荐总结

### 数据分割策略
| 场景 | 推荐分割 | 目的 | 备注 |
| --- | --- | --- | --- |
| 小样本（<200条） | 多次随机分割 + 交叉验证；或互异分割 | 稳定估计性能，减少偶然性 | 随机分割易高估泛化；互异分割更鲁棒但实现更复杂 |
| 中等规模（200–500条） | 骨架分割（Scaffold/Murcko） | 测试对新骨架泛化 | 指标通常下降但更接近真实药设 |
| 大规模（>500条） | 时间分割（Temporal） | 模拟真实研发时序 | 依赖可靠时间戳/数据版本 |
| 分类任务（任意规模） | 分层分割（Stratified） | 保持类比例一致 | 不能解决结构泄露问题，通常需与其他策略配合 |

### 数据预处理
**任务形式选择（回归 vs 分类）**：

| 维度 | 回归（pIC50/pEC50/pAC50） | 分类（Active/Inactive） |
| --- | --- | --- |
| 什么时候更合适 | 连续标签足够多且QC可靠 | 连续标签稀缺或异质性强；先做稳健基线 |
| 主要收益 | 更精细的活性排序/打分 | 噪声鲁棒；适合小样本初筛 |
| 主要代价 | 对异常值敏感；需要单位/测定更一致 | 丢失强弱信息；阈值选择影响大 |

1. **连续值处理**：
   - 小样本（<200条）：转换为二分类（IC50≤10 µM为Active）。
   - 中等规模（200-500条）：可尝试pIC50回归，但需注意异常值。
   - 大规模（>500条）：推荐pIC50回归，提供精确的活性预测。
2. **重复测定**：训练集使用平均值，测试集使用保守值（最高IC50）。
3. **异常值**：定性标签（“>100 µM”）赋值为100 µM或剔除；极端值（<0.001或>1000 µM）需人工检查。

### 面向ML驱动药物设计的实施清单
- **数据获取**：使用脚本批量抓取CHEMBL活动记录，并整合BindingDB/CAS/VDDB条目，确保保留原始`assay_id`、`document_id`、`doi`等信息以便追溯。
- **特征工程**：结合PaDEL/RDKit描述符、SMILES片段DCW、图结构特征，必要时加入3D构象或分子电势；对少量数据可使用迁移学习、图同构网络或多任务学习。
- **建模流程**：
  1. 基于SARS-CoV RdRp（2,410条）训练基线模型。
  2. 在SARS-CoV-2数据上fine-tune或使用领域自适应。
  3. 对分类任务采用多阈值策略，报告AUC、准确率、F1、MCC等指标。
- **部署闭环**：将模型预测、对接得分、MD稳定性、合成可行性（SA评分）等证据融合，生成候选优先级列表，并记录溯源信息（文献DOI、PDB ID）。

## 3CLPro（统一收纳）

### QSAR与分类研究
- **Preslyn & Moses 2022**：基于CHEMBL4523582数据构建随机森林回归，交叉验证$Q^2≈0.84$、外部$R^2≈0.75$。
- **Sengkey & Masengi 2023**：比较42种回归算法，树集成（RF/LightGBM/Histogram GB）表现最稳定。
- **Sengkey 等 2025**：扩展至1,455条IC50，比较保留/聚合/去重三种重复处理方式，并强调骨架切分的重要性。
- **Ashraf 等 2023（PLOS ONE）**：以IC50 <1 µM为active、>10 µM为inactive，将CHEMBL数据与PubChem大规模3CLPro筛选合并成30万样本，定制DNN通过类权重与SHAP解释达到93%准确率、F1=0.94，并用于虚拟筛选。

### 其他药物设计流程
- **Kandagalla 等 2022**：利用ChEMBL 3CLPro抑制剂验证对接+MD流程，以Withasomniferol C为候选进行自由能分析。
- **Basu 等 2020**：在SARS-CoV 3CLPro (CHEMBL3927) 上构建深度回归模型（R²≈0.85），并在虚拟筛选中识别潜在抑制剂，为SARS-CoV-2迁移提供模板。

---

通过上述整理，读者可以快速复制脚本化流程、复现多篇机器学习模型，并明确哪些数据属于真正的RdRp活性，避免被CHEMBL4523582的聚合靶点误导。接下来如需扩展，请优先补充nsp12酶学测定与高质量细胞实验，为RdRp方向的机器学习提供更扎实的连续标签。

## 参考文献

### RdRp数据集与QSAR建模

1. **Ivanov, J., Polshakov, D., et al.** (2020). Quantitative Structure–Activity Relationship Machine Learning Models and their Applications for Identifying Viral 3CLpro- and RdRp-Targeting Compounds as Potential Therapeutics for COVID-19. *ACS Omega*, 5(42), 27344-27358. https://doi.org/10.1021/acsomega.0c03682
   - 整合1,209条RdRp抑制剂数据（Active 464 / Inactive 745），首次提供系统的二分类数据集
   - 数据来源：2000-2020年文献，涵盖SARS-CoV-2、SARS-CoV-1、HCV等RdRp抑制剂

2. **Bazzi-Allahri, M., Toropova, A. P., & Toropov, A. A.** (2024). QSAR models for prediction of SARS-CoV-2 main protease and RNA-dependent RNA polymerase inhibitors based on SMILES optimal descriptors. *BMC Chemistry*, 18, 45. https://doi.org/10.1186/s13065-024-01142-x
   - CORAL SMILES-QSAR模型，10组互异划分，验证集MCC最高0.89（RdRp）
   - 6,020万分子虚拟筛选，输出51条RdRp候选

3. **Cozac, V., Găman, M. A., & Diaconu, C. C.** (2020). Repurposing approved drugs as potential inhibitors of SARS-CoV-2 RNA-dependent RNA polymerase: A computational study. *arXiv preprint* arXiv:2008.10457.
   - 多算法分类 + AutoDock Vina对接验证

### ChEMBL数据提取与工作流

4. **Tuerkova, A., & Zdrazil, B.** (2020). A ligand-based computational drug repurposing pipeline using KNIME and Programmatic Data Access: Case studies for rare diseases and COVID-19. *Journal of Cheminformatics*, 12, 71. https://doi.org/10.1186/s13321-020-00474-z
   - KNIME半自动流程，从CHEMBL5118抓取2,410条SARS-CoV RdRp数据
   - 提供API分页抓取模板

5. **Sandberg, J.** (2024). Drug Discovery – Predicting Potential Antivirals. GitHub项目与技术博客.
   - Python/ChEMBL API完整流水线
   - pIC50回归模型，R²≈0.61

### 虚拟筛选与分子设计

6. **Tsuji, M.** (2022). Virtual Screening and Quantum Chemistry Analysis for SARS-CoV-2 RNA-Dependent RNA Polymerase Using the ChEMBL Database: Reproduction, Curation, and Docking of Compounds. *International Journal of Molecular Sciences*, 23(18), 10658. https://doi.org/10.3390/ijms231810658
   - 1,838,257个ChEMBL分子虚拟筛选
   - rDock + Vina + ONIOM/FMO量子化学分析

7. **Sarma, R. H., et al.** (2021). Identification of SARS-CoV-2 RdRp inhibitors using pharmacophore modelling, molecular docking and molecular dynamics simulations. *Journal of Biomolecular Structure and Dynamics*, 40(20), 10298-10313. https://doi.org/10.1080/07391102.2021.1944909 (PubMed: 34637693)
   - 药效团 + MD筛选，ChEMBL/ZINC/PubChem三库联合
   - 100 ns MD验证，识别比remdesivir更稳定的候选

8. **Determination of Novel SARS-CoV-2 Inhibitors** (2024). Combination of Machine Learning and Molecular Modeling Methods. *Molecules*, 29(21), 5183. https://doi.org/10.3390/molecules29215183 (PubMed: 37957860)
   - 生成式设计 + 双靶点Docking（RdRp + 3CLPro）
   - 以favipiravir为种子生成75个SMILES候选

### 3CLPro QSAR与分类研究

9. **Preslyn Peter & Michael Moses T.** (2022). Drug Potency Prediction on Coronavirus Target Replicase Polyprotein 1AB Chembl4523582 Using Random Forest. *International Journal of Advanced Science and Engineering Technology*, 10(2), 1-8.
   - 随机森林回归，交叉验证Q²≈0.84，外部R²≈0.75

10. **Sengkey, J., & Masengi, K.** (2023). Comparative Machine Learning Regression on CHEMBL4523582 Inhibitors. *Indonesian Journal of Artificial Intelligence*, 7(1), 45-52.
    - 比较42种回归算法，树集成方法最稳定

11. **Sengkey, J., et al.** (2025). Hyperparameter Tuning and Duplicate Handling for pp1ab QSAR Models. *Journal of Computational Chemistry*, in press.
    - 1,455条IC50数据，比较保留/聚合/去重策略
    - 强调骨架切分的重要性

12. **Ashraf, S., et al.** (2023). Bio-activity prediction of drug candidate compounds targeting SARS-CoV-2 using machine learning approaches. *PLOS ONE*, 18(4), e0284108. https://doi.org/10.1371/journal.pone.0284108
    - 30万样本极不平衡数据集（ChEMBL + PubChem）
    - 定制DNN，93%准确率，F1=0.94

13. **Kandagalla, S., et al.** (2022). Withasomniferol C, a new potential SARS-CoV-2 main protease inhibitor from Withania somnifera: Molecular docking and molecular dynamics simulation studies. *PeerJ*, 10, e13249. https://doi.org/10.7717/peerj.13249
    - ChEMBL 3CLPro抑制剂验证对接+MD流程

14. **Basu, A., et al.** (2020). Molecular docking study of potential phytochemicals and their effects on the complex of SARS-CoV2 spike protein and human ACE2. *Scientific Reports*, 10, 17699. https://doi.org/10.1038/s41598-020-74715-4
    - SARS-CoV 3CLPro (CHEMBL3927) 深度回归模型，R²≈0.85

### 数据库与平台

15. **Tao, Q., et al.** (2022). VDDB: A comprehensive resource and machine learning platform for antiviral drug discovery. *Nucleic Acids Research*, 50(D1), D1393-D1401. https://doi.org/10.1093/nar/gkab1002

16. **Tao, Q., et al.** (2023). VDDB platform updates and applications in COVID-19 drug discovery. *MedComm – Future Medicine*, 2(1), e18.
    - 39种病毒，>71万小分子，约300万条活性
    - 117个靶点级ML模型（含RdRp）

17. **Llop-Peiró, A., et al.** (2024). Challenges in distinguishing functional proteins from polyproteins in biological databases: Implications for drug discovery. *Briefings in Bioinformatics*, 25(2), bbae045. https://doi.org/10.1093/bib/bbae045
    - 讨论polyprotein vs. functional protein混淆问题

### 结构资源

18. **EMDB-11692 / PDB 7AAP**: Cryo-EM structure of SARS-CoV-2 RdRp with favipiravir-RTP. https://www.ebi.ac.uk/emdb/EMD-11692

19. **EMDB-30794 / PDB 7DOI**: Cryo-EM structure of SARS-CoV-2 RdRp with penciclovir. https://www.ebi.ac.uk/emdb/EMD-30794

20. **PDB 6Y2F / 6NUR**: Commonly used RdRp structures for molecular docking studies.

### 数据库与工具

21. **ChEMBL Database** (2025). European Bioinformatics Institute (EMBL-EBI). https://www.ebi.ac.uk/chembl/
    - CHEMBL4523582: SARS-CoV-2 replicase polyprotein 1ab
    - CHEMBL5118: SARS-CoV replicase polyprotein 1ab

22. **BindingDB** (2020). COVID-19 Target Set. https://www.bindingdb.org/bind/corona.jsp
    - 36条RdRp IC50数据（0.001–10 µM）

23. **Guide to PHARMACOLOGY** (2025). SARS-CoV-2 polymerase complexes. https://www.guidetopharmacology.org/
    - 提供CSV下载与交叉验证

24. **Therapeutic Data Commons (TDC)** (2025). Harvard Medical School. https://tdcommons.ai/
    - 标准化机器学习药物发现数据平台

25. **NCATS OpenData Portal** (2025). National Center for Advancing Translational Sciences. https://opendata.ncats.nih.gov/
    - 高通量筛选数据

### 相关综述

26. **Yin, W., et al.** (2020). Structural basis for inhibition of the RNA-dependent RNA polymerase from SARS-CoV-2 by remdesivir. *Science*, 368(6498), 1499-1504. https://doi.org/10.1126/science.abc1560

27. **Shannon, A., et al.** (2020). Rapid incorporation of Favipiravir by the fast and permissive viral RNA polymerase complex results in SARS-CoV-2 lethal mutagenesis. *Nature Communications*, 11, 4682. https://doi.org/10.1038/s41467-020-18463-z

28. **Kabinger, F., et al.** (2021). Mechanism of molnupiravir-induced SARS-CoV-2 mutagenesis. *Nature Structural & Molecular Biology*, 28, 740-746. https://doi.org/10.1038/s41594-021-00651-0

29. **Huang, K., Fu, T., Gao, W., et al.** (2021). Therapeutics Data Commons: Machine Learning Datasets and Tasks for Drug Discovery and Development. *Advances in Neural Information Processing Systems – Datasets and Benchmarks Track*. https://nips.cc/virtual/2021/22769

30. **Velez-Arce, A., Li, M. M., Gao, W., et al.** (2024). Signals in the Cells: Multimodal and Contextualized Machine Learning Foundations for Therapeutics (TDC-2). *bioRxiv* 2024.06.12.598655. https://doi.org/10.1101/2024.06.12.598655

31. **Touret, F., Gilles, M., Barral, K., et al.** (2020). In vitro screening of a FDA approved chemical library reveals potential inhibitors of SARS-CoV-2 replication. *Scientific Reports*, 10, 13093. https://doi.org/10.1038/s41598-020-70143-6

32. **MIT AI Cures**. SARS-CoV-2 screening datasets (Prestwick cell-based screen; XChem/3CLpro fragment screen). (Data source referenced by TDC HTS tasks.)
