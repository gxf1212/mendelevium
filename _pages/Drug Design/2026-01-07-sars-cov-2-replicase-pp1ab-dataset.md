---
title: "SARS-CoV-2复制酶多蛋白1ab (CHEMBL4523582) 数据集：机器学习药物设计的真实用例"
date: "2026-01-07"
tags: [sars-cov-2, replicase, pp1ab, rdRp, chembl, dataset, qsar, machine-learning, drug-design]
description: "全面整理CHEMBL4523582（SARS-CoV-2复制酶多蛋白1ab）活性数据集，聚焦RdRp定量数据来源、处理策略与机器学习应用案例，并单列3CLPro研究。"
image: "img/thumbnail_mine/wh-g83qpe.jpg"
thumbnail: "img/thumbnail_mine/wh-g83qpe.jpg"
author: Xufan Gao
lang: zh-CN
---

# SARS-CoV-2复制酶多蛋白1ab (CHEMBL4523582) 数据集：机器学习药物设计的真实用例

## 概要

本文系统整理SARS-CoV-2 RNA依赖的RNA聚合酶（RdRp/nsp12）的连续IC50/EC50数据集现状，并提供机器学习药物设计的完整实施指南。

### CHEMBL和Moonshot不行
- **聚合靶点引入大量"伪RdRp"记录**：CHEMBL4523582将pp1ab所有酶域测定统合在单一靶点ID下，导致RdRp与3CLPro条目严重混杂；直接以靶点ID拉取数据会将3CLPro、PLPro、nsp3/4/13/14等测定全部混在一起，必须依据assay描述、文献DOI与关键字（polymerase、nsp12等）再过滤。
- **数据规模现状**：截至2025年底，CHEMBL4523582包含约1,455条IC50/EC50记录，但真正针对nsp12/RdRp的连续活性数据不足200条，其余多为3CLPro或其他非结构蛋白测定。
- **COVID Moonshot项目聚焦3CLPro**：该开放科学药物发现项目的大规模数据集主要针对主蛋白酶（Mpro/3CLpro），对RdRp研究的直接贡献有限。

### 数据资源
- **8个主要数据源**：ChEMBL（CHEMBL4523582和CHEMBL5118）、BindingDB、CAS汇编、PubChem BioAssay、VDDB、TDC、NCATS OpenData Portal、Guide to PHARMACOLOGY。
- **最系统的数据集**：CAS汇编的1,209条RdRp抑制剂数据（Active 464/Inactive 745，阈值IC50≤10 µM）是目前最完整的二分类数据集，广泛用于QSAR建模。
- **迁移学习资源**：SARS-CoV-1（CHEMBL5118）的完整数据集约2,410条记录，可用于预训练；我们从532行导出数据中筛选出27条nsp12记录，可与SARS-CoV-2数据联合使用。
- **清洗工具**：提供两个自动化筛选脚本：
  - **单病毒版本**：从CHEMBL4523582（17,929行）筛选出101条SARS-CoV-2 RdRp记录（72条IC50 + 9条EC50 + 19条Inhibition）
  - **双病毒版本**：同时处理CHEMBL4523582和CHEMBL5118，输出128条记录（101条SARS-CoV-2 + 27条SARS-CoV-1），自动标注Virus_Type，支持迁移学习

### 方法论贡献
- **5种数据分割策略**：随机分割（适合快速原型）、互异分割（评估稳定性）、骨架分割（测试新骨架泛化）、时间分割（模拟真实场景）、分层分割（保持活性比例）。
- **连续值处理方法**：pIC50转换（使数据分布正态化）、阈值二分类（IC50≤10 µM为Active）、多阈值策略（同时训练多个阈值模型）。
- **重复测定处理**：比较保留所有重复、聚合（平均）、去重（保守值）三种策略，推荐训练集使用聚合、测试集使用保守值。
- **应对数据稀缺**：迁移学习（利用SARS-CoV或HCV RdRp数据）、多任务学习（联合训练多个SARS-CoV-2靶点）、半监督学习（利用无标签分子）、数据增强（分子片段替换、构象采样）。

### 应用案例
- **10个ML/药设案例**：Ivanov等的AutoML QSAR（1,209条，RBF-SVM最优）、Bazzi-Allahri等的CORAL SMILES模型（10组互异划分，验证集MCC最高0.89）、Tsuji的量子化学联动虚拟筛选（1,838,257个分子）、Sarma的药效团+MD筛选（100 ns验证）、生成式分子设计（75个SMILES候选）等。
- **回归vs分类选择标准**：小样本（<200条）优先二分类（RF、SVM、CORAL），中等规模（200-500条）可尝试pIC50回归（树集成方法），大规模（>500条）推荐深度学习（DNN、GNN）。
- **QSAR+药效团+对接三段式虚拟筛选**：可在6,020万分子中缩减至数十个候选，结合合成可行性评分（SA≈3.1）输出优先级。

### 结构资源
- **Cryo-EM结构**：EMD-11692/7AAP（含favipiravir-RTP）、EMD-30794/7DOI（含penciclovir）。
- **对接模板**：PDB 6Y2F/6NUR（常用于分子对接研究）。
- **官方知识库**：Guide to PHARMACOLOGY提供CSV下载与交叉验证。

### 参考文献
- **28篇关键文献**：涵盖数据集构建（Ivanov 2020、Bazzi-Allahri 2024）、QSAR建模（Preslyn & Moses 2022、Sengkey系列）、虚拟筛选（Tsuji 2022、Sarma 2021）、3CLPro研究（Ashraf 2023、Kandagalla 2022）、数据库平台（VDDB、TDC）和核苷类似物机制（remdesivir、favipiravir、molnupiravir）。
- **所有文献均附DOI和详细说明**，便于追溯和复现。
- **3CLPro研究独立章节**：为避免靶点混淆，所有3CLPro相关研究统一收纳至文末独立章节。

## 数据现状与测定来源

### 为什么找不到"5,000条RdRp活性数据"？
1. **聚合式靶点设计**：CHEMBL将pp1ab映射到单一目标（CHEMBL4523582），任何nsp5、nsp12、nsp13、nsp15的测定结果（无论是酶学还是细胞水平）都会归到同一target，导致靶点标签混乱。**关键问题不是测定类型（酶学vs细胞），而是靶点归属不清**。
2. **Guide-to-Pharmacology互链误导**：在Ligand Activity Chart中，JAK3 inhibitor IV、GC-376、SU-3327、PBIT等被标注为"RNA-dependent RNA polymerase"抑制剂，但原始测定是3CLPro荧光FRET，靶点完全错误。
3. **数量级差异的客观存在**：依据Ashraf 等（2023, PLOS ONE），CHEMBL在2020–2021年仅能提供161条pp1ab IC50，其中<1 µM活性仅5条，失活>10 µM有156条；到2024年增至约1,455条，但真正属于RdRp的仍属少数。
4. **COVID Moonshot项目聚焦3CLPro**：COVID Moonshot是开放科学药物发现项目，但其主要靶点是SARS-CoV-2主蛋白酶（Mpro/3CLpro），而非RdRp，因此该项目的大规模数据集对RdRp研究的直接贡献有限。

### RdRp数据源对照表

| 数据源 | 规模/时间 | 连续label数据量 | 采集说明 | 活性类型 |
| --- | --- | --- | --- | --- |
| ChEMBL (CHEMBL4523582) | SARS-CoV-2 pp1ab，持续增长 | **101条**（站内脚本筛选所得） | **数据来源**：DOWNLOAD-TAidibWXc3FFxdRKH17vmM*.csv（17,929行原始数据）<br/>需结合assay描述和关键词"nsp12""polymerase"筛选，从pp1ab聚合靶点中分离出RdRp专属记录；使用`filter_rdrp_csv.py`自动化筛选，输出XLSX格式。 | 酶学IC50、细胞EC50 |
| ChEMBL (CHEMBL5118) | SARS-CoV-1 pp1ab，约2,410条总记录 | **27条**（站内脚本筛选所得）<br/>**2,410条**（完整数据集，可用于迁移学习） | **数据来源**：DOWNLOAD-ZK2SoUFHPJvBJHDgd2lybkLnZbwYEwn1xnXIj_VLyPo*.csv（532行导出数据，筛得27条nsp12记录）<br/>CHEMBL5118完整数据集约2,410条，数据质量相对较高，是SARS-CoV-2数据不足时的重要迁移学习资源[^5]；使用`filter_rdrp_chembl_csv.py`可同时筛选两个靶点并标注Virus_Type。 | IC50、Ki |
| BindingDB COVID-19集合 | 2020年打包发布 | **36条** | 下载“COVID-19 target set”，包含36条RdRp IC50（0.001–10 µM），含部分与ChEMBL重叠的测定，可做交叉验证[^1]。 | IC50 |
| PubChem BioAssay | 零散更新；站内脚本于2026-01-09针对protein name抓取1177个候选AID | **78条**（5个AID，60条IC50 + 18条EC50，均为nM） | 依托`_pages/Drug Design/RdRp_private/pubchem_rdrp_fetcher.py`：禁用全量activity扫描，仅用`assay/target/proteinname/{keyword}`接口拉取“rna-dependent rna polymerase/replicase polyprotein/nsp12”AID，再以关键词过滤与表格解析抽取连续值，输出`pubchem_rdrp_assays.csv`与`pubchem_rdrp_assay_list.csv`。 | IC50、EC50 |
| VDDB（Tao 2022/2023） | 39种病毒、>71万小分子、约300万条活性 | **约300万条**（117靶点合计，未拆分到单靶点） | 自带117个靶点级机器学习模型（含RdRp）；支持下载与在线预测，可用于外部验证或迁移学习[^10]。https://vddb.idruglab.cn/，现在访问不了。。。 | 靶点预测、QSAR输出 |
| Therapeutic Data Commons (TDC) | 核心“TDC-1”提供66个标准化数据集（22个任务），`2024-06-22`发布的“TDC-2”再增加>10种模态、>1,000个多模态数据集（≈8,500万细胞）与Model Hub | **多模态**（小分子、蛋白、单细胞等；可筛选COVID/RdRp相关任务） | Python API在“problem→task→dataset”层级交付ML-ready拆分、指标与评测脚本，另含ADME/安全性预测端点，可用作RdRp候选的外部验证或联合任务。 | 多种活性类型（QSAR、生成、单细胞） |
| NCATS OpenData Portal | 高通量筛选数据 | **待确认** | 美国国家转化科学促进中心（NCATS）的开放数据门户可能包含SARS-CoV-2 RdRp筛选数据，需访问官网查询具体AID。 | 抑制率、IC50 |
| CAS汇编（Ivanov等，2020） | 1,209个化合物 | **1,209条** | Active 464 / Inactive 745（阈值IC50≤10 µM）；覆盖SARS-CoV-2、SARS-CoV-1、HCV等RdRp抑制剂；SMILES与标签在ACS Omega论文附录开放，成为后续QSAR模型的标准训练集[^2]。 | 二分类IC50标签 |

**PubChem增量说明（2026-01-09）**：为了验证PubChem是否存在更多隐含的nsp12条目，我编写了`pubchem_rdrp_fetcher.py`并在2026-01-09运行，核心做法是：
1. 仅通过`https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/proteinname/{keyword}/aids/TXT`接口拉取“rna-dependent rna polymerase | replicase polyprotein | nsp12”三组AID（共1,177个），避免扫描全部32万条activity AID；
2. 批量获取assay description，按关键词组合(`POLYMERASE_KEYWORDS` + `VIRUS_KEYWORDS`)且排除3CL/Mpro等噪声(`EXCLUDE_KEYWORDS`)，命中后再下载CSV；
3. 扫描`Activity/Standard Value|Type|Unit`列抽取连续活性值，并批量拉取SMILES，最终落地两份文件：
   - `pubchem_rdrp_assay_list.csv`：24个匹配assay的摘要（AID、来源、描述等）；
   - `pubchem_rdrp_assays.csv`：5个AID合计78条nM级IC50/EC50记录（34个独立CID，全部来自HEK293T/A549细胞或重组nsp12测定）。
  这些CSV位于`_pages/Drug Design/RdRp_private/`目录，可直接并入后续QSAR流水线；若要扩大范围，只需调整环境变量`RDRP_MAX_AIDS`、`RDRP_TARGET_KEYWORDS`等即可快速复现。

> **TDC补充（截至2026-01-09）**：NeurIPS 2021的TDC-1论文确认平台已整合66个AI-ready数据集、22类学习任务、33个数据函数、23套系统化评测策略、17个分子生成oracle与29条公开排行榜，为药物发现ML提供统一入口（参考文献29）。2024-06-23发布的TDC-2进一步引入10+新模态、>1,000个多模态数据集（累计≈8,500万细胞）、跨尺度上下文化任务，以及可直接调用BBB、CYP3A4等模型的HuggingFace式Model Hub，可将这些ADME/安全性端点并入RdRp候选的后筛或多任务训练流程（参考文献30）。

> 数据对照来自站内PDF《SARS-CoV-2 RdRp抑制剂数据集与机器学习模型综述》以及BMC Chemistry 2024论文附录。

### 数据增长时间线
- **2020年初**：CHEMBL几乎无SARS-CoV-2 RdRp记录，研究者多使用SARS-CoV（CHEMBL5118）或HCV RdRp数据做迁移。
- **2020年10月**：Ivanov 等整合了1,209条RdRp抑制剂，首次提供系统的Active/Inactive标签。
- **2021年**：Eur J. Med. Chem. 报道的HEK293T转染实验提供了十余条nsp12 (±nsp14/nsp10) 细胞系IC50。
- **2022–2023年**：BindingDB、PubChem以及多项分子设计研究陆续上传新的酶学数据，ChEMBL条目累计到千级，但3CLPro/Hydrolase等测定增长更快。
- **2024年**：Bazzi-Allahri 等使用2,377条分子（3CLPro+RdRp）在CORAL中构建SMILES-QSAR模型，并开展6,020万规模的虚拟筛选。

### 结构与表征资源

- **Cryo-EM与PDB**：2020年以来已发布多份含nsp12-nsp8-nsp7-RNA复合物的高分辨率结构，为分子动力学、量子化学与虚拟筛选提供了可靠模板：
  - **EMD-11692 / 7AAP**：掺入favipiravir-RTP的RdRp复合物结构
  - **EMD-30794 / 7DOI**：掺入penciclovir的RdRp复合物结构
  - **PDB 6Y2F / 6NUR**：常用于分子对接的RdRp结构模板
- **官方知识库**：Guide to PHARMACOLOGY将RdRp列为“SARS-CoV-2 polymerase complexes”并提供CSV下载链接，适合作为ChEMBL以外的交叉验证。

### 连续IC50/EC50数据的特殊挑战

与3CLPro相比，SARS-CoV-2 RdRp的连续活性数据面临以下独特挑战：

1. **测定复杂度高**：RdRp是多亚基复合物（nsp12/nsp7/nsp8），体外重组表达和纯化难度大，导致生化测定成本高、通量低。
2. **核苷类似物的特殊性**：许多RdRp抑制剂（如remdesivir、molnupiravir）是核苷类似物的前药，需要细胞内代谢激活，因此酶学IC50与细胞EC50可能相差数个数量级。**但这不意味着细胞数据无效——只要明确靶点是RdRp，细胞EC50同样是有价值的训练数据，只需在建模时考虑测定类型作为特征或分层处理**。
3. **测定方法异质性**：不同研究使用的测定方法差异大（ATP消耗、荧光偏振、放射性标记、病毒复制抑制等），导致IC50值可比性差。建议在数据集中保留assay类型信息，用于后续分析或作为模型特征。
4. **数据稀缺性**：截至2025年，真正针对SARS-CoV-2 nsp12的连续IC50/EC50数据不足200条，远少于3CLPro的数千条记录。

**应对策略**：
- **迁移学习**：利用SARS-CoV（CHEMBL5118，约2,410条）、HCV RdRp等同源靶点数据进行预训练。
- **多任务学习**：联合训练RdRp、3CLPro、PLPro等多个SARS-CoV-2靶点，共享分子表征。
- **半监督学习**：利用大量无标签分子的结构信息（如自监督预训练的图神经网络）。
- **数据增强**：通过分子片段替换、构象采样等方法扩充训练集。

## 数据提取与脚本化清洗

### ChEMBL数据提取流程

1. **API抓取**：分别从两个ChEMBL靶点导出数据
   - **CHEMBL4523582**（SARS-CoV-2 pp1ab）：`https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id=CHEMBL4523582&limit=0`
     - 导出文件：DOWNLOAD-TAidibWXc3FFxdRKH17vmM*.csv（17,929行）
   - **CHEMBL5118**（SARS-CoV-1 pp1ab）：`https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id=CHEMBL5118&limit=0`
     - 导出文件：DOWNLOAD-ZK2SoUFHPJvBJHDgd2lybkLnZbwYEwn1xnXIj_VLyPo*.csv（532行）

2. **关键词过滤**：提供两个自动化筛选脚本，位于`_pages/Drug Design/RdRp_private/`：

   | 脚本 | 病毒范围 | 输出文件 | ChEMBL靶点 | 用途 |
   |------|---------|---------|----------|------|
   | `filter_rdrp_csv.py` | 仅SARS-CoV-2 | 1个文件（101条） | CHEMBL4523582 | COVID-19特异性研究，自动排除SARS-CoV-1 |
   | `filter_rdrp_chembl_csv.py` | SARS-CoV-2 + SARS-CoV-1 | **2个文件**（128条all_sars + 101条sarscov2） | CHEMBL4523582 + CHEMBL5118 | 迁移学习、数据增强，**同时生成双病毒和单病毒版本**，自动添加Virus_Type列标注 |

   **关键词配置**：
   - **INCLUDE_PATTERNS**：`"rdrp"`, `"rna-dependent rna polymerase"`, `"nsp12"`, `"polymerase complex"`, `"recombinant rdrp"`, `"coronavirus rna polymerase"`
   - **EXCLUDE_PATTERNS**：`"3cl"`, `"main protease"`, `"mpro"`, `"papain-like"`, `"plpro"`, `"nsp3/4/5/13/15"`, `"helicase"`, `"methyltransferase"`
   - **SARS_PATTERNS**（双病毒脚本）：自动识别`"sars-cov-2"`、`"covid-19"`、`"severe acute respiratory syndrome-related coronavirus"`等

3. **nsp14测定处理**：
   - **保留**：描述为"nsp12与nsp14/nsp10共转染"的RdRp测定（如CHEMBL4831803，6条记录）
   - **剔除**：`"nsp14 methyltransferase"`, `"nsp14 N7-MTase"`, `"stabilization of sars-cov-2 nsp14"`等以nsp14为主靶的实验

4. **输出文件**：

   **单病毒版本** (`filter_rdrp_csv.py`)：
   - **`rdrp_filtered_sarscov2.xlsx`**：仅处理CHEMBL4523582数据
     - 输入：DOWNLOAD-TAidibWXc3FFxdRKH17vmM*.csv（17,929行原始数据）
     - 输出：101条SARS-CoV-2 RdRp记录
     - 活性类型：IC50 (72条，71.3%) + Inhibition (19条，18.8%) + EC50 (9条，8.9%)
     - Target Organism：100% Severe acute respiratory syndrome coronavirus 2

   **双病毒版本** (`filter_rdrp_chembl_csv.py`)：**同时生成2个输出文件**

   - **文件1: `{prefix}_all_sars.xlsx`**（所有SARS病毒数据，128条记录，2个sheet）
     - **Sheet 1 (CHEMBL4523582)**：101条SARS-CoV-2记录
       - 来源：DOWNLOAD-TAidibWXc3FFxdRKH17vmM*.csv
       - Virus_Type：全部标注为"SARS-CoV-2"
     - **Sheet 2 (CHEMBL5118)**：27条SARS-CoV-1记录
       - 来源：DOWNLOAD-ZK2SoUFHPJvBJHDgd2lybkLnZbwYEwn1xnXIj_VLyPo*.csv（从532行中筛选）
       - Virus_Type：全部标注为"SARS-CoV-1"
       - 与归档文件`rdrp_filtered_zk.csv`100%匹配验证
     - 应用场景：迁移学习（利用SARS-CoV-1预训练）、多任务学习、对比分析

   - **文件2: `{prefix}_sarscov2.xlsx`**（仅SARS-CoV-2数据，101条记录，1个sheet）
     - **Sheet 1 (CHEMBL4523582)**：101条SARS-CoV-2记录
       - 从all_sars文件中自动筛选
       - Virus_Type：全部标注为"SARS-CoV-2"
     - 应用场景：COVID-19特异性研究，与filter_rdrp_csv.py输出等效

   - 新增列：**Virus_Type**（第一列，自动病毒类型标注）

5. **脚本使用示例**：

   ```bash
   # 场景1: 仅筛选SARS-CoV-2数据（COVID-19特异性研究）
   # 输入：CHEMBL4523582导出文件
   # 输出：101条SARS-CoV-2记录
   python3 filter_rdrp_csv.py \
     DOWNLOAD-TAidibWXc3FFxdRKH17vmM-TaM3IfDqBAgTj8yvwaDU_eq_.csv \
     rdrp_filtered_sarscov2.xlsx

   # 场景2: 筛选所有SARS病毒数据（迁移学习、数据增强）
   # 输入：CHEMBL4523582 + CHEMBL5118两个导出文件
   # 输出：自动生成2个文件
   #   - rdrp_filtered_all_sars.xlsx (128条：101条SARS-CoV-2 + 27条SARS-CoV-1)
   #   - rdrp_filtered_sarscov2.xlsx (101条：仅SARS-CoV-2)
   python3 filter_rdrp_chembl_csv.py \
     DOWNLOAD-TAidibWXc3FFxdRKH17vmM-TaM3IfDqBAgTj8yvwaDU_eq_.csv \
     DOWNLOAD-ZK2SoUFHPJvBJHDgd2lybkLnZbwYEwn1xnXIj_VLyPo_eq_.csv \
     rdrp_filtered
   ```

6. **测定类型建议**：保留`assay_type`或`assay_description`字段以供分层建模、特征工程与后处理分析。

7. **CAS Supporting Information**：
   - `quantitative-structure-activity-relationship-machine-learning-models-and-their-applications-for-identifying-viral.pdf`提供Ivanov研究的完整建模流程。
   - `ao0c03682_si_001.xlsx`包含二分类标签、模型AUC与候选列表：Table S4 记录1,209条RdRp Active/Inactive标签及SMILES，可导出为`rdrp_cas_dataset.csv`；其余表格提供3CLpro/RdRp预测命中和累积增益/提升分析（无原始IC50/EC50数值）。

## ML/药设案例详解（RdRp）

### 数据分割策略总结

针对RdRp数据集的小样本特性，不同研究采用了多种数据分割策略：

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

#### 5. 分层分割（Stratified Split）
- **应用**：对于分类任务，确保训练集和测试集中活性/非活性分子的比例一致。
- **Ivanov 等（2020）**：Active（IC50≤10 µM）464条，Inactive（>10 µM）745条，采用分层采样保持比例。

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

#### 异常值处理
- **定性标签**：“>100 µM”、“无活性”等定性描述需要转换为数值（如赋值为100 µM或剔除）。
- **极端值**：IC50 < 0.001 µM或 > 1000 µM的数据需要人工检查，可能是测定错误。

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

### 4. KNIME工作流（Tuerkova & Zdrazil 2020）
- 自动分页抓取CHEMBL5118（SARS-CoV）共2,410条活性，配合RDKit生成Murcko骨架和SMARTS，有助于脚本化数据拆分与结构分析。

### 5. Python流水线（Jennica Sandberg 2024）
- `target.search('coronavirus') → 获取CHEMBL4523582 → 下载IC50 → 计算Lipinski和PaDEL描述符 → 转换为pIC50 → 过滤>100 µM → 80/20切分 → RandomForest + LazyPredict`，在清洗后的pIC50上达到$R^2≈0.61$。

### 6. 虚拟筛选 + 量子化学（Tsuji 2022）
- 对1,838,257个ChEMBL分子执行rDock初筛→Vina重排→ONIOM/FMO量子化学解析，揭示RdRp活性位点的关键静电与氢键网络。

### 7. 药效团 + MD筛选（Sarma 2021）
- 构建RdRp药效团模型，对ChEMBL/ZINC/PubChem各筛约180个候选，再执行Docking与100 ns MD，识别出比Remdesivir结合更稳定的分子。

### 8. 生成式设计 + Docking（2024）
- 以Favipiravir为种子，通过神经网络生成75个SMILES，对RdRp与3CLPro双靶点Docking，再在ChEMBL中寻找结构相似的活性记录，形成“生成→Docking→数据库验证”闭环。

### 9. VDDB平台（Tao 2022/2023）
- 收录39类医学相关病毒、>71万小分子、约300万条活性，并提供117个靶点级机器学习模型，可用于迁移学习或外部验证。

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
- **回归模型**：
  - 优点：保留完整的活性信息，可预测精确的IC50值。
  - 缺点：需要更多高质量数据，对异常值敏感。
  - 适用场景：数据质量高、样本量充足（>500条）、需要精确活性预测。
- **分类模型**：
  - 优点：对噪声和异常值更鲁棒，适合小样本。
  - 缺点：损失活性的精细信息，无法区分强活性和中等活性。
  - 适用场景：数据质量参差不齐、样本量有限（<200条）、主要用于虚拟筛选初筛。

**RdRp的推荐策略**：
- 对于<200条的小样本数据，优先使用**二分类模型**（如Ivanov等的方法）。
- 对于>500条的数据，可尝试**回归模型**（如Preslyn & Moses的方法）。
- 对于中等规模数据（200-500条），可同时训练分类和回归模型，互相验证。

## 数据治理与建模建议
1. **Assay语义过滤**：抓取CHEMBL4523582后，按`assay_description`、`document_id`、`assay_type`过滤含"polymerase""nsp12"等关键字的记录，同时剔除含"main protease""papain-like"等描述。**关键是靶点归属，而非测定类型**。
2. **跨库验证**：将ChEMBL记录与GtoPdb、BindingDB、VDDB等交叉比对，被标注为3CLPro抑制剂的条目标记为"非RdRp"，避免进入训练集。
3. **重复测定处理**：对重复IC50值可采用平均或最保守值；对定性（>100 µM）记录可酌情剔除或赋予高值。
4. **测定类型标注**：保留assay类型信息（酶学vs细胞、IC50 vs EC50），可用于：
   - **分层分析**：评估模型在不同测定类型上的表现
   - **多任务学习**：将测定类型作为辅助任务
   - **特征工程**：作为分类特征或用于数据归一化
5. **结构域标签**：在多结构域复合体（nsp12/nsp7/nsp8）场景下加入结构域或PDB模板标签，供图神经网络或多任务模型使用。
6. **少样本策略**：先在SARS-CoV（CHEMBL5118）数据上训练基线，再迁移到SARS-CoV-2；分类任务建议使用Active <1 µM、Inactive >10 µM，保留中间值用于软标签或回归。
7. **不平衡与主动学习**：采用类权重、分层采样、阈值移动、active learning等方式应对极端失衡，并提供骨架切分、时间切分等多种评估结果。

## 关键发现与推荐总结

### 数据现状
1. **真实RdRp数据稀缺**：CHEMBL4523582虽有1,455条记录，但真正针对nsp12/RdRp的连续IC50/EC50数据不足200条，大部分是3CLPro或其他非结构蛋白的测定。
2. **CAS汇编是最系统的数据源**：Ivanov等（2020）整合的1,209条RdRp抑制剂数据（Active 464 / Inactive 745）是目前最完整的二分类数据集，广泛用于QSAR建模。
3. **SARS-CoV数据可用于迁移**：CHEMBL5118包含约2,410条SARS-CoV RdRp活性记录，可用于预训练或迁移学习。
4. **COVID Moonshot聚焦3CLPro**：该项目的大规模数据集主要针对主蛋白酶，对RdRp研究的直接贡献有限。

### 数据分割策略
1. **小样本场景**（<200条）：推荐**多次随机分割 + 交叉验证**（如Ivanov等的方法），或**互异分割**（如Bazzi-Allahri等的10组划分）。
2. **中等规模**（200-500条）：推荐**骨架分割**，更真实地评估模型对新化学空间的泛化能力。
3. **大规模场景**（>500条）：可尝试**时间分割**，模拟真实药物发现的时间顺序。
4. **分类任务**：必须使用**分层分割**，确保训练集和测试集的活性/非活性比例一致。

### 数据预处理
1. **连续值处理**：
   - 小样本（<200条）：转换为二分类（IC50≤10 µM为Active）。
   - 中等规模（200-500条）：可尝试pIC50回归，但需注意异常值。
   - 大规模（>500条）：推荐pIC50回归，提供精确的活性预测。
2. **重复测定**：训练集使用平均值，测试集使用保守值（最高IC50）。
3. **异常值**：定性标签（“>100 µM”）赋值为100 µM或剔除；极端值（<0.001或>1000 µM）需人工检查。

### 机器学习方法选择
1. **小样本（<200条）**：
   - **推荐**：随机森林、SVM、CORAL SMILES模型。
   - **避免**：深度神经网络（容易过拟合）。
2. **中等规模（200-500条）**：
   - **推荐**：树集成方法（RF、LightGBM、XGBoost）。
   - **可尝试**：浅层神经网络、图神经网络（需结合迁移学习）。
3. **大规模（>500条）**：
   - **推荐**：深度学习（DNN、GNN）、集成学习。
   - **可尝试**：多任务学习（联合训练RdRp、3CLPro等多个靶点）。

### 迁移学习策略
1. **同源靶点迁移**：在SARS-CoV（CHEMBL5118）或HCV RdRp数据上预训练，再fine-tune到SARS-CoV-2。
2. **多任务学习**：联合训练RdRp、3CLPro、PLPro等多个SARS-CoV-2靶点，共享分子表征。
3. **自监督预训练**：在大规模无标签分子库（如ZINC、ChEMBL）上预训练图神经网络，再迁移到RdRp任务。

### 虚拟筛选流程
1. **QSAR初筛**：使用训练好的模型对大规模分子库（如ZINC、ChEMBL）进行初筛，保留预测活性高的分子。
2. **药效团过滤**：基于已知活性分子构建药效团模型，进一步筛选。
3. **分子对接**：对通过前两步的分子进行分子对接，评估结合模式和亲和力。
4. **MD验证**：对对接得分高的分子进行分子动力学模拟，评估结合稳定性。
5. **合成可行性评分**：使用SA评分或其他工具评估合成难度，优先选择易合成的分子。

### 数据获取建议
1. **ChEMBL**：使用API批量抓取CHEMBL4523582，严格过滤assay描述，保留nsp12/RdRp相关记录。
2. **BindingDB**：下载COVID-19 target set，提取RdRp相关数据。
3. **CAS汇编**：从Ivanov等（2020）的ACS Omega论文附录下载1,209条数据。
4. **SARS-CoV数据**：从CHEMBL5118下载约2,410条SARS-CoV RdRp数据，用于迁移学习。
5. **文献挖掘**：系统检索PubMed、Google Scholar等数据库，提取文献中报道的RdRp活性数据。

## 面向ML驱动药物设计的实施清单
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

---

**数据可用性声明**：
- 本文使用的自动化筛选脚本位于`_pages/Drug Design/RdRp_private/`：
  - `filter_rdrp_csv.py`：仅SARS-CoV-2筛选（处理CHEMBL4523582，输出1个文件，101条）
  - `filter_rdrp_chembl_csv.py`：双病毒版本筛选（同时处理CHEMBL4523582 + CHEMBL5118，**自动生成2个文件**，含Virus_Type标注）
  - 配套说明文档：`README_filter.md`和`README_chembl_filter.md`
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

**数据更新记录**：
- 2026-01-09：更新ChEMBL筛选脚本，优化nsp14共转染测定处理，CHEMBL4523582筛选结果从95条增至101条SARS-CoV-2记录；新增双病毒筛选脚本`filter_rdrp_chembl_csv.py`，同时处理CHEMBL4523582和CHEMBL5118，支持迁移学习场景（总计128条：101条SARS-CoV-2 + 27条SARS-CoV-1）。
