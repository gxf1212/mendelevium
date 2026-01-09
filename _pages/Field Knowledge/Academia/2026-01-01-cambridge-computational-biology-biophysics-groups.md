---
title: "剑桥计算生物与生物物理团队全览"
date: "2026-01-01"
tags: [cambridge, computational-biology, biophysics, research-groups]
description: "剑桥大学在计算生物与生物物理方向的核心课题组、研究方法与代表性工作概览"
image: "/assets/img/thumbnail_mine/wh-6o2m9q.jpg"
thumbnail: "/assets/img/thumbnail_mine/wh-6o2m9q.jpg"
author: Xufan Gao
lang: zh-CN
---

# 剑桥计算生物与生物物理团队全览

## 概述
本文汇总了剑桥大学在计算生物与生物物理方向的代表性研究团队，涵盖Yusuf Hamied化学系、遗传学系、生物化学系与MRC分子生物学实验室等单位。每个团队条目都按照研究焦点、关键方法、近期成果与常用工具进行整理，并列出官网等进一步阅读渠道，方便快速对接潜在合作或深度调研。

## 研究团队一览

### Prof. Michele Vendruscolo（Yusuf Hamied化学系）
Vendruscolo课题组发展了多种结合实验约束的分子模拟方法，用于解析疾病相关蛋白的构象与相互作用。团队以全原子分子动力学为主，常常引入NMR约束、Markov模型（如Binding Paths框架）与统计推断来描绘折叠或误折叠通路，并借助增强采样手段计算小分子与无序蛋白的结合轨迹。近期工作集中在淀粉样形成机制以及帕金森病靶点的激酶抑制剂筛选，并在Nature Chemistry 2024年发表了关于MARK激酶的研究。官网：[Vendruscolo Lab](https://www-vendruscolo.ch.cam.ac.uk/)：https://www-vendruscolo.ch.cam.ac.uk/

### Prof. Andreas Bender（Yusuf Hamied化学系）
Bender领导的数据驱动药物发现团队专注于AI、机器学习与化学信息学在化学生物和药物设计中的应用。团队构建深度神经网络、规则模型与大数据分析流程来预测配体性质、毒性和安全窗口，常用数据集包括ToxCast等高通量生物活性库，并结合对接与分子生成技术指导结构优化。2025年发表于Nature Communications的研究展示了深度强化学习在设计高效A₂A受体配体方面的潜力；其他项目也覆盖CNN、Transformer等监督学习模型以及海量化学数据的可解释分析。团队官网：[Bender Group](https://bender.group.ch.cam.ac.uk/)：https://bender.group.ch.cam.ac.uk/

### Dr. Lucy Colwell（Yusuf Hamied化学系）
Colwell实验室以数据科学推动计算结构生物学，善于利用大规模同源序列中的协同演化信号来预测三维结构与功能。通过对同源序列聚类并嵌入AlphaFold2，团队在2023年发表于Nature Biotechnology的工作中展示了预测多种构象的策略。他们结合Potts模型、图神经网络与生成模型，学习残基相关性、蛋白–配体偏好以及酶底物特异性，并打造可解释的机器学习工具（如HMM-logo可视化）以服务蛋白工程。官网：[Colwell Lab](https://www.ch.cam.ac.uk/person/ljc37)：https://www.ch.cam.ac.uk/person/ljc37

### Prof. Jonathan M. Goodman（Yusuf Hamied化学系）
Goodman团队将计算化学、量子化学与AI结合，用于研究小分子的结构与反应性。他们将化学信息学同量子/分子动力学计算整合，发展面向结构验证的机器学习流程，例如把NMR与IR谱图联合输入模型来确认合成产物。团队也尝试语言模型预测反应结果与天然产物骨架，并以DFT、分子对接与定制化ML管线作为日常工具。官网：[Goodman Group](https://www-jmg.ch.cam.ac.uk/)：https://www-jmg.ch.cam.ac.uk/

### Dr. Aleks Reinhardt（Yusuf Hamied化学系）
Reinhardt带领的统计力学团队研究软物质与生物物质的相行为、自组装与凝聚现象。他们在原子级与粗粒度尺度上结合分子动力学、蒙特卡洛模拟，重点关注蛋白/RNA混合物的液–液相分离以及DNA纳米结构装配（J. Chem. Phys. 2023、Biophys. J. 2023）。团队也开发溶剂化与凝聚过程的模拟流程，并通过增强采样计算晶体生长自由能，常用工具包括GROMACS、聚合物与蛋白粗粒模型以及相图统计分析。官网：[Reinhardt Group](https://reinhardt.group.ch.cam.ac.uk/)：https://reinhardt.group.ch.cam.ac.uk/

### Prof. Jonathan Clarke（Yusuf Hamied化学系）
Clarke实验室通过分子动力学与实验结合，研究蛋白折叠动力学与力学性质。他们模拟蛋白结构域的受力解折过程，并与AFM实验匹配；还利用φ值分析约束的MD来解析折叠过渡态。关于titin/FNIII结构域的研究通过蛋白工程和受力MD共同描绘能量景观。团队惯用全原子MD（施加力或约束）、过渡路径采样与自由能计算来连接结构动力学与热力学，目前PI已接近退休。官网：[Clarke Group](https://jclarke.group.ch.cam.ac.uk/computational-studies-protein-folding)：https://jclarke.group.ch.cam.ac.uk/computational-studies-protein-folding

### Prof. Rosana Collepardo（遗传学系）
Collepardo的染色质建模团队构建DNA与染色质的多尺度模型，将粗粒度聚合物模型锚定在全原子分子动力学上，以研究核小体与蛋白如何塑造三维基因组结构。他们模拟染色质纤维来预测接触图，并评估连接组蛋白和转录因子对空间折叠的影响。常用方法包括核小体核心颗粒的全原子MD、千碱基尺度的介观蒙特卡洛模型以及链接序列到空间构象的理论框架。官网：[Collepardo Lab](https://www.gen.cam.ac.uk/research-groups/research-groups/collepardo-group)：https://www.gen.cam.ac.uk/research-groups/research-groups/collepardo-group

### Prof. Laura Itzhaki（药理学系）
Itzhaki实验室专注串联重复蛋白（ankyrin、HEAT、ARM等）的计算设计与功能研究。团队借助原位建模与蛋白工程，绘制重复结构域的折叠能量学并重新设计其结合功能；近期成果包括基于重复结构的抑制剂设计，以及研究内在无序链如何识别结构化重复域。方法涵盖Rosetta等结构建模软件、分子对接、重复框架的共识设计与折叠动力学模拟。官网：[Itzhaki Group](https://www.phar.cam.ac.uk/research/Itzhaki)：https://www.phar.cam.ac.uk/research/Itzhaki

### Prof. Florian Hollfelder（生物化学系）
Hollfelder实验室以实验手段研究酶机制与设计，利用定向进化与微流控来进化多底物酶并探究分子识别原则。某些项目结合X射线晶体学与动力学测试，解析进化后磺酸酯酶的底物结合方向，并与Kamerlin课题组合作开展MD以验证构象变化。团队常用技术包括高通量液滴筛选、突变体的晶体学/NMR以及自建或合作的对比MD模拟。官网：[Hollfelder Lab](https://hollfelder.bioc.cam.ac.uk/)：https://hollfelder.bioc.cam.ac.uk/

### Dr. Joe Greener（MRC分子生物学实验室）
Greener团队开发融合机器学习的分子动力学，训练图神经网络与可微分力场来提升生物大分子模拟精度，目标是让蛋白MD逼近量子化学准确度。他们在Chemical Science 2024年发表的工作展示了面向内在无序蛋白的可微MD力场优化方法，并编写基于Julia/PyTorch的GPU加速MD代码，把ML势能嵌入大规模模拟。官网：[Greener Group](https://www2.mrc-lmb.cam.ac.uk/groups/greener/)：https://www2.mrc-lmb.cam.ac.uk/groups/greener/

## 参考来源
1. [Professor Michele Vendruscolo | Yusuf Hamied Department of Chemistry](https://www.ch.cam.ac.uk/person/mv245)：https://www.ch.cam.ac.uk/person/mv245
2. [The Vendruscolo Laboratory](https://www-vendruscolo.ch.cam.ac.uk/)：https://www-vendruscolo.ch.cam.ac.uk/
3. [Professor Andreas Bender | Data-Driven Drug Discovery and Molecular Informatics](https://bender.group.ch.cam.ac.uk/person/ab454)：https://bender.group.ch.cam.ac.uk/person/ab454
4. [Index | Data-Driven Drug Discovery and Molecular Informatics](https://bender.group.ch.cam.ac.uk/)：https://bender.group.ch.cam.ac.uk/
5. 文中提到的Nat. Commun. 2025 A₂A受体研究，详见Bender团队论文记录。
6. [Dr Lucy Colwell | Yusuf Hamied Department of Chemistry](https://www.ch.cam.ac.uk/person/ljc37)：https://www.ch.cam.ac.uk/person/ljc37
7. [Professor Jonathan Goodman | Yusuf Hamied Department of Chemistry](https://www.ch.cam.ac.uk/person/jmg11)：https://www.ch.cam.ac.uk/person/jmg11
8. [The Goodman Group, Cambridge](https://www-jmg.ch.cam.ac.uk/)：https://www-jmg.ch.cam.ac.uk/
9. Goodman团队关于NMR/IR驱动的结构验证研究，详见其官网出版物。
10. [Dr Aleks Reinhardt | Yusuf Hamied Department of Chemistry](https://www.ch.cam.ac.uk/person/ar732)：https://www.ch.cam.ac.uk/person/ar732
11. [Index | The Reinhardt Group](https://reinhardt.group.ch.cam.ac.uk)：https://reinhardt.group.ch.cam.ac.uk
12. [Computational Studies of Protein Folding | The Clarke Group](https://jclarke.group.ch.cam.ac.uk/computational-studies-protein-folding)：https://jclarke.group.ch.cam.ac.uk/computational-studies-protein-folding
13. [Collepardo Group | Department of Genetics](https://www.gen.cam.ac.uk/research-groups/research-groups/collepardo-group)：https://www.gen.cam.ac.uk/research-groups/research-groups/collepardo-group
14. [Tandem-repeat proteins: Folding, function, role in disease and therapeutic intervention | Department of Pharmacology](https://www.phar.cam.ac.uk/research/Itzhaki)：https://www.phar.cam.ac.uk/research/Itzhaki
15. [Home | Hollfelder Group](https://hollfelder.bioc.cam.ac.uk/)：https://hollfelder.bioc.cam.ac.uk/
16. [Evolutionary repurposing of a promiscuous enzyme | Department of Biochemistry](https://www.bioc.cam.ac.uk/news/archive/2018/evolutionary-repurposing-of-a-promiscuous-enzyme)：https://www.bioc.cam.ac.uk/news/archive/2018/evolutionary-repurposing-of-a-promiscuous-enzyme
17. [Greener Group](https://www2.mrc-lmb.cam.ac.uk/groups/greener/)：https://www2.mrc-lmb.cam.ac.uk/groups/greener/
18. [Publications | Greener Group](https://www2.mrc-lmb.cam.ac.uk/groups/greener/publications/)：https://www2.mrc-lmb.cam.ac.uk/groups/greener/publications/
