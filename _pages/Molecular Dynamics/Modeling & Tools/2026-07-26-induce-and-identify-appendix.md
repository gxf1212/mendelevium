---
title: "MxMD与SiteMap识别隐藏口袋的技术细节与补充分析"
date: "2026-07-25"
last_modified_at: "2026-08-11"
tags: [cryptic-binding-sites, mixed-solvent-md, sitemap, induce-and-identify, drug-discovery, molecular-dynamics, schrodinger, allosteric-pocket]
description: "补充说明MxMD平衡和热点识别协议、SiteMap的逐帧去重流程，以及采样长度、纯水MD和低排序案例的对照结果。"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-rddgwm.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-rddgwm.jpg"
author: Xufan Gao
lang: zh-CN
---

# MxMD与SiteMap识别隐藏口袋的技术细节与补充分析

本页保留正文未展开的模拟参数、去重判据和补充对照。NVT指粒子数、体积和温度固定的系综；NPT指粒子数、压力和温度固定的系综。两者依次用于让溶剂和蛋白在进入生产模拟前稳定下来。

## 本文信息

- **标题**：Identifying Cryptic Binding Sites with Mixed Solvent MD Simulation and SiteMap
- **作者**：Da Shi, Dmitry Lupyan, Steven Jerome, Robert Abel, Yuqi Zhang
- **发表期刊**：Journal of Chemical Information and Modeling
- **发表时间**：2026年7月24日
- **单位**：Schrödinger Inc.（美国加利福尼亚州圣地亚哥、马萨诸塞州剑桥、纽约州纽约）
- **DOI**：https://doi.org/10.1021/acs.jcim.6c01288
- **引用格式**：Shi, D., Lupyan, D., Jerome, S., Abel, R., & Zhang, Y. (2026). Identifying Cryptic Binding Sites with Mixed Solvent MD Simulation and SiteMap. *Journal of Chemical Information and Modeling*. https://doi.org/10.1021/acs.jcim.6c01288
- **代码与数据**：Schrödinger Suite 2025-3（商业软件）：https://www.schrodinger.com/downloads/releases；
  - MxMD + SiteMap工作流的结果文件：https://drive.google.com/drive/folders/1eTCUM-u5QcoS0Mv6Cx_sS8PS0P4IWRha

## 研究方法

### “诱导”阶段：用MxMD采样apo蛋白构象

MxMD的核心思想是把**小分子探针**（Schrödinger默认三种：乙腈、异丙醇、嘧啶）混入水溶剂中跑MD。这三种探针的化学性质互补：

- **乙腈**：极性强，含氮腈基，是典型氢键受体探针
- **异丙醇**：带羟基，兼有氢键给体和疏水异丙基，亲疏水性兼备
- **嘧啶**：芳香杂环，专门探测能形成芳香和π-π相互作用的口袋

探针分子会探索蛋白表面的不同区域，在有利位置形成较高占用。Schrödinger采用两个溶剂化步骤和七个动力学步骤：前7步在10 K或300 K下完成布朗动力学和MD平衡，其中重原子约束持续到第7步；最后两步为无约束生产模拟：

1. 在蛋白外**加一层7 Å厚度的共溶剂层**（乙腈、异丙醇、嘧啶三种）
2. 再加水至**共溶剂与水的体积比为5%**——也就是共溶剂体积占整个体系体积的5%，其余是水和蛋白
3. 24 ps布朗动力学（NVT, 10 K，全原子约束）
4. 24 ps布朗动力学（NVT，10 K，重原子约束）
5. 12 ps MD（NVT，10 K，重原子约束）
6. 12 ps MD（NPT，10 K，重原子约束）
7. 24 ps MD（NPT，300 K，重原子约束）
8. 15 ns生产模拟（NPT，300 K，无约束）
9. 5 ns生产模拟（NPT，300 K，无约束，用于数据采集与分析）

每种探针跑10条重复轨迹，共30条。将每种共溶剂的轨迹末5 ns按**每10帧**采样一次，约得到1000帧结构，三种共溶剂合计约3000帧。

MxMD通过探针在空间中的访问频率筛选候选口袋：

- 把每条轨迹的探针原子坐标按**0.5 Å网格**分箱，得到原始占用计数——探针在哪些位置停留得多
- 把计数转换成**Z-score**（标准化分数），用单连接层次聚类（3 Å截距）把高占用网格点合并成**spot**——排除随机噪声，只保留真正有意义的探针聚集区
- **两个及以上不同探针的spot有重叠时定义为hotspot**，用所有重叠spot的Z-score之和作为**MxMD score**排序。这一条件降低了仅由单一探针产生的局部高占用信号
- 这一步是原MxMD工作流的位点识别方式；在该数据集上，平均每个apo结构识别出约15个hotspot

只依赖探针占用时，短暂开放或低占据的隐藏口袋可能无法形成热点。组合工作流改用SiteMap逐帧评价这些构象。

### “识别”阶段：用SiteMap对所有帧打分排序

MxMD的探针占用率只能反映探针实际到访的开放状态，短暂或低占据的口袋可能不形成热点。因此，工作流把MxMD产生的构象系综交给SiteMap逐帧识别和排序。

SiteMap是Schrödinger的一个基于物理的位点打分工具，沿着蛋白表面做网格扫描，对每个候选位点计算**SiteScore**——综合考虑位点体积、溶剂暴露度、疏水性、亲水性等多个物理特征：

- **体积**（volume）：位点能否容纳配体
- **溶剂暴露度**：位点的可及性
- **疏水性**：与配体疏水相互作用的潜力
- **亲水性**：与极性/带电配体相互作用的潜力

它有两种运行模式，分工如下：

- **site-evaluation模式**：需要把“配体大概在这里”告诉SiteMap，以holo配体位置为中心、6 Å半径的盒子限定扫描区域。这个模式只用于**基准评估**（表3那种apo vs holo对比），**不进工作流**
- **site-detection模式**：不预设任何配体位置，对蛋白表面**全范围扫描**，找遍所有潜在位点。这才是工作流真正用的模式

工作流按以下步骤处理：

1. **提取MxMD帧**：MxMD跑完后，把每个共溶剂的MD轨迹从最后5 ns按**步长10**（每隔10帧取一帧）采样，每种共溶剂约得到1000帧，三种共溶剂总共约**3000帧**蛋白结构
2. **逐帧跑SiteMap**：把全部3000帧都丢给SiteMap用**site-detection模式**打分，每帧平均报出约5个位点，总共得到约**15,000个候选位点**——这15,000就是从3,000个瞬时构象里“识别”出来的候选口袋
3. **用Sphere Exclusion合并去冗**：把所有15,000个位点按SiteScore从高到低排序，**从最高分开始遍历**，选中一个就把它放进最终列表，再把所有“残基组成高度相似”（Tanimoto相似度>0.4）的其他位点标记掉、跳过。这步保证最终列表里每个位点对应蛋白上一个独立的物理区域，不会被同一口袋在不同帧里反复刷出
4. **按SiteScore排序出最终列表**：去冗后剩10到170个位点（取决于蛋白大小），按SiteScore从高到低输出，就是后面表4里的“top-5命中率”——看holo配体是不是排进了前5

配体-位点对应关系不是工作流做的事，而是评估用的判据。**DCC判据**（distance to closest atom）计算位点中心到任一配体原子的最短距离，小于4 Å即视为“配体占据”；4 Å阈值沿用自P2Rank文章。位点周围的蛋白残基定义为距位点中心3.5 Å以内；肽和蛋白配体的位点-配体重叠只统计距受体3 Å以内的配体原子。

> 这一步把物理信号（SiteScore）和位置信号（DCC）分离开，让SiteMap专注于评估位点是否具有配体结合的物理条件，而不需要直接对应到已知配体上。

## 主要结果

### 工作流的进一步验证：采样长度、探针必要性、低排序案例

SI补充了三组对照：

1. **采样长度**：在22个位点子集上，将分析窗口从最后5 ns扩展到最后50 ns，总生产阶段为65 ns，使用的轨迹帧数增加10倍。整体排名反而下降，因为高SiteScore候选位点增加，假阳性增多。3KFL是例外，位点排名从64升至12。
2. **纯水MD对比**：去除探针后，以相同的10条重复模拟运行纯水MD。65个位点按静态SiteMap能否识别分为44个易识别位点和21个难识别位点。纯水MD在易识别位点略好，在难识别位点略差；这支持MxMD对难识别位点提供了更广的构象采样覆盖。
3. **低排序案例**：3KFL、3P53和1FXX/3HL8未进入top-10。下表列出各自的SiteScore、holo结构对照和原因：

| PDB | 工作流SiteScore | holo结构SiteScore | 排名 | 失败原因 |
| --- | --- | --- | --- | --- |
| 3KFL | 0.9 | 1.08 | 64 | MxMD可能未采到足够接近配体结合态的构象；排名第3的工作流位点也在配体附近，只是稍更暴露 |
| 3P53 | 1.0 | 1.07 | 36 | 同上；排名第1的工作流位点也在配体附近 |
| 1FXX/3HL8 | 0.96 | - | 55（默认）/ **9**（浅位点检测） | SSB蛋白与SSB相互作用酶之间的**蛋白-蛋白相互作用**（PPI）位点，默认SiteMap未针对这类浅位点优化；启用浅位点检测后SiteScore升至1.03 |

> 这些案例指向两类限制：MxMD未必采到接近配体结合态的构象，默认SiteMap也未必适合浅的蛋白-蛋白相互作用位点。
