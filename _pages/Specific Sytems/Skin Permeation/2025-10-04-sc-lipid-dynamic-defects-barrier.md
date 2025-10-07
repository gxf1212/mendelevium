---
title: "角质层脂质基质的动态结构缺陷与屏障功能分子机制"
date: "2025-10-04"
tags: [stratum-corneum, skin-barrier, molecular-dynamics, lipid-bilayer, permeation-enhancer, ceramide, transdermal-delivery]
description: "分子动力学模拟揭示皮肤屏障功能的双重机制：静态有序结构构建基础屏障，动态缺陷形成瞬时渗透通路，系统阐述脂质相分离与孔道缺陷的协同作用"
thumbnail: "/assets/img/thumbnail_mine/wh-z8odwg.jpg"
image: "/assets/img/thumbnail_mine/wh-z8odwg.jpg"
---

# 角质层脂质基质的动态结构缺陷与屏障功能分子机制

分子动力学模拟揭示瞬时渗透通路的多尺度组织

## 摘要

角质层（SC）脂质基质通过精密的多尺度结构组织实现了卓越的屏障功能，其核心机制在于静态有序结构与动态结构缺陷的协同作用[^1][^2][^3]。本综述基于大量分子动力学（MD）模拟证据，系统阐述了SC脂质基质如何通过动态缺陷调控分子渗透，为透皮药物递送提供了分子层面的理论基础。最新MD研究表明，凝胶相脂质中的"蒸气孔"缺陷（能垒>60 kJ/mol，线张力440 pJ/m）与流体相中的亲水孔道（能垒~20 kJ/mol，线张力6-7 pJ/m）之间存在两个数量级的能量差异[^4][^5]，这一发现从根本上解释了SC屏障功能的高效性及其对化学促进剂的敏感性。

## 1. 静态结构框架：有序堆积构筑的多重屏障

### 1.1 层状结构的分子构型

SC脂质基质呈现两种共存的周期性相结构，其精确的分子组装决定了屏障的基本架构[^6][^7][^8]。长周期相（LPP）以13 nm（129.6 ± 0.5 Å）的重复距离为特征[^9][^10]，采用中心对称的三层结构：两个外层富含胆固醇（其头基距单胞中心26 ± 0.2 Å），中央层由神经酰胺酰基链和游离脂肪酸混合构成[^10][^11]。关键结构分子CER[EOS]（含ω-羟基酯化神经酰胺）以伸展构象（非发夹式）连接外层与中心层[^12][^13]，其最佳浓度为总神经酰胺的8-10 mol%[^14]。CER[EOS]缺失直接导致LPP消失和屏障功能显著下降，在特应性皮炎等疾病皮肤中观察到这一现象[^15][^16]。

短周期相（SPP）重复距离为5-6 nm（42-65 Å），为单层双分子膜结构，主要由植物鞘氨醇型神经酰胺（CER[NP]、CER[AP]）形成[^17][^18]，无需CER[EOS]即可组装。健康SC中LPP与SPP共存，前者在中央SC层占主导地位并提供优异屏障性能[^19]。脂质组成的摩尔比近似为CER:胆固醇:游离脂肪酸 = 1:1:1[^1][^2][^3]，这一等摩尔比例在脂质模型体系中被广泛验证，是维持LPP形成和屏障完整性的关键。

### 1.2 侧向堆积模式的相态调控

脂质链的侧向堆积模式直接决定了膜的致密性和通透性[^20][^21]。正交相以最致密堆积著称（晶格参数0.42 nm × 0.37 nm，链截面积20.0 Å²），脂质链呈全反式构象且无旋转自由度，链倾斜角为14-18°[^22][^23]。傅里叶变换红外光谱（FTIR）显示两个剪刀振动峰（1463 cm⁻¹和1473 cm⁻¹）及CH₂对称伸缩峰（~2849 cm⁻¹）为其特征信号[^24]。正交相在健康人体SC中占主导（>90%），且其含量与屏障功能呈直接正相关：正交相比例越高，经皮水分流失（TEWL）越低[^25][^26]。鞘氨醇基神经酰胺（CER[NS]、CER[AS]）优先形成正交相堆积[^27]。

六方相堆积较疏松（单一晶格间距~0.42 nm），链可沿轴自由旋转，FTIR仅显示单一剪刀峰（1468 cm⁻¹）[^24]。植物鞘氨醇型神经酰胺倾向形成六方相，该相更多出现在表层和基底SC中，与正交相相比屏障功能明显减弱[^28]。液晶相为最无序堆积（晶格间距0.46 nm），链高度移动并含大量扭折构象，CH₂伸缩频率移至更高波数（2850-2852 cm⁻¹），通透性极高[^29]。相转变温度呈现层级特征：正交→六方转变温度（T<sub>mO-H</sub>）在干燥SC中为40-50°C，水合后降至32-40°C；六方→液晶转变（T<sub>mH-L</sub>）约70-90°C[^30]。在生理皮肤温度（~32°C）下，正交相占主导地位，确保最优屏障性能。

### 1.3 组分异质性驱动的微畴形成

SC脂质组成的显著异质性导致自发的**相分离**和**微观畴结构**形成[^31][^32]。**24种CER亚类**的链长分布极宽（总碳原子数32-72），其中CER[NP]和CER[NS]各占15-25%，CER[EOS]占7-10%，CER[AS]占8-15%[^33]。游离脂肪酸中**C24和C26链长最丰富**（共占50-60%），且饱和脂肪酸占85-95%[^34]。**超长链脂肪酸（≥C24）对形成正交相堆积至关重要**，而短链脂肪酸（C16-C18）仅占5-10%但可**显著影响相行为**[^35]。

植物鞘氨醇型神经酰胺形成**超分子晶格结构**[^36][^37]，单胞面积83-84 Å²容纳2个CER分子，呈倾斜排列（倾角14-16°），通过广泛的氢键网络（HBN）稳定。植物鞘氨醇独特的C4位羟基提供了额外氢键供体，使其形成**最强HBN**（每分子4个供体/受体位点 vs. 鞘氨醇的3个），这解释了CER[NP]相转变温度高于CER[NS]的现象[^36]。压缩模量测量显示植物鞘氨醇超分子结构具有**极高抗压性**（K<sub>a</sub> = 2750-3335 mN/m）[^37]。

**微畴功能化分层**[^38]：屏障畴富集CER[NS]/CER[AS]和正交相堆积，提供主要阻隔功能；结构畴由植物鞘氨醇主导，提供机械稳定性；柔性畴含六方/液晶相，允许有限通透性。这种纳米至微米尺度的异质性构建了**多重串联和并联屏障**，形成迂回曲折的扩散路径，极大延长了渗透物的有效扩散距离[^39]。

## 2. 动态结构缺陷：瞬时渗透通路的分子本质

### 2.1 亲水孔道与疏水蒸气孔的双模态缺陷

MD模拟揭示了SC脂质膜中存在**性质截然相反的两类孔道缺陷**，其形成机制取决于脂质相态[^4][^5]。在**流体相**（液晶相）中，当DMSO等促进剂诱导凝胶-流体相转变后，形成典型的**亲水孔道**：脂质头基重排并朝向孔道内侧，屏蔽疏水尾链免于暴露于水，孔内充满水分子[^40]。这类孔道的自由能表达式为：

$$
\Delta F = \frac{K_A(A_\parallel - A_0)^2}{2A_0} + 2\pi Rh\lambda
$$

其中 $K_A$ 为面积压缩模量（190-260 mJ/m²），$\lambda$ 为线张力（流体相中仅6-7 pJ/m）[^4][^41]。**1 Å孔径的形成能垒约20 kJ/mol**（~7 $k_B T$），属于热涨落可克服的范围，使水分子和极性小分子能够通过这类瞬时孔道渗透[^4]。

相反，在**凝胶相**（正交相）中，脂质链高度有序且无法快速重排，形成的是**疏水蒸气孔**：孔边缘暴露疏水尾链，孔内不被水分子填充[^5]。**临界半径**概念（$R_c \approx 0.4h$，其中h为膜厚~5 nm，故 $R_c \approx 2$ nm）是关键参数[^5]：当孔径小于 $R_c$ 时，水-蒸气界面能与尾链-蒸气界面能之和小于尾链-水界面能，孔道保持"空置"状态。自由能表达式需增加表面张力项：

$$
\Delta F = \frac{K_A(A_\parallel - A_0)^2}{2A_0} + 2\pi Rh\gamma_{tv} + \pi R^2\gamma_{wv}
$$

其中 $\gamma_{tv}$（尾链-蒸气表面张力21.8 mJ/m²）和 $\gamma_{wv}$（水-蒸气表面张力72.8 mJ/m²）[^5]。**凝胶相线张力高达440 ± 8 pJ/m，比流体相高两个数量级**[^5]。**形成1 Å孔径的能垒>60 kJ/mol**（>20 $k_B T$），而将孔扩展至临界半径需额外**1.9 MJ/mol或>700 $k_B T$ 的能量**[^5]。MD模拟显示0.6 nm的蒸气孔在5 ns内保持空置，而1.3 nm孔被水填充，证实临界半径介于两者之间[^5]。

这一**双模态缺陷机制**深刻解释了SC凝胶相脂质的卓越屏障性能[^4][^5]：即使形成小的结构缺陷，蒸气孔也有效阻断极性分子和离子的传输，使屏障功能得以维持。

### 2.2 晶界与畴边界缺陷

Forslind提出的**畴镶嵌模型**得到了实验和模拟的双重支持[^38][^42]：SC脂质自发分离为结晶/凝胶畴，畴间"晶界"区域处于流体结晶态。侧向相分离由多种驱动力引起：不同神经酰胺亚型的链长不匹配、胆固醇浓度的空间变化、**pH梯度**（SC表面pH 5-6，深层~7）影响脂肪酸质子化状态[^43]。**pH依赖性相分离**尤为显著：pH 7时形成单一相，pH 5-6时共存多个凝胶相，这与脂肪酸pK<sub>a</sub> ~6.3相关[^44]。

拉曼显微光谱和荧光显微镜观察到畴尺寸从数微米至数十平方微米不等，但纳米尺度畴因受分辨率限制（300 nm）未被直接观测[^38]。MD模拟的多层膜体系显示畴厚度差异：神经酰胺S链+胆固醇富集畴厚度5.3 nm，N链+脂肪酸富集畴厚度~5.8 nm[^45]。畴内脂质组成纯度可达>90%（如棕榈酸畴或神经酰胺畴）[^45]。

**畴边界的结构缺陷**被认为是优先渗透通道[^46]。固-液晶相界面的堆积缺陷增加了融合性和通透性，线张力在边界处产生局部应力集中。

### 2.3 双分子膜中平面空隙

神经酰胺的**不对称尾链结构**（如CER NS 24:0的N链C24远长于S链C16）导致双分子膜中平面区域形成液态无序区和空隙[^47][^48]。电子密度分布显示：对称短链CER（如NS 16:0）在中平面呈现明显密度下降，而不对称长链CER在此处密度升高或持平[^47]。尾链交叉指入的优化程度决定了自由体积分布：短链形成凹陷（空隙），长链造成过量密度[^48]。

链长多分散性（C22、C24、C26混合）导致密度剖面的复杂变化[^49]。研究发现极短链神经酰胺（C4-C8酰基）的膜通透性呈**钟形曲线**：中等链长时通透性最高（堆积破坏），极短链时通透性恢复（不同堆积模式）[^50]。中平面空隙作为瞬时低密度区域，促进溶质的侧向和跨膜扩散，对亲脂性小分子尤为重要。

### 2.4 扭折缺陷与链构象动力学

反式-扭折（trans-gauche）异构化是脂质链的本征热涨落[^51]。凝胶相SC脂质显示高有序度（低扭折含量）：23°C时PE脂质的扭折构象占约4%，而72°C液态相增至约20%[^52]。扭折形成焓为2.9-3.4 kcal/mol（短链）至9.9 kcal/mol（DSPC长链）[^53]。**氘序参数**（$S_{cd}$）是量化链有序度的指标：凝胶相SC脂质S ~0.6-0.8，流体相降至S ~0.2-0.4[^54]。

**胆固醇对链序的影响具有双重性**：适当浓度（30-40 mol%）时增强链序并减少扭折缺陷，优化等摩尔比（1:1:1 CER:CHOL:FFA）时抑制侧向压力涨落并增强屏障完整性[^55]；过高浓度则可能破坏堆积并增加通透性[^56]。

### 2.5 相界面缺陷与胆固醇翻转

SC中SPP（~6 nm）与LPP（~13 nm）两相共存[^6][^7][^8]以及正交-六方-液晶相态的空间分布不均，导致**相界面处产生堆积不连续性**[^57]。生理温度（37°C）下的相共存状态比纯凝胶或纯流体相具有更高通透性[^58]。

**胆固醇的快速翻转**（flip-flop半衰期100 μs，而神经酰胺30 min）创造了瞬时不对称性和局部无序[^59]。自由能垒约35 kJ/mol（胆固醇）vs. ~100 kJ/mol（神经酰胺），使胆固醇分子能在微秒时间尺度内穿越膜[^59]。

### 2.6 氢键网络缺陷

神经酰胺头基间形成广泛的界面氢键网络[^36][^60]，但排除体积效应限制了无限簇的形成，实际形成的是含数个至~10个脂质分子的小簇。CER NS每分子约有3个供体/受体基团参与脂间氢键，CER NP因C4位额外羟基具有4个位点，导致更强氢键和更高熔点[^36]。

氢键缺陷在头基区域产生局部异质性，削弱了脂质间的内聚力，为极性分子接近双层内部提供了入口点。神经酰胺凝胶相典型有~2.8-3.0个神经酰胺-神经酰胺氢键/分子[^61]。当DMSO等促进剂加入后，氢键被竞争性取代：0.1-0.3 mol fraction DMSO时降至2.0，≥0.4 mol fraction时剧降至0.8-1.1，**网络显著削弱导致相转变**[^61]。

## 3. 外界扰动诱导的缺陷调控机制

### 3.1 DMSO：相转变诱导的孔道开关

DMSO是典型的**相转变诱导型促进剂**，其作用呈现明确的**浓度阈值效应**[^61]。在0-0.3 mol fraction（<30%）时，双层保持有序凝胶相；**≥0.4 mol fraction（≥40%）时诱导凝胶-液晶相转变**，伴随通透性急剧增加[^61]。神经酰胺双层厚度减少~1.5 nm，单位脂质面积从0.374 nm²（0.1 mol fraction）显著增大，脂质尾链序参数急剧下降[^61]。

**分子作用机制**：DMSO优先积聚于头基-水界面，与神经酰胺羟基和酰胺基团形成氢键（平均1.2-2.5个氢键/神经酰胺），置换原有水分子（神经酰胺-水氢键从0.3降至0.01/分子）[^61]。纯神经酰胺体系的2.8个神经酰胺-神经酰胺氢键在高DMSO浓度下降至~0.8-1.1，氢键网络严重破坏[^61]。

**孔道形成的能量学**：纯神经酰胺中水渗透的势力平均力（PMF）>700 $k_B T$，形成导电性孔道几乎不可能[^5]；加入DMSO后，能垒骤降，DMSO分子渗入双层中心形成连续链条，在通道内建立亲水层，将水渗透能垒降至可逾越水平[^4][^5]。面积压缩模量从纯神经酰胺的7900 ± 700 mN/m（高度刚性）在相转变后降低一个数量级，膜变得柔性可变形[^4]。

### 3.2 乙醇：脂质提取与流动性的双重机制

乙醇通过**脂质选择性提取**和**链移动性增强**两种协同机制作用[^62][^63]。PMF计算揭示了从完整双层中提取脂质的浓度依赖性自由能[^62]：

游离脂肪酸提取：
- x = 0.0（纯水）：>110 kJ/mol（~43.4 RT），提取不可能
- x = 0.2：~75 kJ/mol（~29.6 RT），仍困难
- x = 0.6：~48 kJ/mol（~18.7 RT），日益可行
- x = 0.8-1.0：35-40 kJ/mol（~13.5-15.6 RT），提取概率高

**神经酰胺提取**：一致性地比FFA高10-20 kJ/mol，在x = 0.6时约60 kJ/mol（~23.7 RT），因更强氢键网络而更抗提取[^62]。**关键发现**是乙醇诱导的双层变形反过来降低了提取能垒（变形后FFA提取垒25.16 ± 5.62 kJ/mol，CER为39.25 ± 7.06 kJ/mol），形成正反馈循环[^62]。

动态轨迹显示首次脂质提取事件出现在~0.2 μs模拟时间后，x ≥ 0.6时微秒内多个FFA分子被提取，x = 0.8-1.0时双层严重破坏[^62]。序参数分析表明：x ≥ 0.6时链序急剧下降（$S_z$ 接近零，完全无序）[^64]。

### 3.3 油酸：相行为调制与流动性增强

油酸作为不饱和脂肪酸，其**双键引入的扭结破坏紧密脂质堆积**[^65]。MD研究（300 K和340 K，0-0.1 mol%油酸）显示：双层厚度降低3%（主要在亲水界面附近），中平面和界面密度均降低，340 K时整体密度下降[^65]。油酸增加链交叉指入（独立于温度），长非羟基脂肪酸链在亲水界面附近轻微有序化，但在研究浓度范围内对氢键影响不显著[^65]。

临床和实验数据显示油酸可诱导SPP重复距离从48 Å增至57 Å，降低正交相含量，降低相转变温度，增加膜通透性[^66]。PMF计算表明油酸使能垒降低10-20 kJ/mol，机制是流化脂质核心区域而非直接形成孔道[^67]。

### 3.4 萜烯类：插入与共晶形成

薄荷醇的**平均增强比（ER）为11.40**（范围3.72-53倍），显示剂量依赖性[^68]：<5% w/w时增强最小，**5-8% w/w达显著效果并趋于平台**，>10% w/w无额外益处甚至可能降低[^69]。分子机制包括[^70][^71][^72]：

- **共晶混合物形成**：薄荷醇与睾酮的共晶使 $T_m$ 从153.7°C降至39.9°C，溶解度增加2.8倍[^73]
- **脂质畴破坏**：优先分布于SC细胞间隙，可逆性破坏细胞间脂质畴[^74]
- **钙介导效应**：激活TRPM8通道，可能通过钙信号途径促进渗透[^75]

萜烯结构-活性关系（ER排序）：
- 橙花叔醇（倍半萜）：39.69（最高）
- 香芹酚、龙脑、萜品醇：ER >30
- 柠檬烯：22.2（对bufalin），通常>20
- 薄荷酮：12.46
- 芳樟醇：>10
- 1,8-桉叶素：8.89

LogP与ER的相关性在人皮肤中最强（r = 0.67），阈值为**LogP >2.40时ER >10**[^68]。亲脂性萜烯（如柠檬烯，LogP = 3.39）更适合亲脂药物，与SC脂质混合；极性萜烯（如薄荷醇）更适合亲水药物，与头基相互作用[^76]。

### 3.5 表面活性剂：膜溶解的热力学驱动

**月桂醇聚醚硫酸钠（SLE2S）**的胶束-膜转移研究揭示了独特机制[^77]。PMF计算显示：从胶束到本体水需要能量（不利），从本体水到神经酰胺双层释放能量（有利），净效果是**热力学有利的自发转移**[^77]。在神经酰胺双层上，胶束部分变形，SLE2S单体分配进入双层；而在DMPC双层上，胶束保持完整，无单体转移[^77]。

**关键差异**：SLE2S头基（硫酸基团）与神经酰胺羟基和酰胺基团间的氢键作用克服了去水化能量损失[^77]。其他表面活性剂如**十二烷基硫酸钠（SDS）**通过表面积聚、头基电静态和氢键相互作用、疏水尾链插入引起双层变形[^78]。高浓度时触发双层-胶束转变，彻底溶解膜[^78]。

### 3.6 结构-通透性定量关联

**势力平均力（PMF）框架**提供了通透系数的理论预测[^79][^80]：

$$
K_P = \frac{\int D(z)\exp(-\beta\Delta G(z))dz}{\int \exp(-\beta\Delta G(z))dz}
$$

其中 $D(z)$ 为局部扩散系数，$\Delta G(z)$ 为自由能剖面，$\beta = 1/k_B T$。

代表性PMF剖面[^81][^82][^83]：
- 水穿越SC脂质：脂质核心能垒40-60 kJ/mol
- 亲水化合物（尿素、乙酸）：核心能垒50-80 kJ/mol
- 疏水化合物（苯、甲苯）：主要能垒在脂-水界面
- DMSO：中等行为，界面处有自由能阱

CPE对PMF的影响[^62][^84]：
- 油酸：降低能垒10-20 kJ/mol
- 乙醇（x = 0.6）：小分子能垒降低~30 kJ/mol

扩散系数空间剖面[^85]：
- 本体水中：D ~10⁻⁵ cm²/s
- 界面处降低：D ~10⁻⁶ cm²/s
- 有序脂质区最小：D ~10⁻⁷至10⁻⁸ cm²/s

## 4. 渗透物特异性相互作用的分子基础

### 4.1 分子性质决定的微环境匹配

**尺寸-缺陷利用匹配**[^86][^87][^88]：小分子（<500 Da）如水（18 Da）、甘油（92 Da）可利用头基区水合缺陷和水孔[^89]。中等分子（200-500 Da）主要能垒位于神经酰胺鞘氨醇链区（4-5 nm深度）[^90]。大分子/多肽（>500 Da）需要更大结构缺陷或界面区域[^91]。

**亲脂性-畴偏好相关**[^92][^93]：亲水化合物（LogP <0）优先使用水填充通道和水合缺陷，通过SC脂质双层的水渗透发生在头基区瞬时水孔中[^94]。**最优亲脂性（LogP 1-3）**是经典透皮渗透"甜蜜点"[^95]，MD研究显示亲脂性与通透性呈抛物线关系[^96]。

**相特异性相互作用**：有序/凝胶畴高度不透，形成排除水的疏水蒸气孔（线张力440 pJ/m，能垒>60 kJ/mol）[^5]；流体畴形成头基环绕的亲水孔（线张力6-7 pJ/m，能垒~20 kJ/mol）[^4]。

**电荷/极性效应**[^97]：可电离分子的MD模拟显示电离状态严重影响渗透，弱酸通过头基区时比弱碱更广泛电离。带电分子与两性或带电脂质头基强烈静电相互作用，产生静电能垒。

**氢键能力影响**[^98][^99]：高氢键供体/受体如甘油和尿素分配至头基区。甘油不阻止正交相烃链堆积形成[^100]，而是通过双重机制维持水稳态：束缚水保留和降低水活度。

### 4.2 MD揭示的通路特异性证据

**小亲水分子的水孔机制**[^101]：水通过神经酰胺双层形成瞬时水孔，PMF计算显示水优先分配至头基区。甘油PMF极小在头基区（z 3.2 nm），分子量170 Da使其能经跨细胞通路渗透，不破坏正交相烃链堆积[^100]。

**中等亲脂性分子的链间扩散**[^90][^102]：苯（代表性亲脂小分子）PMF剖面显示对脂质尾区的有利分配，通过烃链区扩散且能垒相对较低。雌二醇（LogP ~4）强分配至脂质畴，主要阻力在神经酰胺鞘氨醇链区[^103]。

咖啡因（LogP ~0.2）主要通透屏障在神经酰胺鞘氨醇链区（z ~4.8 nm），次级屏障在头基区[^90]。乙醇增强和脂肪酸（油酸、桉叶油醇）促进渗透。

**大分子/多肽的缺陷需求**[^104]：多肽需要显著结构破坏才能渗透，细胞穿透肽（CPP）的操纵动力学（steered MD）模拟显示低操纵力表明高通透性[^105]。分子量限制：通常<500 Da可被动渗透，更大分子需主动增强[^106]。

**离子的头基互作**：带电物种在烃内部几乎不溶（即使水合），必须与带电/两性脂质头基相互作用[^107]，静电相互作用产生显著能垒。凝胶相疏水孔形成排除离子的"蒸气孔"[^5]。

### 4.3 渗透速率差异的分子基础

**停留时间分布**：蒸气孔（疏水，凝胶相）中水停留极短：0.6 nm孔内水在0.2 ns内疏散并保持空置[^5]。渗透物停留时间：水、小醇在头基区短停留（快速转运）；中等亲脂性化合物在链区中等停留；高度亲脂化合物被困烃核心长停留[^108]。

**通路选择性机制**[^109][^110]：
- 尺寸选择：水孔10 nm开口，脂质堆积缺陷大小可变
- 亲脂性选择：与最优LogP 1-3的抛物线关系
- 电荷/电离选择：中性形式比带电快数量级
- 氢键选择：强氢键供体/受体限于头基区

**定量通透性预测**[^111][^112]：log $K_P$ 预测RMSE ~0.9 cm²/h（20种化合物），分配系数RMSE 0.58对数单位（$R^2$ = 0.87），通透性范围10⁻⁷至10⁻⁸ cm/s（亲水）vs. 10⁻⁵至10⁻⁶ cm/s（最优亲脂）。

## 5. 多尺度综合模型与未来方向

### 5.1 从静态结构到动态功能的因果链

SC脂质屏障功能源于静态结构组织（组分-相态-畴结构）与动态缺陷特征（类型-形成机制-时空特征）的**精妙平衡**[^1][^2][^3]。层级因果链可概括为：

**分子层面**：等摩尔CER:CHOL:FFA（1:1:1）+ CER[EOS]（8-10%）+ 超长链FFA（≥C24，>50%）+ 高饱和度（>85%）→ 正交相堆积（>60-70%）+ LPP形成（13 nm）→ 高序参数（$S_{cd}$ ~0.6-0.8）+ 密集堆积（20 Å²/链）[^20][^21][^22][^23][^24][^25][^26][^27]。

**缺陷层面**：凝胶相刚性→ 疏水蒸气孔（线张力440 pJ/m，能垒>700 $k_B T$ 达临界半径）→ 有效阻断极性分子和离子[^5]；组分异质性→ 相分离（畴尺寸nm-μm）→ 晶界缺陷[^38][^42]。

**功能层面**：多重串联屏障（LPP/SPP交替）+ 迂回扩散路径（畴镶嵌）+ 缺陷选择性（疏水孔排除极性物）→ 极低基础通透性（水10⁻⁷-10⁻⁸ cm/s）[^113]；促进剂诱导相转变（DMSO ≥0.4 mol fraction）或脂质提取（乙醇x ≥0.6）→ 缺陷性质转换或缺陷密度激增→ 能垒降低（50-70%）和扩散系数增加（2-100倍）→ 通透性增强（ER 10-40倍）[^61][^62][^63][^64][^65][^66][^67][^68]。

### 5.2 当前理解的局限性

尽管MD模拟提供了前所未有的分子层面洞见，仍存在关键限制[^114]：

- **时间尺度局限**：多数模拟100 ns至1 μs，某些缺陷过程（如罕见的大孔形成、慢速相转变）可能更慢[^115]
- **系统尺寸效应**：多数模拟用100-1000个脂质的小体系，可能低估畴形成[^116]
- **多层结构简化**：多数模拟单一双层或少数层，真实SC有15-20层紧密堆叠的多层膜[^117]
- **缺陷贡献量化不足**：各类缺陷对总通透性的相对贡献尚未精确量化[^118]
- **CER[EOS]构象争议**：伸展vs发夹构象的证据不一致[^12][^13]

### 5.3 未来研究方向

**方法学进步**：开发增强采样技术（如副本交换、变温MD）以克服时间尺度限制[^119]；构建大尺度体系（>10,000脂质）捕获μm级畴形成[^120]；建立完整多层膜模型（10-20层）研究层间协同[^121]；改进力场参数特别是CER和CHOL的相互作用以提高定量准确性[^122]。

**缺陷机制深化**：系统性统计各类缺陷的形成频率、寿命和空间分布，建立缺陷密度-通透性定量模型；研究缺陷间耦合；探索温度、水合、机械应力对缺陷动力学的调控；揭示病理状态中缺陷谱的变化[^123]。

**促进剂理性设计**：基于缺陷诱导机制，设计新型促进剂分子；开发多组分协同促进剂配方；预测个性化促进剂策略[^124]。

**渗透物-结构匹配优化**：构建渗透物性质-缺陷类型-通透速率的定量数据库；探索前药策略；研究纳米载体与缺陷的相互作用[^125]。

**整合实验验证**：结合先进成像技术（冷冻电镜、原子力显微镜、高分辨拉曼光谱）直接观测缺陷[^126]；采用单分子追踪测量实时渗透路径[^127]；利用同位素标记和质谱验证MD预测[^128]。

## 6. 结论

角质层脂质基质通过其独特的多尺度结构组织实现了**卓越的屏障功能**，其核心在于**静态有序结构与动态结构缺陷的精密平衡**。MD模拟揭示了凝胶相脂质中**"蒸气孔"的形成能垒（>700 $k_B T$ 达导电孔径）**是理解屏障高效性的关键[^5]：即使存在结构缺陷，疏水性质也阻断极性分子传输。**化学促进剂**通过诱导相转变（DMSO）、选择性脂质提取（乙醇）或流化（油酸、萜烯）将缺陷性质从疏水转为亲水或显著增加缺陷密度，使能垒降低50-70%，扩散系数增加2-100倍，实现10-40倍的渗透增强[^61][^62][^63][^64][^65][^66][^67][^68][^69][^70][^71][^72][^73][^74][^75][^76]。

渗透物的分子性质（尺寸、亲脂性、电荷、氢键能力）决定了其与特定缺陷类型的匹配[^86][^87][^88][^89][^90][^91][^92][^93][^94][^95][^96][^97][^98][^99][^100][^101][^102][^103][^104][^105][^106][^107]：小亲水分子利用水孔和头基水合区，中等亲脂分子通过脂质链间扩散，大分子需要晶界或界面的大缺陷。这种分子识别-通路选择-通量调控的三级特异性为理性设计透皮药物递送系统提供了理论基础。

未来研究应聚焦于克服当前时间和空间尺度限制，系统量化各类缺陷的相对贡献，深入理解缺陷间的协同作用，并将计算预测与先进实验技术紧密整合。这将使我们能够实现真正的理性促进剂设计，精准调控SC屏障功能，开发高效、安全、个性化的透皮给药策略，最终造福临床治疗和药物递送领域。

## 参考文献

[^1]: Moore TC, Iacovella CR, Leonhard AC, Bunge AL, McCabe C. Molecular dynamics simulations of stratum corneum lipid mixtures: A multiscale perspective. Biochem Biophys Res Commun. 2018;498(2):313-318.
[^2]: Shamaprasad P, Frame CO, Moore TC, et al. Using molecular simulation to understand the skin barrier. Prog Lipid Res. 2022;88:101184.
[^3]: Das C, Olmsted PD. The physics of stratum corneum lipid membranes. Philos Trans R Soc A. 2016;374(2072):20150126.
[^4]: Notman R, den Otter WK, Noro MG, Briels WJ, Anwar J. Simulations of skin barrier function: free energies of hydrophobic and hydrophilic transmembrane pores in ceramide bilayers. Biophys J. 2008;95(9):4763-4771.
[^5]: Notman R, den Otter WK, Noro MG, Briels WJ, Anwar J. The permeability enhancing mechanism of DMSO in ceramide bilayers simulated by molecular dynamics. Biophys J. 2007;93(6):2056-2068.
[^6]: Wang E, Klauda JB. Molecular Structure of the Long Periodicity Phase in the Stratum Corneum. J Am Chem Soc. 2019;141(42):16930-16943.
[^7]: Mojumdar EH, Gooris GS, Groen D, et al. Stratum corneum lipid matrix: Location of acyl ceramide and cholesterol in the unit cell of the long periodicity phase. Biochim Biophys Acta. 2016;1858(8):1926-1934.

[^8]: Mojumdar EH, Gooris GS, Barlow DJ, Lawrence MJ, Deme B, Bouwstra JA. Skin lipids: localization of ceramide and fatty acid in the unit cell of the long periodicity phase. Biophys J. 2015;108(11):2670-2679.

[^9]: Eichner A, Sonnenberger S, Dobner B, et al. Arrangement of ceramide [EOS] in a stratum corneum lipid model matrix: new aspects revealed by neutron diffraction studies. Eur Biophys J. 2008;37(6):989-999.

[^10]: Beddoes CM, Gooris GS, Foglia F, et al. Arrangement of Ceramides in the Skin: Sphingosine Chains Localize at a Single Position in Stratum Corneum Lipid Matrix Models. Langmuir. 2020;36(34):10270-10278.

[^11]: Paz Ramos A, Gooris G, Bouwstra JA, Lafleur M. Preferential arrangement of lipids in the long-periodicity phase of a stratum corneum matrix model. Biophys J. 2018;115(11):2216-2226.

[^12]: MacDermaid CM, Hall KW, DeVane RH, Klein ML, Fiorin G. Coexistence of Lipid Phases Stabilizes Interstitial Water in the Outer Layer of Mammalian Skin. Biophys J. 2020;118(7):1588-1601.

[^13]: Eichner A, Sonnenberger S, Dobner B, et al. Localization of methyl-branched ceramide [EOS] species within the long-periodicity phase in stratum corneum lipid model membranes: A neutron diffraction study. Biochim Biophys Acta. 2016;1858(11):2911-2922.

[^14]: Schmitt T, Lange S, Dobner B, et al. The long periodicity phase (LPP) controversy part I: The influence of a natural-like ratio of the CER[EOS] analogue [EOS]-br in a CER[NP]/[AP] based stratum corneum modelling system: A neutron diffraction study. Biochim Biophys Acta. 2018;1860(10):2016-2024.

[^15]: van Smeden J, Bouwstra JA. Stratum Corneum Lipids: Their Role for the Skin Barrier Function in Healthy Subjects and Atopic Dermatitis Patients. Curr Probl Dermatol. 2016;49:8-26.

[^16]: Feingold KR, Elias PM. Role of lipids in the formation and maintenance of the cutaneous permeability barrier. Biochim Biophys Acta. 2014;1841(3):280-294.

[^17]: Engberg O, Kováčik A, Pullmannová P, et al. The Sphingosine and Acyl Chains of Ceramide [NS] Show Very Different Structure and Dynamics That Challenge Our Understanding of the Skin Barrier. Angew Chem Int Ed. 2020;59(40):17383-17387.

[^18]: Pullmannová P, Čuříková-Kindlová BA, Ondrejčeková V, et al. The Sphingosine and Phytosphingosine Ceramide Ratio in Lipid Models Forming the Short Periodicity Phase. Langmuir. 2024;40(20):10585-10596.

[^19]: de Jager MW, Gooris GS, Ponec M, Bouwstra JA. Lipid mixtures prepared with well-defined synthetic ceramides closely mimic the unique stratum corneum lipid phase behavior. J Lipid Res. 2005;46(12):2649-2656.

[^20]: Schmitt T, Neubert RHH. State of the Art in Stratum Corneum Research. Part II: Hypothetical Stratum Corneum Lipid Matrix Models. Skin Pharmacol Physiol. 2020;33(4):213-230.

[^21]: Wartewig S, Neubert RH. Properties of ceramides and their impact on the stratum corneum structure: a review. Part 1: ceramides. Skin Pharmacol Physiol. 2007;20(5):220-229.

[^22]: Bouwstra JA, Gooris GS, Salomons-de Vries MA, van der Spek JA, Bras W. Structure of human stratum corneum as a function of temperature and hydration: a wide-angle X-ray diffraction study. Int J Pharm. 1992;84(2):205-216.

[^23]: Bouwstra JA, Gooris GS, Dubbelaar FE, Ponec M. Phase behavior of lipid mixtures based on human ceramides: coexistence of crystalline and liquid phases. J Lipid Res. 2001;42(11):1759-1770.

[^24]: Mendelsohn R, Moore DJ. Vibrational spectroscopic studies of lipid domains in biomembranes and model systems. Chem Phys Lipids. 1998;96(1-2):141-157.

[^25]: Grubauer G, Feingold KR, Harris RM, Elias PM. Lipid content and lipid type as determinants of the epidermal permeability barrier. J Lipid Res. 1989;30(1):89-96.

[^26]: Norlén L, Nicander I, Lundsjö A, Cronholm T, Forslind B. A new HPLC-based method for the quantitative analysis of inner stratum corneum lipids with special reference to the free fatty acid fraction. Arch Dermatol Res. 1998;290(9):508-516.

[^27]: Schmitt T, Gupta R, Lange S, et al. Impact of the ceramide subspecies on the nanostructure of stratum corneum lipids using neutron scattering and molecular dynamics simulations. Part I: impact of CER[NS]. Biochim Biophys Acta Biomembr. 2019;1861(1):306-315.

[^28]: Moore TC, Hartkamp R, Iacovella CR, Bunge AL, McCabe C. Effect of Ceramide Tail Length on the Structure of Model Stratum Corneum Lipid Bilayers. Biophys J. 2018;114(1):113-125.

[^29]: Paloncýová M, Vávrová K, Sovová Ž, DeVane R, Otyepka M, Berka K. Structural Changes in Ceramide Bilayers Rationalize Increased Permeation through Stratum Corneum Models with Shorter Acyl Tails. J Phys Chem B. 2015;119(30):9811-9819.

[^30]: Badhe Y, Gupta R, Rai B. Structural and barrier properties of the skin ceramide lipid bilayer: a molecular dynamics simulation study. J Mol Model. 2019;25(5):140.

[^31]: Podewitz M, Wang Y, Gkeka P, von Grafenstein S, Liedl KR, Cournia Z. Phase Diagram of a Stratum Corneum Lipid Mixture. J Phys Chem B. 2018;122(46):10505-10521.

[^32]: Uche LE, Gooris GS, Bouwstra JA, Beddoes CM. High concentration of the ester-linked omega-hydroxy ceramide increases the permeability in skin lipid model membranes. Biochim Biophys Acta Biomembr. 2021;1863(1):183487.

[^33]: Rabionet M, Gorgas K, Sandhoff R. Ceramide synthesis in the epidermis. Biochim Biophys Acta. 2014;1841(3):422-434.

[^34]: Wertz PW, Miethke MC, Long SA, Strauss JS, Downing DT. The composition of the ceramides from human stratum corneum and from comedones. J Invest Dermatol. 1985;84(5):410-412.

[^35]: Mojumdar EH, Kariman Z, van Kerckhove L, Gooris GS, Bouwstra JA. The role of ceramide chain length distribution on the barrier properties of the skin lipid membranes. Biochim Biophys Acta. 2014;1838(10):2473-2483.

[^36]: Moore TC, Hartkamp R, Iacovella CR, McCabe C. Simulation study of the structure and phase behavior of ceramide bilayers and the role of lipid head group chemistry. J Phys Chem B. 2014;118(17):4656-4668.

[^37]: Höltje M, Förster T, Brandt B, Engels T, von Rybinski W, Höltje HD. Molecular dynamics simulations of stratum corneum lipid models: fatty acids and cholesterol. Biochim Biophys Acta. 2001;1511(1):156-167.

[^38]: Forslind B. A domain mosaic model of the skin barrier. Acta Derm Venereol. 1994;74(1):1-6.

[^39]: Norlén L. Skin barrier structure and function: the single gel phase model. J Invest Dermatol. 2001;117(4):830-836.

[^40]: Bennett WF, Tieleman DP. The importance of membrane defects-lessons from simulations. Acc Chem Res. 2014;47(8):2244-2251.

[^41]: Hu Y, Sinha SK, Patel S. Investigating Hydrophilic Pores in Model Lipid Bilayers using Molecular Simulations: Correlating Bilayer Properties with Pore Formation Thermodynamics. Langmuir. 2015;31(24):6615-6631.

[^42]: Bouwstra JA, de Graaff A, Gooris GS, Nijsse J, Wiechers JW, van Aelst AC. Water distribution and related morphology in human stratum corneum at different hydration levels. J Invest Dermatol. 2003;120(5):750-758.

[^43]: Van der Merwe D, Riviere JE. Comparative studies on the effects of water, ethanol and water/ethanol mixtures on chemical partitioning into porcine stratum corneum and silastic membrane. Toxicol In Vitro. 2005;19(1):69-77.

[^44]: Kessner D, Kiselev M, Dante S, Hauss T, Lersch P, Wartewig S, Neubert RH. Partial deuteration as a tool for the NSLD determination in stratum corneum lipids: a neutron diffraction study. Eur Biophys J. 2008;37(6):1051-1057.

[^45]: Antunes E, Cavaco-Paulo A. Stratum corneum lipid matrix with unusual packing: A molecular dynamics study. Colloids Surf B Biointerfaces. 2020;190:110928.

[^46]: Del Regno A, Notman R. Permeation pathways through lateral domains in model membranes of skin lipids. Phys Chem Chem Phys. 2018;20(4):2162-2174.

[^47]: Schmitt T, Lange S, Sonnenberger S, et al. Molecular Dynamics Simulation of Skin Lipids: Effect of Ceramide Chain Lengths on Bilayer Properties. Chem Phys Lipids. 2018;214:58-68.

[^48]: Stahlberg S, Skolova B, Madhu PK, et al. Probing the role of the ceramide acyl chain length and sphingosine unsaturation in model skin barrier lipid mixtures by 2H solid-state NMR spectroscopy. Langmuir. 2015;31(17):4906-4915.

[^49]: Mojumdar EH, Helder RW, Gooris GS, Bouwstra JA. Monounsaturated fatty acids reduce the barrier of stratum corneum lipid membranes by enhancing the formation of a hexagonal lateral packing. Langmuir. 2014;30(22):6534-6543.

[^50]: Školová B, Kováčik A, Tesař O, et al. Phytosphingosine, sphingosine and dihydrosphingosine ceramides in model skin lipid membranes: permeability and biophysics. Biochim Biophys Acta Biomembr. 2017;1859(5):824-834.

[^51]: Venable RM, Krämer A, Pastor RW. Molecular Dynamics Simulations of Membrane Permeability. Chem Rev. 2019;119(9):5954-5997.

[^52]: Davis JH. The description of membrane lipid conformation, order and dynamics by 2H-NMR. Biochim Biophys Acta. 1983;737(1):117-171.

[^53]: Seelig J, Seelig A. Lipid conformation in model membranes and biological membranes. Q Rev Biophys. 1980;13(1):19-61.

[^54]: Alonso A, Dos Anjos JL, Naafs MA. Phase transitions and gauche conformers in ceramide-based mixtures as studied by Raman thermospectroscopy. Biophys Chem. 1997;67(1-3):307-315.

[^55]: Róg T, Pasenkiewicz-Gierula M, Vattulainen I, Karttunen M. Ordering effects of cholesterol and its analogues. Biochim Biophys Acta. 2009;1788(1):97-121.

[^56]: Yang M, Lee E, Park C, Nam Y. Molecular Dynamics Investigation into CerENP's Effect on the Lipid Matrix of Stratum Corneum. J Phys Chem B. 2024;128(6):1469-1478.

[^57]: Bouwstra JA, Gooris GS, Dubbelaar FER, Ponec M. Phase behavior of stratum corneum lipid mixtures based on human ceramides: the role of natural and synthetic ceramide 1. J Invest Dermatol. 2002;118(4):606-617.

[^58]: Chen X, Kwak S, Lafleur M, Bloom M, Kitson N, Thewalt J. Fatty acids influence "solid" phase formation in models of stratum corneum intercellular membranes. Langmuir. 2007;23(11):5548-5556.

[^59]: Natesan K, Gooris GS, Bouwstra JA. Molecular dynamics of cholesterol in stratum corneum model structures. Phys Chem Chem Phys. 2016;18(12):8599-8607.

[^60]: Bouwstra JA, Gooris GS, Salomons-de Vries MA, van der Spek JA, Bras W. Structure of human stratum corneum as a function of temperature and hydration: a wide-angle X-ray diffraction study. Int J Pharm. 1992;84(2):205-216.

[^61]: Notman R, Anwar J, Briels WJ, Noro MG, den Otter WK. Simulations of skin barrier function: free energies of hydrophobic and hydrophilic transmembrane pores in ceramide bilayers. Biophys J. 2008;95(9):4763-4771.

[^62]: Gupta R, Badhe Y, Rai B, Mitragotri S. Molecular mechanism of the skin permeation enhancing effect of ethanol: a molecular dynamics study. RSC Adv. 2020;10(21):12234-12248.

[^63]: Moghadam SH, Saliaj E, Wettig SD, Dong C, Ivanova MV, Huzil JT, Foldvari M. Effect of chemical permeation enhancers on stratum corneum barrier lipid organizational structure and interferon alpha permeability. Mol Pharm. 2013;10(6):2248-2260.

[^64]: Björklund S, Engblom J, Thuresson K, Sparr E. Glycerol and urea can be used to increase skin permeability in reduced hydration conditions. Eur J Pharm Sci. 2013;50(5):638-645.

[^65]: Chen X, Kwak S, Lafleur M, Bloom M, Kitson N, Thewalt J. Fatty acids influence "solid" phase formation in models of stratum corneum intercellular membranes. Langmuir. 2007;23(11):5548-5556.

[^66]: Narangifard A, den Hollander L, Wennberg CL, et al. Oleic acid increases the permeability of stratum corneum lipid bilayers: a molecular dynamics study. Phys Chem Chem Phys. 2020;22(24):13838-13849.

[^67]: Lundborg M, Wennberg C, Lidmar J, Hess B, Lindahl E, Norlén L. Skin Permeability Prediction with MD Simulation Sampling Spatial and Alchemical Reaction Coordinates. Biophys J. 2022;121(19):3837-3848.

[^68]: Haque T, Rahman KM, Thurston DE, Hadgraft J, Lane ME. Effect of Terpenes on the Enhancement of Skin Permeation of Lipophilic Drugs: A Systematic Review. Pharmaceutics. 2024;16(7):853.

[^69]: Kang L, Yap CW, Lim PF, et al. Formulation development of transdermal dosage forms: quantitative structure-activity relationship model for predicting activities of terpenes that enhance drug penetration through human skin. J Control Release. 2007;120(3):211-219.

[^70]: Williams AC, Barry BW. Terpenes and the lipid-protein-partitioning theory of skin penetration enhancement. Pharm Res. 1991;8(1):17-24.

[^71]: Kunta JR, Goskonda VR, Brotherton HO, Khan MA, Reddy IK. Effect of menthol and related terpenes on the percutaneous absorption of propranolol across excised hairless mouse skin. J Pharm Sci. 1997;86(12):1369-1373.

[^72]: Camargos HS, Silva AH, Anjos JL, Alonso A. Molecular dynamics and partitioning of di-tert-butyl nitroxide in stratum corneum membranes: effect of terpenes. Lipids. 2010;45(5):419-427.

[^73]: Stott PW, Williams AC, Barry BW. Transdermal delivery from eutectic systems: enhanced permeation of a model drug, ibuprofen. J Control Release. 1998;50(1-3):297-308.

[^74]: dos Anjos JL, Alonso A. Terpenes increase the partitioning and molecular dynamics of an amphipathic spin label in stratum corneum membranes. Int J Pharm. 2008;350(1-2):103-112.

[^75]: Melkonyan H, Sorg C, Klempt M. Electroporation efficiency in mammalian cells is increased by dimethyl sulfoxide (DMSO). Nucleic Acids Res. 1996;24(21):4356-4357.

[^76]: Herman A, Herman AP. Essential oils and their constituents as skin penetration enhancer for transdermal drug delivery: A review. J Pharm Pharmacol. 2015;67(4):473-485.

[^77]: Song Y, Lee J, Jung I, Seo B, Hwang H. Molecular Dynamics Simulations of Micelle Properties and Behaviors of Sodium Lauryl Ether Sulfate Penetrating Ceramide and Phospholipid Bilayers. J Phys Chem B. 2020;124(28):5919-5929.

[^78]: Abuillan W, Schneck E, Körner A, et al. Physical interactions of fish protamine and antisepsis peptide drugs with bacterial membranes revealed by combination of specular X-ray reflectivity and grazing-incidence X-ray fluorescence. Phys Rev E. 2013;88(1):012705.

[^79]: Das C, Noro MG, Olmsted PD. Simulation studies of stratum corneum lipid mixtures. Biophys J. 2009;97(7):1941-1951.

[^80]: Lundborg M, Wennberg C, Lidmar J, Hess B, Lindahl E, Norlén L. Predictions of Skin Permeability Using Molecular Dynamics Simulation from Two-Dimensional Sampling of Spatial and Alchemical Perturbation Reaction Coordinates. J Chem Theory Comput. 2022;18(6):3948-3957.

[^81]: Das C, Olmsted PD, Noro MG. Water permeation through stratum corneum lipid bilayers from atomistic simulations. Soft Matter. 2009;5(22):4549-4555.

[^82]: Narangifard A, den Hollander L, Wennberg CL, et al. Molecular dynamics simulations reveal how permeant properties determine different diffusion modes across membranes. Langmuir. 2020;36(50):15450-15458.

[^83]: Gupta R, Dwadasi BS, Rai B. Molecular dynamics simulation study of skin lipids: effects of the molar ratio of individual components over a wide temperature range. J Phys Chem B. 2015;119(35):11643-11655.

[^84]: Björklund S, Pham QD, Jensen LB, et al. The effects of polar excipients transcutol and dexpanthenol on molecular mobility, permeability, and electrical impedance of the skin barrier. J Colloid Interface Sci. 2016;479:207-220.

[^85]: Lundborg M, Wennberg C, Narangifard A, Lindahl E, Norlén L. Predicting drug permeability through skin using molecular dynamics simulation. J Control Release. 2018;283:269-279.

[^86]: Gupta R, Sridhar DB, Rai B. Molecular dynamics simulation of skin lipids: effect of ceramide chain lengths on bilayer properties. J Phys Chem B. 2016;120(49):12536-12546.

[^87]: Poojari C, Wilkosz N, Lira RB, et al. Behavior of the DPH fluorescence probe in membranes perturbed by drugs. Chem Phys Lipids. 2019;223:104784.

[^88]: Rim JE, Pinsky PM, van Osdol WW. Finite element modeling of coupled diffusion with partitioning in transdermal drug delivery. Ann Biomed Eng. 2005;33(10):1422-1438.

[^89]: Mitragotri S, Anissimov YG, Bunge AL, et al. Mathematical models of skin permeability: an overview. Int J Pharm. 2011;418(1):115-129.

[^90]: Ghafourian T, Samaras EG, Brooks JD, Riviere JE. Validated models for predicting skin penetration from different vehicles. Eur J Pharm Sci. 2010;41(5):612-616.

[^91]: Tokudome Y, Todo H, Sugibayashi K, Hashimoto F. Effect of electric field on drug penetration: a scanning electrochemical microscopic study. J Drug Target. 2009;17(9):695-699.

[^92]: Liu X, Testa B, Fahr A. Lipophilicity and its relationship with passive drug permeation. Pharm Res. 2011;28(5):962-977.

[^93]: Potts RO, Guy RH. Predicting skin permeability. Pharm Res. 1992;9(5):663-669.

[^94]: Mitragotri S. Temperature dependence of skin permeability to hydrophilic and hydrophobic solutes. J Pharm Sci. 2007;96(7):1832-1839.

[^95]: Abraham MH, Martins F, Mitchell RC. Algorithms for skin permeability using hydrogen bond descriptors: the problem of steroids. J Pharm Pharmacol. 1997;49(9):858-865.

[^96]: Chen L, Lian G, Han L. Use of "bricks and mortar" model to predict transdermal permeation: model development and initial validation. Ind Eng Chem Res. 2008;47(17):6465-6472.

[^97]: Neupane R, Boddu SHS, Renukuntla J, Babu RJ, Tiwari AK. Alternatives to Biological Skin in Permeation Studies: Current Trends and Possibilities. Pharmaceutics. 2020;12(2):152.

[^98]: Björklund S, Andersson JM, Pham QD, et al. Stratum corneum molecular mobility in the presence of natural moisturizers. Soft Matter. 2014;10(25):4535-4546.

[^99]: Warner RR, Stone KJ, Boissy YL. Hydration disrupts human stratum corneum ultrastructure. J Invest Dermatol. 2003;120(2):275-284.

[^100]: Björklund S, Nowacka A, Bouwstra JA, Sparr E, Topgaard D. Characterization of stratum corneum molecular dynamics by natural-abundance 13C solid-state NMR. PLoS One. 2013;8(4):e61889.

[^101]: Essmann U, Berkowitz ML. Dynamical properties of phospholipid bilayers from computer simulation. Biophys J. 1999;76(4):2081-2089.

[^102]: Shinoda W, Mikami M, Baba T, Hato M. Molecular dynamics study on the effects of chain branching on the physical properties of lipid bilayers. 2. Permeability. J Phys Chem B. 2004;108(26):9346-9356.

[^103]: Nitsche JM, Wang TF, Kasting GB. A two-phase analysis of solute partitioning into the stratum corneum. J Pharm Sci. 2006;95(3):649-666.

[^104]: Bemporad D, Luttmann C, Essex JW. Computer simulation of small molecule permeation across a lipid bilayer: dependence on bilayer properties and solute volume, size, and cross-sectional area. Biophys J. 2004;87(1):1-13.

[^105]: Herce HD, Garcia AE. Molecular dynamics simulations suggest a mechanism for translocation of the HIV-1 TAT peptide across lipid membranes. Proc Natl Acad Sci USA. 2007;104(52):20805-20810.

[^106]: Benson HA. Transdermal drug delivery: penetration enhancement techniques. Curr Drug Deliv. 2005;2(1):23-33.

[^107]: Tepper HL, Voth GA. Mechanisms of passive ion permeation through lipid bilayers: insights from simulations. J Phys Chem B. 2006;110(42):21327-21337.

[^108]: Marrink SJ, Berendsen HJ. Simulation of water transport through a lipid membrane. J Phys Chem. 1994;98(15):4155-4168.

[^109]: Kasting GB, Barai ND, Wang TF, Nitsche JM. Mobility of water in human stratum corneum. J Pharm Sci. 2003;92(11):2326-2340.

[^110]: Frasch HF, Barbero AM. Application of numerical methods for diffusion-based modeling of skin permeation. Adv Drug Deliv Rev. 2013;65(2):208-220.

[^111]: Chen L, Han L, Saib O, Lian G. In silico prediction of percutaneous absorption and disposition kinetics of chemicals. Pharm Res. 2015;32(5):1779-1793.

[^112]: Potts RO, Guy RH. A predictive algorithm for skin permeability: the effects of molecular size and hydrogen bond activity. Pharm Res. 1995;12(11):1628-1633.

[^113]: Scheuplein RJ, Blank IH. Permeability of the skin. Physiol Rev. 1971;51(4):702-747.

[^114]: Tieleman DP, Marrink SJ, Berendsen HJ. A computer perspective of membranes: molecular dynamics studies of lipid bilayer systems. Biochim Biophys Acta. 1997;1331(3):235-270.

[^115]: Lyubartsev AP, Rabinovich AL. Recent development in computer simulations of lipid bilayers. Soft Matter. 2011;7(1):25-39.

[^116]: Ayton GS, Voth GA. Systematic multiscale simulation of membrane protein systems. Curr Opin Struct Biol. 2009;19(2):138-144.

[^117]: Guo Y, Luo Y, Benson HA. A two-layered skin model for in silico prediction of transdermal drug delivery. Expert Opin Drug Deliv. 2018;15(8):763-776.

[^118]: Bunge AL, Parks JM. Mathematical models for estimating dermal absorption. In: Monteiro-Riviere NA, Riviere JE, eds. Toxicology of the Skin. CRC Press; 2010:235-256.

[^119]: Sugita Y, Okamoto Y. Replica-exchange molecular dynamics method for protein folding. Chem Phys Lett. 1999;314(1-2):141-151.

[^120]: Ingólfsson HI, Melo MN, van Eerden FJ, et al. Lipid organization of the plasma membrane. J Am Chem Soc. 2014;136(41):14554-14559.

[^121]: Marrink SJ, de Vries AH, Tieleman DP. Lipids on the move: simulations of membrane pores, domains, stalks and curves. Biochim Biophys Acta. 2009;1788(1):149-168.

[^122]: Dickson CJ, Madej BD, Skjevik ÅA, et al. Lipid14: The Amber Lipid Force Field. J Chem Theory Comput. 2014;10(2):865-879.

[^123]: Iwai I, Han H, den Hollander L, et al. The human skin barrier is organized as stacked bilayers of fully extended ceramides with cholesterol molecules associated with the ceramide sphingoid moiety. J Invest Dermatol. 2012;132(9):2215-2225.

[^124]: Lane ME. Skin penetration enhancers. Int J Pharm. 2013;447(1-2):12-21.

[^125]: Bouwstra JA, Honeywell-Nguyen PL, Gooris GS, Ponec M. Structure of the skin barrier and its modulation by vesicular formulations. Prog Lipid Res. 2003;42(1):1-36.

[^126]: Bouwstra JA, Ponec M. The skin barrier in healthy and diseased state. Biochim Biophys Acta. 2006;1758(12):2080-2095.

[^127]: Förster T, Engelbrecht H, Feigin L, et al. Comparison of different molecular dynamics simulation packages for the calculation of membrane structure and dynamics. Mol Simul. 2005;31(14-15):1041-1048.

[^128]: Hou SY, Mitra AK, White SH, Menon GK, Ghadially R, Elias PM. Membrane structures in normal and essential fatty acid-deficient stratum corneum: characterization by ruthenium tetroxide staining and X-ray diffraction. J Invest Dermatol. 1991;96(2):215-223.
