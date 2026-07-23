---

title: "CDApol极化模型方法论详解：EEM动态电荷平衡的原理与实现"
date: "2026-05-06"
last_modified_at: "2026-05-06"
tags: [metal-ions, force-field, molecular-dynamics, polarization, cationic-dummy-atom, EEM, AMBER, high-valent-ions, methods]
description: "详解CDApol模型的电负性均衡方法（EEM）：动态电荷平衡的能量函数、约束求解、物理图像，以及双步参数化与热力学积分的技术细节"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-g83qpe.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-g83qpe.jpg"
author: Xufan Gao
lang: zh-CN
---
# CDApol极化模型方法详解：EEM动态电荷平衡的原理与实现

对应正文见[固定电荷模型为何难以模拟高价金属离子？关键在于引入动态极化效应](./2026-05-06-polarizable-cationic-dummy-model)。

## 本文信息

- **标题**：A Polarizable Cationic Dummy Metal Ion Model
- **作者**：Ali Rahnamoun, Kurt A. O'Hearn, Mehmet Cagri Kaymak, Zhen Li, Kenneth M. Merz, Jr., Hasan Metin Aktulga
- **发表期刊**：*The Journal of Physical Chemistry Letters*
- **发表时间**：2022年6月8日
- **DOI**：https://doi.org/10.1021/acs.jpclett.2c01279
- **单位**：Michigan State University, USA
- **引用格式**：Rahnamoun, A.; O'Hearn, K. A.; Kaymak, M. C.; Li, Z.; Merz, K. M., Jr.; Aktulga, H. M. (2022). A Polarizable Cationic Dummy Metal Ion Model. *J. Phys. Chem. Lett.*, 13, 5334-5340.
- **相关框架**：Rahnamoun, A.; Kaymak, M. C.; Manathunga, M.; Götz, A. W.; Duin, A. C. T.; Merz, K. M., Jr.; Aktulga, H. M. (2020). ReaxFF/AMBER—A Framework for Hybrid Reactive/Nonreactive Force Field Molecular Dynamics Simulations. *J. Chem. Theory Comput.*, 16, 7645-7654. https://doi.org/10.1021/acs.jctc.0c00874

## 快速结论

- **EEM能量函数**由电负性线性项（驱动力）和硬度矩阵二次项（转移代价）构成，是理解动态电荷平衡的核心
- **总电荷约束**可通过增广线性方程组处理，每步MD只需求解核心区电荷平衡
- **核心7位点是唯一动态电荷未知量**：中心金属离子+6个虚拟原子的电荷每步重排，周围溶剂分子提供瞬时静电环境
- **外层固定电荷如何进入求解**：CDApol主文没有完整展开这套记号；本文采用ReaxFF/AMBER里的mEEM框架来辅助解释
- **两步参数化策略**：第一步训练EEM参数（$\chi_i,\eta_i,\gamma_{ij}$）复现DFT电荷分布，第二步扫描LJ参数（$\varepsilon,R_{\min}/2$）同时匹配实验HFE、IOD和CN

## 方法详解

### EEM能量函数的定义

CDApol模型通过**电负性均衡方法**（Electronegativity Equalization Method，EEM）引入动态极化。首先定义EEM能量函数：

$$
E_{\text{EEM}} = \sum_{i=1}^{N} \chi_i q_i + \dfrac{1}{2} \sum_{i=1}^{N} \sum_{j=1}^{N} q_i J_{ij} q_j
$$

其中$N$是CDApol核心位点数，即7个电荷位点（1个中心金属离子+6个虚拟原子），**不包括周围水分子**。公式中每个符号的含义：

- $q_i$：第$i$个位点的瞬时电荷（可正可负，单位是元电荷$e$）
- $\chi_i$：第$i$个位点的电负性参数（单位是能量，如eV）。在EEM里，$\chi_i$是通过拟合QM电荷分布得到的**可调参数**，不是Mulliken定义的实验量
- $J_{ij}$：位点$i$和$j$之间的相互作用矩阵元——对角项$J_{ii} = \eta_i$是Parr-Pearson硬度参数（防止电荷无限堆积），非对角项$J_{ij}$是带屏蔽的静电耦合（防止短程库仑爆炸）

CDApol文中用$J_{ij}$，ReaxFF/mEEM文中用$H_{ij}$，二者是同一类相互作用核的不同记号。在本文记号体系里，对角项$J_{ii} = H_{ii} = \eta_i$，非对角项$J_{ij} = H_{ij}$。

EEM能量函数也可写成矩阵形式：

$$
E_{\text{EEM}} = \chi^{\mathsf T} q + \dfrac{1}{2} q^{\mathsf T} H q
$$

> **一句话**：EEM不是给整盒水一起「调电荷」，而是只让核心7个位点在总电荷守恒下随环境重排。

### EEM能量函数的物理意义

EEM能量函数的两项分别对应**电荷流动的驱动力**和**电荷重分布的代价**：

#### 第一项：$\chi_i q_i$——电荷流动的驱动力

这一项决定**电荷想往哪里流**。虽然$\chi_i$在EEM中被称为Mulliken电负性参数，但它实际上是一个**可调的拟合参数**，只是借用了电负性的概念。传统的Mulliken电负性定义为 $\chi = \dfrac{I + A}{2}$，其中$I$是电离能，$A$是电子亲和能。

在化学中，电负性越大的原子（如氟、氧）越倾向于吸引电子。但在EEM模型里，$\chi_i$是通过拟合QM电荷分布得到的参数，**可以是正值也可以是负值**，其符号和大小决定了该位点在能量最小化时的电荷分配倾向。

能量项$\chi_i q_i$的物理含义：

- **$\chi_i$越小**：该位点越倾向于失去电荷（带正电）；**$\chi_i$越大**（更负）：越倾向于获得电荷（带负电）
- 如果$\chi_i$较小但仍为正，$q_i > 0$时$\chi_i q_i > 0$，**能量升高**——位点不想要电荷却还带正电，能量当然高；$\chi_i$较大而$q_i < 0$时则势能很低
- 系统会自动调整$q_i$，让总能量$E_{\text{EEM}}$最小——这就是电荷重新分配的**驱动力**

#### 第二项：$\dfrac{1}{2} q_i J_{ij} q_j$——电荷重分布的代价

这一项决定**电荷重分布要付出什么代价**。它包含两部分：

##### 对角项：$J_{ii} = \eta_i$（self energy代价）

对角项对应的是**单个位点上积累电荷的代价**。当$i=j$时，能量项变成：$\dfrac{1}{2} \eta_i q_i^2$。这里$\eta_i$是Parr-Pearson硬度参数，物理上定义为：

$$
\eta_i = \dfrac{I_i - A_i}{2}
$$

也就是电离能和电子亲和能的**差值的一半**。

> 能量项的物理含义：这是一个**二次项**，无论$q_i$是正是负，$q_i^2$总是正的，所以这一项**总是让能量升高**——防止电荷无限制地堆到某一个位点上。**$\eta_i$越大**，电荷积累的代价越高，位点越硬，极化响应越弱；**$\eta_i$越小**，位点越软，极化响应越强

##### 非对角项：$J_{ij}$（位点间相互作用）

非对角项对应的是**两个不同位点之间的静电相互作用**。在CDApol主文里，这部分只强调采用了**electrostatic shielding**来避免近距离的过强排斥；若按ReaxFF/mEEM的写法理解，非对角项对应的是一种**带屏蔽的库仑核**，其强度随位点间距离和屏蔽参数变化。

- 能量项$\dfrac{1}{2} q_i J_{ij} q_j$的物理含义：$q_i$和$q_j$**同号**时相互排斥（能量升高），**异号**时相互吸引（能量降低）。
- 位点越接近、屏蔽越弱，耦合作用通常越强。
- **$\gamma_{ij}$的物理意义**：
  - 如果没有屏蔽项，简单点电荷模型在短程会给出过强排斥
  - 引入屏蔽后，短程相互作用会被软化，用来近似真实电子云不是点电荷这一事实

> **总结**：非对角项$\dfrac{1}{2} q_i J_{ij} q_j$描述**位点间的静电耦合**。它让电荷分布不能随意变化，因为同号电荷会互相排斥，异号电荷会互相吸引。屏蔽参数则用来抑制相邻位点之间的非物理短程排斥。

### 总电荷约束与增广线性方程组求解

EEM真正求解的是一个**带约束的能量最小化问题**：

$$
\min_{\{q_i\}} E_{\text{EEM}}, \quad \sum_{i=1}^{N} q_i = Q_{\text{total}}
$$

在CDApol中，$Q_{\text{total}}$固定为金属离子的形式电荷（$\ce{Zn^{2+}}$的+2、$\ce{Al^{3+}}$的+3或$\ce{Zr^{4+}}$的+4）。电荷可以在中心离子和6个虚拟原子之间自由流动，但7个位点的电荷总和必须守恒。

先构造拉格朗日函数，把约束吸进来：

$$
\mathcal{L}(q_1,\ldots,q_N,\varepsilon) = \sum_i \chi_i q_i + \dfrac{1}{2}\sum_{i,j} q_i H_{ij} q_j + \varepsilon\left(\sum_i q_i - Q_{\text{total}}\right)
$$

对每个位点$i$求偏导并令其为零：

$$
\dfrac{\partial\mathcal{L}}{\partial q_i} = \chi_i + \sum_j H_{ij} q_j + \varepsilon = 0
$$

> 其中$\varepsilon$是拉格朗日乘子（注意这里$\varepsilon$前是负号，从$\varepsilon(\sum_i q_i - Q)$展开后得到$+\varepsilon$，移项后得$-\varepsilon$），它保证在最优解处强制满足总电荷约束——$\varepsilon$本身**不是电荷**，而是核心区平均电化学势的度量，反映系统在坚持$\sum q_i = Q_\text{total}$时付出的代价。

这给出$N$个标量方程，加上约束本身：

$$
\begin{cases}
\chi_i + \sum_{j=1}^N H_{ij} q_j + \varepsilon = 0 & (i=1,\ldots,N) \\
\sum_{j=1}^N q_j = Q_{\text{total}} & (\text{约束})
\end{cases}
$$

写成矩阵形式，就是增广线性方程组：

$$
\begin{bmatrix}
H & \mathbf{1} \\
\mathbf{1}^{\mathsf T} & 0
\end{bmatrix}
\begin{bmatrix}
q \\
\varepsilon
\end{bmatrix}
=
\begin{bmatrix}
-\chi \\
Q_{\text{total}}
\end{bmatrix}
$$

其中$\mathbf{1}$是全1列向量，最后一行对应总电荷约束$\mathbf{1}^{\mathsf T}q = Q_{\text{total}}$。这是一个$8 \times 8$的线性系统，核心7位点每步MD只需一次线性代数求解。其中系数矩阵中的非对角元为 $J_{ij} = F_{ij}$，为了避免极近距离下的库仑发散，SI中明确了其**静电屏蔽参数**（Electrostatic Shielding） $\gamma_{ij}$ 的公式：

$$
F_{ij} = \begin{cases} 
\dfrac{1}{\left( r_{ij}^3 + \gamma_{ij}^{-3} \right)^{1/3}} , & r_{ij} \le r_{\text{nonb}} \\
0, & \text{otherwise}
\end{cases}
$$

其中 $\gamma_{ij} = \sqrt{\gamma_i \cdot \gamma_j}$ 是一对元素相依赖的屏蔽项，确保 $r_{ij} \to 0$ 时静电势保持有限避免模型崩溃。

> **物理图像**：想象一个水池系统，7个水池通过管道连接，水可以在池子之间流动，但总水量不变。每个池子有自己的高度偏好（$\chi_i$）和容量限制（$\eta_i$），池子之间还有流动阻力（$J_{ij}$）。最终水会流到一个平衡状态，让整个系统的势能最低。

### 局部动态极化：外层固定电荷如何驱动核心区

理解EEM时，必须先把「参与方程」和「不作为未知量被优化」分开。CDApol的核心只有7个位点（中心金属离子+6个虚拟原子）是**动态电荷未知量**；周围的水分子和配体是**外层固定电荷**，参与方程但不是未知量。

外层固定电荷对核心区的作用，可以借用ReaxFF/AMBER框架（JCTC 2020）里的mEEM记号来理解。该框架将体系划分为核心区（core）和过渡区/MM区两部分。核心区的未知电荷记为$q_{\text{core}}$，外层固定电荷记为$q_{\text{trans}}$（常数向量，由力场给定，每步MD不重新优化）。

哪些外层原子进入$q_{\text{trans}}$？这由双层筛选机制决定：

- **第一层（分区筛选）**：先按体系划分确定候选身份——transition和MM区原子进入$q_{\text{trans}}$候选池，core区原子进入$q_{\text{core}}$。
- **第二层（距离筛选）**：在每一个MD步，只保留与core区发生有效非键耦合的外层原子，**即与core区原子距离在截断半径$r_{\text{cut}}$以内的那些候选原子**。

$$
\mathcal{S}_{\text{trans}}(t)=\left\{j\in(\text{transition}\cup\text{MM})\mid \exists i\in\text{core},\ r_{ij}(t)<r_{\text{cut}}\right\}
$$

因此，$q_{\text{trans}}$对应的是集合$\mathcal{S}_{\text{trans}}(t)$里这些原子的固定电荷向量。由于水分子和配体都在运动，$\mathcal{S}_{\text{trans}}(t)$会随时间变化，是一个**运行时集合**。

将电荷按core/trans分区后，增广线性方程组可以整理为只含核心区未知量的形式：

$$
\begin{bmatrix}
H_{\text{core}} & \mathbf{1}_c \\
\mathbf{1}_c^{\mathsf T} & 0
\end{bmatrix}
\begin{bmatrix}
q_{\text{core}} \\
\varepsilon
\end{bmatrix}
=
\begin{bmatrix}
-\chi_{\text{core}} \\
Q_{\text{total}}
\end{bmatrix}
-
\begin{bmatrix}
H_{\text{core-trans}} \\
\mathbf{0}^{\mathsf T}
\end{bmatrix}
q_{\text{trans}}
$$

其中$q_{\text{trans}}$是常数向量（AMBER固定电荷），不是新的动态电荷变量。右端第二项$H_{\text{core-trans}}q_{\text{trans}}$是外层固定电荷在核心区产生的静电驱动项，也可以等价写成「有效电负性」形式：

$$
\chi_{\text{core}}^{\text{eff}}=\chi_{\text{core}}+H_{\text{core-trans}}q_{\text{trans}}
$$

> **动态极化的来源**：外层水分子和配体不参与电荷优化，只提供瞬时静电场。随着它们的位置变化，$\chi_{\text{core}}^{\text{eff}}$实时波动，核心7位点重新分配电荷——这就是CDApol中「动态极化」的核心机制：电荷分布随局部构型响应，但7个核心位点的总电荷始终守恒。

## 参数化流程

### 两步串联的参数化策略

CDApol的参数化分两步：**第一步**训练极化力场参数（EEM + dummy骨架几何），**第二步**扫描LJ参数。具体分工：

| 步骤 | 训练目标 | 训练数据 | 参数状态 |
| --- | --- | --- | --- |
| **第一步**：极化力场参数训练 | EEM参数（$\chi_i, \eta_i, \gamma_{ij}$）和dummy骨架几何 | QM能量profile（1-7配位）+ QM电荷分布（1-6配位） | EEM参数和dummy几何参数从无到有；**不涉及任何HFE计算** |
| **第二步**：12-6 LJ参数扫描 | $\varepsilon$和$R_{\min}/2$ | 每个参数组合跑MD+TI，评估HFE、IOD、CN与实验值的偏差 | EEM参数锁定；LJ参数搜索；挑最优组合 |

两步**严格串联**：第一步完全独立于第二步，第一步产出的EEM参数一旦锁定，第二步只动LJ参数。如果同时优化所有参数，EEM的拟合目标（QM电荷）和LJ的拟合目标（实验热力学性质）会互相干扰；分步则各司其职。

> **分步的原因**：EEM的拟合目标是QM电荷分布，LJ的拟合目标是实验热力学性质（HFE/IOD/CN）。两者不在同一个目标空间里，如果同时优化，参数会打架——这也是为什么参数化必须分成两步走。

在每一步MD中，EEM参数固定，EEM通过增广线性方程组计算给定外部环境下的最优电荷；LJ参数则在MD和TI的总体框架中被优化。

![fig2](2026-05-06-polarizable-cationic-dummy-model_figs/fig2.png)

**图2：CDApol参数化管线**。第一步（左）以QM参考训练EEM和dummy几何，第二步（右）用热力学积分在LJ参数空间中搜索最优组合。

### EEM参数训练细节

第一步在指定构象下同时复现**QM能量**和**QM电荷**——电荷和能量一起训练，不是只训练电荷。具体做法：

- **DFT计算**：使用Gaussian 16，在B3LYP/6-311+g(d,p)水平上计算$\ce{Al^{3+}}$与1-7个水分子配位时的势能面，共7个构象。
- **能量基准（Figure S.1）**：图S.1展示了随配位数变化的QM能量曲线，横轴是配位数（1到7），纵轴是相对能量。八面体（6配位）构象能量最低，即全局能量极小点；欠配位或过配位时能量都会升高。

![figS1](cdapol_si_figs/figS1.jpeg)

**图S.1：$\ce{Al^{3+}}$ CDApol模型训练的QM能量曲线**。六配位（Octahedral）构象能量最低，与之偏离的欠配位或过配位构象能量均升高。图中同时标注了各构象的配位类型（Monohydrate至Heptahydrate）。

**电荷基准（Section S.2）**：对1-6配位的每个构象，提取DFT优化的原子电荷作为参考电荷分布。EEM参数（$\chi_i, \eta_i, \gamma_{ij}$）的作用就是让CDApol在给定构象下通过EEM求解得到的电荷分布与QM电荷尽量一致。误差函数同时覆盖能量和电荷两类数据：
$$
e_i = \left(\dfrac{x_{i,\mathrm{QM}} - x_{i,\mathrm{R}}}{w_i}\right)^2
$$

其中$x_{i,\mathrm{QM}}$和$x_{i,\mathrm{R}}$分别是QM参考值和当前ReaxFF计算值，$w_i$是权重参数。参数优化通过最小化该误差函数来完成：对每个训练构象，先**固定几何**（原子坐标取DFT优化后的结构），然后EEM在总电荷约束下求解出7个核心位点的最优电荷分布（与MD中每步的做法相同），再比较与QM电荷的偏差；同时也对整个构象的总能量与QM能量做比较。

权重$w_i$可以按需调节，让电荷项和能量项在总误差中的贡献比例可控。训练数据覆盖1-7配位的水合构象，使CDApol在欠配位（1-5配位）、八面体（6配位）和过配位（7配位）构象中都能复现QM结果，最终在MD模拟中得到正确的配位数。

### LJ参数扫描细节

第二步在$(\varepsilon, R_{\min}/2)$二维参数空间中进行网格搜索：

- **$\varepsilon$扫描范围**：1-3.4 kcal/mol，步长0.2 kcal/mol；**$R_{\min}/2$扫描范围**：0.6-1.0 Å，步长0.1 Å
- 每个$(\varepsilon, R_{\min}/2)$组合都要跑完整的MD+TI计算，评估HFE、IOD和CN三项性质

LJ势函数采用标准AMBER形式：

$$
V_{ij} = \varepsilon_{ij}\left[\left(\dfrac{R_{\min,ij}}{r_{ij}}\right)^{12} - 2\left(\dfrac{R_{\min,ij}}{r_{ij}}\right)^6\right]
$$

结合规则使用Lorentz-Berthelot混合规则，将金属中心的LJ参数与TIP3P水分子的氧原子参数混合，生成成对LJ势。MD模拟在20 Å × 20 Å × 20 Å的TIP3P水盒子中进行，共2736个水分子。0.25 fs是时间步长的保守选择；SI对$\ce{Zn^{2+}}$ CDApol模型测试了0.5 fs、1 fs、1.5 fs和2 fs，结果差异均很小，说明CDApol在较大时间步下仍然稳定：

| 时间步 | IOD (Å) | CN | HFE (kcal/mol) |
| --- | --- | --- | --- |
| 0.5 fs | 2.04 | 6.3 | -464.8 |
| 1.0 fs | 2.04 | 6.3 | -465.6 |
| 1.5 fs | 2.12 | 6.3 | -465.5 |
| 2.0 fs | 2.05 | 6.1 | -465.8 |

#### 热力学积分与三点高斯积分

第二步中每个参数组合的HFE通过**热力学积分**（Thermodynamic Integration，TI）计算。TI的核心思想是沿着一条连接初态和末态的路径，逐步「充电」或「去充电」，然后对路径上的能量导数积分，得到自由能差：

$$
\Delta G = \int_0^1 \left\langle \dfrac{\partial U(\lambda)}{\partial \lambda} \right\rangle_\lambda \mathrm{d}\lambda
$$

其中$\lambda$是耦合参数（$\lambda=0$对应初态，$\lambda=1$对应末态），$U(\lambda)$是$\lambda$状态下的势能，$\langle \cdots \rangle_\lambda$表示在$\lambda$状态下的系综平均。

积分无法解析求解，只能在离散的$\lambda$点上通过MD模拟采样$\langle \partial U/\partial\lambda\rangle_\lambda$，再用数值积分连起来。三点高斯积分（Three-point Gaussian Quadrature）通过精心选择积分点位置和权重，用较少采样点获得较高精度。对于三点高斯积分，$\lambda$点的位置和权重由Legendre多项式的根决定：

$$
\lambda_1 = 0.1127, \quad \lambda_2 = 0.5, \quad \lambda_3 = 0.8873
$$

$$
w_1 = 0.2778, \quad w_2 = 0.4444, \quad w_3 = 0.2778
$$

自由能差近似为：

$$
\Delta G \approx w_1 \left\langle \dfrac{\partial U}{\partial \lambda} \right\rangle_{\lambda_1} + w_2 \left\langle \dfrac{\partial U}{\partial \lambda} \right\rangle_{\lambda_2} + w_3 \left\langle \dfrac{\partial U}{\partial \lambda} \right\rangle_{\lambda_3}
$$

三点高斯积分可以精确积分5阶多项式，对多数较平滑的$\langle \partial U/\partial\lambda\rangle_\lambda$曲线已经够用，常被选作低成本的自由能积分方案。

> **TI在CDApol中的具体作用**：三点Gauss-Legendre积分将连续积分近似为三个加权和，让每个参数组合只需跑三个$\lambda$窗口的MD模拟就能估计HFE——省去了大量中间窗口的采样。

#### TI协议细节

SI中描述的TI协议包含两个独立的自由能变换：

1. **电荷变换**：从$Q=0$到金属离子的形式电荷（+2、+3或+4）
2. **LJ变换**：关闭金属离子与水分子之间的LJ相互作用

每个参数组合在三个$\lambda$窗口内采样（$\lambda = 0.11270, 0.5, 0.88729$）。$\lambda$状态下的势函数采用线性混合：

$$
V(\lambda) = (1 - \lambda)^k \cdot V_0 + \left[1 - (1 - \lambda)^k\right] \cdot V_1
$$

其中$V_0$是初态势能，$V_1$是末态势能。$k=1$时简化为标准线性插值（$V = (1-\lambda)V_0 + \lambda V_1$）。SI测试了不同$k$值，发现超过三个$\lambda$窗口并未显著改善结果，因此采用线性混合（$k=1$）和三点Gauss-Legendre积分即可满足精度需求。真空计算（无水环境）在一个窗口内即可快速收敛。

### CDApol偶极矩计算（SI Section S.3）

CDApol的瞬时偶极矩相对于分子质心计算：

$$
P_x = \sum_{i=1}^n q_i (x_i - x_c),\quad P_y = \sum_{i=1}^n q_i (y_i - y_c),\quad P_z = \sum_{i=1}^n q_i (z_i - z_c)
$$

$$
P = \sqrt{P_x^2 + P_y^2 + P_z^2}
$$

其中$(x_c, y_c, z_c)$是分子质心坐标，$q_i$是原子电荷。质心坐标由原子质量加权平均得到。SI的图S.3展示了50 ps NPT平衡过程中$\ce{Zn^{2+}}$、$\ce{Al^{3+}}$和$\ce{Zr^{4+}}$三种离子的中心离子和虚拟原子的电荷涨落。

偶极矩越大，说明电荷重新分布越明显。大小关系基本是$\ce{Zr^{4+}} > \ce{Zn^{2+}} > \ce{Al^{3+}}$，但并不是简单按价态单调变化：$\ce{Al^{3+}}$的中心离子会出现负电荷补偿，偶极方向也会跟着变。
