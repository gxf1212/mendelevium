---
title: "CDApol极化模型方法论详解：EEM动态电荷平衡的原理与实现"
date: "2026-05-06"
tags: [metal-ions, force-field, molecular-dynamics, polarization, cationic-dummy-atom, EEM, AMBER, high-valent-ions, methods]
description: "详解CDApol模型的电负性均衡方法（EEM）：动态电荷平衡的能量函数、约束求解、物理图像，以及双步参数化与热力学积分的技术细节"
image: "/assets/img/thumbnail_mine/wh-g83qpe.jpg"
author: Xufan Gao
lang: zh-CN
---

# CDApol极化模型方法论详解：EEM动态电荷平衡的原理与实现

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

CDApol模型通过**电负性均衡方法**（Electronegativity Equalization Method，EEM）引入动态极化。它的核心不是把整盒水分子都当成动态电荷一起优化，而是一个**局部动态电荷平衡问题**：每一步MD里，程序对CDApol核心位点（中心金属离子+6个虚拟原子）更新电荷$q_i$，并满足总电荷约束（+2、+3或+4）。

这里要先把“参与方程”和“作为未知量被优化”分开。按照ReaxFF/AMBER框架（JCTC 2020）的mEEM实现，外层原子（如MM水分子与过渡区原子）可以以**固定电荷**形式进入EEM方程，但它们不是动态电荷未知量；真正要求解的是核心区电荷$q_{\text{core}}$。因此，CDApol公式里的$N$仍然指中心金属离子和6个虚拟原子这7个核心位点，**不包括整盒水分子**。这一段是对求解机制的解释性重写，不是2022年CDApol主文逐项给出的原始符号表。

CDApol主文没有逐个定义“哪些水分子”进入$q_{\text{trans}}$。在ReaxFF/AMBER框架里，这个集合由核心区、过渡区、MM区划分和非键截断共同决定。更精确地说，它不是“固定名单”，而是每一步按坐标更新的运行时集合。下面这套记号，主要是为了把这种思路写清楚。

可以把筛选逻辑写成两层：

- **第一层（分区筛选）**：先按体系划分确定候选身份。
  - core 区原子进入动态未知量$q_{\text{core}}$。
  - transition 和 MM 区原子作为外层固定电荷候选，进入$q_{\text{trans}}$候选池。
- **第二层（距离筛选）**：在每一个MD步，只保留与core区发生有效非键耦合的外层原子。

$$
\mathcal{S}_{\text{trans}}(t)=\left\{j\in(\text{transition}\cup\text{MM})\mid \exists i\in\text{core},\ r_{ij}(t)<r_{\text{cut}}\right\}
$$

因此，$q_{\text{trans}}$对应的是集合$\mathcal{S}_{\text{trans}}(t)$里这些原子的固定电荷向量。由于水分子和配体都在运动，$\mathcal{S}_{\text{trans}}(t)$会随时间变化。

这也解释了为什么外层电荷“参与方程”但“不作为未知量”：它们通过$H_{\text{core-trans}}q_{\text{trans}}$给core区提供瞬时静电驱动，core区再在总电荷约束下重排$q_{\text{core}}$。换句话说，变化的是“参与耦合的外层原子集合”和它们相对core的几何关系，而不是把整盒水当作动态电荷一起优化。

这里还要统一符号：CDApol文中常写$J_{ij}$，ReaxFF/mEEM文中常写$H_{ij}$，二者本质上都扮演相互作用核的角色，但具体记号来自不同论文。

- $H_{ij}=\delta_{ij}\eta_i+(1-\delta_{ij})F_{ij}$
- $F_{ij}$是ReaxFF/mEEM里使用的**带屏蔽库仑核**
- 在CDApol记号里可对应理解为：对角项$J_{ii}=H_{ii}=\eta_i$，非对角项$J_{ij}$对应位点间的屏蔽静电耦合

因此，EEM能量可统一写成矩阵形式：

$$
E_{\text{EEM}}=\chi^{\mathsf T}q+\dfrac{1}{2}q^{\mathsf T}Hq
$$

在总电荷约束下，mEEM求解可以写成增广线性方程组的等价形式：

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
q_{\text{net}}
\end{bmatrix}
$$

把电荷分成核心区和外层固定区后，未知量只剩核心区电荷和约束乘子。整理后的mEEM方程可以等价写为两条方程：

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
q_{\text{net}}
\end{bmatrix}
-
\begin{bmatrix}
H_{\text{core-trans}} \\
\mathbf{0}^{\mathsf T}
\end{bmatrix}
q_{\text{trans}}
$$

其中$q_{\text{trans}}$是**常数向量**（AMBER固定电荷），不是新的动态电荷变量。$H_{\text{core-trans}}q_{\text{trans}}$可以理解为外层固定电荷在核心区产生的静电驱动项，也可以等价写成“有效电负性”形式：

$$
\chi_{\text{core}}^{\text{eff}}=\chi_{\text{core}}+H_{\text{core-trans}}q_{\text{trans}}
$$

再解core区的EEM即可。这样就能看清：$H_{ij}$本身不是总能量，而是二次能量项$\dfrac{1}{2}q^{\mathsf T}Hq$里的**系数矩阵**。

$$
E_{\text{EEM}} = \sum_{i=1}^{N} \chi_i q_i + \dfrac{1}{2} \sum_{i=1}^{N} \sum_{j=1}^{N} q_i J_{ij} q_j
$$

> **关键机制**：这里的$N$对应CDApol核心位点（1个中心离子+6个虚拟原子）。周围溶剂与配体的作用通过外部固定电荷产生的静电场耦合进入每步电荷平衡。随着MD中水分子位置变化，这个外场也变化，所以核心7个位点的最优电荷会随时间波动。

这条公式可以直接按**电荷流动方向和电荷重分布代价**来理解。我们先看公式里的每个符号：

- $q_i$：第$i$个位点的瞬时电荷（可正可负，单位是元电荷$e$）
- $\chi_i$：第$i$个位点的电负性参数（单位是能量，如eV）
- $J_{ij}$：位点$i$和$j$之间的相互作用矩阵元（对角项$J_{ii}=\eta_i$是硬度，非对角项是屏蔽库仑相互作用）
- $N$：CDApol核心位点数，即7个电荷位点（1个中心离子+6个虚拟原子）

### 1. 第一项：$\chi_i q_i$——电荷流动的驱动力

这一项决定**电荷想往哪里流**。虽然$\chi_i$在EEM中被称为Mulliken电负性参数，但它实际上是一个**可调的拟合参数**，只是借用了电负性的概念。传统的Mulliken电负性定义为：

$$
\chi = \dfrac{I + A}{2}
$$

其中$I$是电离能，$A$是电子亲和能。在化学中，电负性越大的原子（如氟、氧）越倾向于吸引电子。但在EEM模型里，$\chi_i$是通过拟合QM电荷分布得到的参数，**可以是正值也可以是负值**，其符号和大小决定了该位点在能量最小化时的电荷分配倾向。

能量项$\chi_i q_i$的物理含义：

- **$\chi_i$越小**：该位点越倾向于失去电荷（带正电）；**$\chi_i$越大**（更负）：越倾向于获得电荷（带负电）
- 如果$\chi_i$较小但仍然为正（位点不想要电子），而$q_i > 0$（位点带正电），那么$\chi_i q_i > 0$，**能量升高**——这是合理的，因为位点本来就不想要电荷，现在还带正电，能量当然高；如果$\chi_i$较大而$q_i < 0$（位点带负电），则势能很低
- 系统会自动调整$q_i$，让总能量$E_{\text{EEM}}$最小——这就是电荷重新分配的**驱动力**

### 2. 第二项：$\dfrac{1}{2} q_i J_{ij} q_j$——电荷重分布的代价

这一项决定**电荷重分布要付出什么代价**。它包含两部分：

#### 2.1 对角项：$J_{ii} = \eta_i$（自能代价）

对角项对应的是**单个位点上积累电荷的代价**。当$i=j$时，能量项变成：

$$
\dfrac{1}{2} \eta_i q_i^2
$$

这里$\eta_i$是Parr-Pearson硬度参数，物理上定义为：

$$
\eta_i = \dfrac{I_i - A_i}{2}
$$

也就是电离能和电子亲和能的**差值的一半**。

能量项$\dfrac{1}{2} \eta_i q_i^2$的物理含义：

- 这是一个**二次项**，无论$q_i$是正是负，$q_i^2$总是正的，所以这一项**总是让能量升高**——防止电荷无限制地堆到某一个位点上。**$\eta_i$越大**，电荷积累的代价越高，位点越硬，极化响应越弱；**$\eta_i$越小**，位点越软，极化响应越强

#### 2.2 非对角项：$J_{ij}$（位点间相互作用）

非对角项对应的是**两个不同位点之间的静电相互作用**。在CDApol主文里，这部分只强调采用了**electrostatic shielding**来避免近距离的过强排斥；若按ReaxFF/mEEM的写法理解，非对角项对应的是一种**带屏蔽的库仑核**，其强度随位点间距离和屏蔽参数变化。

- 能量项$\dfrac{1}{2} q_i J_{ij} q_j$的物理含义：$q_i$和$q_j$**同号**时相互排斥（能量升高），**异号**时相互吸引（能量降低）。
- 此外，位点越接近、屏蔽越弱，耦合作用通常越强。
- **$\gamma_{ij}$的物理意义**：
  - 如果没有屏蔽项，简单点电荷模型在短程会给出过强排斥。
  - 引入屏蔽后，短程相互作用会被软化，用来近似真实电子云不是点电荷这一事实。

**总结**：非对角项$\dfrac{1}{2} q_i J_{ij} q_j$描述**位点间的静电耦合**。它让电荷分布不能随意变化，因为同号电荷会互相排斥，异号电荷会互相吸引。屏蔽参数则用来抑制相邻位点之间的非物理短程排斥。

### 3. 总电荷约束：$\sum_i q_i = Q_{\text{total}}$ 与增广线性方程组求解

EEM真正求解的是一个**带约束的最小化问题**：

$$
\min_{\{q_i\}} E_{\text{EEM}}, \quad \sum_{i=1}^{N} q_i = Q_{\text{total}}
$$

在CDApol中，$Q_{\text{total}}$固定为金属离子的形式电荷（$\ce{Zn^{2+}}$的+2、$\ce{Al^{3+}}$的+3或$\ce{Zr^{4+}}$的+4）。这个约束的意思是：

- 电荷可以在中心离子和6个虚拟原子之间**自由流动**
- 但7个位点的电荷**总和必须守恒**
- 系统会自动找到一组$\{q_1, q_2, \dots, q_7\}$，让$E_{\text{EEM}}$最小，同时满足$q_1 + q_2 + \cdots + q_7 = Q_{\text{total}}$

**物理图像**：想象一个水池系统，7个水池通过管道连接，水可以在池子之间流动，但总水量不变。每个池子有自己的高度偏好（$\chi_i$）和容量限制（$\eta_i$），池子之间还有流动阻力（$J_{ij}$）。最终，水会流到一个平衡状态，让整个系统的势能（$E_{\text{EEM}}$）最低。

#### 增广线性方程组：求解动态电荷的机制

在实际计算中，这个带约束的最小化问题等价于求解一个**增广线性方程组**（augmented linear system）。令电负性函数对电荷的导数为零：

$$
\dfrac{\partial E_{\text{EEM}}}{\partial q_i} = \chi_i + \sum_j H_{ij} q_j - \varepsilon = 0
$$

这里的$\varepsilon$是**拉格朗日乘子**（Lagrange multiplier），它的作用是在最优解处**强制满足总电荷约束**。整理后的增广线性方程组可以写成：

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

其中$\mathbf{1}$是全1列向量，约束方程$\mathbf{1}^{\mathsf T}q = Q_{\text{total}}$出现在最后一行。

#### 拉格朗日乘子 $\varepsilon$ 的物理含义

$\varepsilon$ **不是电荷**，而是一个**标量约束惩罚项**。它代表在总电荷守恒约束下，CDApol核心区的**平均电化学势**（electrochemical potential）。具体来说：

- 如果某个位点$i$试图获得额外的电荷，就会增加约束的"成本"，表现为$\varepsilon$增大
- $\varepsilon$的大小反映了**核心区电荷分布的紧张程度**——$\varepsilon$越大，说明系统在坚持满足$\sum q_i = Q_{\text{total}}$这个约束付出的代价越高
- 在MD模拟中，每一步都重新求解这个增广系统，得到新的$q$和新的$\varepsilon$。随着水分子位置变化，外层固定电荷提供的静电驱动（$H_{\text{core-trans}}q_{\text{trans}}$）会变，$\varepsilon$也会波动

> **关键区别**：$\varepsilon$ 是对**约束满足情况的度量**，而不是系统中的第七个电荷位点。核心位点仍然只有7个。

### 4. CDApol中的动态极化

中心离子和虚拟原子分别有自己的$\chi_i$、$\eta_i$和$\gamma_{ij}$参数，这些参数通过QM参考电荷分布训练得到。这样，每一步MD中，水分子或配体的局部构型一变，EEM就会重新分配7个位点上的电荷。

总电荷仍然守恒，但局部电荷分布可以响应环境变化，这就是CDApol里**动态极化**的物理来源。

> CDApol的关键不在于“把电荷做大或做小”，而在于**在总电荷守恒下让电荷分布随局部构型实时重排**。

## 参数化流程

CDApol的参数化分两步，且**严格串联**：

| 步骤 | 具体做法 | 参数状态 |
| --- | --- | --- |
| **第一步**：极化力场参数训练 | 优化中心离子和虚拟原子的EEM参数（$\chi_i, \eta_i, \gamma_{ij}$），以 $\ce{Zn^{2+}}$、$\ce{Al^{3+}}$、$\ce{Zr^{4+}}$ 的7配位水合构象DFT电荷分布为参考 | EEM参数从无到有 |
| **第二步**：12-6 LJ参数扫描 | 用**已冻结的**EEM参数，扫描$\varepsilon$（1-5 kcal/mol，步长0.2）和$R_{\min}/2$（0.1 Å间隔），每个点跑MD+TI计算HFE、IOD、CN | EEM参数锁定；LJ参数扫描；最后挑最优点 |

### 参数冻结时序的关键原理

这个两步设计的核心在于**严格的参数分离**：

1. **第一步完全独立**：只依赖DFT参考电荷，拟合$\chi_i, \eta_i, \gamma_{ij}$。此时不考虑LJ参数，也不计算任何HFE。
2. **第一步结束后**：$\chi_i, \eta_i, \gamma_{ij}$全部锁定，不再改变。
3. **第二步开始**：在固定EEM参数的前提下，进行LJ参数网格搜索。每个$(epsilon, R_{\min}/2)$组合都要：
   - 跑一次完整的MD+TI计算
   - 检查计算出的HFE、IOD、CN与实验值的偏差
4. **第二步的优化目标**：在庞大的LJ参数空间里找一组最优值，使得HFE、IOD、CN同时接近实验值。

#### 为什么要这样做

- 如果第一步和第二步**同时优化**所有参数，系统会陷入过参数化（overparameterization），而且EEM参数的拟合目标（QM电荷）和LJ参数的拟合目标（实验热力学性质）可能互相冲突。
- 通过**参数锁定**，确保EEM参数专一地处理**局部电荷平衡**，LJ参数则处理**非键相互作用**，两者各司其职。
- 这也是为什么你的理解"是不是得先固定参数才能最小化能量"是正确的。在每一步MD中，**EEM参数固定**，EEM只是通过增广线性方程组计算给定外部环境下的最优电荷；LJ参数则在MD和TI的总体框架中被优化。

![fig2](2026-05-06-polarizable-cationic-dummy-model_figs/fig2.png)

**图2：CDApol势函数训练和建模流程图**

- **第一步（左侧）**：优化电负性均衡（EEM）参数，具体包括电负性 $\chi$、硬度 $\eta$ 和静电屏蔽 $\gamma_{ij}$ 三个参数族。参考目标是DFT计算得到的水合构象电荷分布（7配位：6个水分子 + 1个虚拟第七配体）。
- **第二步（右侧）**：扫描Lennard-Jones参数 $\varepsilon$（1-5 kcal/mol，步长0.2）和 $R_{\min}/2$，用热力学积分（TI）计算水合自由能（HFE）。同时验证离子-氧距离（IOD）和配位数（CN）与实验值的一致性。
- **整体逻辑**：两步串联，先定电荷分布（第一步），再调非键参数（第二步），确保动态极化效应（HFE敏感度降低）和结构性质（IOD/CN准确性）同时达标。

图注说明该图描述的是完整的参数化管线，而非单一物理量的拟合。两步之间的箭头表示数据流向：第一步输出的电荷分布约束第二步的LJ参数搜索空间，提高参数优化效率。

### DFT参考数据与计算约束

参数优化以DFT得到的**参考电荷分布**为目标。这里选的是**7配位水合构象**：中心金属离子周围配6个水分子，再加1个虚拟第七配体。这样得到的参考态更接近真实配位壳层，也更适合拿来约束EEM参数。

总电荷约束规定，**电荷可以在中心离子和虚拟原子之间移动，但整套模型的总电荷不能变**。CDApol在每一步MD里都会重新分配电荷，但它们的总和始终等于金属离子的形式电荷，也就是+2、+3或+4。这个约束通常可以通过增广线性方程组或等效的约束处理来实现。

### 最终参数值的文献出处

CDApol论文（Rahnamoun et al. 2022）中，**最终参数值的记载如下**：

- **LJ参数（$\varepsilon, R_{\min}/2$）**：完整列表在论文Table 2中，包括三个金属离子（Zn²⁺、Al³⁺、Zr⁴⁺）在不同参数点下的HFE、IOD、CN计算值与实验值对比。最终被选中的参数组合是使得三个性质同时最接近实验值的那一组。
- **EEM参数（$\chi_i, \eta_i, \gamma_{ij}$）**：详细值在论文Supporting Information中。论文正文第3页指出："Further details regarding the parameter optimization procedures used in this work can be found in the Supporting Information"。SI的标题为"Polarizable and nonpolarizable potential parameter optimizations, CDApol charge fluctuations, dipole moment calculation, and electronegativity equalization method (PDF)"。
- **获取方式**：SI可在ACS Publications官方网站免费下载：https://pubs.acs.org/doi/10.1021/acs.jpclett.2c01279

这一分离提供了**完全的参数化透明度**：LJ参数在主文可查，EEM参数在SI可查。这样的设计便于后续研究者验证、复现或扩展CDApol模型。

### 热力学积分与三点高斯积分

自由能差通过**热力学积分**（Thermodynamic Integration，TI）计算。TI的核心思想是沿着一条连接初态和末态的路径，逐步“充电”或“去充电”，然后对路径上的能量导数积分，得到自由能差：

$$
\Delta G = \int_0^1 \left\langle \dfrac{\partial U(\lambda)}{\partial \lambda} \right\rangle_\lambda \mathrm{d}\lambda
$$

其中$\lambda$是耦合参数（$\lambda=0$对应初态，$\lambda=1$对应末态），$U(\lambda)$是$\lambda$状态下的势能，$\langle \cdots \rangle_\lambda$表示在$\lambda$状态下的系综平均。

**为什么需要数值积分**：上面的积分无法解析求解，因为$\left\langle \dfrac{\partial U}{\partial \lambda} \right\rangle_\lambda$只能通过MD模拟在离散的$\lambda$点上采样得到。因此，我们需要选择若干个$\lambda$点，在每个点上跑MD模拟，然后用数值积分方法把这些点连起来。

**三点高斯积分的优势**：本文采用**三点高斯积分**（Three-point Gaussian Quadrature）来近似上述积分。高斯积分是一种高效的数值积分方法，它通过精心选择积分点的位置和权重，用较少的采样点获得较高精度。

对于三点高斯积分，$\lambda$点的位置和权重由Legendre多项式的根决定。下面这组节点和权重是**标准三点Gauss-Legendre公式**：

$$
\lambda_1 = 0.1127, \quad \lambda_2 = 0.5, \quad \lambda_3 = 0.8873
$$

对应的权重为：

$$
w_1 = 0.2778, \quad w_2 = 0.4444, \quad w_3 = 0.2778
$$

自由能差近似为：

$$
\Delta G \approx w_1 \left\langle \dfrac{\partial U}{\partial \lambda} \right\rangle_{\lambda_1} + w_2 \left\langle \dfrac{\partial U}{\partial \lambda} \right\rangle_{\lambda_2} + w_3 \left\langle \dfrac{\partial U}{\partial \lambda} \right\rangle_{\lambda_3}
$$

**为什么三点就够了**：高斯积分的精度取决于被积函数的光滑程度。对于多项式函数，$n$点高斯积分可以精确积分$2n-1$阶多项式。三点高斯积分可以精确积分5阶多项式，对于多数较平滑的$\left\langle \dfrac{\partial U}{\partial \lambda} \right\rangle_\lambda$曲线通常已经够用；这也是它常被选作低成本自由能积分方案的原因。

**CDApol中的具体实现**：本文在每个$(\varepsilon, R_{\min}/2)$参数组合下，沿着从$Q=0$（无电荷）到目标形式电荷（+2、+3或+4）的路径进行充电，并用三点高斯积分来估计水合自由能（HFE）。把它写成$\lambda$积分形式，是为了更直观地说明这一数值过程。整个过程在TIP3P或OPC水模型中进行，使用周期性边界条件。

## 模型行为验证：瞬时偶极矩评估

偶极矩用来量化电荷分布的波动。本文在50 ps NPT平衡阶段收集**1000个快照**，时间间隔约50 fs，然后统计每个快照上的瞬时偶极矩：

| 离子 | 平均偶极矩 (D) | 物理意义 |
| --- | --- | --- |
| $\ce{Zn^{2+}}$ CDApol | 0.32 | 中等极化，电荷分布相对稳定 |
| $\ce{Al^{3+}}$ CDApol | 0.22 | 最小极化，中心离子负电荷补偿了虚拟原子的极化 |
| $\ce{Zr^{4+}}$ CDApol | 0.53 | 最大极化，高价离子产生强烈的诱导偶极 |

偶极矩越大，说明电荷重新分布越明显。这里整体上是$\ce{Zr^{4+}} > \ce{Zn^{2+}} > \ce{Al^{3+}}$，但它不是简单按价态单调变化，因为$\ce{Al^{3+}}$的中心离子会出现负电荷补偿，偶极方向也会跟着变。
