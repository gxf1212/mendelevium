---
title: "ABCG2电荷模型技术细节：附录"
date: "2025-11-02"
tags: [ABCG2, BCC参数, 电荷模型, ACES方法, 模拟协议]
description: "ABCG2电荷模型的ACES方法、模拟参数和验证协议详细技术附录"
image: "/assets/img/thumbnail_mine/wh-m992d8.jpg"
thumbnail: "/assets/img/La-Mancha.jpg"
author: Xufan Gao
lang: zh-CN
---

# ABCG2技术细节附录

本文档为《优化单一性质≠改善相关性质：ABCG2电荷模型的启示》的技术附录，详细介绍ACES自由能计算方法、模拟参数设置和验证协议。

## 附录A：ACES（Alchemical Enhanced Sampling）自由能计算方法

### A.1 热力学积分框架

ABCG2验证采用**ACES**方法进行高精度自由能计算，这是一种基于哈密顿副本交换分子动力学（HREMD）的热力学积分方法。

**基本原理**：通过λ参数控制初始态和最终态之间的平滑变换，计算自由能差：

$$\Delta G = \int_0^1 \left\langle \frac{\partial H}{\partial \lambda} \right\rangle_\lambda d\lambda$$

其中H为哈密顿量，$\langle \cdot \rangle_\lambda$表示在λ状态下的系综平均。

### A.2 λ状态设置

**炼金术变换参数**：
- **λ状态数量**：11个状态
- **λ值范围**：0.0, 0.1, 0.2, ..., 1.0
- **软核势**：Smooth Step Softcore（用于避免原子碰撞）
- **耦合方案**：VDW和静电相互作用同步耦合

### A.3 HREMD采样策略

**副本交换设置**：
- **交换频率**：每20 MD步尝试一次Hamiltonian交换
- **交换总次数**：每个λ状态进行100,000次交换尝试
- **4次独立运行**：每个系统重复4次相同的模拟

### A.4 模拟协议详细参数

#### 气相系统
- **初始化**：几何最小化（避免立体碰撞）
- **NVT平衡**：0.5 ns at 298 K（Langevin恒温器，衰减系数100 ps^-1）
- **生产阶段**：2.0 ns HREMD
- **总采样深度**：每λ状态等效2,000,000 MD步

#### 液相系统
- **初始化**：几何最小化
- **NVT平衡**：0.5 ns at 298 K
- **NPT平衡**：3.0 ns at 1 atm, 298 K（Monte Carlo压力控制器）
- **溶剂盒设置**：40 Å三斜晶系盒子，与溶质至少2.5 Å间距
- **生产阶段**：2.0 ns HREMD
- **总采样深度**：每λ状态等效2,000,000 MD步

#### 通用MD参数
- **时间步长**：1 fs
- **温度**：298 K
- **压力**：1 atm（仅液相）
- **截断方案**：Particle Mesh Ewald（PME）电磁势，VDW截断12 Å
- **约束条件**：所有含H键约束（SHAKE算法）

## 附录B：数据集详细信息

### B.1 FreeSolv数据库

**数据库特征**：
- **总分子数**：642个中性有机分子
- **分子量范围**：16-499 g/mol
- **官能团覆盖**：30种主要官能团
- **数据来源**：由Dr. J. P. Guthrie精心编制和验证

**分阶段开发**：
- **FreeSolv_p1**：441个单官能团分子
- **FreeSolv_p2**：201个多官能团+含P分子

### B.2 验证数据集

**MNSol数据库**：
- **溶质-溶剂对数**：2068对
- **溶剂种类**：89种有机溶剂
- **用途**：多溶剂环境下的转移自由能验证

**ATB3.0验证集**：
- **分子数**：685个
- **数据要求**：ΔG<sub>exp</sub>误差<1 kcal/mol
- **用途**：高精度基准验证

## 附录C：电荷分配工作流程

### C.1 输入数据处理

**数据来源和格式**：
- FreeSolv：xyz文件
- MNSol：xyz文件
- ATB3.0：xyz文件

**结构检查与修正**：
1. Schrödinger Maestro v11.2进行人工检查
2. 设置正确的键类型和原子参数
3. 转换为统一mol2格式

### C.2 ABCG2电荷分配

**命令行工具**：
```bash
antechamber -i molecule.mol2 -fi mol2 \
  -o molecule.prepi -fo prepi -c abcg2
```

**工作流程**：
1. AM1半经验几何优化（Sqm模块）
2. Mulliken电荷计算
3. BCC参数表查询和应用
4. 最终电荷分配

## 附录D：统计分析方法

### D.1 性能指标定义

**主要指标**：
- **Mean Signed Error (MSE)**：
  $$\text{MSE} = \frac{1}{N}\sum_i (\Delta G_i^{calc} - \Delta G_i^{exp})$$
- **Mean Unsigned Error (MUE)**：
  $$\text{MUE} = \frac{1}{N}\sum_i |\Delta G_i^{calc} - \Delta G_i^{exp}|$$
- **Root Mean Square Error (RMSE)**：
  $$\text{RMSE} = \sqrt{\frac{1}{N}\sum_i (\Delta G_i^{calc} - \Delta G_i^{exp})^2}$$
- **Pearson相关系数 (R)**：线性相关性度量
- **Spearman秩相关系数 (ρ)**：非参数相关性度量

### D.2 统计检验

**配对Student's t检验**：
- 比较三种力场组合的RMSE差异
- 评估差异是否具有统计显著性（p < 0.05）
- 计算95%置信区间

### D.3 误差分析

**误差分布特性**：
- ±1 kcal/mol范围内的数据比例
- ±2 kcal/mol范围内的数据比例
- 离群点（outliers）的识别和分析

## 附录E：相关资源和工具

### 软件工具
- **GROMACS**：分子动力学模拟引擎 (https://www.gromacs.org/)
- **AmberTools**：含ABCG2参数和Antechamber模块 (https://ambermd.org/)
- **pmx**：非平衡炼金术工具 (https://github.com/deGrootLab/pmx)
- **Schrödinger Maestro**：结构准备和验证

### 数据库
- **FreeSolv**：https://github.com/MobleyLab/FreeSolv
- **OpenFE数据集**：https://github.com/OpenFreeEnergy/openfe-data

### 原始论文数据
- **ABCG2原始论文**：He et al., J. Chem. Theory Comput. 2025, 21, 3032–3043
- **评估论文**：Behera et al., J. Chem. Inf. Model. 2025 (Letter)

## 附录F：蛋白-配体RBFE评估的模拟协议

### F.1 数据集来源

**OpenFE蛋白-配体数据集**：
- **来源**：OpenFE协会提供的基准数据集（Ross et al. 2023）
- **规模**：
  - 12个蛋白靶点
  - 273个配体
  - 507个配体微扰（ligand perturbations）
- **覆盖范围**：
  - 'jacs_set'（273个转化）：通用靶点集合
  - 'janssen_bace'（234个转化）：BACE相关靶点（bace_cp, bace_p3等）
- **质量标准**：所有配体均基于临床或实验化合物

### F.2 非平衡炼金术（Nonequilibrium Alchemical Free Energy）协议

**模拟框架**：采用pmx工具进行非平衡FEP（Jarzynski等式和Crooks涨落定理）

#### F.2.1 蛋白系统准备

**结构准备**：
1. 蛋白结构来自PDB数据库或实验提供
2. 质子化状态使用PDB2PQR确定（pH 7.4）
3. 使用Schrödinger Maestro进行配体对接与姿态优化
4. 配体使用GAFF2或GAFF2-ABCG2力场参数化

**力场选择**：
- **配体力场**：GAFF2（基础）+ AM1-BCC或ABCG2电荷
- **蛋白力场**（两种）：
  - AMBER99SB*-ILDN（基准）
  - AMBER14SB（改进版对照）
- **溶剂力场**：TIP3P水（标准）

#### F.2.2 系统构建与平衡

**盒子大小**：
- 蛋白周围距离至少14 Å的水盒子
- 三斜晶系（triclinic）盒子，最小化周期性人工物

**离子补偿**：
- Na⁺/Cl⁻补偿系统电荷
- 最终离子浓度约0.15 M（生理浓度）

**平衡协议**：
1. **几何最小化**：1000步，能量收敛
2. **NVT平衡**（2 ns）：
   - 温度：298 K
   - 恒温器：Langevin，衰减系数100 ps⁻¹
3. **NPT平衡**（3 ns）：
   - 温度：298 K，压力：1 atm
   - 压力控制：Berendsen压力浴
4. **分子约束**：所有含H键约束（SHAKE）

#### F.2.3 非平衡FEP生产阶段

**λ变换参数**：
- **λ状态数量**：5个（0.0, 0.25, 0.5, 0.75, 1.0）
- **变换路径**：VDW和静电相互作用同步耦合（单一λ参数）
- **软核势**：C6/C12软核势用于避免原子碰撞

**模拟参数**：
- **时间步长**：2 fs（使用H-mass repartitioning允许更大时步）
- **运行时间/λ**：1 ns
- **每个转化的总运行时间**：5 ns（5个λ × 1 ns）
- **驱动速度**：λ通常以0.2 ns⁻¹速率驱动（总耗时1 ns）
- **数据采集频率**：每1 ps记录一次配置

**物理常数与截断**：
- **温度控制**：Langevin恒温器（298 K，衰减系数0.1 ps⁻¹）
- **范德华截断**：12 Å
- **静电势**：PME（Particle Mesh Ewald），精度1e⁻6
- **压力控制**：NPT条件下Parrinello-Rahman压力控制器

#### F.2.4 多个独立重复与误差估计

**重复计算**：
- **每个配体微扰**：进行3-5次独立的FEP模拟（不同的初始速度）
- **平衡数据排除**：前100 ps作为平衡期舍弃
- **误差估计**：
  - 使用standard error of the mean（SEM）统计多次运行
  - 使用Jarzynski等式处理不可逆工作
  - 使用动态无偏估计器（BAR, Bennett Acceptance Ratio）整合多条轨迹

### F.3 结果分析与统计

**自由能计算**：
- **相对结合自由能（ΔΔG）**：直接从FEP得到
- **绝对结合自由能（ΔG）**：使用Cinnabar最大似然估计法将ΔΔG累积为ΔG
- **95%置信区间**：基于bootstrap重采样或标准差

**精度评估指标**：
- **RMSE**（Root Mean Square Error）：主要精度指标
- **MUE**（Mean Unsigned Error）：绝对误差平均值
- **Pearson相关系数（r）**：计算与实验的线性相关性
- **Spearman秩相关系数（ρ）**：非参数相关性（化合物排名能力）
- **Kendall's τ**：另一种非参数排名相关性
- **配对Student's t检验**：比较不同力场组合的显著性差异（p值）

### F.4 官能团子分析

**分类标准**：
- 根据配体中改变的官能团分类转化（酮、醚、醇、芳香烃、喹啉等）
- 一个转化可能跨越多个官能团类别（如联苯既属"联苯"也属"芳香烃"）

**统计处理**：
- 仅显示RMSE差异>1 kJ/mol（0.24 kcal/mol）的官能团
- 对所有官能团组进行配对t检验评估显著性
- 补充分析在补充图S16中呈现

### F.5 主要参考配体与案例分析

**两个对比案例**：
1. **叔醇案例**（p38靶点，转化2y→2v）：
   - 实验ΔΔG = 0.81 kcal/mol
   - AM1-BCC预测：2.47 ± 0.26 kcal/mol（偏离）
   - ABCG2预测：0.49 ± 0.20 kcal/mol（接近）
   - ABCG2改进

2. **喹啉案例**（mcl1靶点，转化47→27）：
   - 实验ΔΔG = −0.34 kcal/mol
   - AM1-BCC预测：−0.42 ± 0.52 kcal/mol（接近）
   - ABCG2预测：−3.11 ± 0.23 kcal/mol（严重偏离）
   - ABCG2变差

这两个案例展示了：电荷模型的效能在蛋白环境中具有**化学环境特异性**，同一模型不能保证在所有官能团上都表现一致。

## 附录G：HREMD Reweighting 物理公式总结

### G.1 统计力学基础

**HREMD（Hamiltonian Replica Exchange Molecular Dynamics）**通过在不同 Hamiltonian（lambda 值）间交换构型，实现对复杂自由能面的高效采样。Reweighting 的核心问题是：**如何从多个 lambda replicas 的样本中，准确重构目标 lambda 的系综平均？**

**系综分布关系**：
在温度 $T$ 下，不同 lambda 的系综分布满足：

$$\frac{\rho(\mathbf{r};\lambda_0)}{\rho(\mathbf{r};\lambda_i)} = \frac{Z(\lambda_i)}{Z(\lambda_0)} \exp\left[-\beta\Delta U_{0i}(\mathbf{r})\right]$$

其中：
- $\rho(\mathbf{r};\lambda)$ 是构型 $\mathbf{r}$ 在 lambda $\lambda$ 下的概率密度
- $Z(\lambda)$ 是配分函数
- $\Delta U_{0i}(\mathbf{r}) = U(\mathbf{r};\lambda_0) - U(\mathbf{r};\lambda_i)$ 是势能差
- $\beta = \frac{1}{k_B T}$

### G.2 核心重加权公式

#### 2.1 单 Replica 重加权

对于在目标 lambda $\lambda_0$ 的系综平均，可以从任意 replica $i$ 的样本重加权得到：

$$\langle A \rangle_{\lambda_0} = \frac{\langle A \exp[-\beta\Delta U_{0i}] \rangle_{\lambda_i}}{\langle \exp[-\beta\Delta U_{0i}] \rangle_{\lambda_i}}$$

**通俗解释**：这就像用"汇率"把不同货币的样本转换成目标货币。$\exp[-\beta\Delta U_{0i}]$ 就是转换汇率，把 replica $i$ 的样本值 "折算" 成目标 lambda $\lambda_0$ 的价值。

#### 2.2 多 Replica 综合公式（实际使用）

对于 HREMD 中 $M$ 个 replicas，综合所有样本：

$$\langle A \rangle_{\lambda_0} = \frac{\sum_{i=1}^M \sum_{j=1}^{N_i} A_{i,j} \exp[-\beta\Delta U_{0i}(\mathbf{r}_{i,j})]}{\sum_{i=1}^M \sum_{j=1}^{N_i} \exp[-\beta\Delta U_{0i}(\mathbf{r}_{i,j})]}$$

其中：
- $N_i$ 是 replica $i$ 的样本数
- $A_{i,j}$ 是第 $i$ 个 replica 第 $j$ 个样本的观测值
- $\mathbf{r}_{i,j}$ 是对应的构型
- $\Delta U_{0i}(\mathbf{r}_{i,j}) = U(\mathbf{r}_{i,j};\lambda_0) - U(\mathbf{r}_{i,j};\lambda_i)$

**物理意义**：这是最大似然估计，相当于用所有 replicas 的样本，通过各自的权重，加权平均得到目标 lambda 的期望值。

### G.3 有效样本量和统计质量

#### 3.1 有效样本量计算

由于不同样本的权重不同，实际的有效样本量会减少：

$$
N_{\text{eff}} = \frac{(\sum_{i,j} w_{i,j})^2}{\sum_{i,j} w_{i,j}^2}
$$

其中权重 $w_{i,j} = \exp[-\beta\Delta U_{0i}(\mathbf{r}_{i,j})]$

**重要性**：
- $N_{\text{eff}}/N_{\text{total}} > 0.1$ 通常认为是良好的重叠
- $N_{\text{eff}}$ 太小说明 replica 间重叠不足，误差会很大

#### 3.2 方差估计

重加权估计的方差：

$$
\text{Var}(\langle A \rangle_{\lambda_0}) \approx \frac{1}{N_{\text{eff}}} \frac{\sum_{i,j} w_{i,j} (A_{i,j} - \langle A \rangle_{\lambda_0})^2}{\sum_{i,j} w_{i,j}}
$$

**通俗解释**：有效样本量直接决定了估计的可靠性。如果某些样本的权重特别大（说明它们在目标 lambda 中很重要），但数量很少，那么整个估计就会不稳定。

### G.4 实际应用注意事项

#### 4.1 权重截断策略

**问题**：极端权重会导致数值不稳定和统计偏差
**解决方案**：
- **绝对截断**：设定最大权重 $w_{\max} = \alpha \bar{w}$（通常 $\alpha = 3-5$）
- **相对截断**：使用 $w' = \frac{w}{1 + \epsilon w}$ 进行平滑处理

#### 4.2 交换率优化

**HREMD 交换概率**：
$$P_{\text{acc}}(i \leftrightarrow j) = \min\left[1, \exp\left(-\beta\Delta U_{ji} + \beta\Delta U_{ij}\right)\right]$$

**最优交换率**：一般在 20-40% 之间
- 太低：采样效率不高
- 太高：lambda 间隔太大，重叠不足

#### 4.3 收敛性判断

**收敛标准**：
1. **有效样本量稳定**：$N_{\text{eff}}$ 不再随时间增加
2. **权重分布合理**：避免极端权重（如 $\exp(10)$ 以上）
3. **块平均一致**：不同时间段的平均值应该一致

### G.5 高级方法：WHAM/MBAR

#### 5.1 WHAM（Weighted Histogram Analysis Method）

**基本思想**：同时优化所有 lambda 的配分函数，提高统计效率

**公式**：
$$\hat{F}_i = -\ln \sum_{j=1}^M \sum_{n=1}^{N_j} \frac{\exp(-\beta U_i(\mathbf{x}_{j,n}))}{\sum_{k=1}^M N_k \exp(\hat{F}_k - \beta U_k(\mathbf{x}_{j,n}))}$$

#### 5.2 MBAR（Multistate Bennett Acceptance Ratio）

**优势**：考虑样本间的相关性，理论上更优

**适用场景**：
- 样本数量有限
- 需要 highest precision
- 多个目标态都需要估计

### G.7 常见问题与解决方案

#### 问题1：负权重
**原因**：$\Delta U_{0i} > 0$ 且很大时，$\exp[-\beta\Delta U_{0i}]$ 会很小
**解决**：使用相对权重或截断

#### 问题2：重叠不足
**表现**：$N_{\text{eff}}/N_{\text{total}} < 0.1$
**解决**：增加 lambda 点数，调整 lambda 间隔

#### 问题3：计算成本高
**策略**：
- 使用重要性采样
- 并行化计算
- 预先计算权重

### G.8 物理意义总结

**Reweighting 的本质**：
1. **统计推断**：从容易采样的分布推断难采样的分布
2. **信息利用**：充分利用所有 lambda 的样本信息
3. **误差传播**：样本的统计误差会影响最终结果的精度

**关键洞见**：HREMD reweighting 证明了通过物理定律，我们可以从"不完美"的采样中获得"完美"的统计推断。这就像用散乱的拼图碎片，通过数学方法还原出完整的图像。
