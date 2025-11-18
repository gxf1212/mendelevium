---
title: "分子性质预测：机器学习回归算法详解（一）基础回归模型"
date: "2025-11-10"
tags: [machine-learning, regression, linear-models, svm, knn, molecular-property]
description: "系列第一篇：全面介绍分子性质预测中的基础回归模型，包括线性模型家族、支持向量回归、近邻方法，从简单的线性回归到鲁棒回归、广义线性模型"
thumbnail: "/assets/img/thumbnail_mine/wh-z8odwg.jpg"
image: "/assets/img/thumbnail_mine/wh-z8odwg.jpg"
author: Xufan Gao
lang: zh-CN
---

# 分子性质预测：机器学习回归算法详解（一）基础回归模型

> **系列导航**：
> - **第一篇：基础回归模型**（本文）- 线性模型、支持向量机、近邻方法
> - [第二篇：树模型与梯度提升](2025-11-10-ml-regression-models-part2-trees.md) - 决策树、随机森林、XGBoost/LightGBM等
> - [第三篇：高级模型与应用指南](2025-11-10-ml-regression-models-part3-advanced.md) - 神经网络、概率模型、VAE、模型选择指南

## 导读

在分子性质预测、药物筛选、材料设计等回归任务中，**选对机器学习模型至关重要**。本系列文章将介绍30余种经典和前沿的回归算法，剖析每个模型的原理、公式和适用场景。

**第一篇聚焦基础回归模型**，这些模型是理解现代机器学习的基石：
- **线性模型家族**：从简单的线性回归到鲁棒回归、广义线性模型
- **支持向量回归**：处理非线性关系的经典方法
- **近邻方法**：基于相似性的简单有效算法

所有模型均基于scikit-learn实现，可直接用于实践。

---

## 1. 线性模型家族

### 1.1 核心思想

线性模型假设输入特征与目标值之间存在**线性关系**，通过学习特征权重来进行预测。这是最基础也是最可解释的模型类型。

### 1.2 基础线性模型

#### LinearRegression（线性回归器）

**原理**：最小化预测值与真实值之间的平方误差。

**sklearn实现**：`from sklearn.linear_model import LinearRegression`

**数学公式**：
$$
\hat{y} = w_0 + w_1x_1 + w_2x_2 + \cdots + w_nx_n = \mathbf{w}^T\mathbf{x}
$$

$$
\min_{\mathbf{w}} \sum_{i=1}^{m} (y_i - \mathbf{w}^T\mathbf{x}_i)^2
$$

**特点**：
- ✅ **快速训练**：解析解，无需迭代
- ✅ **高度可解释**：每个特征的权重清晰可见
- ❌ **容易过拟合**：高维数据时权重不稳定
- 📊 **推荐场景**：分子性质预测的baseline模型

---

#### Ridge（岭回归器）

**原理**：在线性回归基础上加入**L2正则化**，防止权重过大。

**sklearn实现**：`from sklearn.linear_model import Ridge`

**数学公式**：
$$
\min_{\mathbf{w}} \sum_{i=1}^{m} (y_i - \mathbf{w}^T\mathbf{x}_i)^2 + \alpha \|\mathbf{w}\|_2^2
$$

**特点**：
- ✅ **缓解共线性**：相关特征的权重更稳定
- ✅ **防止过拟合**：正则化参数 $\alpha$ 控制模型复杂度
- 📊 **推荐场景**：特征数量接近或超过样本数量的高维分子数据

---

#### Lasso（套索回归器）

**原理**：使用**L1正则化**，可将部分特征权重压缩为0，实现**特征选择**。

**sklearn实现**：`from sklearn.linear_model import Lasso`

**数学公式**：
$$
\min_{\mathbf{w}} \sum_{i=1}^{m} (y_i - \mathbf{w}^T\mathbf{x}_i)^2 + \alpha \|\mathbf{w}\|_1
$$

**特点**：
- ✅ **自动特征选择**：不重要的特征权重为0
- ✅ **稀疏解**：结果更简洁易懂
- 📊 **推荐场景**：需要识别关键分子描述符时

---

#### ElasticNet（弹性网络回归器）

**原理**：结合L1和L2正则化，平衡两者优势。

**sklearn实现**：`from sklearn.linear_model import ElasticNet`

**数学公式**：
$$
\min_{\mathbf{w}} \sum_{i=1}^{m} (y_i - \mathbf{w}^T\mathbf{x}_i)^2 + \alpha \rho \|\mathbf{w}\|_1 + \frac{\alpha(1-\rho)}{2} \|\mathbf{w}\|_2^2
$$

其中 $\rho$ 控制L1和L2的比例。

**特点**：
- ✅ **综合优势**：既能特征选择又能处理共线性
- ✅ **灵活调节**：通过 $\rho$ 调整L1/L2权重
- 📊 **推荐场景**：复杂分子数据集的通用首选回归器

---

#### SGDRegressor（随机梯度下降回归器）

**原理**：通过逐样本更新权重实现在线学习，适合超大规模数据。

**sklearn实现**：`from sklearn.linear_model import SGDRegressor`

**特点**：
- ✅ **内存高效**：无需一次加载所有数据
- ✅ **支持在线学习**：数据流式更新模型
- ⚡ **快速收敛**：大数据集训练速度快
- 📊 **推荐场景**：大规模分子数据库的增量学习

---

### 1.3 鲁棒回归家族

当数据中存在**异常值**或**重尾噪声**时，标准最小二乘法会失效。鲁棒回归模型通过特殊的损失函数降低异常值的影响。

#### HuberRegressor（胡伯回归器）

**损失函数**：
$$
L_\delta(r) = \begin{cases}
\frac{1}{2}r^2 & |r| \leq \delta \\
\delta(|r| - \frac{1}{2}\delta) & |r| > \delta
\end{cases}
$$

小误差用平方损失（L2），大误差用绝对值损失（L1），平衡效率和鲁棒性。

**sklearn实现**：`from sklearn.linear_model import HuberRegressor`

**特点**：
- ✅ **对中等异常值鲁棒**：适度降低离群点影响
- ⚙️ **关键参数**：`epsilon`（平方/线性损失转换点，默认1.35）
- 📊 **推荐场景**：包含中等异常值的分子性质回归

---

#### TheilSenRegressor（西尔森回归器）

**核心思想**：基于样本对斜率的**中位数**估计，对异常值具有极强鲁棒性。

**sklearn实现**：`from sklearn.linear_model import TheilSenRegressor`

**特点**：
- ✅ **极强鲁棒性**：可容忍29.3%的异常值
- ❌ **仅适用低维**：特征数 <20
- ❌ **计算复杂度高**：$O(n^2)$
- 📊 **推荐场景**：实验数据质量差，离群点多的分子性质回归

---

#### RANSACRegressor（随机采样一致性回归器）

**核心思想**：**随机采样一致性**（Random Sample Consensus）算法，通过迭代随机采样找到最优内点集。

**sklearn实现**：`from sklearn.linear_model import RANSACRegressor`

**算法流程**：
1. 随机采样最小样本集，拟合模型
2. 计算所有样本的残差
3. 将残差小于阈值的样本标记为内点
4. 重复1-3，选择内点最多的模型
5. 使用所有内点重新拟合最终模型

**特点**：
- ✅ **极强鲁棒性**：可处理>50%异常值
- ❌ **不确定性高**：随机算法，结果有波动
- ⚙️ **关键参数**：
  - `residual_threshold`：内点阈值
  - `max_trials`：迭代次数
- 📊 **推荐场景**：严重污染数据，如高通量分子筛选的批次效应

---

#### QuantileRegressor（分位数回归器）

**核心思想**：不预测均值，而是预测**条件分位数**。

**sklearn实现**：`from sklearn.ensemble import GradientBoostingRegressor`（设置loss='quantile'）

**损失函数**（$\tau$-分位数）：
$$
L_\tau(r) = \begin{cases}
\tau r & r \geq 0 \\
(\tau - 1)r & r < 0
\end{cases}
$$

**特点**：
- ✅ **对异常值不敏感**：中位数回归（$\tau=0.5$）特别鲁棒
- ✅ **不确定性量化**：通过多个分位数（如0.1, 0.5, 0.9）给出预测区间
- ⚙️ **关键参数**：`quantile`（目标分位数，0-1之间）
- 📊 **推荐场景**：关注分布尾部，如毒性阈值预测

---

### 1.4 广义线性模型家族

当响应变量不服从正态分布时，广义线性模型（GLM）通过**链接函数**将线性预测器映射到合适的分布空间。

#### PoissonRegressor（泊松回归器）

**适用场景**：**计数数据**（非负整数）。

**sklearn实现**：`from sklearn.linear_model import PoissonRegressor`

**模型**：
$$
\log(\mu) = \mathbf{w}^T\mathbf{x}
$$
$$
y \sim \text{Poisson}(\mu)
$$

**特点**：
- ✅ **适合计数型目标**：如分子中特定基团数量
- 📊 **推荐场景**：预测计数型分子属性

---

#### GammaRegressor（伽马回归器）

**适用场景**：**正偏态连续数据**（如溶解度、半衰期）。

**sklearn实现**：`from sklearn.linear_model import GammaRegressor`

**模型**：
$$
\log(\mu) = \mathbf{w}^T\mathbf{x}
$$
$$
y \sim \text{Gamma}(\mu, \alpha)
$$

**特点**：
- ✅ **适合右偏数据**：常见于物理化学性质
- 📊 **推荐场景**：药代动力学参数（清除率、分布体积等）

---

#### TweedieRegressor（特威迪回归器）

**适用场景**：**混合分布**（包含零值和正连续值）。

**sklearn实现**：`from sklearn.linear_model import TweedieRegressor`

**特点**：
- ✅ **灵活分布**：通过 `power` 参数调整分布形态
  - `power=0`：正态分布
  - `power=1`：泊松分布
  - `power=2`：伽马分布
  - `1<power<2`：复合泊松-伽马分布
- 📊 **推荐场景**：高通量分子筛选数据

---

### 1.5 概率线性模型

#### BayesianRidge（贝叶斯岭回归器）

**核心思想**：将权重视为随机变量，使用概率分布表示不确定性。

**sklearn实现**：`from sklearn.linear_model import BayesianRidge`

**贝叶斯推断**：
$$
p(\mathbf{w}|\mathcal{D}) \propto p(\mathcal{D}|\mathbf{w})p(\mathbf{w})
$$

**特点**：
- ✅ **自动确定正则化强度**：无需手动调参
- ✅ **提供不确定性估计**：预测值带置信区间
- 📊 **推荐场景**：小样本、需要置信度的药物早期回归研究

---

#### ARDRegressor（自动相关性判定回归器）

**核心思想**：为每个特征赋予独立的精度参数，自动判定特征相关性。

**sklearn实现**：`from sklearn.linear_model import ARDRegression`

**特点**：
- ✅ **极致特征选择**：不相关特征的权重精确为0
- ✅ **贝叶斯框架**：自动正则化
- 📊 **推荐场景**：超高维稀疏分子描述符数据

---

### 1.6 线性模型家族综合对比

| 模型 | sklearn实现 | 核心优势 | 局限性 | 计算复杂度 | 推荐场景 |
|------|-------------|---------|-------|-----------|----------|
| **LinearRegression** | `LinearRegression` | 简单快速，解析解 | 对异常值敏感 | $O(n^3)$ | 基准回归模型 |
| **Ridge** | `Ridge` | 缓解共线性，防止过拟合 | 需要调参$\alpha$ | $O(n^3)$ | 高维数据 |
| **Lasso** | `Lasso` | 自动特征选择，稀疏解 | 特征间相关性高时表现差 | $O(n^3)$ | 特征选择 |
| **ElasticNet** | `ElasticNet` | 兼顾L1/L2优势 | 需要调参$\alpha, \rho$ | $O(n^3)$ | 复杂数据通用 |
| **SGDRegressor** | `SGDRegressor` | 内存高效，支持在线学习 | 需要调参学习率 | $O(n)$ | 大数据流 |
| **BayesianRidge** | `BayesianRidge` | 自动正则化，提供不确定性 | 计算较慢 | $O(n^3)$ | 小样本 |
| **ARDRegressor** | `ARDRegressor` | 极致特征选择 | 仅适用稀疏数据 | $O(n^3)$ | 超高维稀疏 |
| **HuberRegressor** | `HuberRegressor` | 对中等异常值鲁棒 | 需要调参$\epsilon$ | $O(n^3)$ | 含离群点 |
| **TheilSenRegressor** | `TheilSenRegressor` | 极强鲁棒性 | 仅适用低维，计算慢 | $O(n^2)$ | 严重污染数据 |
| **RANSACRegressor** | `RANSACRegressor` | 可处理>50%异常值 | 结果不稳定，随机性 | $O(k \cdot n^2)$ | 严重污染数据 |
| **QuantileRegressor** | `QuantileRegressor` | 预测分位数，不敏感 | 计算慢，需调参$\tau$ | $O(n \cdot \log n)$ | 需要预测区间 |
| **PoissonRegressor** | `PoissonRegressor` | 适合计数数据 | 仅适用非负整数 | $O(n^3)$ | 分子计数属性 |
| **GammaRegressor** | `GammaRegressor` | 适合正偏态数据 | 仅适用正连续值 | $O(n^3)$ | 物理化学性质 |
| **TweedieRegressor** | `TweedieRegressor` | 灵活分布形态 | 需要调参power | $O(n^3)$ | 高通量筛选 |

---

## 2. 支持向量机

### 2.1 核心思想

支持向量回归（SVR）通过**最大间隔回归**来拟合数据，通过核函数可处理非线性回归问题。

### 2.2 SVR（支持向量回归器）

**原理**：容忍预测值在真实值 $\pm \epsilon$ 范围内的误差，只惩罚超出此范围的样本。

**sklearn实现**：`from sklearn.svm import SVR`

**数学公式**：
$$
\min_{\mathbf{w}} \frac{1}{2}\|\mathbf{w}\|^2 + C\sum_{i=1}^{m}\max(0, |y_i - \mathbf{w}^T\mathbf{x}_i| - \epsilon)
$$

**核函数技巧**：
- **Linear Kernel**：$K(\mathbf{x}_i, \mathbf{x}_j) = \mathbf{x}_i^T\mathbf{x}_j$
- **RBF Kernel**：$K(\mathbf{x}_i, \mathbf{x}_j) = \exp(-\gamma\|\mathbf{x}_i - \mathbf{x}_j\|^2)$
- **Polynomial Kernel**：$K(\mathbf{x}_i, \mathbf{x}_j) = (\gamma\mathbf{x}_i^T\mathbf{x}_j + r)^d$

**特点**：
- ✅ **处理非线性**：RBF核可拟合复杂回归关系
- ✅ **鲁棒性强**：对异常值不敏感（$\epsilon$-insensitive）
- ❌ **训练缓慢**：大数据集（>10^4样本）计算成本高
- ⚙️ **关键参数**：
  - `C`：正则化强度（越大越拟合训练数据）
  - `gamma`：RBF核的宽度（越大越关注近邻样本）
  - `epsilon`：容忍误差范围
- 📊 **推荐场景**：复杂非线性分子性质预测

---

## 3. 近邻方法

### 3.1 K-Nearest Neighbors（K近邻回归器）

**核心思想**：预测值由距离最近的 $k$ 个样本的平均值决定。

**sklearn实现**：`from sklearn.neighbors import KNeighborsRegressor`

**数学公式**（回归）：
$$
\hat{y} = \frac{1}{k}\sum_{i \in \mathcal{N}_k(\mathbf{x})} y_i
$$

其中 $\mathcal{N}_k(\mathbf{x})$ 是距离 $\mathbf{x}$ 最近的 $k$ 个样本集合。

**距离度量**：
- **Euclidean Distance**：$d(\mathbf{x}_i, \mathbf{x}_j) = \|\mathbf{x}_i - \mathbf{x}_j\|_2$
- **Manhattan Distance**：$d(\mathbf{x}_i, \mathbf{x}_j) = \|\mathbf{x}_i - \mathbf{x}_j\|_1$

**特点**：
- ✅ **零训练时间**：惰性学习，无需训练过程
- ✅ **天然处理非线性**：基于局部信息
- ❌ **预测缓慢**：需要计算与所有训练样本的距离
- ❌ **对特征缩放敏感**：建议先标准化
- ⚙️ **关键参数**：
  - `n_neighbors`：近邻数量（通常5-15）
  - `weights`：`uniform`（等权）或 `distance`（距离加权）
- 📊 **推荐场景**：小数据集快速baseline，分子相似性搜索

---

## 本篇小结

**第一篇介绍了机器学习回归的基础模型**：

✅ **线性模型家族**：从基础的线性回归、岭回归、Lasso，到鲁棒回归（Huber、TheilSen、RANSAC）、广义线性模型（Poisson、Gamma、Tweedie）和概率模型（BayesianRidge、ARD），形成了完整的线性回归工具箱

✅ **支持向量回归**：通过核函数处理非线性关系，是小样本高维数据的经典选择

✅ **近邻方法**：基于相似性的简单有效算法，零训练时间，适合快速原型

这些基础模型具有以下共同优势：
- **训练速度快**：适合快速建立baseline
- **高度可解释**：特别是线性模型，权重清晰可见
- **理论成熟**：经过数十年验证，稳定可靠

**下一篇**将介绍实战中最常用的**树模型与梯度提升**方法，包括随机森林、XGBoost、LightGBM等竞赛级模型。

