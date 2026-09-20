---
title: "Gen-COMPAS 方法详解：DDPM、VCN、TMD 与 RiteWeight"
date: "2026-09-14"
last_modified_at: "2026-09-14"
tags: [generative-sampling, molecular-dynamics, enhanced-sampling, diffusion-model, committor, variational-committor, reweighting, technical-appendix]
description: "Gen-COMPAS 的技术附录：去噪扩散概率模型、变分 committor 网络、靶向分子动力学和 RiteWeight 重加权的公式与实现要点"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/4K_1080P_compressed/081247t9D67.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/4K_1080P_compressed/081247t9D67.jpg"
author: Xufan Gao
lang: zh-CN
---

# Gen-COMPAS 方法详解：DDPM、VCN、TMD 与 RiteWeight

本文是[《生成式采样打破时间尺度壁垒：无需预设反应坐标的构象转变模拟新框架》](2026-09-14-gen-compas.md)的方法附录，集中说明 Gen-COMPAS 的四个核心模块：去噪扩散概率模型（DDPM）、变分 committor 网络（VCN）、靶向分子动力学（TMD）和 RiteWeight 重加权。正文只保留了算法流程和直观解释，这里给出更完整的公式和实现细节，供想复现或深入理解的读者参考。

## 去噪扩散概率模型（DDPM）

生成组件是一个在**中心化的三维原子坐标**上操作的 DDPM。

通俗地说，扩散模型做的是**洗照片的反过程**。训练时我们拿一张真实的蛋白质结构，一步步往里掺高斯噪声，直到它彻底变成一团无意义的随机坐标，模型要学的就是一个去噪任务，给定任意一个带噪版本，把它还原成干净结构。推理时反过来，从纯噪声出发，模型一步步去噪，最后洗出一张几何上合理的中间结构。这一步只是构造训练目标，模拟本身用不到。

前向过程按照方差调度逐步添加高斯噪声，神经网络被训练从噪声构象和扩散时间预测干净坐标。在推理时，迭代去噪将高斯噪声映射为候选中间构象。

#### 前向扩散过程

给定数据分布 $\mathbf{x}_0 \sim p(\mathbf{x}_0)$，前向马尔可夫过程生成随机变量 $\mathbf{x}_1, \ldots, \mathbf{x}_{t_T}$：

$$
p(\mathbf{x}_1, \ldots, \mathbf{x}_{t_T} | \mathbf{x}_0) = \prod_{t=1}^{t_T} p(\mathbf{x}_t | \mathbf{x}_{t-1}), \quad p(\mathbf{x}_t | \mathbf{x}_{t-1}) = \mathcal{N}(\mathbf{x}_t; \sqrt{1-\beta_t}\mathbf{x}_{t-1}, \beta_t \mathbf{I})
$$

其中 $\beta_t \in (0,1)$ 是方差调度超参数。时间步 $t$ 的噪声状态可直接从 $\mathbf{x}_0$ 采样：

$$
\mathbf{x}_t = \sqrt{\bar{\alpha}_t}\mathbf{x}_0 + \sqrt{1-\bar{\alpha}_t}\boldsymbol{\epsilon}, \quad \boldsymbol{\epsilon} \sim \mathcal{N}(0, \mathbf{I})
$$

其中 $\bar{\alpha}_t := \prod_{s=1}^t \alpha_s$ 是噪声调度 $\alpha_t = 1-\beta_t$ 的累积乘积。为保持平移等变性，所有坐标在加噪前通过减去质心进行中心化。作者采用**余弦方差调度**，训练稳定且高效。

为什么要把坐标中心化？蛋白质在盒子里整体的平移和旋转与构象变化无关，若不减去质心，模型会浪费容量去学这些刚体运动。中心化后，网络只关注分子自身的形状变化。

#### 反向生成过程

从纯噪声 $\mathbf{x}_{t_T} \sim \mathcal{N}(0, \mathbf{I})$ 出发，神经网络迭代去噪。反向过程由另一个马尔可夫过程描述：

$$
p_\theta(\mathbf{x}_{t-1} | \mathbf{x}_t) = \mathcal{N}(\mathbf{x}_{t-1}; \boldsymbol{\mu}_\theta(\mathbf{x}_t, t), \boldsymbol{\Sigma}_\theta(\mathbf{x}_t, t))
$$

其中 $\boldsymbol{\mu}_\theta$ 由网络预测的干净结构 $\hat{\mathbf{x}}_0 = f_\theta(\mathbf{x}_t, t)$ 参数化。训练损失为：

$$
\mathcal{L}(\theta) = \mathbb{E}_{t, \mathbf{x}_0 \sim p(\mathbf{x}_0), \boldsymbol{\epsilon} \sim \mathcal{N}(0,\mathbf{I})} \left[ \| f_\theta(\sqrt{\bar{\alpha}_t}\mathbf{x}_0 + \sqrt{1-\bar{\alpha}_t}\boldsymbol{\epsilon}, t) - \mathbf{x}_0 \|_2^2 \right]
$$

#### 去噪网络架构

核心是一个**分层图神经网络（GNN）**，包含原子级和残基级处理。原子级使用 SchNet 风格的连续滤波卷积层，残基级使用多头自注意力层捕获长程相互作用。动态 k-NN 图基于当前坐标构建，与共价键静态图结合定义完整的相互作用网络。蛋白质和配体坐标联合去噪，而不是分别生成再对接。

为什么用图网络而不是普通卷积？蛋白质不是规则网格，原子之间的关系由化学键和空间邻近决定。图网络天然适配这种不规则拓扑，动态 k-NN 还能随构象变化更新谁和谁相互作用。联合去噪则保证配体在生成时就处在结合口袋的相对位置，避免生成完再对接带来的错位。

#### 采样实现

在迭代精修循环中，噪声尺度逐步降低。初始轮次使用高噪声尺度鼓励结构多样性，后续轮次逐步降低至零，产生越来越确定的高精度结构。

## 变分 committor 网络（VCN）

#### committor 函数

committor 函数 $q(\mathbf{x}_0)$ 定义为从构象 $\mathbf{x}_0$ 出发的轨迹在到达反应物 A 之前先到达产物 B 的概率：

$$
q(\mathbf{x}_0) = \mathbb{P}(\tau_B(\mathbf{x}_0) < \tau_A(\mathbf{x}_0))
$$

其中 $\tau_S(\mathbf{x}_0) = \inf\{t > 0 | \mathbf{x}(t; \mathbf{x}_0) \in S\}$ 是首次击中时间。$q=0$ 对应确定返回 A，$q=1$ 对应确定到达 B，$q=1/2$ 的超曲面（separatrix）是最优动力学分界面，定义了过渡态系综（TSE）。

#### VCN 的变分原理

VCN 基于变分原理，将标量通量 $J_{AB}$ 与 $q$ 联系起来：

$$
J_{AB}[q; \tau] = \frac{C[q; \tau]}{\tau}, \quad C[q; \tau] = \frac{1}{2} \langle (q(\tau) - q(0))^2 \rangle
$$

VCN 寻找最小化标量通量（或等价地，相关泛函 $C$）的函数。网络输出通过 sigmoid 激活函数约束在 $[0,1]$ 范围内，并添加边界损失项强制在 $(A \cup B)^c$ 与盆地边界处的连续映射。

这背后的物理很漂亮，实现上却很直接。理论告诉我们，真正的最优 committor 就是让短时间内的 $q$ 变化尽量小的那个函数，VCN 用一个神经网络去逼近它，训练目标就是最小化 $C[q;\tau]$。sigmoid 把输出压在概率区间，边界损失则保证在 A、B 盆地内部 $q$ 分别锁死为 $0$ 和 $1$，中间平滑过渡。

#### Z-matrix 内部坐标

这是 Gen-COMPAS 的一个关键设计。传统方法需要预先选择少量物理动机的 CVs，但这可能遗漏慢自由度。Gen-COMPAS 采用完整的 3N-6 个内部坐标（Z-matrix 坐标）作为 VCN 的输入：

- **完备性**：完整内部坐标表示与原始构象空间（模刚体运动）同构，保留了所有本质自由度。
- **无预设 CV 选择**：网络自行学习哪些内部自由度对转变有信息量。
- **物理可解释性**：距离、角度和二面角自然编码局部结构变化。

实际中，Z-matrix 表示以层次化方式指定键长、价角和扭转角，每个原子相对于少量参考原子定义。

为什么不用笛卡尔坐标直接喂给网络？笛卡尔坐标里混着整体的旋转和平移，这些都是与构象变化无关的噪音，网络会浪费容量去学它们。内部坐标（键长、键角、二面角）直接描述分子自身的形状，把刚体运动自然剔除在外，网络只需关注真正的构象自由度。

## 靶向分子动力学（TMD）

TMD 将系统从初始构象 $\mathbf{x}(0)$ 驱动向参考目标构象 $\mathbf{x}_{\text{ref}}$，通过偏置势强制 RMSD 按预定调度降低：

$$
V_{\text{TMD}}(t) = \frac{k_{\text{TMD}}}{2N} (\text{RMSD}(t) - \text{RMSD}^*(t))^2
$$

其中 $\text{RMSD}^*(t)$ 从初始 RMSD 线性降至零。Gen-COMPAS 中，从 A 和 B 分别向同一生成中间结构运行 TMD，系统性识别从转变两端都可及的物理合理区域。

## RiteWeight 重加权

Gen-COMPAS 生成的轨迹来自非平衡系综（许多从 TMD 端点或 separatrix 附近启动），不直接服从 Boltzmann 分布。为恢复热力学可观测量，作者使用 RiteWeight 算法：

- 基于固定时间滞后 $\tau$ 的转变对，迭代重加权轨迹片段。
- 每轮迭代使用随机聚类，构建离散状态转移矩阵 $T^{(k)}$，求解平稳分布 $\pi^{(k)}$。
- 平稳概率用于更新轨迹片段权重，收敛后得到准连续的稳态分布估计。

这一步解决的是样本权重的公平性问题。TMD 把系统摆到了特定构象附近，无偏 MD 从那里出发，所以这些轨迹的起点并不是按真实平衡分布取的，而是被 TMD 人为摆位过的。如果直接把它们当平衡样本去统计自由能景观，盆地高度会被系统性偏置。RiteWeight 的思路是不去改动轨迹，而是给每条轨迹片段重新算权重，利用这些片段之间按时间滞后 $\tau$ 互相转移的频率，反推出稳态分布，把非平衡的样本纠回平衡。没有了这个重加权，自由能景观的盆地高低就不可信。

一旦获得收敛权重 $w(x)$，沿任意 CV $\xi(x)$ 的概率密度和自由能可直接估计：

$$
P(\xi) \propto \sum_{x \in \mathcal{D}} w(x) \delta(\xi(x) - \xi), \quad F(\xi) = -k_B T \ln P(\xi) + C
$$

## 关键提醒

这些公式和实现细节直接决定了 Gen-COMPAS 能否从端点结构走到可靠的过渡态。但单纯看公式容易忽略一点：**生成的结构只是提议，必须放回真实分子哈密顿量中精修和验证**。TMD 和无偏 MD 这一步不是后处理，而是把扩散模型给出的“合理构象”锚定在物理力场上的核心环节。如果跳过这一步，生成模型可能给出漂亮但动力学上无意义的中间态，committor 和 FEL 都会跟着失真。