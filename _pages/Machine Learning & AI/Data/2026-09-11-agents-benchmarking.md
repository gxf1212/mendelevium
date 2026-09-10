---
title: "让 Agent 替你跑评测：三个化学基础模型在 ADMET 数据集上的头对头对比"
date: "2026-09-11"
last_modified_at: "2026-09-11"
tags: [agent, benchmark, admet, foundation-model, chemprop, tabular-foundation-model, machine-learning, virtual-screening]
description: "借助 Claude 与一组 agent，无需手动 clone 任何仓库，即可在 OpenADMET ExpansionRx 与 Biogen 两个 ADMET 数据集上头对头评测 MEGA-CL、Monroe、Mol-JEPA 三个化学基础模型，并用 Tukey HSD 与数据泄露检验严谨对比。"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/Wallpaper_compressed/nature-3616194_1920.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/Wallpaper_compressed/nature-3616194_1920.jpg"
author: Xufan Gao
lang: zh-CN
---

> 本文翻译自 Pat Walters 的博客 [Let the Agents Do the Benchmarking](https://patwalters.github.io/Let-the-Agents-Do-the-Benchmarking/)，原文发表于 2026-08-29。以下为直译，版权归原作者所有。

## 一个全新的世界

在计算工作中，很少有事情比读了一篇精彩的论文、从 GitHub 上拽下代码、再在可信数据集上跑一遍验证更让人满足。然而，那份最初的热情常常被现实的麻烦冲散：大把时间消耗在调试 CUDA、PyTorch、PyTorch Lightning 之间此起彼伏的版本冲突，以及一张没完没了的依赖网上。即便代码真跑通了，论文的价值也常被对有缺陷数据集（比如 TDC 或 MoleculeNet）的依赖，以及缺乏严格的统计比较所削弱。

幸好，那些令人沮丧的日子已经彻底成为过去，也希望对你而言同样如此。借助我的朋友 Claude 以及它手下的一组 agent，我最近做 benchmark 连一个 Git 仓库都没手动 clone 过。这种行云流水般的方式与整体成效，着实让我吃惊。

为了启动 benchmark 流程，我先指定目标数据集。鉴于我聚焦 ADMET 建模，我选了近期 OpenADMET ExpansionRx Blind Challenge 的数据集，以及 Cheng Fang 及其同事在 2023 年论文里提供的 Biogen 数据集。几个关键因素让它们成为出类拔萃的评测选择。

这些数据并非像 Franken-dataset（弗兰肯数据集）那样从几十篇互不相干的文献拼凑而来，而是由同一实验室的同一批科学家统一产出。我对那些整理早期文献集、采集原始测量值的先驱绝无冒犯之意；他们的贡献对让这个领域起步至关重要。但我们现在完全有能力达到更高的标准。都 2026 年了，还依赖过时又有缺陷的基准数据，完全说不过去。想和更老的论文保持一致也不是正当理由——就像我妈以前常说的：「要是你朋友跳桥，你也跟着跳吗？」

这些化合物覆盖的化学空间，忠实地反映了标准的药物发现场景。同时纳入 ExpansionRx 与 Biogen 两套基准集，让我们能评估两种截然不同的真实世界情境：ExpansionRx 数据集来自活跃的发现项目、高度同系物（congeneric），是 hit-to-lead 与 lead 优化阶段的准确写照；而 Biogen 集多半由商业筛选分子构成，捕捉了早期发现阶段的特征。

这些数据集拥有真实的动态范围。文献 benchmark 里最让我长期不爽的一点，就是它们夸张得离谱的人为量程；那些声称能预测跨越十几个数量级的水溶性的论文，至今还在发表，实在让我费解（说真的，这到底是谁在审稿？）。在如此夸张的数据上声称强性能，比在浴缸里钓鱼还假。

## 让评测跑起来

为了确保 benchmark 流程遵循严格的统计规范，在选定基线数据集之后，我让 Claude 参考了我们 2025 年的论文《Practically significant method comparison protocols for machine learning in small molecule drug discovery》，以及几篇相关博客。我请 Claude 搭建一个包含四个候选模型的初始套件：用 LightGBM 的 Morgan 指纹基线、ChemProp 单任务、ChemProp 多任务，以及接入 CheMeleon 基础模型的 ChemProp。Claude 自行定位了所需的代码仓库与文献。虽然它起初拉到的是 CheMeleon 预印本，但很快更新为近期发表的 JCIM 论文。

我没有依赖自己的笔记本，而是让 Claude 把 benchmark 的计算负载卸载到地下室的 Linux 服务器上以加速。为了做会话编排、并在合上笔记本时不中断任务，我依赖 herdr——一个为 agent 驱动的工作流量身打造、受 tmux 启发的现代工具。如果你还没试过 herdr，我强烈推荐你去看看（我的朋友们已经听我反复安利听到烦了！）。

我们用论文里写的 5×5 交叉验证协议评估了每个方法。每个数据集都有固定的留出测试集，25 个复现模型来自对训练分子做 5 次五折交叉验证、按聚类分组，于是每个方法在每一折都看到完全相同的训练分子，并在同一个未被触碰的测试集上打分。划分设定如下：

- **ExpansionRx 数据集**：用挑战赛自带的 train/test 划分，5,326 训练 / 2,282 测试，70/30。
- **Biogen 数据集**：数据本身不带划分，于是把整个 BitBIRCH 聚类留出，直到测试集达到同样的 30%。
- **基础模型**：CheMeleon 与 MEGA-CL 在训练折上微调；Monroe 与 Mol-JEPA 冻结编码器、做上下文内（in-context）预测，完全不做下游训练。

到此为止做的这些都不错，但算不上特别激动人心。过去几周，我们见证了三个新的化学基础模型登场。

### MEGA-CL（arXiv:2607.24314）

走图对比学习路线。它把一个增强版 GCN+ 消息传递骨干（加了残差连接和层归一化）与多头图外部注意力模块配对。模型在约 1 亿分子上用 NT-Xent 对比学习目标预训练，再用作者发布的 checkpoint 对每个 endpoint 微调。由于该架构为单目标任务设计，每个 endpoint 都需要单独训练一个模型。

### Monroe（arXiv:2608.18982）

是一个 58.5M 参数的图 transformer，建立在 GRIT 架构之上。它在 1,152 个同步任务上做了广泛预训练：为 8,100 万分子预测 62 个量子化学性质（通过 PM6）、1,089 个来自 PCBA 的二值生物 assay，外加一个构象去噪目标。它的图表示很特别——额外加边来编码立体化学构型，因而能区分立体异构体。下游任务里编码器保持冻结：每个分子被映射成一个 720 维向量，TabPFN 用单次前向直接从这些表示做上下文内（in-context）预测。

### Mol-JEPA（arXiv:2608.22642）

约 50M 参数，来自 Boehringer Ingelheim、Tübingen 大学、Brown 大学与 UT Austin。它放弃了常见的「扰动结构」自监督路线（作者认为这不适合化学），改用联合嵌入预测架构（JEPA, joint-embedding predictive architecture）：在 14 类分子数据（图、ECFP/MOE 描述符、xTB/DFT 计算、各类实验标签）上掩盖整个模态，用一个 transformer 从其余模态预测缺失的潜在表示，训练用 469 万分子。推理时模型只需要 SMILES 字符串；冻结的 512 维 CLS token（[CLS] 令牌）随后交给 TabICL 做下游任务。

Monroe 与 Mol-JEPA 属于同一方法论类别，与 MEGA-CL 区分开来。前二者冻结编码器、靠上下文内表格模型做适应，而 MEGA-CL 对整个网络做针对每个 endpoint 的完整微调。

在三个模型的所有 case 里，我都只是提供了预印本 PDF，让 Claude 去下载、测试代码，然后跑和之前方法相同的 benchmark。整个过程竟如此轻松，让我惊讶。过去，benchmark 一个方法要花我好几天。我通常得花几小时只为解决代码跑通，再花更多时间搞清楚输入、跑代码、写分析脚本。而这次，我给了几条 prompt、上床睡觉，醒来就收到一份漂亮的报告。说实话，中间出了一个小故障，但那不是 Claude 的错。下载 MEGA-CL 的 GitHub 仓库后，Claude 指出模型权重缺失。我翻了仓库和 Zenodo 档案，果然确实缺失。我邮件联系了 MEGA-CL 的作者，他们很快把模型权重补进了仓库。谢谢 Jinfeng Liu！

## 真章在报告里

我让 Claude 把结果写成一份 artifact。我不得不承认，Claude 的报告比我凭自己能力写出的任何东西都好。它包含我钟爱的 Tukey 诚实显著性差异（HSD, honestly significant difference）图，Claude 还做了一版更好的「可怕的加粗表」——不只是把「最佳」加粗，而是标出哪些结果在统计上等价。更妙的是，全部分析代码与数据都在独立的 GitHub 仓库里保持更新，我凭一个 URL 就能分享报告。

https://claude.ai/code/artifact/fa406d3a-3b09-4e17-9e1c-ecd2d92903d2

此时此刻，你一定迫不及待想知道哪个基础模型表现最好；跨两套数据集，都有一个清晰的赢家。报告很长，我力劝你去看看。为了简洁起见，我只分享几张关键图。先看 ExpansionRx 数据集中 9 个 endpoint 的平均绝对误差（MAE, mean absolute error）。既然看的是 MAE，越低越好。对不太熟悉这些图的人说明一下：最佳方法标蓝，其他统计等价的方法标灰。若某方法显著差于最佳，它的点标红。评测 ExpansionRx 时，两个明显的赢家浮现：Monroe + TabPFN 拿下 LogD、两个微粒体稳定性 endpoint、两个 Caco-2 endpoint；而 ChemProp + CheMeleon 拿下全部三个蛋白/组织结合 endpoint 加水溶性（LogS）。胜负按 assay 切分，而非按数据量多少。LogS 的训练数据比任何 endpoint 都多，5,128 条，ChemProp + CheMeleon 仍拿下它，Monroe + TabPFN 紧随其后。同样值得注意的是榜上缺席者：MEGA-CL 与 Mol-JEPA 在全部九个 endpoint 上，既无一最佳、也无并列最佳。

![ExpansionRx 九个 endpoint 的 MAE Tukey HSD 图](https://patwalters.github.io/assets/images/2026_08_AGENTS/tukey_mae_expansion.png)

Biogen 数据集上，冠军更没有悬念。每一个 case，Monroe + TabPFN 都压过其他模型。ExpansionRx 那种按 assay 切分的格局在这里不复现。两个血浆蛋白结合 endpoint 本该是 ChemProp + CheMeleon 的强项，却恰恰是它在整页里输得最惨的地方。

![Biogen 六个 endpoint 的 MAE Tukey HSD 图](https://patwalters.github.io/assets/images/2026_08_AGENTS/tukey_mae_biogen.png)

我也很喜欢 Claude 对「可怕的加粗表」的替代方案：不是只把「最佳」方法加粗，而是把统计等价的方法标为「并列（tied）」。

![九个 ExpansionRx endpoint 的 MAE，标注最佳与并列方法](https://patwalters.github.io/assets/images/2026_08_AGENTS/mae_table_expansion.png)

有一个 caveat 值得专门 flag，而且是 Claude 在报告里主动提出的，不是我。Monroe 的作者本人在 ExpansionRx 上也跑过自己的模型，那自然要问：模型是不是见过这些测量值？答案是没有，而且论文写得足够具体，能定下这个结论。Monroe 在 1,152 个具名任务上预训练：62 个来自 PM6 的半经验量子性质、1,089 个来自 PCBA 的二值生物 assay、一个构象去噪目标——没有一个是 ADME endpoint。折、变换、测试集都是我们的，也没有任何 Monroe 超参在它们上做过 tuning。

Monroe 作者 Blazej Banaszewski 走得更远，直接查了分子本身，发现零标签重叠。恰好有一个 ExpansionRx 测试分子（带 LogD 与 LogS 值）出现在 PM6 里。Biogen 集重叠重得多：约 53% 的测试分子出现在 PM6，约 8% 出现在 PCBA。对一份 8,100 万分子的语料来说，覆盖商业来源 ADME 集的一半并不意外，关键是这些结构挂了什么标签——PM6 给的是与这些 assay 无关的量子化学描述符；PCBA 里少数生物 assay 与 Biogen endpoint 生物相邻，但是不同的实验、报不同的标签。Monroe 见过不少这些分子，但从没见过它们测的是什么。

Mol-JEPA 是唯一外人能独立核查的 case，因为它的作者公开了完整的预训练表——466 万行，每个分子带一个 InChIKey。问题在这里也最尖锐：14 个模态中有两个是来自 ChEMBL、PCBA、Therapeutic Data Commons（TDC）的实验标签向量，而 TDC 确实含 ADME 任务。但 7,608 个 ExpansionRx 分子按精确 InChIKey 无一出现在该表里。Biogen 是公开集，其中 59 / 3,521 个分子出现（1.7%），经 nabla-DFT、PubChem BioAssay、TDC 进入，另有 132 个共享连接块。InChIKey 命中不等于 Mol-JEPA 见过 Biogen 数据——事实是它没有：这些行没有一行带 Biogen 测量值，为这些 endpoint 命名的 6 列在全部 466 万行里都是空的。结构略有重叠，标签从不重叠。

## 按下「快捷键」

直到今年春天，我还坚定地站在 LLM 怀疑论阵营。随着 Claude Opus 4.8 与 GPT-5.5 的出现，这一立场剧烈改变。我对药物化学问题的回答，从拿到一堆胡言乱语，变成了被那些让我心想「我怎么没想到？」的回答真正震撼。如今，搭建新协议、跑 benchmark，简单得就像写一条 prompt。我不再花好几天评估一个新模型才把它加进工作流，而是设好 prompt、上床过夜，醒来就有一份全面的报告。

这项工作，加上近期的 OpenADMET PXR Blind Challenge，带来一个关键洞见：表格基础模型（TFM, Tabular Foundation Model）的威力惊人。把 TFM 与来自 Monroe、CheMeleon 这类基础模型的预训练表示配对，能产出真正亮眼的表现。想快速上手这些技术的人，我强烈推荐去看 Karim Ben Hicham 的 webinar，它对 TFM 及其在分子性质预测中的角色做了极友好的概览。

如我所愿地展示的那样，严谨地搭建并运行 benchmark，如今已变得不可思议地简单。既然已经这么简单，就看大家会不会去做了。如果他们不做，我的 agent 和我来做。小心了！

本次分析用到的全部代码与 prompt 都发布在 GitHub 上。

https://github.com/PatWalters/expansion-ml-comparison
