---
title: "把小活交给AI：一份可委托的科研与事务任务清单"
date: "2026-11-24"
last_modified_at: "2026-11-24"
tags: [writing, ai-workflow, automation, research-productivity, task-delegation]
description: "整理适合交给AI独立完成、人工只需短暂验收的科研与事务小活，涵盖代码交叉验证、文献核查、集群管理和个人事务"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/4K_1080P_compressed/074934fn3HH.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/4K_1080P_compressed/074934fn3HH.jpg"
author: Xufan Gao
lang: zh-CN
---

我觉得还可以换个思路：找一些**平时值得做、但又不值得自己专门花半天时间做的小活**。最好不需要你持续参与决策，也不涉及新的科研方向，Codex做完后花几分钟验收即可。

除了前面说过的代码审计、论文整理和套磁准备，以下这些可能更适合你。

## **一、科研相关的小活**

| 任务 | 具体做法|
| ---------- | ------------------------------------------------------------------- |
| MD分析脚本交叉验证 | 用MDAnalysis、VMD、CPPTRAJ分别计算同一轨迹的RMSD、Rg、SASA、氢键等，检查自己的脚本是否存在定义或单位偏差 |
| 建立标准测试体系 | 整理几个小型蛋白、小分子、聚合物和膜体系，作为以后开发程序时通用的测试输入 |
| 力场参数对照表| 整理CHARMM、AMBER、GAFF2、OPLS、Martini中常见官能团的参数、适用范围和注意事项|
| MD软件功能对照 | 系统整理GROMACS、NAMD、Amber、PLUMED支持的算法、GPU加速范围及兼容性|
| 科研绘图工具箱| 把过去的自由能面、热图、RDF、PMF、相关性散点图等整理成独立可调用的绘图函数|
| 统计方法小工具| 实现block averaging、自相关时间、bootstrap置信区间等常用分析方法，并编写测试和示例 |
| 文献引用核查 | 遍历文章和笔记中的参考文献，通过Crossref检查DOI、题目、作者、年份是否正确|
| 文献追踪 | 检查以前收藏的预印本是否正式发表，是否更新了版本，是否已经公开代码 |
| 已有数据的可视化 | 将以前没有充分利用的模拟分析结果做成可交互HTML，方便观察数据和寻找规律 |
| 小型算法复现 | 复现一些以前读过、感兴趣但没有亲自实现过的分子表征、轨迹分析或统计力学算法 |

其中我尤其推荐第一项。你自己写过大量MD分析程序，而科研代码中最难发现的往往是程序能够正常运行、结果却存在细微的科学定义偏差。科学软件开发指南也特别强调数值容差、独立实现对照和回归测试。

文献核查也很适合自动化，[Crossref REST API](https://www.crossref.org/documentation/retrieve-metadata/rest-api/?utm_source=chatgpt.com)可以直接根据DOI获取论文元数据，不需要让AI逐篇手动搜索。

## **二、结合你最近的服务器管理工作**

这些工作比较琐碎，但做完以后可能真的有用。

| 任务| 具体做法 |
| ----------- | ------------------------------------------------------------ |
| 集群硬件档案| 整理各节点的CPU、GPU、NUMA、内存、硬盘、驱动、CUDA版本，生成统一网页|
| Slurm历史任务统计 | 分析过去一年的GPU使用率、任务失败率、排队时间和不同用户的资源使用情况 |
| 软件环境盘点| 检查各服务器上重复安装的Conda环境、Python包、编译软件和占用空间|
| 历史存储用量分析| 利用你以前的`monitor_disk_size`记录，分析个人及集群存储增长趋势|
| 磁盘冗余文件盘点| 只读检查重复轨迹、重复压缩包、过期checkpoint等，生成清单，不自动删除|
| SSH配置整理 | 生成节点连接配置、跳板机关系图、常用连接命令和故障排查文档|
| 集群新人手册| 将软件安装、Slurm提交、GPU使用、存储规范、常见报错写成组内教程|
| 危险命令检查工具| 对涉及`rm`、`dd`、`mkfs`、`badblocks`、`rsync --delete`的操作设计检查和确认机制 |
| Shell脚本体检 | 利用ShellCheck检查变量引用、路径处理、异常退出、误覆盖等潜在问题|

考虑到你最近经历的硬盘数据恢复问题，**开发一个Linux高风险操作预检工具**其实非常适合你。可以检查设备UUID、挂载状态、源目录和目标目录、文件系统状态，并提示命令是否可能覆盖现有数据。默认只读，不实际执行危险命令。

这个工具规模不必很大，但以后管理服务器可能经常用到。

## **三、个人效率和日常事务**

| 任务 | 具体做法|
| ---------- | ------------------------------------------------------------- |
| 半月报自动归档| 将历年半月报按时间和课题归类，生成各项目的进展时间线|
| 科研工作量统计| 根据Git历史、半月报和项目日志，整理每年实际开展过哪些工作|
| 博士期间成果时间轴| 论文、专利、软著、会议、合作研究、学生工作统一整理 |
| 个人文件搜索工具 | 给本地科研文件建立全文索引，支持按关键词、课题和时间检索|
| 微信公众号文章归档| 整理自己写过的文章、图片、参考文献、发布日期和阅读数据 |
| PDF批量整理| 根据论文元数据统一命名，识别重复文献、缺失SI和错误归类|
| 个人网站完善 | 增加科研项目页面、论文列表、软件演示页面和英文介绍 |
| 学院英文新闻写作工具 | 从过去修改过的英文新闻稿中提炼术语表、写作规范和事实核查流程|
| 科研汇报素材库| 把历年PPT中的流程图、技术示意图、代表性结果整理出来，方便复用|
| AI生成文字检查工具 | 根据你的anti-defensive-writing和humanizer-zh要求，检查Markdown中的常见AI式表达 |

半月报归档尤其值得做。你已经持续进行了几年科研汇报，这些材料实际上是很好的科研历史记录。将来写毕业论文、个人陈述、申请材料和年度总结时，可以直接检索，而不需要重新回忆。

## **四、开发一些自己真正会用的小软件**

这类工作比较有趣，也不必承担发表论文的压力。

| 小工具| 功能|
| ------------------ | ------------------------------------------- |
| MD Inspector | 自动识别模拟目录中的软件、力场、体系规模、模拟时长、完成进度|
| MD Compare | 比较两个MD输入目录的参数差异，解释哪些差异可能影响结果|
| Trajectory Preview | 自动生成轨迹预览动画及基本统计图|
| Slurm Doctor | 输入Job ID，自动分析作业失败原因 |
| Research Diff| 比较同一篇文章不同版本中数据、结论和图表的变化 |
| Paper Watch| 根据已有文献列表检查正式发表、代码发布及更正情况|
| Figure Inspector | 检查论文图片中的字体、DPI、尺寸、裁切、单位和图例一致性 |
| Markdown Doctor| 检查博客中的公式渲染、图片路径、死链、Markdown格式及中英文间距 |
| Prompt Library | 整理自己常用的AI提示词，支持分类、检索和版本管理 |
| Environment Doctor | 检查CUDA、PyTorch、MPI、编译器、GROMACS、PLUMED的环境兼容性 |

## **五、把重复工作做成Codex Skills**

我觉得这是一个你之前可能还没有充分利用的方向。与其每次重新写一大段prompt，不如让Codex将你已经成熟的工作流程封装成可重复调用的Skill。

例如：

- `review-scientific-writing`：科研文章专业性、逻辑和AI味检查。
- `review-md-setup`：MD体系及参数检查。
- `review-ml-benchmark`：数据泄漏、划分、随机种子和评价指标检查。
- `review-paper-blog`：对照论文原文审查公众号文章。
- `prepare-software-copyright`：软件著作权材料准备。
- `prepare-biweekly-report`：从项目日志提取进展和困难。
- `diagnose-hpc`：集群故障的只读诊断。
- `prepare-handoff`：整理长时间agent工作的上下文和交接文档。

OpenAI官方也建议将反复出现的工作流程封装成Skills，并将项目特定的约束、测试方法和验收要求放入`AGENTS.md`。

[image](https://www.google.com/s2/favicons?domain=https://learn.chatgpt.com\&sz=32)

官方参考：Codex Best Practices、Iterating Development Workflows with Codex。

如果今天只是想消耗掉剩余的额度，我个人最推荐五件事：

1. 归档历年半月报，生成科研项目时间轴。
2. 开发Linux高风险命令预检工具。
3. 整理自己的MD分析脚本，并开展跨软件数值验证。
4. 建立科研绘图工具箱和历史PPT素材库。
5. 将常用的科研、写作、编程流程封装为Codex Skills。

这些事情做完就能放下，不会再给你增加一个需要长期维护的科研课题。

我会先做4个技能：**聊天归档、论文审校、科研图表核对、机器学习评估审计**。仓库里反复出现的是同一套核对工作：确认材料版本，追到原始数据或代码，最后检查交付物。这些步骤跨项目复用的价值最高。

**优先做**

1. **`prompt-archive` 聊天记录归档**。给一份导出文件，结合时间、回复关系和上下文判断项目归属，把有用的原句插入现有主题，核对遗漏、重复和Markdown结构。交付是整理后的项目笔记，以及可复用要求的简短归纳。可以直接沿用本仓库的归类规则和`markdown_formatter.py`。 [整理规范](/mnt/e/github-repo/prompts/instructions.md:35) · [格式脚本说明](/mnt/e/github-repo/prompts/tools/README.md)
2. **`paper-check` 论文与投稿前审校**。沿用笔记里已经提到的`paper-check`，优先补上三件事：逐项核对正文、图表、SI和原始结果，提取导师DOCX里的批注与修订并对照Markdown，按目标期刊作者指南检查投稿材料。交付按严重程度排序的修改清单，每条附位置和依据。 [OP基元投稿记录](/mnt/e/github-repo/prompts/04_OP基元.md:765) · [通用投稿要求](/mnt/e/github-repo/prompts/34_通用研究协作与项目理解.md:235)
3. **`research-figure` 科研图表生成与核对**。“重画图3”“把几个子图拼到SI”“核对图6至图8”都能触发。先确认每个子图的数据、分析脚本和统计定义，再处理配色、单位、误差条、组图与图注。交付可重跑的脚本或配置，以及打开检查过的图。不同课题的配色和命名从各自项目读取。 [HA-OP论文图表](/mnt/e/github-repo/prompts/03_HA-OP分子模拟和分析.md:2277) · [BioMISTRA图表记录](/mnt/e/github-repo/prompts/09_metal.md:446)
4. **`ml-evaluation-audit` 机器学习评估审计**。沿着数据划分、特征预处理、CV、HPO和最终测试逐层核对，从数据库或预测文件复算指标。重点查特征筛选是否用了测试集信息、跨种子结果是否误复用、Stage 1和Stage 2是否用了同一批组合比较。交付含样本数、划分、种子、指标和异常来源的审计表。 [MolBlender Benchmark](/mnt/e/github-repo/prompts/02_Polyglotmol机器学习.md:3995) · [FIND-QSAR泄漏检查](/mnt/e/github-repo/prompts/10_FIND-QSAR.md:282)

**第二批**

- **`research-run-monitor` 长任务巡检与续跑**。跨MD、FEP和ML复用。用作业状态、资源利用率、日志和产物增长判断任务是在推进还是卡住，再依据已有输出和缓存判断哪些可以续跑。交付实际完成数、失败原因和可接续的位置。 [技术运维](/mnt/e/github-repo/prompts/33_技术运维相关.md:662) · [PRISM运行记录](/mnt/e/github-repo/prompts/05_PRISM建模和分析.md:1927)
- **`software-doc-audit` 软件手册与软著核对**。从源码、CLI帮助和实际运行结果核对手册、截图与代码清单，复用`softcopy-toolkit`编译并检查PDF。源码摘录在提交副本中处理。 [OP基元软著](/mnt/e/github-repo/prompts/04_OP基元.md:755) · [软著文档规范](/mnt/e/github-repo/prompts/34_通用研究协作与项目理解.md:86)
- **`paper-to-post` 论文转中文科研推送**。从正文、关键图和SI组织一篇可发布的简讯，逐项核对数字和结论，并按现有栏目模板处理图注与排版。它和论文审校的区别在于读者和交付物，值得独立触发。 [推送写作记录](/mnt/e/github-repo/prompts/01_推送写作相关.md:957)

Beamer已经有内容规划、布局、编译和质量检查技能，我会直接复用，并把“逐图比对原文、逐页检查PDF”补进现有流程。中文标点、加粗和空行规则适合留在`AGENTS.md`，字数统计、DOCX导出和组图适合复用脚本。会话恢复与多会话协作在笔记里已有`session-hygiene`和`claude-task-runner`的记录，优先沿用这些入口。 [Beamer现有技能](/mnt/e/github-repo/IQB-JC-master/.claude/skills) · [会话协作记录](/mnt/e/github-repo/prompts/34_通用研究协作与项目理解.md:173)