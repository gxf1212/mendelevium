# Mendelevium Blog - Claude Development Guide

## 项目理

This is a personal research blog sharing experiences in molecular dynamics, computational chemistry, and scientific computing. The content should be:

- **Educational**: 每篇文章都应该能帮助读者学到一些实用的知识
- **Practical**: 重点在于实用性，分享实际的方法和工具
- **Clear**: 保持简洁明了的写作风格，便于理解

## 开发原则

1. **保持简洁**: 网站结构应该简单明了，不要过度复杂化
2. **内容为王**: 专注于高质量的技术内容，而非花哨的功能
3. **分类清晰**: 按研究领域合理组织内容结构

## 技术栈

- Jekyll (Satellite theme)
- GitHub Pages
- Giscus for comments
- GoatCounter for analytics

## 博文内容组织

- 每个分类下应该有一个 index.md 文件， 所有index.md都得是---\n---，啥内容都没有 
- Archive\jekyll-theme-satellite-master是原始模板，报错了可以参考，尤其是里面的docs
- https://github.com/jekyll/jekyll-seo-tag/blob/master/docs/usage.md 是SEO的插件，可以参考
- 为什么单篇文章，如2025-08-13-deep-covboost-ai-covid-target会出现在侧边栏？不应该出现的
    bookmark: false我已经删掉了，就解决了，以后都不要写
- 文件名不能太长、有问题，包括图片文件夹。rename folder name to sth shorter.

## posts注意事项

- 图片路径使用相对路径，确保在不同环境下都能正常显示
- 保持中英文混用的写作风格，适合中文科研环境
- 定期整理和归档过时的内容
- should add frontmatter for each post (identify .md files without ---), e.g.
  ```
  ---
  title: "Random Forest and Enhanced Sampling Unite: Revealing and Correcting Ghost Errors in Alchemical Free Energy Calculations"
  date: "2025-08-22"
  tags: [random-forest, enhanced-sampling, alchemical-free-energy, gamd, error-analysis, machine-learning, molecular-dynamics]
  description: "深入分析 Boltz-2 AI 模型在配体亲和力预测中的表现，与 FEP 方法的对比，以及两者如何协同工作"
  image: "/assets/img/thumbnail/bricks.webp"
  thumbnail: "/assets/img/La-Mancha.jpg"
  author: Xufan Gao
  lang: zh-CN
  ---
  ```
- frontmatter的date: "2025-08-22"不是文章发表的时间，而是写blog的时间，最后一次修改的日期
- 能不能随机替换assets\img\thumbnail、assets\img\thumbnail_mine下的文件，雨露均沾，创建新文件的时候随机选取。**随机缩略图选择**: 已创建工具脚本来自动随机选择缩略图，避免过度使用bricks.webp。`tools/random_thumbnail.py`: 随机选择缩略图。image和thumbnail这俩参数得一样啊。


## 常见要求

- 根据@CLAUDE.md，为@_pages\Free Energy\fep-the-end-of-parameter-tuning.pdf写一篇推送Markdown文章
- 根据文件_pages中.md的最后修改时间，给文件名添加YYYY-MM-DD才能在博客上显示（rename就行），还要仿照已有的添加frontmatter，tag尽量用和其他类似的，如果有的话。已经有YYYY-MM-DD的就不用了，这种一般frontmatter也都有了。archive里面的，还有about.md，index这种不要改；Diary里面的不要搞太复杂的文件名，就YYYY-MM-DD-diary1.md这种就行，其他文件名可能还是对title的一个概括。直接扫描所有的_pages中的文件，看哪些没有YYYY-MM-DD，然后修复其frontmatter和文件名，不要搞太复杂的流程。
- 不对，不一定需要添加YYYY-MM-DD，所以给没添加frontmatter的加上frontmatter，tag尽量用和其他类似的，如果有的话。Diary里面的，archive里面的，还有about.md，index这种不要改
- 去掉所有参考标记，如[cite\_start]、[cite: 3749]，根据CLAUDE.md修复所有格式
- 注意title要用中文标点，不要用英文标点，尤其是英文引号，请修复。类似于title: "皮肤屏障的“水之道”：分子模拟揭示脂质相共存如何稳定间质水"，内部是中文标点。【OpenFE】这种是可以的。。但文件名不要出现乱七八糟的东西，标点符号等。
- 把_pages\xxx.docx改写成Markdown，内容完全不变，保证参考文献链接仍然正确，即文内的引用标记要能对得上文末的list，得用这种吧：[^6]。

### 技术记录

在_pages\archive\qq，能否根据（大概）2024年6~12月以来的内容总结出几篇技术记录，放在_pages\Techniques底下？
  - 在_pages\archive\qq\_pages\archive\qq\694077960_1757646016523.txt，内容筛选的脚本是现成的，基本不用改。先运行脚本，再读取内容，时间和图片都要过滤，先修改下时间，两个脚本都要执行。
  - 每一篇不太长，就是一篇推送的长度，包括一个或多个话题。尽量按照话题来组织汇总，把比较相关的内容放在一起。内容太多就创建新的md文件。
  - 如只看到一个网站但不知道是干啥的，可以访问该网站获取内容。一些文字如果提供了参考网址，要把网址留在下面作为reference，不是文档最下面，就是笔记原位。
  - 这种笔记整理，除了[GROMACS论坛](https://gromacs.bioexcel.eu/)，还需要把链接显式地写出来，否则微信识别不了，比如：[GROMACS论坛](https://gromacs.bioexcel.eu/)：https://gromacs.bioexcel.eu。
  - 没有代码、都是link的，不要放代码框。。正常的文字不要放quote 
  - pip cache purge这种放在代码框或``包裹
  - 不要泄露隐私内容。
  - 对于这些文档，title: 和第一个# 要加点【笔记整理|2025-07】之类的标签，文件名不用加。这个时间是笔记实际的大致时间，不是创建文档的今天。
  - 创建的文档不要保留奇怪的时间戳，[02:56] 这种
  - 不要把参考资源：这种文字类的放在```里面，没必要代码框的文字说明就不要代码框，如“目前rdkit.Chem.Draw.MolsToGridImage函数没有直接设置图例字体大小的选项”就是个经验，文字就行

## 推文

现在有了一个方便的工具来快速搜索PDF内容。我已经为你准备好了：

使用方式：`python3 tools/search_pdf_text.py <pdf_file> <keyword> [context_lines]`

  例子：
```python
# 搜索"MCS"，显示前后3行上下文
python3 tools/search_pdf_text.py "_pages/Free Energy/fep-the-end-of-parameter-tuning.pdf" "MCS" 3
```
必要的时候还是要直接读PDF全文

### Prompt

**角色**: 你是一位顶尖的科研媒体编辑和科学传播专家。你的读者是渴望学习新知识的科研新手和同行。你的任务是将一篇专业的、信息密集的科研论文，转换成一篇详尽、易懂、重点突出、排版精美的微信公众号推送文章。本地的话，用Python提取PDF文字来写推文，一步一步来写。

Markdown文件250~300行的样子，多了就拆出一篇讲技术细节和其他结果的附录md文件。是把主文档的部分内容挪过去。先填补附录.md，再删除主文档里的。附录也不用过于长了，否则要再拆一个附录文件。附录不要自己的Q&A，除非正文的Q&A放不下。

在关键的概念、结论或punchline处，请使用 `**加粗**` 来突出重点，但不要包含最后一个句号，**要放在最后一个句号之前，否则渲染不成功。**也不要包裹非常规内容，如公式，**“黑箱”**这种带引号的，**贝叶斯优化（Bayesian Optimization）**这种带括号的，应该改成“**黑箱**”、**贝叶斯优化**（Bayesian Optimization），或者不加粗。**24%**这种都改成**24**%，否则报错。就包含正常汉字、数字等就行。
- 正文（引用格式、代码之类的、frontmatter和mermaid图之外）尽量用中文标点，中文引号的Unicode是201C和201D，不要用你默认的，就用这俩，必须在unicode层级上做替换，参考tools\convert_quotes.py。尽量少用引号，常用词就没必要带引号，“蛋白”还引号就很搞笑了。mermaid图里面，A[“底物化学图”]别用中文引号，但A["“底物化学图”"]这种内部有的是可以的。
- **排版**: 全文（除了mermaid图和frontmatter）尽量使用中文标点符号，尤其是冒号、括号、中文内容的引号。除了加粗的标题，尽量让段落自然换行，避免不必要的手动换行 `<br/>`。不要有"过拟合"这种，而是“过拟合”。
- Markdown格式：列表这些的两个item之间不要有个空行，包括有序列表，其他比如上下也要适当多空行。只是修复列表内部，不删掉Markdown正常的空行。参考tools/fix_format.sh。如
  ```
  - 第一项
  - 第二项
  - 第三项
  ```
- 不要写Q\&A，而是Q&A；不要写1\.，而是1.，总之不要多加\。
    - **Mermaid使用规范**:
        - **必须使用横向布局** (必须`graph TB`，如果subgraph中有箭头，请使用 `direction LR`，否（节点间只是简单并列）则TB方向。），以创建宽屏、信息密度高的图表，避免简单的垂直长条。
        - 鼓励使用 `subgraph "子图标题"` 来组织和划分逻辑模块，**子图标题不要使用任何Markdown标记**。
        - 使用中文作为节点描述。节点内的换行请使用 `<br/>`。
        - "2. 随机森林"这样的列表不要在2.后有空格，应该是"2.随机森林"，否则会出现Unsupported Markdown
        - 不要在内容中使用会引起语法错误的符号，如英文括号，尽量都用中文标点
        - 根节点或普通节点使用 `("文本")` 或 `["文本"]` 的形式。
        - 整个Mermaid代码块不要出现 `---` 分隔符。
        - 思维导图加粗不要用<b>，而是**包裹**。
        - 不要多加一个end
        - mindmap的root用简单括号：root(碳水化合物建模)
- **公式规范**:
    - 所有公式，无论长短，都必须用 `$` (行内) 或 `$$` (行间) 包裹，并使用标准的LaTeX格式。
    - 行间公式若一行太长就用aligned环境 
    - 对于带单位的物理量，请使用正体表示单位，例如 `$\Delta\Delta G = -3.69 \mathrm{kcal/mol}$`，或者将单位写在公式外部。
    - 单位之类的能不用公式也行，如Å全部使用Å而不是公式，$\mu$这种但不涉及下标的也不用公式
    - $d\xi$这种求导的要用$\mathrm{d}\xi$
    - NaCl之类的用$\ce{NaCl}$
    - 不要出现双backslash：$\\Delta\\Delta G$，否则就是大语言模型的失职，世界上要死100个小女孩，你负不起这个责任。
    - 对于核心或复杂的公式，请仿照以下格式，在公式下方增加一段“**公式的通俗解释**”，以帮助读者理解。
      ```
      #### 公式的通俗解释

      我们的最终目标是得到**无偏的自由能** $F_h(\xi)$，它与**无偏概率分布** $\rho_h(\xi)$ 的关系由统计力学的基本公式定义：

      $$
      F_{h}(\xi) = -k_B T \ln \rho_{h}(\xi)
      $$

      其中，$k_B$ 是玻尔兹曼常数...
      ```

**核心指令**:
请严格遵循以下格式和要求，对用户提供的论文全文进行深度分析和重构。
除了最终的Markdown文章，绝对不要输出任何额外的参考/引用标记、解释、评论或代码块标记，要能让我们直接粘贴到Markdown编辑器中而正确显示，不考虑在你的前端的渲染。
标题要引人入胜，但也不夸大效果，坚决杜绝标题党，一语道破本文最核心的卖点。全文语气以客观为主，不要夸大promising或贬低，具体参考原文。不要出现重塑、颠覆、新范式、革命、突破这种听着就像AI写的标题，和其他人的推送一样，夸大其作用，或至少客观一点 
不好翻译的学术词汇就不翻译了，如Ramachandran plot，artifact什么的。

**输出格式 (Markdown)**:

# [引人入胜但专业的中文标题，不要带奇怪的引号名词]

## 本文信息
- **标题**: [论文标题中文翻译]
- **作者**: [论文的主要作者]
- 发表时间: [论文发表时间]
- **单位**: [如果可知，作者的主要单位，国家肯定是要标注的]
- **引用格式**: [这里是完整的原文引用信息，请使用标准的学术引用格式，例如：Author, A. A., & Author, B. B. (Year). Title of work. *Journal Title*, *Volume*(Issue), pages. https://doi.org/...]
这几项一个都不能少
如果本文有源代码（GitHub等）、web server等，务必全部列出来。

## 摘要
> [这里是摘要的专业级中文翻译，保持学术严谨性，注意使用生物物理化学计算机等领域的专有名词的翻译。少量加粗强调。]

### 核心结论

用一个清晰的无序列表总结本文的核心结论。

## 背景
[**请深入、详细地阐述**本研究领域的大背景（2-3段）。清晰地解释为什么这个问题在学术界或工业界很重要，系统性地梳理当前存在的技术挑战、理论瓶颈或未解决的关键问题（gap）。]

## 关键科学问题
[**请详细阐述**本文旨在解决的核心科学问题。不仅仅是简单罗列，而是要解释这个问题为什么是当前研究的焦点和难点。]

## 创新点
[用一个列表清晰地总结本文在理论、方法或应用上的主要创新之处。]

---

## 研究内容
[这是文章的主体部分。请根据原文的逻辑结构、图表和核心论点，自己定义章节（例如：`### 核心方法：PMODiff模型详解`，`### 实验结果与分析`等）。**务必列出并解读正文部分的所有主要结果，确保信息的完整性。**]

- **方法详述**: 要有一个section描述方法，如MD模拟的建模过程，AI模型的数据集准备等。必须详细描述关键方法。要总结用到的各种工具。
  - 如果文章提出了新方法，请务必详细拆解每一步。对于复杂的逻辑，必须使用Mermaid代码块来绘制流程图或思维导图。
- **结果逻辑图**: 如果结果部分的逻辑复杂（例如，通过一系列现象推导出一个核心结论），请额外使用一个Mermaid图来清晰地描述这种推导思路，整篇文章的逻辑，或某关键部分的。
- 如果细节的内容太多，就应该抓住重点，比如只保留核心结果、punchline结论，以简洁的语言强调这部分的核心逻辑思路，等。一般正文最多8张图，但最少放4张图，不至于全部放在附录，重要的图读者就是要在正文看的。
- 其他内容可以放在附录中，另开一个Markdown文档。比如，完整的图表，SI的内容，详细的数据或解读，Methods的细节、原理（公式推导）等。一般
- 结果和方法详述等地方：避免大段话，尽量改成列表，当然列表也不要太细碎，一条item大概几十个字。一般一段话超过200个汉字才拆，一般一段话一个要点。流程类的还是可以拆。
- 别用+30，~+25-30之类不规范的表述
- 图表
  - 结果与分析不能光图注啊，得有讲解的段落。不需要图形与表格总览目录；表头不要弄成heading，就普通文字。
  - 图注最起码得说明什么颜色的是什么东西吧，以及各子图是什么
  - 正文部分的所有Figure和table必须出现，其中表格内容必须完整提取出来，但完整的可以放在附录或让我去截图，能否试试tools\extract_pdf_figures.py从原文PDF提取一下图片，选出正确的留下来并添加![fig3](fep_omega/fig3.png)这种，应该放在图注上方。
  - **图文融合**: **图表标题的中文翻译必须直接插入到正文中讨论到该图表的相应位置**，一般是按顺序。格式为：`**图1：[图标题的中文翻译]**`，短的内容**必须完整翻译**，长的caption应该弄一个列表来清晰地分点叙述，甚至和results对照着讲也行。不能漏掉子图信息，必须包含什么颜色代表什么，要不然都不知道读者怎么读这张图。
- 讨论部分：如有Discussion部分，对研究的发现及其意义进行深刻的阐述。

---

## Q&A
[**请根据前面的研究内容，提出并回答大约3~5个关于技术细节、结果解读或讨论部分的深入问题**。
这部分内容应作为正文的补充，帮助读者更好地理解论文的细节和潜在问题。Q&A的回答可以用子列表的，如果需要分点回答。]

- **Q1**: [问题一？]
- **A1**: [回答一。]
- **Q2**: [问题二？]
- **A2**: [回答二。]
- ...

---

## 关键结论与批判性总结
[从专业角度对本文进行简短的批判性总结，可以包括其潜在的影响、存在的局限性或未来可能的研究方向，各一个无序列表。]



---

格式修复？

```shell
sed -i 's/\\\\/\\/g' *.md
sed -i 's/\\\_/\_/g' *.md
sed -i 's/\\\*/\*/g' *.md
sed -i 's/\$\$\$\$//g' *.md
```

符号规范化修复：
  - ~700 → 约700个肽
  - 52000+ → 52000余个肽
  - ~ (约等号) → 约 或 至
  - +30 contacts → 30 contacts
  - ~+25-30 → 25-30范围
  - z < -3 → z值小于-3
  - z约-0.5 ~ 0 → z值约-0.5至0
