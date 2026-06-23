# Mendelevium Blog - Claude Development Guide

本文件只做索引，不重复规则正文。执行任务时先按用户要求定位到对应 rule，再读取具体文件。

## Skills可用性
本项目配置了多个skills用于特定任务，当用户调用相关skill时会自动加载对应规则：
- `blog` - 根据PDF论文生成科研博客文章，会自动引用格式规则
- `format` - 检查和修复推送文章格式问题，会自动引用格式规范  
- `literature-deepread` - 文献精读流程，会自动引用阅读和分析规则
- `fix-ai-writing` - 修复AI写作特征，会自动引用 [`.claude/skills/fix-ai-writing.md`](.claude/skills/fix-ai-writing.md) 去AI味规则
- `adversarial-paper-reading` - 对抗式论文阅读，会自动引用审稿和批判性分析规则
- 其他skills按需引用相关规则文件

Skills会根据任务类型自动加载对应的rules文件，确保格式规范的一致性。

## Rule目录总览

### 01-project：项目结构和博客组织
处理项目理念、技术栈、栏目、frontmatter、文件命名、新文章元数据时：
- [01-project/01-project-basics.md](.claude/rules/01-project/01-project-basics.md)
- [01-project/02-posts-and-common-requests.md](.claude/rules/01-project/02-posts-and-common-requests.md)

### 02-tools：工具、PDF和技术记录  
用工具、查PDF、整理技术记录、处理PMC/PubMed下载时：
- [02-tools/01-tools-pdf-and-technical-notes.md](.claude/rules/02-tools/01-tools-pdf-and-technical-notes.md)

### 03-article-workflow：推文角色、工作流和核心指令
"根据PDF写推文""写微信公众号文章""读论文写文章"时：
- [03-article-workflow/01-role-workflow-and-core-instructions.md](.claude/rules/03-article-workflow/01-role-workflow-and-core-instructions.md)

### 04-article-formatting：文章格式细则
"修格式""中文标点""加粗""列表""去AI味""检查全文格式"时：
- [04-article-formatting/01-bold.md](.claude/rules/04-article-formatting/01-bold.md)
- [04-article-formatting/02-punctuation-headings-abbreviations.md](.claude/rules/04-article-formatting/02-punctuation-headings-abbreviations.md)
- [04-article-formatting/03-markdown-lists.md](.claude/rules/04-article-formatting/03-markdown-lists.md)
- [04-article-formatting/04-ai-writing-and-quote-blocks.md](.claude/rules/04-article-formatting/04-ai-writing-and-quote-blocks.md)

### 05-visual-elements：可视化元素
处理Mermaid、图、表、图片提取、图注、图文融合时：
- [05-visual-elements/01-mermaid.md](.claude/rules/05-visual-elements/01-mermaid.md)
- [05-visual-elements/02-figures-tables.md](.claude/rules/05-visual-elements/02-figures-tables.md)

### 06-math-chemistry：数学和化学规范
处理公式、单位、化学式、SMILES、数学比较符号时：
- [06-math-chemistry/01-formulas.md](.claude/rules/06-math-chemistry/01-formulas.md)

### 07-article-structure：推文结构模板
写单篇研究论文解读、检查推文结构、补齐章节时：
- [07-article-structure/01-title-and-paper-info.md](.claude/rules/07-article-structure/01-title-and-paper-info.md)
- [07-article-structure/02-abstract-background.md](.claude/rules/07-article-structure/02-abstract-background.md)
- [07-article-structure/03-research-content-and-final-check.md](.claude/rules/07-article-structure/03-research-content-and-final-check.md)

### 08-maintenance-and-interaction：维护、符号规范和交互风格
"格式修复""符号规范化""回答简洁点""按项目交互风格反馈"时：
- [08-maintenance-and-interaction/01-format-symbols-and-interaction.md](.claude/rules/08-maintenance-and-interaction/01-format-symbols-and-interaction.md)

## 常见任务索引

### "根据PDF写推文/写一篇论文解读"
1. [03-article-workflow/01-role-workflow-and-core-instructions.md](.claude/rules/03-article-workflow/01-role-workflow-and-core-instructions.md)
2. [02-tools/01-tools-pdf-and-technical-notes.md](.claude/rules/02-tools/01-tools-pdf-and-technical-notes.md)
3. [07-article-structure/01-title-and-paper-info.md](.claude/rules/07-article-structure/01-title-and-paper-info.md)
4. [07-article-structure/02-abstract-background.md](.claude/rules/07-article-structure/02-abstract-background.md)
5. [07-article-structure/03-research-content-and-final-check.md](.claude/rules/07-article-structure/03-research-content-and-final-check.md)
6. [05-visual-elements/02-figures-tables.md](.claude/rules/05-visual-elements/02-figures-tables.md)
7. [06-math-chemistry/01-formulas.md](.claude/rules/06-math-chemistry/01-formulas.md)

### "检查/修复文章格式、中文标点、加粗、列表"
1. [04-article-formatting/01-bold.md](.claude/rules/04-article-formatting/01-bold.md)
2. [04-article-formatting/02-punctuation-headings-abbreviations.md](.claude/rules/04-article-formatting/02-punctuation-headings-abbreviations.md)
3. [04-article-formatting/03-markdown-lists.md](.claude/rules/04-article-formatting/03-markdown-lists.md)
4. [04-article-formatting/04-ai-writing-and-quote-blocks.md](.claude/rules/04-article-formatting/04-ai-writing-and-quote-blocks.md)
5. [08-maintenance-and-interaction/01-format-symbols-and-interaction.md](.claude/rules/08-maintenance-and-interaction/01-format-symbols-and-interaction.md)

### "创建新文章/补frontmatter/改文件名/选缩略图"
1. [01-project/02-posts-and-common-requests.md](.claude/rules/01-project/02-posts-and-common-requests.md)
2. [01-project/01-project-basics.md](.claude/rules/01-project/01-project-basics.md)

### "整理技术记录/QQ记录/笔记整理"
1. [02-tools/01-tools-pdf-and-technical-notes.md](.claude/rules/02-tools/01-tools-pdf-and-technical-notes.md)
2. [01-project/02-posts-and-common-requests.md](.claude/rules/01-project/02-posts-and-common-requests.md)
3. [04-article-formatting/03-markdown-lists.md](.claude/rules/04-article-formatting/03-markdown-lists.md)

### "处理图表/重新提取fig/scheme/检查图片编号"
1. [05-visual-elements/02-figures-tables.md](.claude/rules/05-visual-elements/02-figures-tables.md)
2. [02-tools/01-tools-pdf-and-technical-notes.md](.claude/rules/02-tools/01-tools-pdf-and-technical-notes.md)

### "修公式/单位/化学式/SMILES"
1. [06-math-chemistry/01-formulas.md](.claude/rules/06-math-chemistry/01-formulas.md)
2. [04-article-formatting/01-bold.md](.claude/rules/04-article-formatting/01-bold.md)
