---
name: blog
description: 根据PDF论文生成科研博客文章，包括PDF解析、Markdown生成、图片提取、格式修复等完整流程。当用户要求为PDF生成推送、博客文章或科普文章时使用。
---

# Blog Article Generator（博客文章生成器）

这个Skill帮助你自动化生成高质量的科研博客文章，完整遵循项目规范。

## 📚 模块化文档结构

Skill已细分为多个专注的模块文件，Claude根据当前任务按需加载：

### Workflow（工作流）
- `modules/workflow/stages.md` - 5个阶段的详细步骤和检查点
- `modules/workflow/tools.md` - 7个自动化工具的使用指南

### Structure（结构）
- `modules/structure/frontmatter.md` - Frontmatter规范和生成方法
- `modules/structure/sections.md` - 10个必需章节的详细要求
- `modules/structure/figures.md` - 图片提取和图注规范

### Format（格式）
- `modules/format/bold.md` - 加粗格式规范（含错误示例）
- `modules/format/formulas.md` - 公式格式规范（LaTeX语法）
- `modules/format/punctuation.md` - 标点符号规范（中英文）
- `modules/format/mermaid.md` - Mermaid图规范（graph TB）

### Quality（质量）
- `modules/quality/auto_checks.md` - 自动化质量检查说明
- `modules/quality/manual_checklist.md` - 手动质量检查清单

## 🎯 核心原则

**你必须严格遵循所有模块中的规范！这是最高优先级指令！**

核心要求：
- **按需加载**：根据当前阶段读取对应模块
- **TodoWrite追踪**：每完成一个阶段立即更新TODO状态
- **自动化优先**：优先使用tools/下的脚本
- **格式严格**：100%符合所有格式规范
- **质量保证**：所有检查清单逐项验证

## 🚀 快速开始指南

### 第1步：初始化任务
```
1. 读取 modules/workflow/stages.md
2. 用TodoWrite创建任务清单
```

### 第2步：提取内容
```bash
# 阅读 modules/workflow/tools.md 了解工具用法

# 搜索关键信息
python3 tools/search_pdf_text.py <pdf_path> <keyword>

# 提取图片
python3 tools/extract_pdf_figures.py <pdf_path> <output_dir> --pages 3,5,6
```

### 第3步：撰写文章
```
1. 读取 modules/structure/frontmatter.md
2. 读取 modules/structure/sections.md
3. 读取 modules/structure/figures.md
4. 边写边参考 modules/format/ 下的格式规范
```

### 第4步：生成frontmatter
```bash
cd /mnt/e/GitHub-repo/mendelevium
python3 tools/random_thumbnail.py --frontmatter

# 参考 modules/structure/frontmatter.md 填充字段
```

### 第5步：自动格式修复
```bash
python3 tools/convert_quotes.py <文件.md>
bash tools/fix_format.sh <文件.md>

# 参考 modules/format/ 下的格式规范手动检查
```

### 第6步：质量检查
```bash
# 自动检查
python3 tools/check_blog_quality.py <文件.md>

# 手动检查
# 1. 读取 modules/quality/auto_checks.md
# 2. 读取 modules/quality/manual_checklist.md
# 3. 逐项验证
```

## 📖 详细文档索引

### 工作流阶段划分
详见 `modules/workflow/stages.md`：
- 阶段0：初始化（TodoWrite）
- 阶段1：内容提取（PDF + 图片）
- 阶段2：文章撰写（主文档 + 附录）
- 阶段3：Frontmatter生成
- 阶段4：自动化格式修复
- 阶段5：质量检查

### 自动化工具使用
详见 `modules/workflow/tools.md`：
1. search_pdf_text.py - PDF文本搜索
2. extract_pdf_figures.py - 图片提取
3. random_thumbnail.py - 随机缩略图
4. convert_quotes.py - 引号修复
5. fix_format.sh - 格式修复
6. check_blog_quality.py - 质量检查

### 文章结构要求
- `modules/structure/frontmatter.md` - Frontmatter各字段详解
- `modules/structure/sections.md` - 10个必需章节
- `modules/structure/figures.md` - 图片插入和图注规范

### 格式规范详解
- `modules/format/bold.md` - 加粗规范（标点、百分号、引号、括号）
- `modules/format/formulas.md` - 公式规范（双backslash、微分、化学式、单位）
- `modules/format/punctuation.md` - 标点规范（中文标点、引号、转义符）
- `modules/format/mermaid.md` - Mermaid规范（graph TB、---禁止、列表格式）

### 质量检查
- `modules/quality/auto_checks.md` - 自动化检查说明
- `modules/quality/manual_checklist.md` - 手动检查清单

## 🚨 常见高危错误速查

### 最严重（导致渲染失败）
- 公式双backslash: `$\\Delta$` → `$\Delta$`
- 加粗包含标点: **重要.** → **重要**.
- Mermaid包含`---`分隔符

### 严重（破坏格式）
- frontmatter的date不是今天
- image和thumbnail不一致
- 英文引号代替中文引号
- 英文括号代替中文括号（新增！）

**详细错误列表和修复方法请查阅对应模块**

## 📁 输出位置

- 主文档：`_pages/<分类>/YYYY-MM-DD-<简短标题>.md`
- 附录（如有）：`_pages/<分类>/YYYY-MM-DD-<简短标题>-appendix.md`
- 图片文件夹：`_pages/<分类>/YYYY-MM-DD-<简短标题>/`

## 💡 使用提示

### 对于Claude

**执行每个阶段前，必须先读取对应的模块文档！**

| 当前任务 | 必读模块 |
|---------|---------|
| 开始任务 | workflow/stages.md |
| 使用工具 | workflow/tools.md |
| 撰写文章 | structure/frontmatter.md, sections.md, figures.md |
| 检查格式 | format/bold.md, formulas.md, punctuation.md, mermaid.md |
| 质量验证 | quality/auto_checks.md, manual_checklist.md |

**每一步自问**：
1. 我是否读取了对应模块的规范？
2. 我是否避免了所有高危错误？
3. 我是否用TodoWrite更新了进度？

### 对于用户

- **信任自动化**：工具已高度优化，覆盖常见错误
- **关注TODO**：通过TODO了解进度
- **最终验证**：运行`check_blog_quality.py`确认无误

## 📚 相关文件

- **项目规范**: `/mnt/e/GitHub-repo/mendelevium/CLAUDE.md`
- **原始完整指南**: `workflows/article_guide.md`（保留作为参考）
- **工具源码**: `/mnt/e/GitHub-repo/mendelevium/tools/`
- **详细README**: `README.md`（模块化说明）

## 触发场景

当用户说以下内容时，自动激活此Skill：
- "请为xxx.pdf做一篇推送"
- "根据xxx.pdf写一篇博客"
- "为xxx.pdf生成科普文章"
- "根据@CLAUDE.md为xxx.pdf写推送"
- "根据SI也写一下"
