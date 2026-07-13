# Blog Skill 模块目录

## 工作流（Workflow）
- **stages.md**：分阶段说明推送的完整流程，包含TodoWrite模板和每个阶段的检查点。
- **tools.md**：所有辅助脚本的功能、参数与示例，覆盖PDF搜索、图片提取、缩略图生成、格式修复和质量验证。

## 结构（Structure）
- **frontmatter.md**：Frontmatter字段、命名规则与缩略图生成要求（含last_modified_at必填、title中文标点）。
- **sections.md**：标准文章必须包含的9个章节及写作规范。
- **figures.md**：图片数量、命名、图注和图文融合的详细指南。

## 格式（Format）
- **bold.md**：加粗使用准则、常见错误、2026新增约束（禁包LaTeX/化学式、短句限制等）。
- **formulas.md**：公式、化学式、单位、微分符号及"公式的通俗解释"模板（含连续公式规范）。
- **punctuation.md**：中文标点、引号、括号与转义符规范。
- **mermaid.md**：Mermaid布局、子图、列表格式和mindmap写法。
- **lists.md**：列表与段落/表格的选择准则、长度规范、间距要求与典型错误示例。

## 质量（Quality）
- **auto_checks.md**：`tools/check_blog_quality.py` 的检测项与输出示例。
- **manual_checklist.md**：提交前的手动核对清单，覆盖结构、格式、内容与frontmatter。

## 使用建议
1. 初始化任务前阅读workflow模块并利用TodoWrite跟踪进度。
2. 撰写阶段同步参考structure与format模块，确保结构、图文与符号一次到位。
3. 完稿后先运行自动化检查，再依manual_checklist逐项确认。
4. 如遇图像、公式或Mermaid问题，回看对应模块即可快速定位规范。
