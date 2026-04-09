# Blog Tools 工具集

本目录包含自动化工具，用于博客文章生成和质量控制。

## 📋 工具列表

### 格式修复工具

| 工具 | 功能 | 主要用途 |
|-----|------|---------|
| `fix_markdown.py` | **综合格式修复** | 一键修复所有格式问题（推荐使用） |
| `convert_quotes.py` | 引号修复 | 英文引号 → 中文引号 |
| `fix_parentheses.py` | 括号修复 | 英文括号 → 中文括号（智能识别） |
| `fix_format.sh` | 快速修复 | Shell脚本快速修复 |

### 内容生成工具

| 工具 | 功能 | 主要用途 |
|-----|------|---------|
| `search_pdf_text.py` | PDF文本搜索 | 快速定位论文关键词和上下文 |
| `extract_pdf_figures.py` | PDF图片提取 | 自动提取论文中的Figure |
| `random_thumbnail.py` | 随机缩略图 | 为文章生成随机缩略图路径 |

### 质量检查工具

| 工具 | 功能 | 主要用途 |
|-----|------|---------|
| `check_blog_quality.py` | 质量检查 | 11项自动化质量检查 |
| `verify_blog_figures.py` | 图片验证 | **最终验证**：检查图片编号、文件和图注 |
| `verify_scheme_in_image.py` | Scheme验证 | OCR验证图片内是否包含Scheme文字 |
| `update_frontmatter_dates.py` | 日期更新 | 批量更新frontmatter日期 |

## 🚀 快速开始

### 标准博客生成工作流

```bash
# 1. 提取PDF内容
python3 tools/search_pdf_text.py paper.pdf "关键词" 3
python3 tools/extract_pdf_figures.py paper.pdf output_dir/ --pages 3,5,6

# 2. 生成frontmatter
cd /mnt/e/GitHub-repo/mendelevium
python3 tools/random_thumbnail.py --frontmatter

# 3. 格式修复（按顺序）
python3 tools/convert_quotes.py article.md
python3 tools/fix_parentheses.py article.md
bash tools/fix_format.sh article.md

# 4. 质量检查
python3 tools/check_blog_quality.py article.md

# 5. 最终验证（必须执行）
python3 tools/verify_blog_figures.py article.md
```

### 批量处理多个文件

```bash
for file in main.md appendix.md; do
    python3 tools/convert_quotes.py "$file"
    python3 tools/fix_parentheses.py "$file"
    bash tools/fix_format.sh "$file"
    python3 tools/check_blog_quality.py "$file"
done
```

## 📖 详细文档

完整的使用指南和规范请参考：
- `.claude/skills/blog/modules/workflow/tools.md` - 工具使用详细说明
- `.claude/skills/blog/skill.md` - 完整工作流程

## 🔧 工具依赖

```bash
# Python依赖
pip install pymupdf pillow

# 系统依赖
# (所有工具均使用Python标准库或已安装的包)
```

## 📝 使用示例

### 1. PDF搜索
```bash
python3 tools/search_pdf_text.py "_pages/paper.pdf" "Author" 3
```

### 2. 提取图片
```bash
# 按Figure编号提取
python3 tools/extract_pdf_figures.py paper.pdf output/ --figures "1:3,2:5"

# 按页码提取
python3 tools/extract_pdf_figures.py paper.pdf output/ --pages 3,5,6
```

### 3. 生成缩略图
```bash
cd /mnt/e/GitHub-repo/mendelevium
python3 tools/random_thumbnail.py --frontmatter
```

### 4. 修复标点符号
```bash
# 引号
python3 tools/convert_quotes.py article.md

# 括号
python3 tools/fix_parentheses.py article.md
```

### 5. 质量检查
```bash
python3 tools/check_blog_quality.py article.md
# 返回值: 0=无错误, 1=有错误
```

### 6. 图片验证（最终步骤）
```bash
# 验证所有图片编号、文件和图注
python3 tools/verify_blog_figures.py article.md

# 可选：验证Scheme图片内容
python3 tools/verify_scheme_in_image.py paper.pdf
```

## ⚠️ 注意事项

1. **工作目录**: `random_thumbnail.py` 必须在项目根目录运行
2. **执行顺序**: 格式修复工具应按 引号 → 括号 → 其他格式 的顺序执行
3. **质量检查**: 在格式修复完成后再运行质量检查
4. **最终验证**: **必须**在完成文章后运行 `verify_blog_figures.py` 验证所有图片
5. **文件编码**: 所有工具均使用 UTF-8 编码

## 🔗 相关资源

- **Skills文档**: `.claude/skills/blog/`
- **项目规范**: `CLAUDE.md`
- **Git仓库**: [mendelevium](https://github.com/your-repo)
