# Tools 使用指南

本目录包含若干Markdown和博客编辑工具。

## 1. convert_quotes.py - 引号转换工具

**功能**: 将英文双引号转换为中文引号（符合CLAUDE.md规范）

**用法**:
```bash
python3 tools/convert_quotes.py <file> [file2] [file3] ...
```

**示例**:
```bash
# 转换单个文件
python3 tools/convert_quotes.py _pages/article.md

# 转换多个文件
python3 tools/convert_quotes.py _pages/a.md _pages/b.md _pages/c.md

# 转换目录中所有markdown文件（在bash中执行）
python3 tools/convert_quotes.py _pages/**/*.md
```

**说明**:
- 将英文双引号 `"` (U+0022) 转换为中文引号对
- 第1、3、5...个引号转为左引号 `"` (U+201C)
- 第2、4、6...个引号转为右引号 `"` (U+201D)
- 自动配对，避免手动修改

---

## 2. random_thumbnail.py - 随机缩略图选择

**功能**: 为新建文档随机选择缩略图路径

**用法**:
```bash
python3 tools/random_thumbnail.py
```

---

## 3. search_pdf_text.py - PDF文本搜索

**功能**: 搜索PDF中的特定关键词，显示上下文

**用法**:
```bash
python3 tools/search_pdf_text.py <pdf_file> <keyword> [context_lines]
```

**示例**:
```bash
python3 tools/search_pdf_text.py "_pages/paper.pdf" "算法" 3
```

---

## 4. blog_tools.py - 博客工具集

**功能**: 综合的博客编辑工具

**用法**: 根据脚本内文档使用

---

## 快速参考

| 工具 | 用途 | 主要命令 |
|-----|------|--------|
| convert_quotes.py | 引号规范化 | `python3 tools/convert_quotes.py file.md` |
| random_thumbnail.py | 缩略图选择 | `python3 tools/random_thumbnail.py` |
| search_pdf_text.py | PDF搜索 | `python3 tools/search_pdf_text.py file.pdf keyword` |
| blog_tools.py | 综合工具 | 见脚本文档 |
