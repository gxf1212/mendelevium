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

## 4. extract_pdf_figures.py - PDF图片提取工具

**功能**: 从PDF中自动或交互式提取Figure图片

**用法**:

### 自动模式（按Figure编号）
```bash
python3 tools/extract_pdf_figures.py <pdf_file> <output_dir> --figures "1:3,2:5,3:6"
```
- 提取Figure 1从第3页 → 保存为`fig1.png`
- 提取Figure 2从第5页 → 保存为`fig2.png`
- 提取Figure 3从第6页 → 保存为`fig3.png`

### 自动模式（按页码）
```bash
python3 tools/extract_pdf_figures.py <pdf_file> <output_dir> --pages 3,5,6
```

### 交互模式
```bash
python3 tools/extract_pdf_figures.py <pdf_file> <output_dir> --interactive
```
- 显示所有图片列表（按面积排序）
- 手动选择需要的图片
- 自定义文件名

### 调整尺寸过滤
```bash
python3 tools/extract_pdf_figures.py <pdf_file> <output_dir> --pages 3,5,6 --min-size 100000
```
- 默认最小面积：50000px²（过滤小图标）
- 提高阈值只提取大图，降低阈值包含较小图片

**实际案例**:
```bash
# ML/MM论文图片提取
python3 tools/extract_pdf_figures.py \
  "_pages/Free Energy/paper.pdf" \
  "_pages/Free Energy/figures" \
  --figures "1:3,2:5,3:6"
```

**在Markdown中引用**:
```markdown
![fig2](figures/fig2.png)

**图2：实验结果对比**
```

**参数说明**:
- `pdf_path`: PDF文件路径（必需）
- `output_dir`: 输出目录（必需）
- `--pages`: 页码列表，如`3,5,6`
- `--figures`: Figure映射，如`"1:3,2:5,3:6"`
- `--interactive` / `-i`: 启用交互模式
- `--min-size`: 最小图片面积（px²），默认50000

**工作原理**:
1. 扫描PDF中所有图片
2. 按面积排序（主图通常最大）
3. 过滤掉小于min-size的图片
4. 提取每页最大的图片
5. 保存为原格式（png/jpg等）

---

## 5. blog_tools.py - 博客工具集

**功能**: 综合的博客编辑工具

**用法**: 根据脚本内文档使用

---

## 快速参考

| 工具 | 用途 | 主要命令 |
|-----|------|--------|
| convert_quotes.py | 引号规范化 | `python3 tools/convert_quotes.py file.md` |
| random_thumbnail.py | 缩略图选择 | `python3 tools/random_thumbnail.py` |
| search_pdf_text.py | PDF搜索 | `python3 tools/search_pdf_text.py file.pdf keyword` |
| extract_pdf_figures.py | PDF图片提取 | `python3 tools/extract_pdf_figures.py pdf out --figures "1:3,2:5"` |
| blog_tools.py | 综合工具 | 见脚本文档 |
