# 工具使用指南

所有自动化脚本的详细说明和使用方法。

## 1. PDF文本搜索

### tools/search_pdf_text.py

**功能**：快速搜索PDF中的关键词，获取上下文

**用法**：
```bash
python3 tools/search_pdf_text.py <pdf_path> <keyword> [context_lines]
```

**参数**：
- `pdf_path`: PDF文件路径
- `keyword`: 要搜索的关键词
- `context_lines`: 上下文行数（默认3行）

**示例**：
```bash
# 搜索作者信息
python3 tools/search_pdf_text.py paper.pdf "Author" 3

# 搜索关键结论
python3 tools/search_pdf_text.py paper.pdf "conclusion" 5
```

---

## 2. 图片提取

### tools/extract_pdf_figures.py

**功能**：从PDF中提取指定的图片

**用法**：
```bash
# 自动模式：按页码提取
python3 tools/extract_pdf_figures.py <pdf_path> <output_dir> --pages 3,5,6

# 指定Figure编号
python3 tools/extract_pdf_figures.py <pdf_path> <output_dir> --figures "1:3,2:5,3:6"

# 交互模式
python3 tools/extract_pdf_figures.py <pdf_path> <output_dir> --interactive
```

**参数**：
- `pdf_path`: PDF文件路径
- `output_dir`: 输出目录（建议使用文章同名文件夹）
- `--pages`: 指定页码（逗号分隔）
- `--figures`: 指定Figure编号和页码的映射
- `--interactive`: 交互模式

**示例**：
```bash
# 提取第3、5、6页的图片
python3 tools/extract_pdf_figures.py "paper.pdf" "article_folder/" --pages 3,5,6

# Figure 1在第3页，Figure 2在第5页
python3 tools/extract_pdf_figures.py "paper.pdf" "article_folder/" --figures "1:3,2:5"
```

**建议**：
- 优先使用自动模式
- 图片保存到与markdown文件同名的文件夹

---

## 3. 随机缩略图生成

### tools/random_thumbnail.py

**功能**：从两个缩略图目录中随机选择一张图片

**用法**：
```bash
cd /mnt/e/GitHub-repo/mendelevium
python3 tools/random_thumbnail.py --frontmatter
```

**参数**：
- `--frontmatter`: 输出frontmatter格式（可直接复制粘贴）

**输出示例**：
```
  thumbnail: "/assets/img/thumbnail_mine/wh-z8p9rj.jpg"
  image: "/assets/img/thumbnail_mine/wh-z8p9rj.jpg"
```

**注意**：
- 必须在项目根目录运行（/mnt/e/GitHub-repo/mendelevium）
- image和thumbnail会自动一致
- 从`assets/img/thumbnail/`和`assets/img/thumbnail_mine/`中随机选择

---

## 4. 中文引号修复

### tools/convert_quotes.py

**功能**：将英文引号转换为中文引号（Unicode 201C/201D）

**用法**：
```bash
python3 tools/convert_quotes.py <markdown_file>
```

**参数**：
- `markdown_file`: 要修复的Markdown文件路径

**示例**：
```bash
python3 tools/convert_quotes.py "article.md"
```

**修复内容**：
- 英文引号 `"` → 中文引号 `"`（U+201C）
- 英文引号 `"` → 中文引号 `"`（U+201D）

**排除区域**：
- frontmatter
- 代码块（\`\`\`）
- 行内代码（\`\`）
- mermaid图

---

## 5. 中文括号修复

### tools/fix_parentheses.py

**功能**：自动修复英文括号为中文括号

**用法**：
```bash
python3 tools/fix_parentheses.py <markdown_file>
```

**参数**：
- `markdown_file`: 要修复的Markdown文件路径

**示例**：
```bash
python3 tools/fix_parentheses.py "article.md"
```

**修复规则**：
- 中文后的括号：`塑料(PET)` → `塑料（PET）`
- 括号内包含中文：`(北京大学)` → `（北京大学）`
- 混合内容智能判断：`方法(QM/MM)` → `方法（QM/MM）`

**保留英文括号的情况**：
- 公式中：`$f(x)$`
- 代码中：`print("hello")`
- 坐标：`(1, 2, 3)` 或 `(x, y, z)`
- 纯英文ID：`(PDB ID: 5XG0)`

**排除区域**：
- frontmatter
- 代码块
- 行内代码
- 公式（行内和行间）
- mermaid图

**注意**：此工具使用上下文感知算法，自动判断是否需要转换括号

---

## 6. 格式修复

### tools/fix_format.sh

**功能**：修复常见的Markdown格式问题

**用法**：
```bash
bash tools/fix_format.sh <markdown_file>
```

**参数**：
- `markdown_file`: 要修复的Markdown文件路径

**示例**：
```bash
bash tools/fix_format.sh "article.md"
```

**修复内容**：
1. 删除列表项之间的空行
2. 修复加粗包含百分号（**24%** → **24**%）
3. 其他常见格式问题

---

## 7. 质量自动检查

### tools/check_blog_quality.py

**功能**：全面的文章质量自动检查

**用法**：
```bash
python3 tools/check_blog_quality.py <markdown_file>
```

**参数**：
- `markdown_file`: 要检查的Markdown文件路径

**示例**：
```bash
python3 tools/check_blog_quality.py "article.md"
```

**检查项目**：
- ✅ Frontmatter完整性
- ✅ 文章结构完整性
- ✅ 公式格式（双backslash检测）
- ✅ 加粗格式（标点、百分号、引号、括号）
- ✅ 中文标点（逗号、分号、引号、**括号**）
- ✅ 转义符
- ✅ 不规范符号
- ✅ Mermaid语法
- ✅ 图片数量和图注
- ✅ 文章长度
- ✅ 标题党词汇

**输出**：
- ❌ 错误（必须修复）
- ⚠️ 警告（建议修复）
- ✅ 所有检查通过

**返回值**：
- 0: 无错误（可能有警告）
- 1: 有错误

---

## 8. 图片验证（最终验证）⭐

### tools/verify_blog_figures.py

**功能**：验证博客文章中所有Figure/Scheme的编号、文件和图注

**用法**：
```bash
python3 tools/verify_blog_figures.py <markdown_file> [pdf_file]
```

**参数**：
- `markdown_file`: 要验证的Markdown文件路径
- `pdf_file`: 原始PDF文件路径（可选，用于对照验证）

**示例**：
```bash
python3 tools/verify_blog_figures.py "article.md"
python3 tools/verify_blog_figures.py "article.md" "original.pdf"
```

**验证项目**：
1. **图片引用提取**：提取所有 `![figX]` 和 `![schemeX]` 引用
2. **文件存在性**：检查所有引用的图片文件是否存在
3. **编号连续性**：检查Figure和Scheme的编号是否连续（1,2,3...）
4. **图注匹配**：检查图注（**图X：**）和图片引用是否一一对应

**输出示例**：
```
🔍 验证博客文章图片: article.md

✓ 找到 8 个图片引用:
  - fig1: petase_mechanism/fig1.png
  - fig2: petase_mechanism/fig2.png
  - scheme1: petase_mechanism/scheme1.png

📁 检查图片文件是否存在:
  ✓ fig1: /path/to/fig1.png (125.3 KB)
  ✗ fig2: 文件不存在

🔢 检查编号连续性:
  ✓ Figure编号连续: 1-5
  ⚠️ Scheme编号不连续，缺少: [2]

📝 找到 5 个图注:
  - 图1: PETase的晶体结构...

🔗 检查图片引用和图注匹配:
  ⚠️ 有图片引用但无图注: [6]

✅ 验证通过 / ❌ 验证失败
```

**返回值**：
- 0: 验证通过
- 1: 验证失败（有缺失文件）

---

## 9. Scheme验证（OCR辅助）

### tools/verify_scheme_in_image.py

**功能**：使用OCR验证图片中是否包含"Scheme"文字

**用法**：
```bash
# 扫描PDF中所有包含Scheme的图片
python3 tools/verify_scheme_in_image.py <pdf_path>

# 检查指定图片
python3 tools/verify_scheme_in_image.py <pdf_path> <page_num> <img_index>
```

**参数**：
- `pdf_path`: PDF文件路径
- `page_num`: 页码（1-based，可选）
- `img_index`: 图片索引（0-based，可选）

**示例**：
```bash
# 扫描所有Scheme图
python3 tools/verify_scheme_in_image.py "paper.pdf"

# 检查第3页的第1张图
python3 tools/verify_scheme_in_image.py "paper.pdf" 3 0
```

**输出示例**：
```
Scanning paper.pdf for images containing 'Scheme'...

✓ Page 2, Image 1: Scheme 1 (1250x556)
✓ Page 2, Image 2: Scheme 2 (667x556)
✓ Page 3, Image 2: Scheme 3 (667x641)

=== Summary ===
Found 3 images containing 'Scheme':
  Scheme 1: Page 2, Image 1 (1250x556)
  Scheme 2: Page 2, Image 2 (667x556)
  Scheme 3: Page 3, Image 2 (667x641)
```

**注意**：
- 需要安装依赖：`pip install pdf2image pytesseract poppler-utils`
- 此工具为**辅助验证**，不是必须步骤
- 主要判断标准仍然是**图注位置验证**（图注与图片在同一页）

---

## 工具使用顺序

**标准工作流**：
```bash
# 1. 提取内容
python3 tools/search_pdf_text.py paper.pdf "keyword"
python3 tools/extract_pdf_figures.py paper.pdf output/ --pages 3,5,6

# 2. 生成frontmatter
cd /mnt/e/GitHub-repo/mendelevium
python3 tools/random_thumbnail.py --frontmatter

# 3. 格式修复（按顺序执行）
python3 tools/convert_quotes.py article.md        # 修复引号
python3 tools/fix_parentheses.py article.md       # 修复括号
bash tools/fix_format.sh article.md               # 其他格式

# 4. 质量检查
python3 tools/check_blog_quality.py article.md

# 5. 图片最终验证（必须执行）⭐
python3 tools/verify_blog_figures.py article.md
```

---

## 批量处理

如果有主文档和附录，可以批量处理：

```bash
# 格式修复
for file in main.md appendix.md; do
    python3 tools/convert_quotes.py "$file"
    python3 tools/fix_parentheses.py "$file"  # 🆕
    bash tools/fix_format.sh "$file"
done

# 质量检查
for file in main.md appendix.md; do
    python3 tools/check_blog_quality.py "$file"
done
```
