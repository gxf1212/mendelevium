# 工作流阶段划分

每个阶段独立完成，完成后立即用TodoWrite标记。

## 阶段0：初始化

### TodoWrite任务模板
```
[ ] 提取PDF全文和关键信息
[ ] 提取所有需要的图片
[ ] 撰写主文档初稿
[ ] 拆分附录（如果主文档超长）
[ ] 生成frontmatter和缩略图
[ ] 自动格式修复
[ ] 质量检查清单验证
[ ] 图片最终验证（必须执行）
```

---

## 阶段1：内容提取

### 步骤1.1：提取PDF全文
```bash
# 方法1：快速搜索关键信息
python3 tools/search_pdf_text.py <pdf_path> <keyword> [context_lines]

# 方法2：读取完整PDF
Read工具直接读取PDF文件
```

### 步骤1.2：提取图片
```bash
# 自动模式（推荐）
python3 tools/extract_pdf_figures.py <pdf_path> <output_dir> --pages 3,5,6

# 指定Figure编号
python3 tools/extract_pdf_figures.py <pdf_path> <output_dir> --figures "1:3,2:5,3:6"
```

### ✅ Checkpoint
- [ ] PDF全文已提取
- [ ] 关键图片已提取（4-8张）
- [ ] 图片保存在正确的文件夹
- [ ] TodoWrite标记为completed

---

## 阶段2：文章撰写

### 步骤2.1：撰写主文档
严格遵循 `structure/template.md` 的结构要求

### 步骤2.2：检查长度
主文档必须**严格250-300行**

### 步骤2.3：拆分附录（如需要）
如果主文档超过300行，创建附录文件：
- 文件名：`YYYY-MM-DD-<标题>-appendix.md`
- 附录不需要自己的Q&A

### ✅ Checkpoint
- [ ] 主文档结构完整
- [ ] 主文档长度250-300行
- [ ] 附录已拆分（如需要）
- [ ] frontmatter的date是**今天的日期**
- [ ] TodoWrite标记为completed

---

## 阶段3：Frontmatter生成

### 步骤3.1：生成随机缩略图
```bash
cd /mnt/e/GitHub-repo/mendelevium
python3 tools/random_thumbnail.py --frontmatter
```

### 步骤3.2：填充frontmatter
```yaml
---
title: "文章标题（中文，使用中文标点）"
date: "YYYY-MM-DD"  # 今天的日期！
tags: [tag1, tag2, tag3]
description: "简短描述"
image: "/assets/img/thumbnail/xxx.jpg"
thumbnail: "/assets/img/thumbnail/xxx.jpg"  # 必须与image一致！
author: Xufan Gao
lang: zh-CN
---
```

### ✅ Checkpoint
- [ ] 缩略图已随机生成
- [ ] image和thumbnail参数**完全一致**
- [ ] date是**今天的日期**
- [ ] 未使用empty.jpg
- [ ] TodoWrite标记为completed

---

## 阶段4：自动化格式修复

### 步骤4.1：中文标点修复
```bash
python3 tools/convert_quotes.py <主文档.md>
python3 tools/convert_quotes.py <附录.md>  # 如有附录
```

### 步骤4.2：列表格式修复
```bash
bash tools/fix_format.sh <主文档.md>
bash tools/fix_format.sh <附录.md>  # 如有附录
```

### 步骤4.3：符号规范化检查
手动检查并修复：
- `~700` → `约700个`
- `52000+` → `52000余个`
- `+30` → `30`
- `z < -3` → `z值小于-3`

### ✅ Checkpoint
- [ ] 中文引号已修复（201C/201D）
- [ ] 列表格式已修复
- [ ] 不规范符号已修复
- [ ] TodoWrite标记为completed

---

## 阶段5：质量检查

### 步骤5.1：运行自动化检查
```bash
python3 tools/check_blog_quality.py <主文档.md>
python3 tools/check_blog_quality.py <附录.md>  # 如有附录
```

### 步骤5.2：手动检查清单
参考 `quality/manual_checklist.md` 逐项验证

### ✅ Checkpoint
- [ ] 自动化检查无错误
- [ ] 手动检查清单全部通过
- [ ] TodoWrite标记为completed

---

## 阶段6：图片最终验证（必须执行）⭐

### 步骤6.1：自动化图片验证
```bash
# 验证所有图片编号、文件和图注
python3 tools/verify_blog_figures.py <主文档.md>
python3 tools/verify_blog_figures.py <附录.md>  # 如有附录
```

### 步骤6.2：手动对照原文PDF
逐一确认所有图片：

1. **编号正确性**
   ```bash
   # 使用search_pdf_text.py确认图注在同一页面
   python3 tools/search_pdf_text.py paper.pdf "Figure 1." 3
   python3 tools/search_pdf_text.py paper.pdf "Scheme 1." 3
   ```

2. **内容一致性**
   - 检查提取的图片内容是否与原文图注描述一致
   - 确认没有提取错误的图片（如Figure和Scheme混淆）
   - 确认没有遗漏重要图片

3. **图注准确性**
   - 验证中文翻译是否准确
   - 检查子图说明是否完整
   - 检查颜色说明是否完整

4. **文件完整性**
   ```bash
   # 检查所有图片引用
   grep -E "!\[fig|!\[scheme" <主文档.md>
   ```

### 步骤6.3：常见错误检查

- [ ] **图注位置验证**：图注文字必须与图片在同一页面（唯一判断标准）
- [ ] **Scheme特殊性**：Scheme图注文字必须与图片在同一页
- [ ] **混排情况**：同一页可能有Figure和Scheme，需根据图注位置判断
- [ ] **不存在的图**：如原文只有文字描述但无大图，不要创建该图片
- [ ] **编号连续性**：Figure/Scheme各自编号应连续（1,2,3...）

### ✅ Checkpoint
- [ ] `verify_blog_figures.py` 无错误
- [ ] 所有图片编号与原文一致
- [ ] 所有图片文件存在且内容正确
- [ ] 所有图注完整准确
- [ ] TodoWrite标记为completed

---

## 最终验证

- [ ] 所有阶段TODO已标记为completed
- [ ] 文章可以正常渲染
- [ ] 图片显示正常
- [ ] 公式渲染正确
- [ ] Mermaid图显示正确
