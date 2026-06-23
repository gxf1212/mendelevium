# posts注意事项和frontmatter

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
- frontmatter的date: "2025-08-22"不是文章发表的时间，而是初次写blog的时间，last_modified_at是最后一次修改的日期
- **必须同时添加last_modified_at字段**：由于jekyll-sitemap插件对html_pages只识别`last_modified_at`，所有博客文章必须同时包含`date`和`last_modified_at`字段，否则sitemap会显示错误的1900-01-01日期。新文章创建时必须同时添加这两个字段。
  ```
  ---
  date: "2026-06-20"
  last_modified_at: 2026-06-20
  ---
  ```
- 能不能随机替换assets\img\thumbnail、assets\img\thumbnail_mine下的文件，雨露均沾，创建新文件的时候随机选取。**随机缩略图选择**: 已创建工具脚本来自动随机选择缩略图，避免过度使用bricks.webp。`tools/random_thumbnail.py`: 随机选择缩略图。image和thumbnail这俩参数得一样啊。

# 常见要求

## 常见要求

- 根据@CLAUDE.md，为@_pages\Free Energy\fep-the-end-of-parameter-tuning.pdf写一篇推送Markdown文章
- 根据文件_pages中.md的最后修改时间，给文件名添加YYYY-MM-DD才能在博客上显示（rename就行），还要仿照已有的添加frontmatter，tag尽量用和其他类似的，如果有的话。已经有YYYY-MM-DD的就不用了，这种一般frontmatter也都有了。archive里面的，还有about.md，index这种不要改；Diary里面的不要搞太复杂的文件名，就YYYY-MM-DD-diary1.md这种就行，其他文件名可能还是对title的一个概括。直接扫描所有的_pages中的文件，看哪些没有YYYY-MM-DD，然后修复其frontmatter和文件名，不要搞太复杂的流程。
- 不对，不一定需要添加YYYY-MM-DD，所以给没添加frontmatter的加上frontmatter，tag尽量用和其他类似的，如果有的话。Diary里面的，archive里面的，还有about.md，index这种不要改
- 去掉所有参考标记，如[cite\_start]、[cite: 3749]，根据CLAUDE.md修复所有格式
- 注意title要用中文标点，不要用英文标点，尤其是英文引号，请修复。类似于title: "皮肤屏障的“水之道”：分子模拟揭示脂质相共存如何稳定间质水"，内部是中文标点。【OpenFE】这种是可以的。。但文件名不要出现乱七八糟的东西，标点符号等。
- 把_pages\xxx.docx改写成Markdown，内容完全不变，保证参考文献链接仍然正确，即文内的引用标记要能对得上文末的list，得用这种吧：[^6]。
- 写一篇【科研简讯】栏目的，200行左右，框架类似，内容精简，每一行内容可以多点。挑重点、关键方法和结果来写，太不符合逻辑主线可以不放了
- 写一篇【科研快讯】栏目的，125行以内，超级精简版本，只保留最核心的内容
