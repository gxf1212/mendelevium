# Mendelevium Blog - Claude Development Guide

## 项目理念 / Project Philosophy

This is a personal research blog sharing experiences in molecular dynamics, computational chemistry, and scientific computing. The content should be:

- **Educational**: 每篇文章都应该能帮助读者学到一些实用的知识
- **Practical**: 重点在于实用性，分享实际的方法和工具
- **Clear**: 保持简洁明了的写作风格，便于理解

## 开发原则 / Development Principles

1. **保持简洁**: 网站结构应该简单明了，不要过度复杂化
2. **内容为王**: 专注于高质量的技术内容，而非花哨的功能
3. **分类清晰**: 按研究领域合理组织内容结构

## 技术栈 / Tech Stack

- Jekyll (Satellite theme)
- GitHub Pages
- Giscus for comments
- GoatCounter for analytics

## 内容组织 / Content Organization

```
_pages/
├── Drug Design/          # 药物设计相关研究
├── Machine Learning & AI/    # 机器学习在分子科学中的应用
├── Molecular Dynamics/       # 分子动力学模拟研究
├── MD Tools/                # MD软件工具和实用程序  
├── Techniques/              # 通用技术和方法
├── Field Knowledge/         # 领域基础知识
└── Nano Polymers/          # 纳米聚合物研究
```

- Archive\jekyll-theme-satellite-master是原始模板，报错了可以参考，尤其是里面的docs
- 为什么单篇文章，如2025-08-13-deep-covboost-ai-covid-target会出现在侧边栏？不应该出现的
    bookmark: false我已经删掉了，就解决了，以后都不要写
- 常见要求：根据文件_pages中.md的最后修改时间，给文件名添加YYYY-MM-DD才能在博客上显示（rename就行），还要仿照已有的添加frontmatter，tag尽量用和其他类似的，如果有的话。已经有YYYY-MM-DD的就不用了，这种一般frontmatter也都有了。about.md，index这种不要改。


## 注意事项 / Notes

- 图片路径使用相对路径，确保在不同环境下都能正常显示
- 保持中英文混用的写作风格，适合中文科研环境
- 定期整理和归档过时的内容
- 公式格式大概如下，$$的上下要和文字空一行：

    我们的最终目标是得到**无偏的自由能** $F_h(\xi)$，它与**无偏概率分布** $\rho_h(\xi)$ 的关系由统计力学的基本公式定义：

    $$
    F_{h}(\xi) = -k_B T \ln \rho_{h}(\xi)
    $$

    其中，$k_B$ 是玻尔兹曼常数...

- should add frontmatter for each post (identify .md files without ---), e.g.
---
title: "Random Forest and Enhanced Sampling Unite: Revealing and Correcting Ghost Errors in Alchemical Free Energy Calculations"
date: "2025-08-22"
tags: [random-forest, enhanced-sampling, alchemical-free-energy, gamd, error-analysis, machine-learning, molecular-dynamics]
---
  