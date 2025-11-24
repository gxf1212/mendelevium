# Frontmatter规范

## 标准模板

```yaml
---
title: "文章标题（中文，使用中文标点）"
date: "YYYY-MM-DD"  # 博客创建/修改日期，非论文发表日期
tags: [tag1, tag2, tag3]  # 尽量复用已有tags
description: "简短描述"
image: "/assets/img/thumbnail/xxx.jpg"
thumbnail: "/assets/img/thumbnail/xxx.jpg"  # 必须与image一致
author: Xufan Gao
lang: zh-CN
---
```

## 各字段要求

### title
- **必须使用中文标点**（包括引号）
- 引人入胜但不标题党
- 一语道破核心卖点

**✅ 好的标题**：
```yaml
title: "DFT/MM揭示PETase酶降解PET塑料的反应机理：从原子到动力学的完整图景"
```

**❌ 差的标题**：
```yaml
title: "颠覆性突破！革命性PETase研究震惊学术界"  # 标题党
title: "PETase酶的研究"  # 太笼统
```

### date
- **必须是今天的日期**（博客创建日期）
- **不是论文发表日期**
- 格式：YYYY-MM-DD

**✅ 正确**：
```yaml
date: "2025-11-22"  # 今天写博客的日期
```

**❌ 错误**：
```yaml
date: "2021-09-03"  # 论文发表日期
```

### tags
- 尽量复用已有标签
- 3-10个标签
- 使用小写和连字符

**✅ 好的标签**：
```yaml
tags: [petase, dft-mm, qm-mm, umbrella-sampling, serine-hydrolase, pet-degradation]
```

**❌ 差的标签**：
```yaml
tags: [PETase, DFT/MM, 伞形采样]  # 大写、斜杠、中文
```

### description
- 1-2句话，简洁明了
- 概括文章核心内容

### image & thumbnail
- **必须完全一致**
- 使用`tools/random_thumbnail.py`生成
- **禁止使用empty.jpg**

**✅ 正确**：
```yaml
image: "/assets/img/thumbnail_mine/wh-z8p9rj.jpg"
thumbnail: "/assets/img/thumbnail_mine/wh-z8p9rj.jpg"
```

**❌ 错误**：
```yaml
image: "/assets/img/thumbnail/empty.jpg"  # 使用empty.jpg
thumbnail: "/assets/img/thumbnail_mine/wh-z8p9rj.jpg"  # 不一致
```

### author
固定值：`Xufan Gao`

### lang
固定值：`zh-CN`

## 生成流程

### 1. 生成随机缩略图
```bash
cd /mnt/e/GitHub-repo/mendelevium
python3 tools/random_thumbnail.py --frontmatter
```

### 2. 复制输出并填充其他字段
工具会输出：
```
  thumbnail: "/assets/img/thumbnail_mine/wh-z8p9rj.jpg"
  image: "/assets/img/thumbnail_mine/wh-z8p9rj.jpg"
```

### 3. 手动填充其他字段
- title: 根据文章内容起标题
- date: **今天的日期**
- tags: 根据文章主题选择
- description: 写1-2句话概括

## 检查清单

- [ ] title使用中文标点
- [ ] date是今天的日期（不是论文日期）
- [ ] image和thumbnail完全一致
- [ ] 未使用empty.jpg
- [ ] tags全部小写，使用连字符
- [ ] description简洁明了
- [ ] author和lang正确
