# wechat_html

把博客 Markdown 转成**可直接粘进微信公众号编辑器**的 HTML。所有样式以
内联方式写在微信会保留的标签上（`<h1>~<h6>`、`<p>`、`<span>`、`<strong>`、
`<blockquote>`、`<pre>`），绕开 ueditor 对 `<div>`/`<ul>`/`<li>` 的容器重写，
并修复「加粗词：正文」在微信里被孤行断开的问题。

## 运行环境

脚本依赖 Pillow，**必须用 venv 解释器**运行（裸 managed python 没有 PIL）：

```bash
PY="C:/Users/Lenovo/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
```

（仅首次需要建环境：`python -m venv ...\envs\default` 后 `pip install pillow`）

## 基本用法

```bash
cd E:/GitHub-repo/mendelevium/tools/wechat_html
"$PY" build_wechat_html.py "E:/GitHub-repo/mendelevium/_pages/.../2026-07-22-xxx.md"
```

- 默认输出：`tools/wechat_html/wechat_preview/xxx.wechat.html`（脚本同级目录，已 git 忽略，不误提交）
- 用 `--out 路径` 指定输出位置（如临时对比多版时）
- 生成后用浏览器打开，全选 → 复制 → 粘进公众号编辑器即可

### 常用开关

| 参数 | 作用 |
|---|---|
| `--banner` | 渲染顶部横版 banner 图（默认关闭） |
| `--subtitle "文字"` | banner 副标题（默认「东山月光下 · 科研精读」） |
| `--no-footer` | 不渲染文末「关注」引导卡 |
| `--math svg` | 公式渲染方式（默认 `svg` 矢量内联，清晰不糊；`--math png` 则栅格图兜底） |
| `--quality 82` | 图片 JPEG 质量（默认 82） |
| `--max-img-w 920` | 图片最大宽度（默认 920） |
| `--diagram a.png b.png` | 按顺序替换文中的 ```mermaid 块（无需此参数时自动渲染） |

## 主题（配色 / 字号与代码分离）

主题是一份 JSON（`themes/*.json`），含 `colors` / `headings` / `body` /
`blockquote` / `list` / `table` / `code` / `image` / `layout` / `styles` 字段。
用 `--theme` 选择：

```bash
"$PY" build_wechat_html.py --theme themes/dreamy.json "路径.md"
```

内置主题（共 10 套，均浅色）：

| 文件 | 名称 |
|---|---|
| `default.json` | 经典蓝（原样式，默认） |
| `dreamy.json` | 柔光月夜（奶白底、蓝紫粉） |
| `shanchui.json` | 山吹（暖黄） |
| `rose.json` | 蔷薇紫 |
| `fullstack-blue.json` | 全栈蓝 |
| `cute-green.json` | 萌绿 |
| `frontend.json` | 前端之巅（青蓝） |
| `orange-blue.json` | 橙蓝风 |
| `cute-pink.json` | 萌粉 |
| `black-white.json` | 极简（黑白灰） |

加新主题：往 `gen_themes.py` 的 `PALETTES` 加一个色板 dict，重跑
`python gen_themes.py` 即可批量生成，不用手写整份 JSON。

## 组件样式变体

每个组件都有「简单版 + 花哨版」，可显式指定，或 `--random` 一键随机
（显式参数优先于随机）。

| 组件 | 参数 | 可选值 |
|---|---|---|
| 标题 h2/h3 | `--heading-style` | `plain` · `bar` · `pill` · `number`（h2 与 h3 各自渲染、层级已拉开：bar 下 h2 粗块+满宽底 / h3 细条无底；pill 下 h2 深底白字 / h3 浅底小标签；number 下 h2 深底大徽章 / h3 浅底小徽章） |
| 有序列表 | `--list-style` | `plain`（裸 `N.`）· `circle`（圆形徽章） · `square`（圆角方块徽章） |
| 引用块 | `--quote-style` | `bar`（左竖线+浅底） · `quote`（大引号「"」装饰） |

```bash
# 全指定
"$PY" build_wechat_html.py --theme themes/dreamy.json \
    --heading-style pill --list-style square --quote-style quote "路径.md"

# 全部随机
"$PY" build_wechat_html.py --theme themes/dreamy.json --random "路径.md"

# 序号固定圆形，其余随机
"$PY" build_wechat_html.py --theme themes/dreamy.json --list-style circle --random "路径.md"
```

> 所有新样式只用 `background`+`border-radius`+`padding` 内联实现，
> 刻意避开 div 定位 / 伪元素 / 渐变 / 阴影（会被微信剥掉）。
> 构建结束会打印实际选中的样式组合，方便确认随机结果。

## 文件说明

- `build_wechat_html.py` —— 主脚本（MD → 公众号 HTML）
- `gen_themes.py` —— 批量生成主题 JSON（改色板后重跑）
- `render_math*.py / *.js` —— 公式渲染（LaTeX → PNG/SVG）
- `render_mermaid.py` —— mermaid 流程图渲染
- `poems.txt` —— 文末随机诗词库（带作者《标题》）
- `themes/` —— 主题 JSON 目录
