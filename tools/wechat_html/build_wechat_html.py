# -*- coding: utf-8 -*-
"""
build_wechat_html.py — 博客 Markdown 转「公众号粘贴版」HTML 生成器
==================================================================

把任意一篇博客 Markdown 转换成【单文件、全 inline-CSS、图片 base64 内嵌】的
HTML，粘进公众号草稿编辑器后排版与图片零修改即可。

设计系统（dreamy / 月夜蓝）：
  - 顶部 banner 默认关闭（加 --banner 开启，由 frontmatter title 自动绘制为图片）
  - 正文 16px / 行高 1.75 / 两端缩进
  - H2 左蓝边、H3/H4 蓝色标题、引用圆角蓝边、代码深底
  - 表格：表头背景 background-color + bgcolor 双写，规避 WeChat 粘贴丢背景
  - 文末「关注」引导卡：绘制为图片（避免纯文字复制到草稿里变回文字）

依赖（隔离 venv）：
  pip install Pillow markdown
字体：C:/Windows/Fonts/msyh.ttc（微软雅黑）

用法：
  python build_wechat_html.py <文章.md> [--out out.html]
                               [--subtitle "东山月光下 · 科研精读"]
                               [--diagram diag1.png [diag2.png ...]]
                               [--banner] [--no-footer]
                               [--quality 82] [--max-img-w 920]

说明：
  - 相对路径图片（如 jcim6c01288_figs/fig1.png）自动以 base64 内嵌；
    http(s) 图片与 data: 图片保持原样。
  - ```mermaid 代码块：优先用 mmdc 渲染；若不可用则调用内置 Pillow 渲染器解析
    graph/flowchart + subgraph + 箭头链并转成图片；均失败时生成占位图。公众号不
    渲染 Mermaid/LaTeX，公式建议预先转成图片再引用。
  - 微信对「个人主体未认证」公众号已收回全部 raw API（建草稿/发布/群发/素材），
    但手动粘贴与第三方 OAuth 工具（秀米/壹伴）不受影响。
"""
import argparse
import base64
import html as html_mod
import io
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile

from PIL import Image, ImageDraw, ImageFont
import markdown

# 引入同目录下的 Mermaid 解析渲染器（无需 mmdc/Chromium）
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)
try:
    from render_mermaid import render_to_bytes as render_mermaid_png
except Exception:
    render_mermaid_png = None
try:
    from render_math import render_math, inline_img_style
except Exception:
    render_math, inline_img_style = None, None
try:
    from render_math_svg import render_math_svg_batch
except Exception:
    render_math_svg_batch = None

# ===== 路径 / 字体 =====
FONT_R = "C:/Windows/Fonts/msyh.ttc"
FONT_B = "C:/Windows/Fonts/msyhbd.ttc"

# ===== 设计系统配色（dreamy / 月夜蓝）=====
C = {
    "blue":     "#5b8def",   # 主色：月夜蓝
    "blue_d":   "#3a6bc4",   # 深蓝（标签/边框）
    "coral":    "#e07a5f",   # 暖珊瑚（强调）
    "blue_bg":  "#eef4fb",   # 浅蓝底（交替行/引用）
    "ink":      "#3f3f3f",   # 正文
    "ink_d":    "#2c3e50",   # 标题
    "ink_h3":   "#3a5a8c",   # H3
    "ink_h4":   "#5a6b85",   # H4
    "gray":     "#888888",   # 图注
    "border":   "#d6e4f7",   # 浅蓝边框
    "quote_bg": "#f4f8fd",   # 引用底
    "code_bg":  "#2c3e50",   # 代码底
    "code_fg":  "#e6edf3",   # 代码字
}

FONT_STACK = "-apple-system,'Segoe UI','Microsoft YaHei',sans-serif"
MONO_STACK = "Consolas,'Courier New',monospace"


# ===== 字体 / 绘图工具 =====
def font(size, bold=False):
    try:
        return ImageFont.truetype(FONT_B if bold else FONT_R, size)
    except Exception:
        return ImageFont.load_default()


def wrap(draw, text, f, max_w):
    lines, cur = [], ""
    for ch in text:
        if draw.textlength(cur + ch, font=f) <= max_w:
            cur += ch
        else:
            lines.append(cur)
            cur = ch
    if cur:
        lines.append(cur)
    return lines


# ===== 图片 → base64（JPEG 压缩控制体积）=====
def img_to_b64(path, max_w=1000, quality=82):
    im = Image.open(path).convert("RGB")
    if im.width > max_w:
        h = int(im.height * max_w / im.width)
        im = im.resize((max_w, h))
    buf = io.BytesIO()
    im.save(buf, format="JPEG", quality=quality, optimize=True)
    return base64.b64encode(buf.getvalue()).decode()


# ===== 横版 dreamy banner（由标题绘制）=====
def render_banner(title, subtitle=""):
    W, H = 1080, 420
    base = Image.new("RGB", (W, H), (33, 46, 84))
    d = ImageDraw.Draw(base)
    top, bot = (33, 46, 84), (91, 141, 239)
    for y in range(H):
        t = y / (H - 1)
        d.line([(0, y), (W, y)],
               fill=tuple(int(top[i] + (bot[i] - top[i]) * t) for i in range(3)))
    ov = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    od = ImageDraw.Draw(ov)
    mx, my = 880, 110
    for rr in range(130, 44, -6):
        a = int(45 * (1 - (130 - rr) / 86))
        od.ellipse([mx - rr, my - rr, mx + rr, my + rr], fill=(250, 246, 227, a))
    random.seed(7)
    for _ in range(30):
        bx, by = random.randint(0, W), random.randint(0, H)
        br = random.randint(6, 30)
        od.ellipse([bx - br, by - br, bx + br, by + br],
                  fill=(255, 255, 255, random.randint(18, 60)))
    for _ in range(50):
        sx, sy = random.randint(0, W), random.randint(0, H // 2)
        od.ellipse([sx - 1, sy - 1, sx + 1, sy + 1],
                  fill=(255, 255, 255, random.randint(80, 200)))
    base = Image.alpha_composite(base.convert("RGBA"), ov).convert("RGB")
    d = ImageDraw.Draw(base)
    # 副标题置于顶部
    if subtitle:
        d.text((60, 66), subtitle, font=font(26, True), fill=(214, 226, 250))
    # 标题自适应字号：优先单行；过长则均衡换行（避免长标题被丑陋地折行）
    fs = 42
    while fs > 26 and d.textlength(title, font=font(fs, True)) > W - 120:
        fs -= 2
    if d.textlength(title, font=font(fs, True)) <= W - 120:
        lines = [title]
    else:
        lines = wrap(d, title, font(fs, True), W - 120)
    lh = fs + 12
    block_h = len(lines) * lh
    ty = (H - block_h) / 2 + 8
    for i, line in enumerate(lines):
        d.text((60, ty + i * lh), line, font=font(fs, True), fill="white")
    buf = io.BytesIO()
    base.save(buf, format="JPEG", quality=85, optimize=True)
    return base64.b64encode(buf.getvalue()).decode()


# ===== 文末「关注」引导卡（绘制为图片）=====
def render_footer(title="东山月光下", poems_file=None, fallback_poems=(
        "醉后不知天在水，满船清梦压星河。——唐温如《题龙阳县青草湖》",
        "采菊东篱下，悠然见南山。——陶渊明《饮酒·其五》",
        "行到水穷处，坐看云起时。——王维《终南别业》")):
    if poems_file is None:
        poems_file = os.path.join(_SCRIPT_DIR, "poems.txt")
    poems = list(fallback_poems)
    try:
        with open(poems_file, encoding="utf-8") as fh:
            for line in fh.read().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    poems.append(line)
    except FileNotFoundError:
        pass
    W, H = 1080, 300
    slogan = random.choice(poems)
    img = Image.new("RGB", (W, H), (33, 46, 84))
    d = ImageDraw.Draw(img)
    top, bot = (33, 46, 84), (91, 141, 239)
    for y in range(H):
        t = y / (H - 1)
        d.line([(0, y), (W, y)],
               fill=tuple(int(top[i] + (bot[i] - top[i]) * t) for i in range(3)))
    ov = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    od = ImageDraw.Draw(ov)
    mx, my = 930, 70
    for rr in range(80, 34, -5):
        a = int(50 * (1 - (80 - rr) / 46))
        od.ellipse([mx - rr, my - rr, mx + rr, my + rr], fill=(250, 246, 227, a))
    random.seed(11)
    for _ in range(26):
        bx, by = random.randint(0, W), random.randint(0, H)
        br = random.randint(5, 26)
        od.ellipse([bx - br, by - br, bx + br, by + br],
                  fill=(255, 255, 255, random.randint(15, 55)))
    img = Image.alpha_composite(img.convert("RGBA"), ov).convert("RGB")
    d = ImageDraw.Draw(img)
    d.text((70, 70), f"关注「{title}」", font=font(40, True), fill="white")
    for i, line in enumerate(wrap(d, slogan, font(24, False), W - 200)):
        d.text((72, 140 + i * 38), line, font=font(24, False), fill=(220, 232, 255))
    buf = io.BytesIO()
    img.save(buf, format="JPEG", quality=86, optimize=True)
    return base64.b64encode(buf.getvalue()).decode()


# ===== mermaid 占位图（无 mmdc 且无 --diagram 时）=====
def render_placeholder(text):
    W, H = 720, 160
    img = Image.new("RGB", (W, H), (244, 248, 253))
    d = ImageDraw.Draw(img)
    d.rectangle([0, 0, W - 1, H - 1], outline=C["border"], width=2)
    for i, line in enumerate(wrap(d, text, font(18, True), W - 60)):
        d.text((30, 50 + i * 30), line, font=font(18, True), fill=C["blue_d"])
    buf = io.BytesIO()
    img.save(buf, format="JPEG", quality=85)
    return base64.b64encode(buf.getvalue()).decode()


def try_mmdc(content, workdir):
    """尝试用 mermaid-cli 渲染 mermaid 文本为 PNG，返回路径或 None。"""
    mmdc = shutil.which("mmdc")
    if not mmdc:
        return None
    try:
        fd_in = tempfile.NamedTemporaryFile("w", suffix=".mmd", dir=workdir,
                                            delete=False, encoding="utf-8")
        fd_out = tempfile.NamedTemporaryFile(suffix=".png", dir=workdir,
                                             delete=False)
        fd_in.write(content)
        fd_in.close()
        fd_out.close()
        subprocess.run([mmdc, "-i", fd_in.name, "-o", fd_out.name],
                       check=True, capture_output=True, timeout=120)
        if os.path.getsize(fd_out.name) > 0:
            return fd_out.name
    except Exception:
        return None
    return None


# ===== frontmatter 解析 =====
def parse_frontmatter(text):
    if not text.startswith("---"):
        return {}, text
    parts = text.split("---", 2)
    if len(parts) < 3:
        return {}, text
    fm_raw, body = parts[1], parts[2]
    fm = {}
    for line in fm_raw.splitlines():
        m = re.match(r'\s*(\w+)\s*:\s*"?([^"]*)"?\s*$', line)
        if m:
            fm[m.group(1)] = m.group(2)
    return fm, body


# ===== 公式：提取 / 还原 + 代码块修复 =====
def stash_math(body):
    """markdown 转换前提取 $...$ / $$...$$ 为占位 token，避免被 markdown 破坏。
    代码围栏内的 $ 不当公式。返回 (stashed_body, token_map)。
    token_map: token -> (tex, display, raw)
    """
    code_blocks = []
    def stash_code(m):
        code_blocks.append(m.group(0))
        return "@@CODE%d@@" % (len(code_blocks) - 1)
    b = re.sub(r"```.*?```", stash_code, body, flags=re.S)
    b = re.sub(r"(?s)~~~.*?~~~", stash_code, b)

    token_map = {}
    # 块公式 $$...$$（可跨行）
    def blk(m):
        key = "@@MATHB%d@@" % len(token_map)
        token_map[key] = (m.group(1).strip(), True, m.group(0).strip())
        return key
    b = re.sub(r"\$\$([\s\S]+?)\$\$", blk, b)
    # 行内公式 $...$（排除 $$ 已处理）。博客中 $ 专用于数学；
    # 仅当内容像「纯货币/数字」（如 $5、$1,200）时才跳过后保留原文，
    # 其余一律当作数学公式渲染（含 $>0.4$、单一变量 $x$ 等）。
    def inl(m):
        tex = m.group(1).strip()
        # 明显是货币/纯数字（如 $5、$1,200）才跳过，避免把金额当公式
        if re.fullmatch(r"[0-9][0-9,\. ]*", tex):
            return m.group(0)
        key = "@@MATHI%d@@" % len(token_map)
        token_map[key] = (tex, False, m.group(0).strip())
        return key
    b = re.sub(r"(?<!\$)\$(?!\$)([^$\n]+?)\$(?!\$)", inl, b)
    for i, c in enumerate(code_blocks):
        b = b.replace("@@CODE%d@@" % i, c)
    return b, token_map


def restore_math(html, token_map, math_mode="png"):
    """把占位 token 还原为公式。math_mode:
      - "png": 渲染为 base64 PNG <img>（栅格，字号靠缩放系数凑）
      - "svg": 渲染为内联矢量 <svg>（与正文同字号、可缩放、不糊），
               单条失败时回退 PNG，再失败保留原文。
    """
    if not token_map:
        return html
    # 收集顺序，便于批量调 SVG
    toks = list(token_map.items())
    if math_mode == "svg" and render_math_svg_batch is not None:
        items = [(tex, display) for _, (tex, display, _) in toks]
        svgs = render_math_svg_batch(items)
    else:
        svgs = [None] * len(toks)

    for (tok, (tex, display, raw)), svg in zip(toks, svgs):
        if svg is not None:
            repl = svg
        elif render_math is not None:
            res = render_math(tex, display)
            if res is None:
                repl = raw
            elif display:
                repl = (f'<img src="data:image/png;base64,{res[0]}" '
                        f'style="display:block;margin:16px auto;'
                        f'max-width:96%;vertical-align:middle;">')
            else:
                repl = (f'<img src="data:image/png;base64,{res[0]}" '
                        f'style="{inline_img_style(res[2])}">')
        else:
            repl = raw
        html = html.replace(tok, repl)
    return html


def fix_code_blocks(html):
    """代码块 <pre><code> 内部不应套行内 code 的浅蓝底，改为透明底继承 pre 深色。"""
    def repl(m):
        block = m.group(0)
        block = re.sub(
            r'^<pre\b([^>]*)>',
            r'<pre\1 style="background-color:#2c3e50;color:#e6edf3;'
            r'padding:12px 14px;border-radius:8px;overflow-x:auto;'
            r'font-size:13px;line-height:1.6;margin:12px 0;">',
            block, count=1)
        block = re.sub(
            r'<code\b([^>]*)>',
            r'<code\1 style="background:transparent;color:inherit;padding:0;'
            r'font-family:Consolas,\'Courier New\',monospace;">',
            block)
        return block
    return re.sub(r'<pre\b[^>]*>.*?</pre>', repl, html, flags=re.S)


# ===== 图片内嵌 =====
def embed_images(html, base_dir, quality):
    def repl(m):
        tag = m.group(0)
        alt_m = re.search(r'alt="([^"]*)"', tag)
        src_m = re.search(r'src="([^"]*)"', tag)
        if not src_m:
            return tag
        alt = alt_m.group(1) if alt_m else ""
        src = src_m.group(1)
        if src.startswith(("http://", "https://", "data:")):
            return tag
        path = os.path.join(base_dir, src)
        if not os.path.exists(path):
            return tag
        b64 = img_to_b64(path, quality=quality)
        return (f'<img alt="{alt}" src="data:image/jpeg;base64,{b64}" '
                f'style="max-width:96%;display:block;margin:0 auto;'
                f'border:1px solid #eee;border-radius:6px;">')
    return re.sub(r'<img\s+[^>]*?/?>', repl, html)


# ===== 表格：背景双写（规避 WeChat 丢背景）=====
def style_tables(html):
    def table_repl(m):
        tbl = m.group(0)
        # 注意负向先行断言 (?![a-zA-Z])：只匹配 <th> / <th ...>，
        # 不能匹配 <thead>（否则 <thead> 会被错误替换成游离的 <th>，
        # 在表头上方多出一行空单元格）
        tbl = re.sub(
            r'<th(?![a-zA-Z])([^>]*)>',
            (f'<th style="background-color:{C["blue"]};color:#ffffff;'
             f'border:1px solid {C["border"]};padding:9px 8px;'
             f'font-weight:bold;text-align:left;" bgcolor="{C["blue"]}">'),
            tbl)
        # 数据行斑马纹
        def tbody_repl(tm):
            body = tm.group(0)
            rows = re.findall(r'<tr>(.*?)</tr>', body, re.S)
            out = body
            for i, row in enumerate(rows):
                bg = C["blue_bg"] if i % 2 == 1 else "#ffffff"
                new = re.sub(
                    r'<td([^>]*)>',
                    (f'<td style="border:1px solid {C["border"]};'
                     f'padding:8px 8px;color:{C["ink"]};'
                     f'background-color:{bg};" bgcolor="{bg}">'),
                    row)
                out = out.replace(row, new, 1)
            return out
        tbl = re.sub(r'<tbody>.*?</tbody>', tbody_repl, tbl, flags=re.S)
        tbl = re.sub(r'<table([^>]*)>',
                     r'<table style="border-collapse:collapse;'
                     r'width:100%;font-size:16px;margin:14px 0;">', tbl)
        return tbl
    return re.sub(r'<table>.*?</table>', table_repl, html, flags=re.S)


# ===== mermaid 代码块渲染 =====
def render_mermaid(html, base_dir, diagrams, quality):
    blocks = re.findall(r'<pre><code class="language-mermaid">(.*?)</code></pre>',
                        html, re.S)
    for i, raw_content in enumerate(blocks):
        # markdown 库会把代码块里的 " & < > 等转义成 HTML 实体，
        # 渲染前必须还原；但替换回 HTML 时要用原始的转义字符串，否则 replace 匹配不上。
        content = html_mod.unescape(raw_content)
        b64 = None
        note = ""
        # 1) 显式 --diagram 覆盖
        if i < len(diagrams) and diagrams[i] and os.path.exists(diagrams[i]):
            b64 = img_to_b64(diagrams[i], quality=quality)
        # 2) mermaid-cli
        if b64 is None:
            png = try_mmdc(content, base_dir)
            if png:
                b64 = img_to_b64(png, quality=quality)
        # 3) 内置 Pillow 渲染器
        if b64 is None and render_mermaid_png:
            b64, _note = render_mermaid_png(content, quality=quality)
        # 4) 占位图
        if b64 is None:
            b64 = render_placeholder("流程图（Mermaid）请预渲染为图片后插入")
        repl = (f'<p style="text-align:center;margin:6px 0 2px;">'
                f'<img src="data:image/jpeg;base64,{b64}" '
                f'style="max-width:96%;display:block;margin:0 auto;'
                f'border:1px solid #eee;border-radius:6px;"></p>'
                f'<p style="text-align:center;color:{C["gray"]};'
                f'font-size:13px;margin:2px 0 14px;">{note}</p>')
        html = html.replace(f'<pre><code class="language-mermaid">{raw_content}</code></pre>',
                            repl, 1)
    return html


# ===== 行内样式注入 =====
def apply_inline_styles(html):
    subs = [
        (r'<p>', '<p style="margin:12px 0;color:#3f3f3f;">'),
        (r'<h1>', '<h1 style="font-size:22px;font-weight:bold;color:#2c3e50;'
                  'border-bottom:3px solid #5b8def;padding-bottom:10px;'
                  'margin:18px 0 14px;">'),
        (r'<h2>', '<h2 style="font-size:19px;font-weight:bold;color:#2c3e50;'
                  'border-left:4px solid #5b8def;padding-left:10px;'
                  'margin:26px 0 12px;">'),
        (r'<h3>', '<h3 style="font-size:17px;font-weight:bold;color:#3a5a8c;'
                  'margin:22px 0 10px;">'),
        (r'<h4>', '<h4 style="font-size:16px;font-weight:bold;color:#5a6b85;'
                  'margin:16px 0 8px;">'),
        (r'<h5>', '<h5 style="font-size:16px;font-weight:bold;color:#7a8aa0;'
                  'margin:14px 0 6px;">'),
        (r'<h6>', '<h6 style="font-size:16px;font-weight:bold;color:#9aa7ba;'
                  'margin:12px 0 6px;">'),
        (r'<ul>', '<ul style="padding-left:22px;margin:8px 0;">'),
        (r'<ol>', '<ol style="padding-left:22px;margin:8px 0;">'),
        (r'<li>', '<li style="margin:8px 0;">'),
        (r'<blockquote>',
         '<blockquote style="border-left:4px solid #5b8def;'
         'background-color:#f4f8fd;padding:12px 16px;color:#5a6b85;'
         'margin:14px 0;border-radius:0 8px 8px 0;">'),
        (r'<pre>',
         f'<pre style="background-color:{C["code_bg"]};color:{C["code_fg"]};'
         'padding:12px 14px;border-radius:8px;overflow-x:auto;font-size:13px;'
         'line-height:1.6;margin:12px 0;">'),
        (r'<code>',  # 仅行内 code（无属性）；代码块内的 code 颜色继承自 pre
         f'<code style="font-family:{MONO_STACK};background-color:#eef4fb;'
         f'color:#3a6bc4;padding:1px 5px;border-radius:4px;font-size:13px;">'),
        (r'<strong>', '<strong style="font-weight:bold;color:#2c3e50;">'),
        (r'<em>', '<em style="font-style:italic;color:#5a6b85;">'),
        (r'<a ', '<a style="color:#5b8def;text-decoration:none;'
                  'border-bottom:1px solid #bcd0f0;" '),
        (r'<sup\b', '<sup style="font-size:0.72em;color:#5b8def;'
                    'font-weight:bold;" '),
        (r'<div class="footnote[s]?"',
         r'<div class="footnote" style="border-top:1px solid #d6e4f7;'
         r'margin-top:24px;padding-top:12px;font-size:13px;color:#888;'
         r'line-height:1.7;">'),
        (r'<hr\s*/?>', '<hr style="border:none;border-top:1px solid #d6e4f7;'
                       'margin:24px 0;">'),
    ]
    for pat, rep in subs:
        html = re.sub(pat, rep, html)
    return html


def flatten_lists(html):
    """微信编辑器基于 ueditor：粘贴时会把 <ul>/<ol> 强制转成 <p>，且 ueditor
    的「进入编辑器的 li 要套 p 标签」逻辑（ueditor.all.js 约 15112 行）会遍历
    <li> 的 children 重新打包。对 <li><strong>词</strong>句子</li> 这种
    「加粗元素 + 裸文本」混合结构，会把 strong 与后面的文本拆成两个独立块，
    导致「加粗词：句子」在粗体结束后孤行换行（与字数无关）。
    之前所有「在 li 内打补丁」的方案（冒号并进 strong / span 替代 / 零宽空格）
    都无效，因为 ueditor 的 li→p 转换是结构性重写，会吞掉补丁。
    彻底绕开：在生成阶段就把列表展开为一组 <p>，无序用「• 」、有序用「N. 」前缀；
    <p> 是微信主力标签，<p> 内的 <strong> 安全，不会再被拆块。
    嵌套列表递归展开，子项用 &nbsp; 缩进区分层级（子项作为独立段落提升到父项之后）。"""
    _ULOL_OPEN = re.compile(r'<(ul|ol)\b[^>]*>', re.I)
    _ULOL_CLOSE = re.compile(r'</(ul|ol)>', re.I)
    _LI_OPEN = re.compile(r'<li\b[^>]*>', re.I)
    _LI_CLOSE = re.compile(r'</li>', re.I)

    def match_close(s, i, open_re, close_re):
        """从 i 起平衡匹配，返回与起始标签配对的闭标签索引。"""
        depth = 1
        while True:
            cm = close_re.search(s, i)
            if cm is None:
                return len(s)
            om = open_re.search(s, i)
            if om and om.start() < cm.start():
                depth += 1
                i = om.end()
            else:
                depth -= 1
                if depth == 0:
                    return cm.start()
                i = cm.end()

    def expand_li_content(content, depth):
        """把单个 <li> 内容拆成「直接文本」+「嵌套列表展开的独立段落」。
        嵌套 <ul>/<ol> 先用平衡匹配保护起来，避免被 loose-list 的 <p> 剥离误伤。"""
        prot = []  # 嵌套 <ul>/<ol> 原始块
        parts = []
        pos = 0
        while True:
            m = _ULOL_OPEN.search(content, pos)
            if not m:
                parts.append(content[pos:])
                break
            parts.append(content[pos:m.start()])
            tag = m.group(1)
            close_start = match_close(content, m.end(), _ULOL_OPEN, _ULOL_CLOSE)
            block = content[m.start():close_start + len('</' + tag + '>')]
            prot.append(block)
            pos = close_start + len('</' + tag + '>')
        text = re.sub(r'<p[^>]*>(.*?)</p>', r'\1', ''.join(parts), flags=re.S).strip()
        nested = [expand_region(b, depth + 1) for b in prot]
        return text, nested

    def expand_list(inner, ordered, depth):
        out = []
        idx = 0
        i = 0
        n = len(inner)
        while i < n:
            lm = _LI_OPEN.search(inner, i)
            if not lm:
                break
            content_start = lm.end()
            content_end = match_close(inner, content_start, _LI_OPEN, _LI_CLOSE)
            content = inner[content_start:content_end]
            i = content_end + len('</li>')
            idx += 1
            text, nested = expand_li_content(content, depth)
            bullet = (f'<span style="font-size:18px;line-height:1;'
                      f'font-weight:bold;color:#33312e;">{idx}.</span> '
                      if ordered else
                      f'<span style="font-size:20px;line-height:1;'
                      f'font-weight:bold;color:#33312e;">•</span> ')
            indent = '&nbsp;' * (4 * depth)
            out.append(f'<p style="margin:6px 0;">{indent}{bullet}{text}</p>')
            out.extend(nested)
        return '\n'.join(out)

    def expand_region(s, depth):
        out = []
        pos = 0
        while True:
            m = _ULOL_OPEN.search(s, pos)
            if not m:
                out.append(s[pos:])
                break
            out.append(s[pos:m.start()])
            tag = m.group(1)
            close_start = match_close(s, m.end(), _ULOL_OPEN, _ULOL_CLOSE)
            inner = s[m.end():close_start]
            ordered = (tag == 'ol')
            out.append(expand_list(inner, ordered, depth))
            pos = close_start + len('</' + tag + '>')
        return ''.join(out)

    return expand_region(html, 0)


def fix_bold_colon(html):
    """公众号编辑器会把「**短标签**：正文」里的短标签孤行断开（粗体 run 与
    普通文本间的 run 边界在紧跟全角冒号时易断）。把紧跟在 </strong> 后的
    全角冒号并进粗体 run（</strong>： → ：</strong>），让「标签：」成为一个
    不可分割的整体，消除孤行。仅对 ≤12 字符的短标签生效，避免长句中的冒号被加粗。
    """
    return re.sub(
        r'<strong([^>]*)>([^<]{1,12}?)</strong>：',
        r'<strong\1>\2：</strong>', html)


def fix_bold_break(html):
    """微信 contenteditable 在粘贴时，会把 <strong> 闭合后紧跟正文的「内联边界」
    重写成换行（列表里「加粗词：句子」在粗体一结束就孤行断开）。
    修复：在 </strong> 与后续文本之间插入零宽空格(U+200B, &#8203;) 作为桥接，
    让编辑器把正文视为粗体 run 的延续而非新块，从而避免孤行换行。
    仅对「闭合后紧跟非空白、非标签」的边界插入（如 </strong>句子 / </strong>：）。"""
    return re.sub(r'(</strong>)(?=[^\s<])', r'\1&#8203;', html)


def suppress_dup_title(html, title):
    if not title:
        return html
    m = re.search(r'<h1[^>]*>(.*?)</h1>', html, re.S)
    if m and m.group(1).strip() == title.strip():
        html = html.replace(m.group(0), "", 1)
    return html


def build_doc(banner_b64, body_html, footer_b64, no_banner, no_footer):
    parts = []
    if not no_banner:
        parts.append(
            '<p style="margin:0 0 6px;">'
            f'<img src="data:image/jpeg;base64,{banner_b64}" '
            'style="width:100%;display:block;border-radius:10px;"></p>')
    parts.append(
        f'<div style="font-size:16px;line-height:1.75;color:{C["ink"]};'
        f'font-family:{FONT_STACK};max-width:680px;margin:0 auto;padding:0 8px;">')
    parts.append(body_html)
    parts.append('</div>')
    if not no_footer:
        parts.append(
            '<p style="text-align:center;margin:28px 0 10px;">'
            f'<img src="data:image/jpeg;base64,{footer_b64}" '
            'style="width:100%;max-width:680px;display:block;margin:0 auto;'
            'border-radius:12px;"></p>')
    return "\n".join(parts)


def main():
    ap = argparse.ArgumentParser(description="博客 MD → 公众号粘贴版 HTML")
    ap.add_argument("md", help="博客 Markdown 文件路径")
    ap.add_argument("--out", help="输出 HTML 路径（默认与 md 同目录同名 .wechat.html）")
    ap.add_argument("--subtitle", default="东山月光下 · 科研精读", help="banner 副标题")
    ap.add_argument("--diagram", nargs="*", default=[],
                    help="预渲染的 mermaid 图 PNG，按顺序替换 ```mermaid 块")
    ap.add_argument("--banner", action="store_true",
                    help="渲染顶部 banner 图片（默认关闭）")
    ap.add_argument("--no-footer", action="store_true",
                    help="不渲染文末「关注」引导卡")
    ap.add_argument("--quality", type=int, default=82, help="图片 JPEG 质量")
    ap.add_argument("--max-img-w", type=int, default=920, help="图片最大宽度")
    ap.add_argument("--math", choices=["png", "svg"], default="png",
                    help="公式渲染方式：png=栅格图（默认，稳妥）；"
                         "svg=矢量内联 SVG（与正文同字号、可缩放不糊，微信可能再次包裹）")
    args = ap.parse_args()

    md_path = os.path.abspath(args.md)
    if not os.path.exists(md_path):
        sys.exit(f"找不到文件：{md_path}")
    base_dir = os.path.dirname(md_path)
    raw = open(md_path, encoding="utf-8").read()
    fm, body = parse_frontmatter(raw)
    title = fm.get("title", "")

    # 0) 提取公式，避免被 markdown 破坏（代码围栏内的 $ 不当公式）
    body_stashed, math_map = stash_math(body)

    # 1) MD → 原始 HTML（footnotes 支持 [^n] 参考文献）
    html = markdown.markdown(body_stashed, extensions=[
        "tables", "fenced_code", "sane_lists", "footnotes"])
    # 2) mermaid → 图片
    html = render_mermaid(html, base_dir, args.diagram, args.quality)
    # 2.5) 修复代码块内 code 被误套浅蓝底
    html = fix_code_blocks(html)
    # 3) 图片内嵌
    html = embed_images(html, base_dir, args.quality)
    # 4) 表格背景双写 + 斑马纹
    html = style_tables(html)
    # 5) 行内样式
    html = apply_inline_styles(html)
    # 5.4) 把 <ul>/<ol>/<li> 展开为 <p> + 项目符号（· / N.），
    #      彻底绕开 ueditor 把 li 内「加粗+裸文本」拆块导致的孤行换行
    html = flatten_lists(html)
    # 5.5) 短粗体标签后的全角冒号并进粗体，保持「标签：」整体加粗
    html = fix_bold_colon(html)
    # 5.6) 在 </strong> 与后续正文间插入零宽空格(U+200B)桥接，
    #      作为 p 内粗体边界的双保险（p 内 strong 本就安全）
    html = fix_bold_break(html)
    # 6) 抑制与标题重复的 h1
    html = suppress_dup_title(html, title)
    # 7) 还原公式（svg=矢量内联 / png=栅格图；失败回退原文）
    html = restore_math(html, math_map, args.math)

    # 资源
    banner_b64 = render_banner(title, args.subtitle) if args.banner else ""
    footer_b64 = render_footer() if not args.no_footer else ""
    doc = build_doc(banner_b64, html, footer_b64, not args.banner, args.no_footer)

    out = args.out or os.path.join(base_dir,
                                   os.path.splitext(os.path.basename(md_path))[0]
                                   + ".wechat.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(doc)
    print(f"written: {out}")
    print(f"size(KB): {round(os.path.getsize(out)/1024, 1)}")


if __name__ == "__main__":
    main()
