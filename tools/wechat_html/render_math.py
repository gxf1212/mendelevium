# -*- coding: utf-8 -*-
"""
render_math.py — 用本机 LaTeX 工具链把公式渲染为透明背景 PNG（base64）。

公众号不渲染 LaTeX，故在生成 HTML 前把 $...$ / $$...$$ 预先转成图片。
优先用本机 latex + dvipng，质量最高，支持化学式、分式、求和、希腊字母等；
若工具链不可用则返回 None，由调用方 fallback。

依赖：系统 PATH 中需有 latex 与 dvipng（MiKTeX / TeX Live 均可）。
"""
import os
import re
import sys
import shutil
import base64
import subprocess
import tempfile
from io import BytesIO

from PIL import Image

# 正文基准字号（px），用于将公式图高度换算成 em，实现与文字基线对齐
BASE_FONT_PX = 16
DPI = 200

_PREAMBLE = (
    r"\documentclass[preview,border=1pt]{standalone}"
    r"\usepackage{amsmath,amssymb,amstext,mhchem}"
    r"\begin{document}"
    # 用与公众号正文一致的 16pt 字号渲染，使公式视觉大小≈正文（16px = 1em）
    r"\fontsize{16}{19}\selectfont"
)
_POSTAMBLE = r"\end{document}"


def _engine_ok():
    return bool(shutil.which("latex")) and bool(shutil.which("dvipng"))


def render_math(tex, display=False):
    """渲染一条 LaTeX 公式为透明背景 PNG。

    返回 (b64_png, width_px, height_px)；失败返回 None。
    display=True 渲染为块级公式（用 \\[ ... \\]），否则行内（$ ... $）。
    """
    if not _engine_ok():
        return None
    env = (r"\[" + tex + r"\]") if display else ("$" + tex + "$")
    doc = _PREAMBLE + env + _POSTAMBLE
    d = tempfile.mkdtemp()
    texpath = os.path.join(d, "f.tex")
    try:
        with open(texpath, "w", encoding="utf-8") as f:
            f.write(doc)
        r = subprocess.run(
            ["latex", "-interaction=nonstopmode", "-halt-on-error",
             "-output-directory=" + d, texpath],
            capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            return None
        dvi = os.path.join(d, "f.dvi")
        if not os.path.exists(dvi):
            return None
        png = os.path.join(d, "f.png")
        r2 = subprocess.run(
            ["dvipng", "-D", str(DPI), "-bg", "Transparent", "-o", png, dvi],
            capture_output=True, text=True, timeout=60)
        if r2.returncode != 0 or not os.path.exists(png):
            return None
        im = Image.open(png).convert("RGBA")
        buf = BytesIO()
        im.save(buf, format="PNG")
        return (base64.b64encode(buf.getvalue()).decode(),
                im.width, im.height)
    except Exception:
        return None


def inline_img_style(h_px):
    """公式已用 16pt 渲染（与正文 16px 同字号），按真实高度 1:1 显示，
    不再缩小，使其视觉大小与正文一致（≈1em）。含上下标的高公式略限高。"""
    em = round(h_px * 72.0 / DPI / BASE_FONT_PX, 2)
    em = min(em, 1.25)
    return f"vertical-align:middle;height:{em}em;"
