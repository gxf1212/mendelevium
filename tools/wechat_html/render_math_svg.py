# -*- coding: utf-8 -*-
"""
render_math_svg.py — 用 MathJax(tex-svg, 经 Node) 把公式渲染为微信友好的内联 SVG。

与 render_math.py（PNG 栅格图）的区别：
- 输出是矢量 <svg>（每个字形为 <path>），不依赖外部字体，不依赖 <defs>；
- 字号用 ex 相对单位，缩放与正文一致、任意放大不糊；
- 微信保存时会丢 currentColor 继承，故在 Node 侧已把颜色写死为深炭灰。

依赖：Node（环境内置 22）+ mathjax-full（装于 node 托管工作区）。
失败（如某公式 MathJax 解析不了）返回 None，由调用方回退 PNG 或原文。
"""
import os
import sys
import json
import shutil
import subprocess

_HERE = os.path.dirname(os.path.abspath(__file__))
NODE = shutil.which("node")
if not NODE:
    NODE = r"C:\Users\Lenovo\.workbuddy\binaries\node\versions\22.22.2\node.exe"
NODE_PATH = r"C:\Users\Lenovo\.workbuddy\binaries\node\workspace\node_modules"
_SCRIPT = os.path.join(_HERE, "render_math_svg.js")


def _available():
    return os.path.exists(NODE) and os.path.exists(_SCRIPT) \
        and os.path.exists(os.path.join(NODE_PATH, "mathjax-full"))


def render_math_svg_batch(items):
    """批量渲染。items: list of (tex, display)。
    返回与 items 等长的 list：svg 字符串 或 None（失败）。"""
    if not _available():
        return [None] * len(items)
    payload = {"items": [{"id": str(i), "tex": t, "display": d}
                        for i, (t, d) in enumerate(items)]}
    try:
        env = dict(os.environ)
        env["NODE_PATH"] = NODE_PATH
        p = subprocess.run([NODE, _SCRIPT],
                           input=json.dumps(payload),
                           capture_output=True, text=True,
                           timeout=180, env=env)
        if p.returncode != 0:
            sys.stderr.write("[render_math_svg] node error: %s\n" % p.stderr[:500])
            return [None] * len(items)
        out = json.loads(p.stdout)
        return [out.get(str(i)) for i in range(len(items))]
    except Exception as e:
        sys.stderr.write("[render_math_svg] %s\n" % e)
        return [None] * len(items)
