#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""批量生成 wechat_html 浅色主题 JSON。

结构以 default.json 为准，这里只定义每个主题的「色板 palette」，其余
（字号/边距/布局/圆角等非颜色字段）由 build_theme() 固定填充。
想加新主题：往 PALETTES 里加一个 dict，重跑本脚本即可。

用法：
  python gen_themes.py          # 在 themes/ 目录生成/覆盖所有主题 JSON
"""
import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.join(_HERE, "themes")

# 固定非颜色结构（与 default.json 一致）
FONT_STACK = "-apple-system,'Segoe UI','Microsoft YaHei',sans-serif"
MONO = "Consolas,'Courier New',monospace"


def build_theme(p):
    """p: 色板 dict。返回完整主题 dict。"""
    return {
        "meta": {"name": p["name"], "desc": p["desc"]},
        "styles": {
            "heading": {"default": "plain", "choices": ["plain", "pill", "number", "bar"]},
            "list":    {"default": "plain", "choices": ["plain", "circle", "square"]},
            "quote":   {"default": "bar",   "choices": ["bar", "quote"]},
        },
        "font": {"stack": FONT_STACK, "mono": MONO,
                 "base_size": 16, "line_height": 1.75},
        "colors": {
            "blue": p["blue"], "blue_d": p["blue_d"], "coral": p["coral"],
            "blue_bg": p["blue_bg"], "ink": p["ink"], "ink_d": p["ink_d"],
            "ink_h3": p["ink_h3"], "ink_h4": p["ink_h4"], "ink_h5": p["ink_h5"],
            "ink_h6": p["ink_h6"], "gray": p["gray"], "border": p["border"],
            "quote_bg": p["quote_bg"], "code_bg": p["code_bg"], "code_fg": p["code_fg"],
            "code_inline_bg": p["code_inline_bg"], "code_inline_fg": p["code_inline_fg"],
            "link": p["blue"], "link_border": p["link_border"], "hr": p["border"],
            "img_border": p["img_border"], "strong": p["ink_d"], "em": p["ink_h4"],
            "sup": p["blue"], "bullet": p["blue_d"], "footnote": p["gray"],
            "footer_top": p["footer_top"], "footer_bottom": p["footer_bottom"],
            "footer_text": "#ffffff", "footer_sub": p["footer_sub"],
        },
        "headings": {
            "h1": {"size": 22, "color": p["ink_d"], "weight": "bold",
                   "border_bottom": f'3px solid {p["blue"]}', "padding_bottom": "10px",
                   "margin": "18px 0 14px"},
            "h2": {"size": 19, "color": p["ink_d"], "weight": "bold",
                   "border_left": f'4px solid {p["blue"]}', "padding_left": "10px",
                   "margin": "26px 0 12px"},
            "h3": {"size": 17, "color": p["ink_h3"], "weight": "bold",
                   "margin": "22px 0 10px"},
            "h4": {"size": 16, "color": p["ink_h4"], "weight": "bold",
                   "margin": "16px 0 8px"},
            "h5": {"size": 16, "color": p["ink_h5"], "weight": "bold",
                   "margin": "14px 0 6px"},
            "h6": {"size": 16, "color": p["ink_h6"], "weight": "bold",
                   "margin": "12px 0 6px"},
        },
        "body": {"color": p["ink"], "size": 16, "line_height": 1.75,
                 "margin": "12px 0"},
        "blockquote": {
            "bg": p["quote_bg"], "border_left": f'4px solid {p["blue"]}',
            "color": p["ink_h4"], "padding": "12px 16px", "radius": "0 8px 8px 0",
            "margin": "14px 0", "size": 16, "line_height": 1.75},
        "list": {"margin": "6px 0", "indent_per_level": 4, "bullet_symbol": "•",
                 "bullet_size": 20, "ordered_size": 18, "bullet_color": p["blue_d"],
                 "bullet_weight": "bold"},
        "table": {"header_bg": p["blue"], "header_fg": "#ffffff",
                  "border": p["border"], "zebra_bg": p["blue_bg"],
                  "cell_padding": "8px 8px", "header_padding": "9px 8px",
                  "cell_size": 16, "cell_line_height": 1.6, "margin": "14px 0"},
        "code": {"block_bg": p["code_bg"], "block_fg": p["code_fg"],
                 "block_size": 13, "block_line_height": 1.6,
                 "block_padding": "12px 14px", "block_radius": "8px",
                 "block_margin": "12px 0", "inline_bg": p["code_inline_bg"],
                 "inline_fg": p["code_inline_fg"], "inline_size": 13,
                 "inline_padding": "1px 5px", "inline_radius": "4px"},
        "image": {"border": f'1px solid {p["img_border"]}', "radius": "6px",
                  "max_width": "96%", "caption_size": 13, "caption_color": p["gray"]},
        "layout": {"max_width": 680, "outer_padding": "0 8px"},
    }


# ===== 浅色主题色板 =====
PALETTES = {
    "shanchui": {
        "name": "山吹", "desc": "温暖金黄 · 山吹色，温馨内容、读书笔记",
        "blue": "#D99A2B", "blue_d": "#B87E1E", "coral": "#E07A5F",
        "blue_bg": "#FBF3E4", "ink": "#4A4234", "ink_d": "#3A3226",
        "ink_h3": "#B87E1E", "ink_h4": "#8C6A2F", "ink_h5": "#A88E5F",
        "ink_h6": "#BFAE8A", "gray": "#8C8070", "border": "#EFE3C8",
        "quote_bg": "#FBF6EA", "code_bg": "#3A3226", "code_fg": "#F5EDDB",
        "code_inline_bg": "#FBF3E4", "code_inline_fg": "#B87E1E",
        "link_border": "#EAD9B0", "img_border": "#EFE8DA",
        "footer_top": "#3A3226", "footer_bottom": "#D99A2B", "footer_sub": "#F3E7C8",
    },
    "rose": {
        "name": "蔷薇紫", "desc": "优雅紫 · 蔷薇色调，优质文章、深度分享",
        "blue": "#9B6DB6", "blue_d": "#7E549A", "coral": "#E58AB0",
        "blue_bg": "#F5EFFA", "ink": "#4A4156", "ink_d": "#3A3050",
        "ink_h3": "#7E549A", "ink_h4": "#9A7FB0", "ink_h5": "#B49CC8",
        "ink_h6": "#C6B3D6", "gray": "#9A8FA8", "border": "#E8DDF0",
        "quote_bg": "#F8F2FB", "code_bg": "#3A3050", "code_fg": "#F3EDF8",
        "code_inline_bg": "#F5EFFA", "code_inline_fg": "#7E549A",
        "link_border": "#E0D3EC", "img_border": "#ECE4F2",
        "footer_top": "#3A3050", "footer_bottom": "#9B6DB6", "footer_sub": "#EADEF2",
    },
    "fullstack-blue": {
        "name": "全栈蓝", "desc": "专业蓝 · 全栈技术风，技术文章、教程",
        "blue": "#2F6FED", "blue_d": "#1E55C8", "coral": "#E07A5F",
        "blue_bg": "#EBF2FE", "ink": "#3F4858", "ink_d": "#2C3A55",
        "ink_h3": "#1E55C8", "ink_h4": "#4A6FA8", "ink_h5": "#7A92BC",
        "ink_h6": "#A0B2D0", "gray": "#8A93A3", "border": "#D8E4F8",
        "quote_bg": "#F2F7FE", "code_bg": "#2C3A55", "code_fg": "#EDF2FB",
        "code_inline_bg": "#EBF2FE", "code_inline_fg": "#1E55C8",
        "link_border": "#C9D8F5", "img_border": "#E3EAF5",
        "footer_top": "#2C3A55", "footer_bottom": "#2F6FED", "footer_sub": "#DCE8FB",
    },
    "cute-green": {
        "name": "萌绿", "desc": "清新绿 · 自然系，轻松阅读、科普",
        "blue": "#45A97A", "blue_d": "#2F8A5F", "coral": "#E8935B",
        "blue_bg": "#EAF6F0", "ink": "#40504A", "ink_d": "#2E4038",
        "ink_h3": "#2F8A5F", "ink_h4": "#5A9A7E", "ink_h5": "#85B29E",
        "ink_h6": "#A8C6B8", "gray": "#8A9892", "border": "#D8EAE1",
        "quote_bg": "#F0F8F4", "code_bg": "#2E4038", "code_fg": "#EDF5F1",
        "code_inline_bg": "#EAF6F0", "code_inline_fg": "#2F8A5F",
        "link_border": "#C9E3D6", "img_border": "#E3EEE8",
        "footer_top": "#2E4038", "footer_bottom": "#45A97A", "footer_sub": "#DCEFE5",
    },
    "frontend": {
        "name": "前端之巅", "desc": "青蓝科技 · 前端之巅同款气质，技术分享",
        "blue": "#1890C0", "blue_d": "#0E7490", "coral": "#E07A5F",
        "blue_bg": "#E8F5FA", "ink": "#3F4A52", "ink_d": "#2A383F",
        "ink_h3": "#0E7490", "ink_h4": "#4A8899", "ink_h5": "#7AA6B2",
        "ink_h6": "#A0C0C8", "gray": "#8A969B", "border": "#D6E8EE",
        "quote_bg": "#F0F8FA", "code_bg": "#2A383F", "code_fg": "#EDF5F7",
        "code_inline_bg": "#E8F5FA", "code_inline_fg": "#0E7490",
        "link_border": "#C5E0E8", "img_border": "#E2EDF0",
        "footer_top": "#2A383F", "footer_bottom": "#1890C0", "footer_sub": "#D9EEF3",
    },
    "orange-blue": {
        "name": "橙蓝风", "desc": "橙蓝撞色 · 活力明快，清单、盘点、评测",
        "blue": "#E8823A", "blue_d": "#C96A24", "coral": "#2F6FED",
        "blue_bg": "#FCF0E6", "ink": "#4A423C", "ink_d": "#3A3228",
        "ink_h3": "#C96A24", "ink_h4": "#B0804E", "ink_h5": "#C89E7C",
        "ink_h6": "#D8BCA2", "gray": "#9A8F86", "border": "#F0DFCE",
        "quote_bg": "#FCF5EC", "code_bg": "#3A3228", "code_fg": "#F5EDE4",
        "code_inline_bg": "#FCF0E6", "code_inline_fg": "#C96A24",
        "link_border": "#F0D9C0", "img_border": "#F0E6DA",
        "footer_top": "#3A3228", "footer_bottom": "#E8823A", "footer_sub": "#F5E0CA",
    },
    "cute-pink": {
        "name": "萌粉", "desc": "温柔粉 · 樱花系，生活分享、随笔",
        "blue": "#D977A8", "blue_d": "#B85A8A", "coral": "#9B6DB6",
        "blue_bg": "#FCEFF5", "ink": "#4A4148", "ink_d": "#3A3038",
        "ink_h3": "#B85A8A", "ink_h4": "#C0809E", "ink_h5": "#D0A2B8",
        "ink_h6": "#DDBEC8", "gray": "#9A8F94", "border": "#F0DDE6",
        "quote_bg": "#FDF3F8", "code_bg": "#3A3038", "code_fg": "#F5EDF1",
        "code_inline_bg": "#FCEFF5", "code_inline_fg": "#B85A8A",
        "link_border": "#F0D2E0", "img_border": "#F2E3EA",
        "footer_top": "#3A3038", "footer_bottom": "#D977A8", "footer_sub": "#F5DEE8",
    },
    "black-white": {
        "name": "极简", "desc": "黑白灰 · 极简主义，白底黑字，克制留白",
        "blue": "#4A4A4A", "blue_d": "#2E2E2E", "coral": "#E07A5F",
        "blue_bg": "#F5F5F5", "ink": "#4A4A4A", "ink_d": "#2E2E2E",
        "ink_h3": "#2E2E2E", "ink_h4": "#5A5A5A", "ink_h5": "#7A7A7A",
        "ink_h6": "#9A9A9A", "gray": "#8A8A8A", "border": "#E0E0E0",
        "quote_bg": "#F7F7F7", "code_bg": "#2E2E2E", "code_fg": "#F0F0F0",
        "code_inline_bg": "#F5F5F5", "code_inline_fg": "#2E2E2E",
        "link_border": "#D0D0D0", "img_border": "#E8E8E8",
        "footer_top": "#2E2E2E", "footer_bottom": "#4A4A4A", "footer_sub": "#DEDEDE",
    },
}


def main():
    os.makedirs(_OUT, exist_ok=True)
    for key, palette in PALETTES.items():
        path = os.path.join(_OUT, f"{key}.json")
        with open(path, "w", encoding="utf-8") as f:
            json.dump(build_theme(palette), f, ensure_ascii=False, indent=2)
            f.write("\n")
        print(f"written: {path}")
    print(f"共 {len(PALETTES)} 个主题")


if __name__ == "__main__":
    main()
