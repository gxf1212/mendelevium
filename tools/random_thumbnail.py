#!/usr/bin/env python3
"""
随机选择缩略图脚本
用于为新的博客文章随机选择缩略图，输出 GitHub raw CDN URL
"""

import random
import sys
from pathlib import Path

# GitHub raw CDN 基础 URL（仓库公开，可直接访问）
RAW_BASE = "https://raw.githubusercontent.com/gxf1212/mendelevium/main"

# 缩略图目录（相对于仓库根目录）
THUMBNAIL_DIRS = [
    "assets/img/thumbnail",
    "assets/img/thumbnail_mine",
    "assets/img/Wallpaper_compressed",
    "assets/img/4K_1080P_compressed"
]

IMAGE_EXTS = {'.jpg', '.jpeg', '.png', '.webp'}
EXCLUDED = {"empty.jpg"}


def get_random_thumbnail():
    """
    从所有缩略图目录中随机选择一个图片

    Returns:
        str: GitHub raw CDN URL（可直接用于 frontmatter 的 image/thumbnail 字段）
    """
    all_images = []
    for dir_path in THUMBNAIL_DIRS:
        p = Path(dir_path)
        if p.exists():
            for f in p.iterdir():
                if f.is_file() and f.suffix.lower() in IMAGE_EXTS and f.name not in EXCLUDED:
                    all_images.append(str(f))  # 保持完整路径如 assets/img/xxx/yyy.jpg

    if not all_images:
        print("错误：没有找到可用的缩略图文件", file=sys.stderr)
        sys.exit(1)

    chosen = random.choice(all_images)
    return f"{RAW_BASE}/{chosen}"


def print_random_thumbnail():
    """打印随机选择的缩略图（GitHub raw CDN URL）"""
    print(get_random_thumbnail())


def generate_frontmatter_thumbnail():
    """生成用于 frontmatter 的缩略图配置"""
    url = get_random_thumbnail()
    print(f'  thumbnail: "{url}"')
    print(f'  image: "{url}"')


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "--frontmatter":
            generate_frontmatter_thumbnail()
        elif sys.argv[1] == "--help":
            print("用法:")
            print("  python random_thumbnail.py              # 随机缩略图URL")
            print("  python random_thumbnail.py --frontmatter  # frontmatter格式")
            print("  python random_thumbnail.py --help       # 帮助")
        else:
            print("未知参数。使用 --help 查看帮助。")
    else:
        print_random_thumbnail()
