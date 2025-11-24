#!/usr/bin/env python3
"""
随机选择缩略图脚本
用于为新的博客文章随机选择缩略图
"""

import os
import random
import sys
from pathlib import Path

def get_random_thumbnail():
    """
    从assets/img/thumbnail和assets/img/thumbnail_mine目录中随机选择一个缩略图

    Returns:
        str: 随机选择的缩略图路径（相对于assets目录）
    """
    # 定义缩略图目录
    thumbnail_dirs = [
        "assets/img/thumbnail",
        "assets/img/thumbnail_mine"
    ]

    # 支持的图片扩展名
    image_extensions = {'.jpg', '.jpeg', '.png', '.webp'}

    all_images = []

    # 遍历所有目录，收集图片文件
    for dir_path in thumbnail_dirs:
        full_path = Path(dir_path)
        if full_path.exists():
            for file in full_path.iterdir():
                if file.is_file() and file.suffix.lower() in image_extensions:
                    # 返回相对于assets目录的路径
                    relative_path = str(file).replace("assets/", "")
                    all_images.append(relative_path)

    if not all_images:
        # 如果没有找到图片，返回默认的empty.jpg
        return "img/thumbnail/empty.jpg"

    # 随机选择一个图片
    return random.choice(all_images)

def print_random_thumbnail():
    """打印随机选择的缩略图路径"""
    thumbnail = get_random_thumbnail()
    print(thumbnail)

def generate_frontmatter_thumbnail():
    """生成用于frontmatter的缩略图配置"""
    thumbnail = get_random_thumbnail()
    full_path = f"/assets/{thumbnail}"

    print(f'  thumbnail: "{full_path}"')
    print(f'  image: "{full_path}"')

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "--frontmatter":
            generate_frontmatter_thumbnail()
        elif sys.argv[1] == "--help":
            print("用法:")
            print("  python random_thumbnail.py              # 随机选择缩略图路径")
            print("  python random_thumbnail.py --frontmatter  # 生成frontmatter格式")
            print("  python random_thumbnail.py --help       # 显示帮助")
        else:
            print("未知参数。使用 --help 查看帮助。")
    else:
        print_random_thumbnail()