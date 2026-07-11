#!/usr/bin/env python3
"""
随机选择缩略图脚本
用于为新的博客文章随机选择缩略图
"""

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
        "assets/img/thumbnail_mine",
        "assets/img/Wallpaper_compressed",
        "assets/img/4K_1080P_compressed"
    ]

    # 支持的图片扩展名
    image_extensions = {'.jpg', '.jpeg', '.png', '.webp'}

    all_images = []

    # 遍历所有目录，收集图片文件
    excluded = {"empty.jpg"}
    for dir_path in thumbnail_dirs:
        full_path = Path(dir_path)
        if full_path.exists():
            for file in full_path.iterdir():
                if file.is_file() and file.suffix.lower() in image_extensions:
                    if file.name in excluded:
                        continue
                    # 返回相对于assets目录的路径
                    relative_path = str(file).replace("assets/", "")
                    all_images.append(relative_path)

    if not all_images:
        # 如果没有找到图片，报错并退出
        print("错误：没有找到可用的缩略图文件", file=sys.stderr)
        sys.exit(1)

    # 随机选择一个图片
    return random.choice(all_images)

def print_random_thumbnail():
    """打印随机选择的缩略图路径"""
    thumbnail = get_random_thumbnail()
    # 输出完整路径，与frontmatter格式一致
    full_path = f"/assets/{thumbnail}"
    print(full_path)

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