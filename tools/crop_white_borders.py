#!/usr/bin/env python3
"""
裁剪图片的白色边距
"""
from PIL import Image
from pathlib import Path

def crop_white_borders(image_path, output_path=None, threshold=250):
    """
    裁剪图片四周的白边

    Args:
        image_path: 输入图片路径
        output_path: 输出图片路径（None则覆盖原文件）
        threshold: 白色阈值（0-255，默认250）
    """
    img = Image.open(image_path)

    # 转为灰度图
    gray = img.convert('L')

    # 将接近白色的像素设为255，其他设为0
    # 然后用getbbox()找到非白色区域的边界框
    bbox = gray.point(lambda x: 0 if x > threshold else 255).getbbox()

    if bbox:
        # 裁剪
        img_cropped = img.crop(bbox)

        # 保存
        save_path = output_path if output_path else image_path
        img_cropped.save(save_path, 'PNG', optimize=True)

        print(f"✓ {Path(image_path).name}: {img.width}x{img.height} -> {img_cropped.width}x{img_cropped.height}")
        return True
    else:
        print(f"✗ {Path(image_path).name}: 无法检测到边界")
        return False

if __name__ == "__main__":
    trio_figs_dir = Path("trio_figs")

    print("开始裁剪白边...")
    print("=" * 60)

    for fig_file in sorted(trio_figs_dir.glob("fig*.png")):
        crop_white_borders(fig_file)

    print("=" * 60)
    print("完成！")
