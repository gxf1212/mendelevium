#!/usr/bin/env python3
"""
验证图片中是否包含Scheme标题文字
使用OCR检测图片中的文字
"""

import pymupdf
from pdf2image import convert_from_path
import pytesseract
from PIL import Image
import io
import sys

def check_image_contains_scheme(pdf_path, page_num, img_index=0):
    """
    检查指定页面的指定图片是否包含'Scheme'文字

    Args:
        pdf_path: PDF文件路径
        page_num: 页码（1-based）
        img_index: 图片索引（0-based）

    Returns:
        bool: 是否包含Scheme文字
    """
    doc = pymupdf.open(pdf_path)
    page = doc[page_num - 1]
    image_list = page.get_images()

    if img_index >= len(image_list):
        print(f"Error: Page {page_num} only has {len(image_list)} images")
        doc.close()
        return False

    # 提取图片
    xref = image_list[img_index][0]
    base_image = doc.extract_image(xref)
    image_bytes = base_image["image"]

    # 转换为PIL Image
    image = Image.open(io.BytesIO(image_bytes))

    # 使用OCR提取文字
    try:
        text = pytesseract.image_to_string(image)
        contains_scheme = 'Scheme' in text or 'scheme' in text

        if contains_scheme:
            print(f"✓ Page {page_num}, Image {img_index+1}: Contains 'Scheme'")
            # 打印包含Scheme的行
            for line in text.split('\n'):
                if 'Scheme' in line or 'scheme' in line:
                    print(f"  Found: {line.strip()}")
        else:
            print(f"✗ Page {page_num}, Image {img_index+1}: No 'Scheme' found")

        doc.close()
        return contains_scheme
    except Exception as e:
        print(f"Error during OCR: {e}")
        doc.close()
        return False


def scan_all_images_for_schemes(pdf_path):
    """扫描PDF中所有图片，查找包含Scheme的图片"""
    doc = pymupdf.open(pdf_path)

    print(f"Scanning {pdf_path} for images containing 'Scheme'...\n")

    scheme_images = []

    for page_num in range(len(doc)):
        page = doc[page_num]
        image_list = page.get_images()

        if not image_list:
            continue

        for img_index, img in enumerate(image_list):
            xref = img[0]
            try:
                base_image = doc.extract_image(xref)
                image_bytes = base_image["image"]
                width = base_image.get("width", 0)
                height = base_image.get("height", 0)

                # 转换为PIL Image
                image = Image.open(io.BytesIO(image_bytes))

                # OCR
                text = pytesseract.image_to_string(image)

                if 'Scheme' in text:
                    # 提取Scheme编号
                    scheme_num = None
                    for line in text.split('\n'):
                        if 'Scheme' in line:
                            if 'Scheme 1' in line:
                                scheme_num = 1
                            elif 'Scheme 2' in line:
                                scheme_num = 2
                            elif 'Scheme 3' in line:
                                scheme_num = 3
                            elif 'Scheme 4' in line:
                                scheme_num = 4

                            if scheme_num:
                                break

                    scheme_images.append({
                        'page': page_num + 1,
                        'img_index': img_index,
                        'scheme_num': scheme_num,
                        'width': width,
                        'height': height,
                        'area': width * height
                    })

                    print(f"✓ Page {page_num+1}, Image {img_index+1}: Scheme {scheme_num} ({width}x{height})")

            except Exception as e:
                continue

    doc.close()

    print(f"\n=== Summary ===")
    print(f"Found {len(scheme_images)} images containing 'Scheme':")
    for info in sorted(scheme_images, key=lambda x: x.get('scheme_num', 999)):
        print(f"  Scheme {info['scheme_num']}: Page {info['page']}, Image {info['img_index']+1} ({info['width']}x{info['height']})")

    return scheme_images


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python3 verify_scheme_in_image.py <pdf_path>")
        print("  python3 verify_scheme_in_image.py <pdf_path> <page_num> <img_index>")
        sys.exit(1)

    pdf_path = sys.argv[1]

    if len(sys.argv) == 4:
        # 检查单个图片
        page_num = int(sys.argv[2])
        img_index = int(sys.argv[3])
        check_image_contains_scheme(pdf_path, page_num, img_index)
    else:
        # 扫描所有图片
        scan_all_images_for_schemes(pdf_path)
