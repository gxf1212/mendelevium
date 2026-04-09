#!/usr/bin/env python3
"""
PDF Page Extraction Tool for Claude Code Debugging
用于快速提取PDF单页为临时图片供Claude视觉模型分析

Usage:
    python tools/extract_pdf_page.py <pdf_path> <page_number> [output_path]

Example:
    python tools/extract_pdf_page.py examples/output/membrane-pore-jc.pdf 1
    python tools/extract_pdf_page.py examples/output/membrane-pore-jc.pdf 5 /tmp/page5.png
"""

import sys
import os
from pathlib import Path

def extract_pdf_page(pdf_path: str, page_number: int, output_path: str = None) -> str:
    """
    提取PDF的指定页面为PNG图片

    Args:
        pdf_path: PDF文件路径
        page_number: 页码（从1开始）
        output_path: 输出路径（可选），默认为 /tmp/pdf_page_{page}.png

    Returns:
        生成的图片路径
    """
    try:
        from pdf2image import convert_from_path
    except ImportError:
        print("❌ Error: pdf2image not installed")
        print("📦 Install with: pip install pdf2image")
        print("📦 Also requires poppler-utils:")
        print("   Ubuntu/Debian: sudo apt install poppler-utils")
        print("   macOS: brew install poppler")
        sys.exit(1)

    # 验证PDF文件存在
    if not os.path.exists(pdf_path):
        print(f"❌ Error: PDF file not found: {pdf_path}")
        sys.exit(1)

    # 设置默认输出路径
    if output_path is None:
        output_path = f"/tmp/pdf_page_{page_number}.png"

    try:
        # 转换指定页面 (first_page和last_page从1开始)
        print(f"📄 Extracting page {page_number} from {pdf_path}...")
        images = convert_from_path(
            pdf_path,
            first_page=page_number,
            last_page=page_number,
            dpi=300  # 高质量分辨率（300 DPI标准打印质量）
        )

        if not images:
            print(f"❌ Error: Could not extract page {page_number}")
            sys.exit(1)

        # 保存图片
        images[0].save(output_path, 'PNG')

        # 获取文件大小
        file_size = os.path.getsize(output_path)
        size_kb = file_size / 1024

        print(f"✅ Success! Page {page_number} extracted to: {output_path}")
        print(f"📊 Image size: {images[0].width}x{images[0].height} pixels ({size_kb:.1f} KB)")
        print(f"\n🔍 Use with Claude Code:")
        print(f"   Please view {output_path}")

        return output_path

    except Exception as e:
        print(f"❌ Error during extraction: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    pdf_path = sys.argv[1]

    try:
        page_number = int(sys.argv[2])
        if page_number < 1:
            print("❌ Error: Page number must be >= 1")
            sys.exit(1)
    except ValueError:
        print("❌ Error: Page number must be an integer")
        sys.exit(1)

    output_path = sys.argv[3] if len(sys.argv) > 3 else None

    extract_pdf_page(pdf_path, page_number, output_path)


if __name__ == "__main__":
    main()
