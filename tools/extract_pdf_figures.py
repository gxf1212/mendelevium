#!/usr/bin/env python3
"""
通用PDF图片提取工具
支持自动模式和交互模式
"""

import pymupdf
import os
import sys
import argparse

def extract_figures_auto(pdf_path, output_dir, figure_pages=None, min_size=50000):
    """
    自动模式：提取指定页面的最大图片

    Args:
        pdf_path: PDF文件路径
        output_dir: 输出目录
        figure_pages: 字典，如 {1: 3, 2: 5, 3: 6} 表示Figure 1在第3页，Figure 2在第5页
        min_size: 最小图片面积（像素），用于过滤小图标
    """
    os.makedirs(output_dir, exist_ok=True)
    doc = pymupdf.open(pdf_path)
    print(f"PDF总页数: {len(doc)}")

    extracted_count = 0

    for fig_num, page_num in sorted(figure_pages.items()):
        page_idx = page_num - 1  # 转为0-based索引
        print(f"\n提取 Figure {fig_num} (页面 {page_num})...")

        try:
            page = doc[page_idx]
            image_list = page.get_images()

            if not image_list:
                print(f"  警告: 页面 {page_num} 没有图片")
                continue

            # 找最大的图片（面积大于min_size）
            max_area = 0
            best_img = None

            for img in image_list:
                xref = img[0]
                try:
                    base_image = doc.extract_image(xref)
                    width = base_image.get("width", 0)
                    height = base_image.get("height", 0)
                    area = width * height

                    if area > max_area and area >= min_size:
                        max_area = area
                        best_img = base_image
                except:
                    continue

            if best_img:
                output_name = f"fig{fig_num}.{best_img['ext']}"
                output_path = os.path.join(output_dir, output_name)
                with open(output_path, "wb") as f:
                    f.write(best_img["image"])
                print(f"  ✓ 已保存: {output_name} ({best_img['width']}x{best_img['height']}px)")
                extracted_count += 1
            else:
                print(f"  ✗ 未找到符合条件的图片（最小面积: {min_size}px²）")

        except Exception as e:
            print(f"  ✗ 错误: {e}")

    doc.close()
    print(f"\n完成! 成功提取 {extracted_count}/{len(figure_pages)} 张图片")
    print(f"输出目录: {output_dir}")


def extract_figures_interactive(pdf_path, output_dir):
    """
    交互模式：显示所有图片并让用户选择
    """
    os.makedirs(output_dir, exist_ok=True)
    doc = pymupdf.open(pdf_path)
    print(f"PDF总页数: {len(doc)}")

    # 提取所有图片信息
    all_images = []
    for page_num in range(len(doc)):
        page = doc[page_num]
        image_list = page.get_images()

        if image_list:
            for img_index, img in enumerate(image_list):
                xref = img[0]
                try:
                    base_image = doc.extract_image(xref)
                    all_images.append({
                        'page': page_num + 1,
                        'index': img_index,
                        'xref': xref,
                        'width': base_image.get("width", 0),
                        'height': base_image.get("height", 0),
                        'ext': base_image.get("ext", "png"),
                        'bytes': base_image["image"]
                    })
                except:
                    pass

    # 按面积排序
    all_images.sort(key=lambda x: x['width'] * x['height'], reverse=True)

    print(f"\n总共找到 {len(all_images)} 张图片")
    print(f"\n{'序号':<6} {'页码':<6} {'尺寸':<20} {'格式':<8} {'面积':<12}")
    print("-" * 60)

    for i, img in enumerate(all_images[:20]):  # 只显示前20张
        area = img['width'] * img['height']
        print(f"{i+1:<6} {img['page']:<6} {img['width']}x{img['height']:<15} {img['ext']:<8} {area:<12}")

    if len(all_images) > 20:
        print(f"... (还有 {len(all_images)-20} 张图片未显示)")

    # 交互选择
    while True:
        choice = input("\n请输入要保存的图片序号（逗号分隔，如: 1,3,5 或 1-5），'q'退出: ").strip()

        if choice.lower() == 'q':
            break

        selected_indices = []
        try:
            parts = choice.split(',')
            for part in parts:
                part = part.strip()
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    selected_indices.extend(range(start - 1, end))
                else:
                    selected_indices.append(int(part) - 1)
        except ValueError:
            print("输入格式错误，请重试")
            continue

        selected_indices = sorted(set(selected_indices))

        # 保存选中的图片
        for idx in selected_indices:
            if 0 <= idx < len(all_images):
                img = all_images[idx]
                default_name = f"fig{idx+1}.{img['ext']}"
                custom_name = input(f"图片 {idx+1} (页{img['page']}, {img['width']}x{img['height']}): 文件名 (回车={default_name}): ").strip()

                filename = custom_name if custom_name else default_name
                output_path = os.path.join(output_dir, filename)

                with open(output_path, "wb") as f:
                    f.write(img['bytes'])
                print(f"  ✓ 已保存: {filename}")

        cont = input("\n继续选择? (y/n): ").strip().lower()
        if cont != 'y':
            break

    doc.close()
    print("\n完成!")


def main():
    parser = argparse.ArgumentParser(
        description='从PDF中提取图片',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 自动模式：提取指定页面的图片
  python3 tools/extract_pdf_figures.py input.pdf output_dir --pages 3,5,6

  # 指定Figure编号和页码
  python3 tools/extract_pdf_figures.py input.pdf output_dir --figures "1:3,2:5,3:6"

  # 交互模式：手动选择要提取的图片
  python3 tools/extract_pdf_figures.py input.pdf output_dir --interactive

  # 设置最小图片尺寸（过滤小图标）
  python3 tools/extract_pdf_figures.py input.pdf output_dir --pages 3,5,6 --min-size 100000
        """
    )

    parser.add_argument('pdf_path', help='PDF文件路径')
    parser.add_argument('output_dir', help='输出目录')
    parser.add_argument('--pages', help='要提取图片的页码，逗号分隔，如: 3,5,6')
    parser.add_argument('--figures', help='Figure编号和页码映射，如: "1:3,2:5,3:6"')
    parser.add_argument('--interactive', '-i', action='store_true', help='交互模式')
    parser.add_argument('--min-size', type=int, default=50000, help='最小图片面积（像素²），默认50000')

    args = parser.parse_args()

    if not os.path.exists(args.pdf_path):
        print(f"错误: PDF文件不存在: {args.pdf_path}")
        sys.exit(1)

    print(f"PDF文件: {args.pdf_path}")
    print(f"输出目录: {args.output_dir}")
    print("=" * 60 + "\n")

    if args.interactive:
        # 交互模式
        extract_figures_interactive(args.pdf_path, args.output_dir)
    elif args.figures:
        # 使用Figure映射
        figure_pages = {}
        for item in args.figures.split(','):
            fig_num, page_num = map(int, item.split(':'))
            figure_pages[fig_num] = page_num
        extract_figures_auto(args.pdf_path, args.output_dir, figure_pages, args.min_size)
    elif args.pages:
        # 使用页码列表
        pages = [int(p.strip()) for p in args.pages.split(',')]
        figure_pages = {i+1: p for i, p in enumerate(pages)}
        extract_figures_auto(args.pdf_path, args.output_dir, figure_pages, args.min_size)
    else:
        print("错误: 请指定 --pages、--figures 或 --interactive 模式")
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
