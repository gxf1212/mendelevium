#!/usr/bin/env python3
"""
博客文章图片验证工具
在完成文章撰写后，验证所有Figure/Scheme的编号、文件和图注是否正确
"""

import os
import re
import sys
from pathlib import Path


def extract_figure_references(md_file):
    """从Markdown文件中提取所有图片引用"""
    with open(md_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # 匹配 ![figX](path) 和 ![schemeX](path)
    pattern = r'!\[(fig|scheme)(\d+)\]\(([^)]+)\)'
    matches = re.findall(pattern, content)

    references = []
    for match in matches:
        fig_type = match[0]  # 'fig' or 'scheme'
        fig_num = int(match[1])
        fig_path = match[2]
        references.append({
            'type': fig_type,
            'number': fig_num,
            'path': fig_path,
            'full_name': f'{fig_type}{fig_num}'
        })

    return references


def extract_figure_captions(md_file):
    """从Markdown文件中提取所有图注"""
    with open(md_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    captions = []
    for i, line in enumerate(lines):
        # 匹配 **图X：** 或 **Scheme X：**
        match = re.match(r'\*\*(?:图|Scheme)\s*(\d+)：(.+?)\*\*', line.strip())
        if match:
            fig_num = int(match.group(1))
            caption_text = match.group(2)
            captions.append({
                'number': fig_num,
                'text': caption_text,
                'line': i + 1
            })

    return captions


def verify_blog_figures(md_file, pdf_file=None):
    """
    验证博客文章中的图片

    Args:
        md_file: Markdown文件路径
        pdf_file: 原始PDF文件路径（可选，用于对照验证）
    """
    print(f"🔍 验证博客文章图片: {md_file}\n")
    print("=" * 80)

    # 1. 提取所有图片引用
    references = extract_figure_references(md_file)
    print(f"\n✓ 找到 {len(references)} 个图片引用:")
    for ref in references:
        print(f"  - {ref['full_name']}: {ref['path']}")

    # 2. 检查文件是否存在
    print(f"\n📁 检查图片文件是否存在:")
    md_dir = os.path.dirname(os.path.abspath(md_file))
    missing_files = []

    for ref in references:
        file_path = os.path.join(md_dir, ref['path'])
        if os.path.exists(file_path):
            file_size = os.path.getsize(file_path) / 1024  # KB
            print(f"  ✓ {ref['full_name']}: {file_path} ({file_size:.1f} KB)")
        else:
            print(f"  ✗ {ref['full_name']}: 文件不存在 - {file_path}")
            missing_files.append(ref['full_name'])

    # 3. 检查编号连续性
    print(f"\n🔢 检查编号连续性:")

    # 分别检查fig和scheme
    figs = sorted([r for r in references if r['type'] == 'fig'], key=lambda x: x['number'])
    schemes = sorted([r for r in references if r['type'] == 'scheme'], key=lambda x: x['number'])

    # 检查Figure编号
    if figs:
        fig_numbers = [f['number'] for f in figs]
        expected_figs = list(range(1, max(fig_numbers) + 1))
        missing_figs = set(expected_figs) - set(fig_numbers)

        if missing_figs:
            print(f"  ⚠️  Figure编号不连续，缺少: {sorted(missing_figs)}")
        else:
            print(f"  ✓ Figure编号连续: 1-{max(fig_numbers)}")

    # 检查Scheme编号
    if schemes:
        scheme_numbers = [s['number'] for s in schemes]
        expected_schemes = list(range(1, max(scheme_numbers) + 1))
        missing_schemes = set(expected_schemes) - set(scheme_numbers)

        if missing_schemes:
            print(f"  ⚠️  Scheme编号不连续，缺少: {sorted(missing_schemes)}")
        else:
            print(f"  ✓ Scheme编号连续: 1-{max(scheme_numbers)}")

    # 4. 提取图注
    captions = extract_figure_captions(md_file)
    print(f"\n📝 找到 {len(captions)} 个图注:")
    for cap in captions:
        print(f"  - 图{cap['number']}: {cap['text'][:50]}...")

    # 5. 检查图片引用和图注是否匹配
    print(f"\n🔗 检查图片引用和图注匹配:")
    ref_numbers = set([r['number'] for r in references if r['type'] == 'fig'])
    caption_numbers = set([c['number'] for c in captions])

    only_in_refs = ref_numbers - caption_numbers
    only_in_captions = caption_numbers - ref_numbers

    if only_in_refs:
        print(f"  ⚠️  有图片引用但无图注: {sorted(only_in_refs)}")
    if only_in_captions:
        print(f"  ⚠️  有图注但无图片引用: {sorted(only_in_captions)}")
    if not only_in_refs and not only_in_captions:
        print(f"  ✓ 所有图片都有对应的图注")

    # 6. 总结
    print(f"\n" + "=" * 80)
    print(f"📊 验证总结:")
    print(f"  - 图片引用总数: {len(references)}")
    print(f"    - Figures: {len(figs)}")
    print(f"    - Schemes: {len(schemes)}")
    print(f"  - 图注总数: {len(captions)}")
    print(f"  - 缺失文件: {len(missing_files)}")

    if missing_files:
        print(f"\n❌ 验证失败: {len(missing_files)} 个文件缺失")
        return False
    elif only_in_refs or only_in_captions:
        print(f"\n⚠️  验证通过但有警告")
        return True
    else:
        print(f"\n✅ 验证通过: 所有图片文件和图注都正确")
        return True


def main():
    if len(sys.argv) < 2:
        print("用法: python3 verify_blog_figures.py <markdown_file> [pdf_file]")
        print("\n示例:")
        print("  python3 verify_blog_figures.py article.md")
        print("  python3 verify_blog_figures.py article.md original.pdf")
        sys.exit(1)

    md_file = sys.argv[1]
    pdf_file = sys.argv[2] if len(sys.argv) > 2 else None

    if not os.path.exists(md_file):
        print(f"错误: Markdown文件不存在: {md_file}")
        sys.exit(1)

    success = verify_blog_figures(md_file, pdf_file)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
