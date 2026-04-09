#!/usr/bin/env python3
"""
Blog文章图片验证工具
验证图片引用、图注和文件完整性
"""

import re
import os
import hashlib
from pathlib import Path
import argparse


def calculate_file_hash(filepath):
    """计算文件的MD5哈希值"""
    try:
        with open(filepath, 'rb') as f:
            return hashlib.md5(f.read()).hexdigest()
    except:
        return None


def verify_blog_figures(md_file, verbose=True):
    """
    验证blog文章中的图片引用和图注

    Args:
        md_file: markdown文件路径
        verbose: 是否输出详细信息

    Returns:
        dict: 验证结果
    """
    if not os.path.exists(md_file):
        print(f"❌ 错误: 文件不存在: {md_file}")
        return None

    with open(md_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # 提取图片引用
    figure_refs = re.findall(r'!\[fig(\d+)\]', content)
    scheme_refs = re.findall(r'!\[scheme(\d+)\]', content)

    # 提取图注（支持中文和英文图注）
    figure_captions = re.findall(r'\*\*图(\d+)：', content) + re.findall(r'\*\*Figure (\d+):', content)
    scheme_captions = re.findall(r'\*\*Scheme (\d+)：', content) + re.findall(r'\*\*Scheme (\d+):', content)

    # 提取图片路径
    figure_paths = re.findall(r'!\[fig\d+\]\(([^)]+)\)', content)
    scheme_paths = re.findall(r'!\[scheme\d+\]\(([^)]+)\)', content)

    # 验证编号连续性
    figures = sorted(set(int(f) for f in figure_refs))
    schemes = sorted(set(int(s) for s in scheme_refs))

    expected_figures = list(range(1, len(figures) + 1)) if figures else []
    expected_schemes = list(range(1, len(schemes) + 1)) if schemes else []

    # 验证文件存在性
    base_dir = Path(md_file).parent
    missing_files = []
    existing_files = []
    file_hashes = {}

    for path in figure_paths + scheme_paths:
        full_path = base_dir / path
        if full_path.exists():
            existing_files.append(path)
            file_hash = calculate_file_hash(full_path)
            if file_hash:
                file_hashes[path] = file_hash
        else:
            missing_files.append(path)

    # 检测重复的哈希值
    hash_to_files = {}
    for path, hash_val in file_hashes.items():
        if hash_val not in hash_to_files:
            hash_to_files[hash_val] = []
        hash_to_files[hash_val].append(path)

    duplicate_hashes = {h: files for h, files in hash_to_files.items() if len(files) > 1}

    # 输出结果
    if verbose:
        print(f"\n{'='*60}")
        print(f"验证文件: {md_file}")
        print(f"{'='*60}\n")

        # 图片引用统计
        print(f"📊 图片引用统计:")
        print(f"  Figure引用: {len(figure_refs)}个 (编号: {', '.join(map(str, figures))})")
        print(f"  Scheme引用: {len(scheme_refs)}个 (编号: {', '.join(map(str, schemes))})")

        # 图注统计
        print(f"\n📝 图注统计:")
        print(f"  Figure图注: {len(figure_captions)}个")
        print(f"  Scheme图注: {len(scheme_captions)}个")

        # 编号连续性检查
        print(f"\n🔢 编号连续性检查:")
        if figures == expected_figures:
            print(f"  ✓ Figure编号连续: {', '.join(map(str, figures))}")
        else:
            print(f"  ❌ Figure编号不连续:")
            print(f"     实际: {', '.join(map(str, figures))}")
            print(f"     期望: {', '.join(map(str, expected_figures))}")

        if schemes == expected_schemes:
            print(f"  ✓ Scheme编号连续: {', '.join(map(str, schemes))}")
        elif schemes:
            print(f"  ❌ Scheme编号不连续:")
            print(f"     实际: {', '.join(map(str, schemes))}")
            print(f"     期望: {', '.join(map(str, expected_schemes))}")

        # 图注完整性检查
        print(f"\n📋 图注完整性检查:")
        if len(figure_refs) == len(figure_captions):
            print(f"  ✓ Figure引用和图注数量一致 ({len(figure_refs)}个)")
        else:
            print(f"  ⚠️  Figure引用({len(figure_refs)})和图注({len(figure_captions)})数量不匹配")

        if len(scheme_refs) == len(scheme_captions):
            print(f"  ✓ Scheme引用和图注数量一致 ({len(scheme_refs)}个)")
        elif scheme_refs:
            print(f"  ⚠️  Scheme引用({len(scheme_refs)})和图注({len(scheme_captions)})数量不匹配")

        # 文件存在性检查
        print(f"\n📁 文件存在性检查:")
        if existing_files:
            print(f"  ✓ 存在的文件: {len(existing_files)}个")
        if missing_files:
            print(f"  ❌ 缺失的文件: {len(missing_files)}个")
            for path in missing_files:
                print(f"     - {path}")

        # 重复哈希检查
        if duplicate_hashes:
            print(f"\n⚠️  检测到重复的图片（MD5相同）:")
            for hash_val, files in duplicate_hashes.items():
                print(f"  哈希值: {hash_val}")
                for path in files:
                    file_path = base_dir / path
                    size = os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    print(f"    - {path} ({size:,} bytes)")
            print(f"\n  💡 请手动确认这些图片是否为重复或同一页面的子图")
        else:
            print(f"\n✓ 未检测到重复的图片")

        print(f"\n{'='*60}")
        print(f"验证完成!")
        print(f"{'='*60}\n")

    return {
        'figures': figures,
        'schemes': schemes,
        'figure_refs': figure_refs,
        'scheme_refs': scheme_refs,
        'figure_captions': figure_captions,
        'scheme_captions': scheme_captions,
        'figure_paths': figure_paths,
        'scheme_paths': scheme_paths,
        'existing_files': existing_files,
        'missing_files': missing_files,
        'duplicate_hashes': duplicate_hashes,
        'file_hashes': file_hashes
    }


def main():
    parser = argparse.ArgumentParser(
        description='验证blog文章中的图片引用和图注',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 验证主文档
  python3 tools/verify_blog_figures.py article.md

  # 验证附录
  python3 tools/verify_blog_figures.py article-appendix.md

  # 简洁模式（仅输出错误）
  python3 tools/verify_blog_figures.py article.md --quiet
        """
    )

    parser.add_argument('md_file', help='markdown文件路径')
    parser.add_argument('--quiet', '-q', action='store_true', help='简洁模式，仅输出错误和警告')

    args = parser.parse_args()

    result = verify_blog_figures(args.md_file, verbose=not args.quiet)

    if result is None:
        exit(1)

    # 如果有错误，返回非零退出码
    has_errors = (
        result['missing_files'] or
        (result['figures'] != list(range(1, len(result['figures']) + 1)) if result['figures'] else False) or
        len(result['figure_refs']) != len(result['figure_captions'])
    )

    if has_errors:
        exit(1)
    else:
        exit(0)


if __name__ == "__main__":
    main()
