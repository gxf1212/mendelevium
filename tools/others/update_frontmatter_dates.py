#!/usr/bin/env python3
"""
更新所有Markdown文件frontmatter中的date字段为文件的最后修改时间。

使用方法:
    python3 tools/update_frontmatter_dates.py --dry-run  # 预览更改
    python3 tools/update_frontmatter_dates.py            # 实际更新
"""

import os
import re
from pathlib import Path
from datetime import datetime
import argparse


def get_file_mtime(file_path: Path) -> str:
    """获取文件的最后修改时间，格式化为YYYY-MM-DD"""
    mtime = os.path.getmtime(file_path)
    return datetime.fromtimestamp(mtime).strftime("%Y-%m-%d")


def extract_frontmatter(content: str) -> tuple[str, str, str]:
    """
    提取frontmatter、正文内容

    返回:
        (frontmatter, body, delimiter) 或 (None, content, None) 如果没有frontmatter
    """
    # 匹配frontmatter: 开头的 --- ... ---
    pattern = r'^---\n(.*?)\n---\n(.*)$'
    match = re.match(pattern, content, re.DOTALL)

    if match:
        return match.group(1), match.group(2), '---'
    else:
        return None, content, None


def update_date_in_frontmatter(frontmatter: str, new_date: str) -> tuple[str, bool]:
    """
    更新frontmatter中的date字段

    返回:
        (updated_frontmatter, changed)
    """
    # 匹配 date: "YYYY-MM-DD" 或 date: YYYY-MM-DD
    pattern = r'date:\s*["\']?(\d{4}-\d{2}-\d{2})["\']?'
    match = re.search(pattern, frontmatter)

    if match:
        old_date = match.group(1)
        if old_date == new_date:
            return frontmatter, False
        # 替换日期
        updated = re.sub(pattern, f'date: "{new_date}"', frontmatter)
        return updated, True
    else:
        # 没有date字段，不做修改
        return frontmatter, False


def should_skip_file(file_path: Path) -> bool:
    """判断是否应该跳过这个文件"""
    skip_patterns = [
        '/archive/',
        '/about.md',
        'index.md',
        '/Diary/'
    ]

    file_str = str(file_path)
    for pattern in skip_patterns:
        if pattern in file_str:
            return True

    return False


def process_file(file_path: Path, dry_run: bool = False) -> dict:
    """
    处理单个Markdown文件

    返回:
        {
            'file': str,
            'status': 'skipped'|'no_frontmatter'|'no_date_field'|'unchanged'|'updated',
            'old_date': str or None,
            'new_date': str or None
        }
    """
    result = {
        'file': str(file_path.relative_to('/mnt/e/GitHub-repo/mendelevium')),
        'status': 'unknown',
        'old_date': None,
        'new_date': None
    }

    # 检查是否应该跳过
    if should_skip_file(file_path):
        result['status'] = 'skipped'
        return result

    # 读取文件
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        result['status'] = f'error_reading: {e}'
        return result

    # 提取frontmatter
    frontmatter, body, delimiter = extract_frontmatter(content)

    if frontmatter is None:
        result['status'] = 'no_frontmatter'
        return result

    # 获取文件修改时间
    new_date = get_file_mtime(file_path)
    result['new_date'] = new_date

    # 更新date字段
    updated_frontmatter, changed = update_date_in_frontmatter(frontmatter, new_date)

    if not changed:
        # 检查是否有date字段
        if 'date:' not in frontmatter:
            result['status'] = 'no_date_field'
        else:
            result['status'] = 'unchanged'
        return result

    # 提取旧日期
    date_match = re.search(r'date:\s*["\']?(\d{4}-\d{2}-\d{2})["\']?', frontmatter)
    if date_match:
        result['old_date'] = date_match.group(1)

    result['status'] = 'updated'

    # 如果不是dry-run，写入文件
    if not dry_run:
        new_content = f"{delimiter}\n{updated_frontmatter}\n{delimiter}\n{body}"
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(new_content)
        except Exception as e:
            result['status'] = f'error_writing: {e}'

    return result


def main():
    parser = argparse.ArgumentParser(
        description='更新Markdown文件frontmatter中的date字段为文件最后修改时间'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='预览更改，不实际修改文件'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='显示详细信息'
    )

    args = parser.parse_args()

    # 查找所有Markdown文件
    pages_dir = Path('/mnt/e/GitHub-repo/mendelevium/_pages')
    md_files = list(pages_dir.rglob('*.md'))

    print(f"找到 {len(md_files)} 个Markdown文件")
    if args.dry_run:
        print("【预览模式】不会实际修改文件\n")
    else:
        print("【实际修改模式】\n")

    # 统计结果
    stats = {
        'skipped': 0,
        'no_frontmatter': 0,
        'no_date_field': 0,
        'unchanged': 0,
        'updated': 0,
        'error': 0
    }

    updated_files = []

    # 处理每个文件
    for md_file in sorted(md_files):
        result = process_file(md_file, dry_run=args.dry_run)

        status = result['status']
        if status.startswith('error'):
            stats['error'] += 1
        elif status in stats:
            stats[status] += 1

        if status == 'updated':
            updated_files.append(result)
            if args.verbose:
                print(f"✓ {result['file']}")
                print(f"  {result['old_date']} → {result['new_date']}")
        elif args.verbose and status not in ['skipped', 'unchanged']:
            print(f"- {result['file']}: {status}")

    # 打印汇总
    print("\n" + "="*60)
    print("统计结果:")
    print(f"  跳过的文件（archive/Diary/index等）: {stats['skipped']}")
    print(f"  无frontmatter: {stats['no_frontmatter']}")
    print(f"  无date字段: {stats['no_date_field']}")
    print(f"  日期未变: {stats['unchanged']}")
    print(f"  已更新: {stats['updated']}")
    print(f"  错误: {stats['error']}")
    print("="*60)

    # 显示更新的文件列表
    if updated_files:
        print(f"\n需要更新的 {len(updated_files)} 个文件:")
        for item in updated_files[:20]:  # 只显示前20个
            print(f"  {item['file']}")
            print(f"    {item['old_date']} → {item['new_date']}")

        if len(updated_files) > 20:
            print(f"  ... 还有 {len(updated_files) - 20} 个文件")

    if args.dry_run and updated_files:
        print("\n提示: 使用不带 --dry-run 参数运行以实际更新文件")


if __name__ == '__main__':
    main()
