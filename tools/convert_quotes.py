#!/usr/bin/env python3
"""
将Markdown文件中的英文双引号转换为中文引号。

用法:
    python3 convert_quotes.py <file_path>

功能:
    - 将英文双引号 (U+0022) 转换为中文左右引号对 (U+201C/U+201D)
    - 支持单个文件或多个文件
    - 自动配对：第1、3、5...个引号为左引号，第2、4、6...个为右引号
    - 智能排除特殊情况：
      * Frontmatter（YAML元数据）
      * 代码块（```code``` 和 `inline_code`）
      * Mermaid图表（```mermaid...```）
      * HTML属性（image路径等）

示例:
    # 转换单个文件
    python3 convert_quotes.py _pages/Specific\ Sytems/Enzyme/enzy_control_正文.md

    # 转换多个文件
    python3 convert_quotes.py file1.md file2.md file3.md

    # 使用glob模式（需要在bash中展开）
    python3 convert_quotes.py _pages/**/*.md
"""

import sys
import re
from pathlib import Path


def convert_quotes_in_file(file_path: str) -> tuple[int, int]:
    """
    将文件中的英文双引号转换为中文引号，智能排除特殊情况。

    Args:
        file_path: 要处理的文件路径

    Returns:
        (转换前的英文引号数, 转换后的引号对数)
    """
    file_path = Path(file_path)

    if not file_path.exists():
        print(f"❌ 文件不存在: {file_path}")
        return 0, 0

    if not file_path.is_file():
        print(f"❌ 不是文件: {file_path}")
        return 0, 0

    try:
        # 读取文件
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
    except Exception as e:
        print(f"❌ 读取文件失败: {file_path}")
        print(f"   错误: {e}")
        return 0, 0

    # 计算转换前的英文双引号数
    english_quote_count = content.count(chr(0x0022))

    if english_quote_count == 0:
        print(f"✓ {file_path.name}: 无需转换（已全部是中文引号）")
        return 0, 0

    # 智能转换：排除特殊情况
    new_content = smart_convert_quotes(content)
    quote_pair_count = count_converted_quotes(new_content)

    # 写回文件
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(new_content)
    except Exception as e:
        print(f"❌ 写入文件失败: {file_path}")
        print(f"   错误: {e}")
        return english_quote_count, 0

    quotes_converted = english_quote_count - (new_content.count(chr(0x0022)))
    quote_pair_count = quotes_converted // 2
    print(f"✓ {file_path.name}: {quotes_converted} 个英文引号 → {quote_pair_count} 对中文引号")

    return quotes_converted, quote_pair_count


def smart_convert_quotes(content: str) -> str:
    """智能转换英文引号为中文引号，排除特殊情况"""

    # 定义需要排除的区域模式
    patterns = [
        # Frontmatter: --- ... ---
        (r'^---.*?^---', re.MULTILINE | re.DOTALL),
        # 代码块: ```lang ... ```
        (r'```.*?```', re.MULTILINE | re.DOTALL),
        # 内联代码: `code`
        (r'`[^`]+`', re.MULTILINE),
        # Mermaid图表: ```mermaid ... ```
        (r'```mermaid.*?```', re.MULTILINE | re.DOTALL),
        # HTML标签属性: <tag attr="value"> 或 [text](url"param")
        (r'<[^>]*"[^"]*"[^>]*>', re.MULTILINE),
        (r'\[[^\]]*\]\([^)]*"[^)]*\)', re.MULTILINE),
        # 图片路径: ![alt](path"title")
        (r'!\[[^\]]*\]\([^)]*"[^)]*\)', re.MULTILINE),
    ]

    # 找到所有需要排除的区域
    exclude_regions = []
    for pattern, flags in patterns:
        for match in re.finditer(pattern, content, flags):
            exclude_regions.append((match.start(), match.end()))

    # 合并重叠的区域
    exclude_regions = merge_regions(exclude_regions)

    # 逐个字符处理
    new_content = []
    i = 0
    quote_count = 0
    n = len(content)

    while i < n:
        # 检查当前位置是否在排除区域内
        in_exclude = False
        for start, end in exclude_regions:
            if start <= i < end:
                in_exclude = True
                # 直接添加整个排除区域
                new_content.append(content[i:end])
                i = end
                break

        if in_exclude:
            continue

        if i >= n:
            break

        # 不在排除区域内，检查是否是英文双引号
        if content[i] == chr(0x0022):  # 英文双引号 U+0022
            if quote_count % 2 == 0:
                new_content.append(chr(0x201C))  # 左引号 U+201C "
            else:
                new_content.append(chr(0x201D))  # 右引号 U+201D "
            quote_count += 1
        else:
            new_content.append(content[i])

        i += 1

    return ''.join(new_content)


def merge_regions(regions: list) -> list:
    """合并重叠或相邻的区域"""
    if not regions:
        return []

    # 按起始位置排序
    regions = sorted(regions, key=lambda x: x[0])

    merged = [regions[0]]
    for current in regions[1:]:
        last = merged[-1]
        # 如果重叠或相邻，合并
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)

    return merged


def count_converted_quotes(content: str) -> int:
    """统计转换后的中文引号对数"""
    left_quotes = content.count(chr(0x201C))
    right_quotes = content.count(chr(0x201D))
    return min(left_quotes, right_quotes)


def main():
    if len(sys.argv) < 2:
        print("用法: python3 convert_quotes.py <file_path> [file_path2] ...")
        print("\n示例:")
        print("  python3 convert_quotes.py _pages/article.md")
        print("  python3 convert_quotes.py file1.md file2.md file3.md")
        sys.exit(1)

    file_paths = sys.argv[1:]

    total_quotes = 0
    total_pairs = 0

    print(f"\n开始转换 {len(file_paths)} 个文件...\n")

    for file_path in file_paths:
        quotes, pairs = convert_quotes_in_file(file_path)
        total_quotes += quotes
        total_pairs += pairs

    print(f"\n总计: {total_quotes} 个英文引号 → {total_pairs} 对中文引号")
    print("✓ 转换完成！\n")


if __name__ == "__main__":
    main()
