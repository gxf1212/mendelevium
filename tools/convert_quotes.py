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


def fix_unpaired_quotes_by_position(content: str) -> str:
    """
    修复未配对的中文引号：奇数位置的引号改为左引号，偶数位置的改为右引号。
    只处理有配对问题的行（左引号数≠右引号数）。

    Args:
        content: 文件内容（字符串）

    Returns:
        修复后的文件内容
    """
    LEFT_QUOTE = chr(0x201C)  # 左引号 "
    RIGHT_QUOTE = chr(0x201D)  # 右引号 "

    lines = content.split('\n')
    fixed_lines = []

    for line_num, line in enumerate(lines, 1):
        # 跳过frontmatter（前15行）
        if line_num <= 15 and line.startswith('---'):
            fixed_lines.append(line)
            continue

        # 找出所有中文引号的位置和类型（用ord判断避免字面量问题）
        left_positions = []
        right_positions = []
        for i, char in enumerate(line):
            code = ord(char)
            if code == 0x201C:  # 左引号
                left_positions.append(i)
            elif code == 0x201D:  # 右引号
                right_positions.append(i)

        total = len(left_positions) + len(right_positions)
        if total <= 1:
            # 引号少于2个，无需修复
            fixed_lines.append(line)
            continue

        # 检查是否有配对问题
        if len(left_positions) == len(right_positions):
            # 配对正确，无需修复
            fixed_lines.append(line)
            continue

        # 按出现顺序合并所有引号位置（保持原始顺序）
        all_positions = []
        for pos in left_positions:
            all_positions.append((pos, LEFT_QUOTE))
        for pos in right_positions:
            all_positions.append((pos, RIGHT_QUOTE))
        all_positions.sort(key=lambda x: x[0])

        # 按奇偶位置重置引号：奇数位置为左，偶数位置为右
        fixed_chars = list(line)
        for idx, (pos, _) in enumerate(all_positions):
            if idx % 2 == 0:  # 奇数位置
                fixed_chars[pos] = LEFT_QUOTE
            else:  # 偶数位置
                fixed_chars[pos] = RIGHT_QUOTE
        fixed_lines.append(''.join(fixed_chars))

    return '\n'.join(fixed_lines)


def check_quote_pairing(content: str, content_bytes: bytes) -> list:
    """
    检查中文引号是否配对，从未配对的引号中找出问题位置。

    Args:
        content: 文件内容（字符串）
        content_bytes: 文件内容（字节，用于unicode级别检查）

    Returns:
        问题位置列表 [(line_num, context, issue_type), ...]
    """
    issues = []
    lines = content.split('\n')

    # 检查左引号 U+201C (0xE2 0x80 0x9C in UTF-8)
    # 检查右引号 U+201D (0xE2 0x80 0x9D in UTF-8)
    left_quote_bytes = b'\xe2\x80\x9c'
    right_quote_bytes = b'\xe2\x80\x9d'

    for line_num, line in enumerate(lines, 1):
        # 跳过frontmatter（前10行通常是frontmatter）
        if line_num <= 10:
            continue

        line_bytes = line.encode('utf-8')
        left_count = line_bytes.count(left_quote_bytes)
        right_count = line_bytes.count(right_quote_bytes)

        # 检查是否未配对
        if left_count != right_count:
            # 生成问题上下文
            preview = line.strip()[:100]
            if left_count > right_count:
                issue_type = f"多{left_count - right_count}个左引号"
            else:
                issue_type = f"多{right_count - left_count}个右引号"
            issues.append((line_num, preview, issue_type))

    return issues


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
        # 以二进制模式读取，确保准确识别所有ASCII引号
        with open(file_path, "rb") as f:
            content_bytes = f.read()
    except Exception as e:
        print(f"❌ 读取文件失败: {file_path}")
        print(f"   错误: {e}")
        return 0, 0

    # 计算转换前的ASCII双引号数（0x22）
    english_quote_count = content_bytes.count(b'\x22')

    # DEBUG
    import os
    if os.environ.get('DEBUG_QUOTES'):
        print(f"DEBUG: 文件大小={len(content_bytes)}, ASCII引号数={english_quote_count}")

    if english_quote_count == 0:
        print(f"✓ {file_path.name}: 无需转换（无ASCII引号）")
        return 0, 0

    # 解码为字符串进行处理
    try:
        content = content_bytes.decode('utf-8')
    except Exception as e:
        print(f"❌ 解码文件失败: {file_path}")
        print(f"   错误: {e}")
        return 0, 0

    # 智能转换：排除特殊情况
    new_content = smart_convert_quotes(content)

    # 修复未配对的引号（如果启用FIX_QUOTES环境变量）
    import os
    if os.environ.get('FIX_QUOTES'):
        new_content = fix_unpaired_quotes_by_position(new_content)
        print(f"✓ {file_path.name}: 已修复未配对的引号")

    # 检查中文引号配对（仅在CHECK_QUOTES环境变量启用时）
    import os
    if os.environ.get('CHECK_QUOTES'):
        new_bytes = new_content.encode('utf-8')
        pairing_issues = check_quote_pairing(new_content, new_bytes)
        if pairing_issues:
            print(f"⚠️  发现 {len(pairing_issues)} 处引号配对问题：")
            for line_num, preview, issue in pairing_issues:
                print(f"   第{line_num}行: {preview}")
                print(f"   问题: {issue}")
            # 不阻止转换，只报告问题

    # DEBUG
    import os
    if os.environ.get('DEBUG_QUOTES'):
        new_bytes = new_content.encode('utf-8')
        ascii_q = b'\x22'
        left_q = b'\xe2\x80\x9c'
        right_q = b'\xe2\x80\x9d'
        print(f"DEBUG: 转换后文件大小={len(new_bytes)}")
        print(f"DEBUG: 转换后ASCII引号数={new_bytes.count(ascii_q)}")
        print(f"DEBUG: 转换后左引号数={new_bytes.count(left_q)}")
        print(f"DEBUG: 转换后右引号数={new_bytes.count(right_q)}")

    # 写回文件
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(new_content)
    except Exception as e:
        print(f"❌ 写入文件失败: {file_path}")
        print(f"   错误: {e}")
        return english_quote_count, 0

    # 重新统计ASCII引号数量
    new_content_bytes = new_content.encode('utf-8')
    remaining_quotes = new_content_bytes.count(b'\x22')
    quotes_converted = english_quote_count - remaining_quotes
    quote_pair_count = quotes_converted // 2

    print(f"✓ {file_path.name}: {quotes_converted} 个ASCII引号 → {quote_pair_count} 对中文引号")

    return quotes_converted, quote_pair_count


def smart_convert_quotes(content: str) -> str:
    """智能转换英文引号为中文引号，排除特殊情况"""

    # DEBUG
    import os
    debug = os.environ.get('DEBUG_QUOTES')

    # 定义需要排除的区域模式
    patterns = [
        # 代码块: ```lang ... ``` (先匹配，避免和其他模式冲突)
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

    # 特殊处理：Frontmatter只匹配文件开头的
    if content.startswith('---\n') or content.startswith('---\r\n'):
        # 找到第二个 ---
        second_dash = content.find('\n---\n', 4)
        if second_dash == -1:
            second_dash = content.find('\n---\r\n', 4)
        if second_dash != -1:
            exclude_regions = [(0, second_dash + 5)]  # 包括 \n---\n
            if debug:
                print(f"DEBUG: Frontmatter排除区域 0-{second_dash + 5}")
        else:
            exclude_regions = []
    else:
        exclude_regions = []

    # 找到其他需要排除的区域
    for pattern, flags in patterns:
        for match in re.finditer(pattern, content, flags):
            exclude_regions.append((match.start(), match.end()))
            if debug:
                print(f"DEBUG: 排除区域 {match.start()}-{match.end()}, 长度{match.end()-match.start()}")

    # 合并重叠的区域
    exclude_regions = merge_regions(exclude_regions)

    if debug:
        print(f"DEBUG: 合并后排除区域数={len(exclude_regions)}")
        if exclude_regions:
            print(f"DEBUG: 第一个排除区域={exclude_regions[0]}")

    # 逐个字符处理
    new_content = []
    i = 0
    quote_count = 0
    n = len(content)
    total_converted = 0

    while i < n:
        # 检查当前位置是否在排除区域内
        in_exclude = False
        for start, end in exclude_regions:
            if start <= i < end:
                in_exclude = True
                # 直接添加整个排除区域
                chunk = content[i:end]
                if debug and '"' in chunk:
                    quote_count_in_chunk = chunk.count('"')
                    print(f"DEBUG: Position {i}-{end}: SKIPPING {quote_count_in_chunk} quotes in exclude region")
                new_content.append(chunk)
                i = end
                break

        if in_exclude:
            continue

        if i >= n:
            break

        # 不在排除区域内，检查是否是ASCII双引号
        # 注意：'"' 在Python中就是 chr(0x22)，无论如何输入都是同一个字符
        if content[i] == '"':  # ASCII双引号 (U+0022 或 0x22)
            if quote_count % 2 == 0:
                converted_char = '\u201C'  # 左引号 " (U+201C)
            else:
                converted_char = '\u201D'  # 右引号 " (U+201D)
            new_content.append(converted_char)
            if debug:
                print(f"DEBUG: Position {i}: Quote #{quote_count+1} → {repr(converted_char)} (U+{ord(converted_char):04X})")
            quote_count += 1
            total_converted += 1
        else:
            new_content.append(content[i])

        i += 1

    if debug:
        print(f"DEBUG: smart_convert_quotes转换了 {total_converted} 个引号")
        print(f"DEBUG: new_content列表长度={len(new_content)}")

    result = ''.join(new_content)
    if debug:
        print(f"DEBUG: 返回字符串长度={len(result)}")

    return result


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
        print("\n环境变量:")
        print("  CHECK_QUOTES=1  检查引号配对问题（不修复）")
        print("  FIX_QUOTES=1    修复未配对的引号（按奇偶位置重置）")
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
