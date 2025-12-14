#!/usr/bin/env python3
"""
全面修复Markdown文件中的中文标点符号。

转换规则：
- 英文逗号 , → 中文逗号 ，
- 英文句号 . → 中文句号 。（仅在中文句子末尾）
- 英文分号 ; → 中文分号 ；
- 英文感叹号 ! → 中文感叹号 ！
- 英文问号 ? → 中文问号 ？

排除区域：
- Frontmatter
- 代码块
- 公式
- URL
- 数字（如 3.14, 1.5）
- 英文句子
"""

import sys
import re
from pathlib import Path


def fix_all_punctuation(file_path: str) -> dict:
    """修复文件中的所有标点符号"""
    file_path = Path(file_path)

    if not file_path.exists():
        print(f"❌ 文件不存在： {file_path}")
        return {}

    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    new_content = smart_convert_punctuation(content)

    # 统计转换数量
    stats = {
        '逗号': content.count(',') - new_content.count(','),
        '句号': content.count('.') - new_content.count('.'),
        '分号': content.count(';') - new_content.count(';'),
        '感叹号': content.count('!') - new_content.count('!'),
        '问号': content.count('?') - new_content.count('?'),
    }

    with open(file_path, "w", encoding="utf-8") as f:
        f.write(new_content)

    print(f"✓ {file_path.name}：")
    for punct, count in stats.items():
        if count > 0:
            print(f"  - {punct}： {count} 个")

    return stats


def smart_convert_punctuation(content: str) -> str:
    """智能转换标点符号"""

    # 找到排除区域
    exclude_regions = []

    # 1. Frontmatter
    if content.startswith('---\n') or content.startswith('---\r\n'):
        second_dash = content.find('\n---\n', 4)
        if second_dash == -1:
            second_dash = content.find('\n---\r\n', 4)
        if second_dash != -1:
            exclude_regions.append((0, second_dash + 5))

    # 2. 代码块
    for match in re.finditer(r'```.*?```', content, re.DOTALL):
        exclude_regions.append((match.start(), match.end()))

    # 3. 内联代码
    for match in re.finditer(r'`[^`]+`', content):
        exclude_regions.append((match.start(), match.end()))

    # 4. URL
    for match in re.finditer(r'https?://[^\s\)）]+', content):
        exclude_regions.append((match.start(), match.end()))

    # 5. 公式
    for match in re.finditer(r'\$[^\$]+\$', content):
        exclude_regions.append((match.start(), match.end()))
    for match in re.finditer(r'\$\$.*?\$\$', content, re.DOTALL):
        exclude_regions.append((match.start(), match.end()))

    # 6. HTML标签
    for match in re.finditer(r'<[^>]+>', content):
        exclude_regions.append((match.start(), match.end()))

    # 7. 表格分隔符行
    for match in re.finditer(r'\|[\-\s\|]+\|', content):
        exclude_regions.append((match.start(), match.end()))

    exclude_regions = merge_regions(exclude_regions)

    # 逐字符处理
    new_content = []
    i = 0
    n = len(content)

    while i < n:
        # 检查是否在排除区域
        in_exclude = False
        for start, end in exclude_regions:
            if start <= i < end:
                in_exclude = True
                new_content.append(content[i:end])
                i = end
                break

        if in_exclude:
            continue

        if i >= n:
            break

        char = content[i]

        # 转换标点符号
        if char == ',':
            if should_convert_comma(content, i):
                new_content.append('，')
            else:
                new_content.append(',')
        elif char == '.':
            if should_convert_period(content, i):
                new_content.append('。')
            else:
                new_content.append('.')
        elif char == ';':
            if should_convert_semicolon(content, i):
                new_content.append('；')
            else:
                new_content.append(';')
        elif char == '!':
            if should_convert_exclamation(content, i):
                new_content.append('！')
            else:
                new_content.append('!')
        elif char == '?':
            if should_convert_question(content, i):
                new_content.append('？')
            else:
                new_content.append('?')
        else:
            new_content.append(char)

        i += 1

    return ''.join(new_content)


def should_convert_comma(text: str, pos: int) -> bool:
    """判断逗号是否应该转换"""
    # 检查前后是否有中文字符
    before = get_context_before(text, pos, 5)
    after = get_context_after(text, pos, 5)

    # 数字中的逗号不转换（如 1,000）
    if re.search(r'\d$', before) and re.search(r'^\d', after):
        return False

    # 如果前后有中文字符，转换
    if has_chinese(before) or has_chinese(after):
        return True

    return False


def should_convert_period(text: str, pos: int) -> bool:
    """判断句号是否应该转换"""
    before = get_context_before(text, pos, 10)
    after = get_context_after(text, pos, 5)

    # 数字中的小数点不转换（如 3.14）
    if re.search(r'\d$', before) and re.search(r'^\d', after):
        return False

    # 文件扩展名不转换（如 .md, .py）
    if re.search(r'\w$', before) and re.search(r'^[a-z]{2,4}(\s|$)', after):
        return False

    # 英文缩写不转换（如 Dr., Mr.）
    if re.search(r'\b[A-Z][a-z]?$', before):
        return False

    # 如果前面有中文且后面是换行或空格，转换
    if has_chinese(before) and re.match(r'^(\s|$|\n)', after):
        return True

    return False


def should_convert_semicolon(text: str, pos: int) -> bool:
    """判断分号是否应该转换"""
    before = get_context_before(text, pos, 10)
    after = get_context_after(text, pos, 10)

    # 如果前后有中文字符，转换
    if has_chinese(before) or has_chinese(after):
        return True

    return False


def should_convert_exclamation(text: str, pos: int) -> bool:
    """判断感叹号是否应该转换"""
    before = get_context_before(text, pos, 10)

    # 如果前面有中文字符，转换
    if has_chinese(before):
        return True

    return False


def should_convert_question(text: str, pos: int) -> bool:
    """判断问号是否应该转换"""
    before = get_context_before(text, pos, 10)

    # 如果前面有中文字符，转换
    if has_chinese(before):
        return True

    return False


def get_context_before(text: str, pos: int, length: int) -> str:
    """获取前面的上下文"""
    start = max(0, pos - length)
    return text[start:pos]


def get_context_after(text: str, pos: int, length: int) -> str:
    """获取后面的上下文"""
    end = min(len(text), pos + length + 1)
    return text[pos+1:end]


def has_chinese(text: str) -> bool:
    """检查是否包含中文字符"""
    return bool(re.search(r'[\u4e00-\u9fff]', text))


def merge_regions(regions: list) -> list:
    """合并重叠区域"""
    if not regions:
        return []
    regions = sorted(regions, key=lambda x: x[0])
    merged = [regions[0]]
    for current in regions[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法： python3 fix_all_punctuation.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    fix_all_punctuation(file_path)
