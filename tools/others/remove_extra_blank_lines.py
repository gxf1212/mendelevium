#!/usr/bin/env python3
"""
删除文件中多余的空行，只保留单个空行。

用法：
    python3 tools/remove_extra_blank_lines.py <file_path>

示例：
    python3 tools/remove_extra_blank_lines.py _pages/test.md
"""

import sys
import re
from pathlib import Path


def remove_extra_blank_lines(file_path: str) -> None:
    """
    删除文件中的连续空行，只保留一个空行。

    Args:
        file_path: 要处理的文件路径
    """
    path = Path(file_path)

    if not path.exists():
        print(f"错误：文件不存在: {file_path}")
        sys.exit(1)

    # 读取文件内容
    content = path.read_text(encoding='utf-8')

    # 使用正则表达式将连续的空行替换为单个空行
    # \n\n+ 匹配两个或更多连续换行符（即一个或多个空行）
    # 替换为 \n\n（即一个空行）
    new_content = re.sub(r'\n\n+', '\n\n', content)

    # 如果内容有变化，写回文件
    if new_content != content:
        path.write_text(new_content, encoding='utf-8')
        print(f"✓ 已处理: {file_path}")
    else:
        print(f"- 无需处理: {file_path}")


def main():
    if len(sys.argv) < 2:
        print("用法: python3 remove_extra_blank_lines.py <file1.md> [file2.md] ...")
        print("\n示例:")
        print("  python3 remove_extra_blank_lines.py _pages/test.md")
        print("  python3 remove_extra_blank_lines.py _pages/*.md")
        sys.exit(1)

    # 处理所有输入的文件
    for file_path in sys.argv[1:]:
        try:
            remove_extra_blank_lines(file_path)
        except Exception as e:
            print(f"✗ 处理失败 {file_path}: {e}")


if __name__ == "__main__":
    main()
