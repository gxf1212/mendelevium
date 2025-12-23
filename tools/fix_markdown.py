#!/usr/bin/env python3
"""
Markdown格式综合修复工具 - 合并所有格式修复功能

功能：
1. 加粗格式修复
   - **text。** → **text**。
   - **"text"** → "**text**"
   - **text（note）** → **text**（note）
   - **X%** → **X**%

2. 中文标点符号转换
   - 英文引号 " → 中文引号 " "
   - 英文冒号 : → 中文冒号 ：
   - 英文逗号 , → 中文逗号 ，
   - 英文句号 . → 中文句号 。
   - 英文分号 ; → 中文分号 ；
   - 英文感叹号 ! → 中文感叹号 ！
   - 英文问号 ? → 中文问号 ？
   - 英文括号 () → 中文括号 （）

3. 列表格式
   - 删除列表项之间的空行

4. 智能排除特殊区域
   - Frontmatter (YAML)
   - 代码块和行内代码
   - 公式 ($...$, $$...$$)
   - Mermaid图表
   - URL和链接
   - 引用的英文原文

用法：
    python3 fix_markdown.py <file.md> [file2.md...]
    python3 fix_markdown.py --check <file.md>  # 仅检查不修改
"""

import sys
import re
from pathlib import Path
from typing import List, Tuple


class MarkdownFixer:
    def __init__(self):
        self.stats = {
            '加粗格式': 0,
            '引号': 0,
            '冒号': 0,
            '逗号': 0,
            '句号': 0,
            '括号': 0,
            '列表': 0,
            '其他标点': 0,
        }

    def find_excluded_regions(self, content: str) -> List[Tuple[int, int]]:
        """找到所有需要排除的区域"""
        regions = []

        # 1. Frontmatter
        if content.startswith('---\n') or content.startswith('---\r\n'):
            end = content.find('\n---\n', 4)
            if end == -1:
                end = content.find('\n---\r\n', 4)
            if end != -1:
                regions.append((0, end + 5))

        # 2. 代码块
        for m in re.finditer(r'```.*?```', content, re.DOTALL):
            regions.append((m.start(), m.end()))

        # 3. 行内代码
        for m in re.finditer(r'`[^`]+`', content):
            regions.append((m.start(), m.end()))

        # 4. 公式
        for m in re.finditer(r'\$\$.*?\$\$', content, re.DOTALL):
            regions.append((m.start(), m.end()))
        for m in re.finditer(r'\$[^$\n]+\$', content):
            regions.append((m.start(), m.end()))

        # 5. URL和链接
        for m in re.finditer(r'\[([^\]]+)\]\(([^)]+)\)', content):
            regions.append((m.start(), m.end()))
        for m in re.finditer(r'https?://[^\s\)]+', content):
            regions.append((m.start(), m.end()))
        for m in re.finditer(r'DOI:\s*[^\s]+', content, re.IGNORECASE):
            regions.append((m.start(), m.end()))

        # 6. HTML标签
        for m in re.finditer(r'<[^>]+>', content):
            regions.append((m.start(), m.end()))

        # 7. 表格分隔符
        for m in re.finditer(r'\|[\-\s\|]+\|', content):
            regions.append((m.start(), m.end()))

        # 8. 引用的英文原文
        for m in re.finditer(r'"[A-Za-z0-9\s,.\-:;()\']+?"', content):
            if not re.search(r'[\u4e00-\u9fff]', m.group()):
                regions.append((m.start(), m.end()))

        # 合并重叠区域
        regions.sort()
        merged = []
        for start, end in regions:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))

        return merged

    def is_excluded(self, pos: int, regions: List[Tuple[int, int]]) -> bool:
        for start, end in regions:
            if start <= pos < end:
                return True
        return False

    def fix_bold_formatting(self, content: str) -> str:
        """修复加粗格式"""
        count = 0

        # 1. 移除加粗内的标点
        new_content = re.sub(r'\*\*([^*]+?)([。！？，；])\*\*', r'**\1**\2', content)
        count += len(re.findall(r'\*\*[^*]+?[。！？，；]\*\*', content))

        # 2. 引号移到外面
        new_content = re.sub(r'\*\*"([^"]+?)"\*\*', r'"**\1**"', new_content)
        new_content = re.sub(r'\*\*"([^"]+?)"\*\*', r'"**\1**"', new_content)

        # 3. 括号移到外面
        new_content = re.sub(r'\*\*([^*（）]+?)（([^）]+?)）\*\*', r'**\1**（\2）', new_content)
        new_content = re.sub(r'\*\*([^*\(\)]+?)\(([^\)]+?)\)\*\*', r'**\1**(\2)', new_content)

        # 4. 百分号移到外面
        new_content = re.sub(r'\*\*(\d+)%\*\*', r'**\1**%', new_content)

        if new_content != content:
            self.stats['加粗格式'] += count

        return new_content

    def fix_list_spacing(self, content: str) -> str:
        """删除列表项之间的空行"""
        lines = content.split('\n')
        new_lines = []
        i = 0
        count = 0

        while i < len(lines):
            # 检查是否是两个列表项之间有空行
            if (i < len(lines) - 2 and
                re.match(r'^- ', lines[i]) and
                lines[i + 1].strip() == '' and
                re.match(r'^- ', lines[i + 2])):
                new_lines.append(lines[i])
                # 跳过空行
                i += 2
                count += 1
            else:
                new_lines.append(lines[i])
                i += 1

        if count > 0:
            self.stats['列表'] = count

        return '\n'.join(new_lines)

    def fix_punctuation(self, content: str) -> str:
        """修复中文标点符号"""
        regions = self.find_excluded_regions(content)
        result = []
        i = 0
        quote_count = 0

        while i < len(content):
            # 跳过排除区域
            if self.is_excluded(i, regions):
                for start, end in regions:
                    if start <= i < end:
                        result.append(content[i:end])
                        i = end
                        break
                continue

            char = content[i]
            prev = content[i-1] if i > 0 else ''
            next_char = content[i+1] if i < len(content)-1 else ''

            # 1. 英文引号 → 中文引号
            if char == '"':
                result.append('"' if quote_count % 2 == 0 else '"')
                quote_count += 1
                self.stats['引号'] += 1
                i += 1
                continue

            # 2. 英文冒号 → 中文冒号（中文/引号后）
            if char == ':' and re.match(r'[\u4e00-\u9fff"""）]', prev):
                result.append('：')
                self.stats['冒号'] += 1
                i += 1
                continue

            # 3. 英文逗号 → 中文逗号（非数字间）
            if char == ',':
                if not (prev.isdigit() and next_char.isdigit()):
                    # 检查附近是否有中文
                    context = content[max(0, i-3):min(len(content), i+4)]
                    if re.search(r'[\u4e00-\u9fff]', context):
                        result.append('，')
                        self.stats['逗号'] += 1
                        i += 1
                        continue

            # 4. 英文句号 → 中文句号（句末）
            if char == '.' and re.match(r'[\u4e00-\u9fff]', prev):
                if next_char in ['', '\n', ' ', '"', '"', '）', ')']:
                    result.append('。')
                    self.stats['句号'] += 1
                    i += 1
                    continue

            # 5. 英文括号 → 中文括号（含中文）
            if char == '(':
                depth, j = 1, i + 1
                has_chinese = False
                while j < len(content) and depth > 0:
                    if content[j] == '(':
                        depth += 1
                    elif content[j] == ')':
                        depth -= 1
                    elif re.match(r'[\u4e00-\u9fff]', content[j]):
                        has_chinese = True
                    j += 1

                if has_chinese and depth == 0:
                    result.append('（')
                    self.stats['括号'] += 1
                    i += 1
                    continue

            if char == ')':
                depth, j = 1, i - 1
                has_chinese = False
                while j >= 0 and depth > 0:
                    if content[j] == ')':
                        depth += 1
                    elif content[j] == '(':
                        depth -= 1
                    elif re.match(r'[\u4e00-\u9fff]', content[j]):
                        has_chinese = True
                    j -= 1

                if has_chinese and depth == 0:
                    result.append('）')
                    i += 1
                    continue

            # 6. 其他标点（中文后）
            if char in [';', '!', '?'] and re.match(r'[\u4e00-\u9fff]', prev):
                punct_map = {';': '；', '!': '！', '?': '？'}
                result.append(punct_map[char])
                self.stats['其他标点'] += 1
                i += 1
                continue

            result.append(char)
            i += 1

        return ''.join(result)

    def process_file(self, file_path: Path, check_only: bool = False) -> bool:
        """处理单个文件"""
        if not file_path.exists():
            print(f"❌ 文件不存在: {file_path}")
            return False

        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
        except Exception as e:
            print(f"❌ 读取失败: {e}")
            return False

        # 重置统计
        self.stats = {k: 0 for k in self.stats}

        # 应用修复
        original = content
        content = self.fix_bold_formatting(content)
        content = self.fix_list_spacing(content)
        content = self.fix_punctuation(content)

        # 显示结果
        total = sum(self.stats.values())
        if total > 0:
            print(f"\n{'[CHECK] ' if check_only else ''}✓ {file_path.name}")
            for key, count in self.stats.items():
                if count > 0:
                    print(f"  - {key}: {count} 处")

            if not check_only:
                try:
                    with open(file_path, 'w', encoding='utf-8') as f:
                        f.write(content)
                except Exception as e:
                    print(f"❌ 写入失败: {e}")
                    return False
        else:
            print(f"✓ {file_path.name}: 无需修复")

        return True


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    check_only = '--check' in sys.argv
    files = [arg for arg in sys.argv[1:] if not arg.startswith('--')]

    if not files:
        print("错误: 请指定至少一个文件")
        sys.exit(1)

    fixer = MarkdownFixer()

    print("\n" + "="*60)
    print("Markdown格式综合修复工具")
    print("="*60)

    if check_only:
        print("⚠️  检查模式 - 不会修改文件\n")

    success = 0
    for file_path in files:
        if fixer.process_file(Path(file_path), check_only):
            success += 1

    print(f"\n{'='*60}")
    print(f"完成! 处理 {success}/{len(files)} 个文件")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
