#!/usr/bin/env python3
"""
博客文章质量自动检查脚本
根据 .claude/skills/blog/skill.md 的质量检查清单进行自动化验证
"""

import re
import sys
from pathlib import Path
from typing import List, Tuple


class BlogQualityChecker:
    def __init__(self, file_path: str):
        self.file_path = Path(file_path)
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.content = ""
        self.frontmatter = ""
        self.body = ""

    def load_file(self):
        """读取文件并分离frontmatter和正文"""
        with open(self.file_path, 'r', encoding='utf-8') as f:
            self.content = f.read()

        parts = self.content.split('---', 2)
        if len(parts) >= 3:
            self.frontmatter = parts[1]
            self.body = parts[2]
        else:
            self.errors.append("❌ 缺少frontmatter")
            self.body = self.content

    def check_frontmatter(self):
        """检查frontmatter完整性"""
        required_fields = ['title', 'date', 'tags', 'description', 'image', 'thumbnail', 'author', 'lang']

        for field in required_fields:
            if f'{field}:' not in self.frontmatter:
                self.errors.append(f"❌ Frontmatter缺少字段: {field}")

        # 检查image和thumbnail是否一致
        image_match = re.search(r'image:\s*"([^"]+)"', self.frontmatter)
        thumbnail_match = re.search(r'thumbnail:\s*"([^"]+)"', self.frontmatter)

        if image_match and thumbnail_match:
            if image_match.group(1) != thumbnail_match.group(1):
                self.errors.append(f"❌ image和thumbnail不一致: {image_match.group(1)} vs {thumbnail_match.group(1)}")

        # 检查是否使用了empty.jpg
        if 'empty.jpg' in self.frontmatter:
            self.errors.append("❌ 使用了empty.jpg作为缩略图，应该使用随机缩略图")

        # 检查date格式
        date_match = re.search(r'date:\s*"(\d{4}-\d{2}-\d{2})"', self.frontmatter)
        if not date_match:
            self.errors.append("❌ date格式不正确，应为YYYY-MM-DD")

    def check_structure(self):
        """检查文章结构"""
        required_sections = [
            '## 本文信息',
            '## 摘要',
            '### 核心结论',
            '## 背景',
            '### 关键科学问题',
            '### 创新点',
            '## 研究内容',
            '## Q&A',
            '## 关键结论与批判性总结'
        ]

        for section in required_sections:
            if section not in self.body:
                self.warnings.append(f"⚠️ 可能缺少章节: {section}")

        # 检查本文信息的完整性
        if '## 本文信息' in self.body:
            info_section = self.body.split('## 本文信息', 1)[1].split('##', 1)[0]
            required_info = ['标题', '作者', '发表时间', '单位', '引用格式']
            for info in required_info:
                if info not in info_section:
                    self.warnings.append(f"⚠️ 本文信息可能缺少: {info}")

    def check_formulas(self):
        """检查公式格式"""
        # 检查双backslash（高危错误）
        if re.search(r'\$[^\$]*\\\\[^\$]*\$', self.body):
            self.errors.append("❌ 🚨 发现双backslash公式（高危错误）: $\\\\xxx$")

        # 检查是否有未转义的微分符号
        if re.search(r'\$[^\$]*\bd[A-Za-z]', self.body):
            self.warnings.append("⚠️ 可能存在未使用\\mathrm的微分符号，应为 $\\mathrm{d}\\xi$")

    def check_bold_format(self):
        """检查加粗格式"""
        # 更精确的检测：逐行检查，避免跨越多个加粗标记

        # 检查加粗包含标点 - 逐行检查，避免跨行匹配
        punct_in_bold_errors = []
        for line in self.body.split('\n'):
            # 在单行内查找所有**...**模式
            # 排除列表项（以-或>开头，后面有**xxx**：格式）
            if line.strip().startswith(('- **', '> **')):
                # 这是列表标题，检查是否有：紧跟在**之后
                if re.search(r'\*\*[^*]*[。，；：！?][^*]*\*:', line):
                    # 冒号在**外，这是正确的
                    continue

            # 查找行内所有**...**模式
            bold_matches = re.finditer(r'\*\*([^*]+?)\*\*', line)
            for match in bold_matches:
                content = match.group(1)
                # 检查内容结尾是否有标点（不是数字或字母的标点）
                if content and content[-1] in '。，；：！？.':
                    # 排除特殊情况：如果内容是纯数字或单个字符，标点可以在内
                    if not (content[:-1].isdigit() or len(content) <= 2):
                        punct_in_bold_errors.append(f"❌ 发现加粗包含标点: **{content[:30]}...**")
                        break  # 每行只报告一次

        self.errors.extend(punct_in_bold_errors[:5])  # 最多显示5个错误

        # 检查加粗包含百分号
        pct_bold = re.findall(r'\*\*(\d+%)\*\*', self.body)
        if pct_bold:
            for match in pct_bold[:3]:
                clean_num = match[:-1]  # 去掉百分号
                self.errors.append(f"❌ 发现加粗包含百分号: **{match}**，应为 **{clean_num}**%")

        # 检查加粗包含引号
        quote_bold = re.findall(r'\*\*("[^"]+")\*\*', self.body)
        if quote_bold:
            for match in quote_bold[:3]:
                clean_text = match.replace('"', '')
                self.errors.append(f"❌ 发现加粗包含引号: **{match}**，应为 \"{clean_text}\"")

        # 检查加粗包含括号（括号必须在**外面）
        # 逐行检查 **xxx（yyy）** 这种模式
        paren_bold_errors = []
        for line in self.body.split('\n'):
            # 查找行内所有**xxx（yyy）**模式
            bold_matches = re.finditer(r'\*\*([^*]+?（[^）]*？）)\*\*', line)
            for match in bold_matches:
                content = match.group(1)
                # 检查是否是 **xxx（yyy）** 格式
                if '（' in content and content.endswith('）'):
                    parts = content.split('（', 1)
                    if len(parts) == 2:
                        before_paren = parts[0]
                        after_paren = parts[1][:-1]  # 去掉最后的）
                        paren_bold_errors.append(f"❌ 发现加粗包含括号: **{content[:30]}...**，应为 **{before_paren}**（{after_paren}）")
                        break  # 每行只报告一次

        self.errors.extend(paren_bold_errors[:3])  # 最多显示3个错误

    def check_punctuation(self):
        """检查中文标点"""
        # 保护代码块和公式
        protected_body = re.sub(r'```[\s\S]*?```', '', self.body)
        protected_body = re.sub(r'`[^`]+`', '', protected_body)
        protected_body = re.sub(r'\$[^\$]+\$', '', protected_body)

        # 检查中文后的英文逗号
        if re.search(r'[\u4e00-\u9fff],', protected_body):
            self.errors.append("❌ 发现中文后的英文逗号，应使用中文逗号，")

        # 🆕 检查数字+单位+英文逗号（如 15 ps, 或 3.5 Å,）
        if re.search(r'\d+\.?\d*\s*(Å|ps|ns|kcal|mol),', protected_body):
            self.errors.append("❌ 发现数字单位后的英文逗号，如 15 ps, 应为 15 ps，")

        # 🆕 检查**xxx**,的情况
        if re.search(r'\*\*[^\*]+\*\*,', protected_body):
            self.errors.append("❌ 发现加粗后紧跟英文逗号，如 **xxx**,应为 **xxx**，")

        # 检查中文后的英文分号
        if re.search(r'[\u4e00-\u9fff];', protected_body):
            self.errors.append("❌ 发现中文后的英文分号，应使用中文分号；")

        # 检查英文引号（除了代码）
        if re.search(r'[\u4e00-\u9fff]"[^"]*"', protected_body):
            self.warnings.append("⚠️ 可能存在英文引号，应使用中文引号""")

        # 检查中文内容中的英文括号
        # 匹配：中文文字 + ( + 内容 + ) 或 中文后紧跟(
        if re.search(r'[\u4e00-\u9fff]\([^)]*\)', protected_body):
            self.errors.append("❌ 发现中文内容中的英文括号，如(PET)，应使用中文括号（）")

        # 额外检查：纯中文在括号内
        if re.search(r'\([\u4e00-\u9fff]+[，、。；：！？]*[\u4e00-\u9fff]*\)', protected_body):
            self.errors.append("❌ 发现括号内纯中文内容使用英文括号，应使用中文括号（）")

    def check_escape_characters(self):
        """检查多余的转义符"""
        if r'\&' in self.body:
            self.errors.append(r"❌ 发现 \& 转义符，应为 &")

        if re.search(r'\d+\\.', self.body):
            self.warnings.append(r"⚠️ 可能存在 1\. 这样的转义符，应为 1.")

    def check_irregular_symbols(self):
        """检查不规范符号"""
        irregular_patterns = [
            (r'~\d+', '~700应改为约700个'),
            (r'\d+\+(?!\d)', '52000+应改为52000余个'),
            (r'\+\d+(?!\d)', '+30应改为30'),
        ]

        for pattern, msg in irregular_patterns:
            if re.search(pattern, self.body):
                self.warnings.append(f"⚠️ 发现不规范符号: {msg}")

    def check_mermaid(self):
        """检查Mermaid图语法"""
        mermaid_blocks = re.findall(r'```mermaid([\s\S]*?)```', self.body)

        for i, block in enumerate(mermaid_blocks):
            # 检查是否使用graph TB
            if 'graph TB' not in block and 'mindmap' not in block:
                self.warnings.append(f"⚠️ Mermaid图{i+1}可能未使用graph TB布局")

            # 检查是否有---分隔符
            if '---' in block:
                self.errors.append(f"❌ Mermaid图{i+1}包含---分隔符，会被误认为frontmatter分隔符")

            # 检查列表格式（不应有空格）
            if re.search(r'\d+\.\s+[^\s]', block):
                self.warnings.append(f"⚠️ Mermaid图{i+1}列表格式可能有误，应为2.xx而非2. xx")

    def check_figures(self):
        """检查图片和图注"""
        # 查找所有图片
        figures = re.findall(r'!\[([^\]]*)\]\(([^)]+)\)', self.body)

        if len(figures) < 4:
            self.warnings.append(f"⚠️ 正文图片数量较少（{len(figures)}张），建议至少4张")
        elif len(figures) > 8:
            self.warnings.append(f"⚠️ 正文图片数量较多（{len(figures)}张），建议最多8张")

        # 检查图注完整性
        figure_captions = re.findall(r'\*\*图\d+[：:][^*]+\*\*', self.body)
        if len(figure_captions) != len(figures):
            self.warnings.append(f"⚠️ 图片数量（{len(figures)}）与图注数量（{len(figure_captions)}）不匹配")

    def check_length(self):
        """检查文章长度"""
        lines = self.content.split('\n')
        total_lines = len(lines)

        # 主文档应该250-300行
        if 'appendix' not in str(self.file_path):
            if total_lines < 250:
                self.warnings.append(f"⚠️ 主文档行数偏少（{total_lines}行），建议250-300行")
            elif total_lines > 300:
                self.warnings.append(f"⚠️ 主文档行数偏多（{total_lines}行），建议拆分附录")

    def check_title_words(self):
        """检查标题党词汇"""
        title_spam_words = ['重塑', '颠覆', '革命', '突破', '新范式']

        for word in title_spam_words:
            if word in self.content:
                self.warnings.append(f"⚠️ 发现标题党词汇: {word}")

    def run_all_checks(self) -> Tuple[List[str], List[str]]:
        """运行所有检查"""
        self.load_file()

        print(f"\n🔍 正在检查: {self.file_path.name}")
        print("=" * 60)

        self.check_frontmatter()
        self.check_structure()
        self.check_formulas()
        self.check_bold_format()
        self.check_punctuation()
        self.check_escape_characters()
        self.check_irregular_symbols()
        self.check_mermaid()
        self.check_figures()
        self.check_length()
        self.check_title_words()

        return self.errors, self.warnings

    def print_results(self):
        """打印检查结果"""
        if self.errors:
            print("\n❌ 错误（必须修复）:")
            for error in self.errors:
                print(f"  {error}")

        if self.warnings:
            print("\n⚠️ 警告（建议修复）:")
            for warning in self.warnings:
                print(f"  {warning}")

        if not self.errors and not self.warnings:
            print("\n✅ 所有检查通过！")

        print("\n" + "=" * 60)
        print(f"总计: {len(self.errors)}个错误, {len(self.warnings)}个警告")


def main():
    if len(sys.argv) < 2:
        print("用法: python3 check_blog_quality.py <markdown文件路径>")
        print("示例: python3 check_blog_quality.py '_pages/xxx/2025-11-22-xxx.md'")
        sys.exit(1)

    file_path = sys.argv[1]

    if not Path(file_path).exists():
        print(f"❌ 文件不存在: {file_path}")
        sys.exit(1)

    checker = BlogQualityChecker(file_path)
    errors, warnings = checker.run_all_checks()
    checker.print_results()

    # 返回错误码：有错误返回1，仅有警告返回0
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
