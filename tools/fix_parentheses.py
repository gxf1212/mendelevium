#!/usr/bin/env python3
"""
自动修复Markdown文件中的英文括号为中文括号

功能：
1. 将中文内容中的英文括号 () 转换为中文括号 （）
2. 保护特殊区域：frontmatter、代码块、行内代码、公式、mermaid图
3. 保留特殊场景：纯英文内容、坐标、代码中的括号

使用方法：
    python3 tools/fix_parentheses.py <markdown_file>

规则：
- 中文后的括号：塑料(PET) → 塑料（PET）
- 括号内纯中文：(北京大学) → （北京大学）
- 混合内容需判断：使用了方法(QM/MM) → 使用了方法（QM/MM）

保留：
- 公式中：$f(x)$
- 代码中：print("hello")
- 坐标：(1, 2, 3)
- 纯英文：(PDB ID: 5XG0)
- Frontmatter: tags: [tag1, tag2]
"""

import re
import sys
from pathlib import Path


def is_pure_english_content(content: str) -> bool:
    """判断括号内是否为纯英文内容（可保留英文括号）"""
    # 坐标格式：(x, y, z) 或 (1, 2, 3)
    if re.match(r'^[a-zA-Z0-9\s,\.\-]+$', content):
        return True

    # PDB ID等纯英文缩写
    if re.match(r'^[A-Z\s:]+[0-9A-Z]+$', content):
        return True

    return False


def fix_parentheses(content: str) -> str:
    """修复括号：英文 () → 中文 （）"""

    # 提取并保护 frontmatter
    frontmatter_match = re.match(r'^---\n(.*?\n)---\n', content, re.DOTALL)
    if frontmatter_match:
        frontmatter = frontmatter_match.group(0)
        body = content[len(frontmatter):]
    else:
        frontmatter = ""
        body = content

    # 保护区域：代码块、行内代码、公式、mermaid
    protected_regions = []

    # 1. 保护代码块
    def protect_code_block(match):
        idx = len(protected_regions)
        protected_regions.append(match.group(0))
        return f"__PROTECTED_BLOCK_{idx}__"

    body = re.sub(r'```[\s\S]*?```', protect_code_block, body)

    # 2. 保护行内代码
    body = re.sub(r'`[^`]+`', protect_code_block, body)

    # 3. 保护行间公式
    body = re.sub(r'\$\$[\s\S]*?\$\$', protect_code_block, body)

    # 4. 保护行内公式
    body = re.sub(r'\$[^\$]+\$', protect_code_block, body)

    # ========== 核心修复逻辑 ==========

    # 统一处理所有括号：遍历所有 (xxx) 结构，根据上下文判断是否需要转换
    def fix_all_parentheses(match):
        before = match.group(1) if match.group(1) else ""  # 括号前的字符
        content = match.group(2)  # 括号内容
        after = match.group(3) if match.group(3) else ""   # 括号后的字符

        # 保留规则1: 纯英文坐标或ID (x, y, z) 或 (PDB ID: 5XG0)
        if is_pure_english_content(content):
            # 但如果前面是中文且不是全大写缩写，仍需转换
            # 例外：PDB ID(5XG0) 保留英文括号
            if before and re.search(r'[\u4e00-\u9fff]$', before):
                # 检查是否是 "ID(xxx)" 这种缩写格式
                if not re.search(r'[A-Z]{2,}$', before):
                    return f"{before}（{content}）{after}"
            return match.group(0)

        # 转换规则1: 前面有中文
        if before and re.search(r'[\u4e00-\u9fff]', before):
            return f"{before}（{content}）{after}"

        # 转换规则2: 内容是纯中文或中文为主
        if re.search(r'[\u4e00-\u9fff]', content):
            # 如果包含中文，转为中文括号
            return f"{before}（{content}）{after}"

        # 转换规则3: 内容是英文缩写，但前面紧邻中文
        # 通过lookahead检查（已在before中捕获）

        # 默认保留
        return match.group(0)

    # 匹配模式：(可选的前导字符)(括号内容)(可选的后续字符)
    # 使用更宽的匹配来捕获上下文
    body = re.sub(
        r'([\u4e00-\u9fff]?)\(([^)]+)\)([\u4e00-\u9fff]?)',
        fix_all_parentheses,
        body
    )

    # 补充修复：处理任何残留的半转换情况
    body = re.sub(r'（([^）]+)\)', r'（\1）', body)  # 左中右英 → 全中
    body = re.sub(r'\(([^)]+)）', r'（\1）', body)   # 左英右中 → 全中

    # ========== 恢复保护区域 ==========

    for idx, region in enumerate(protected_regions):
        body = body.replace(f"__PROTECTED_BLOCK_{idx}__", region)

    return frontmatter + body


def main():
    if len(sys.argv) < 2:
        print("用法: python3 tools/fix_parentheses.py <markdown_file>")
        print("\n示例:")
        print("  python3 tools/fix_parentheses.py _pages/blog.md")
        sys.exit(1)

    file_path = Path(sys.argv[1])

    if not file_path.exists():
        print(f"❌ 文件不存在: {file_path}")
        sys.exit(1)

    if not file_path.suffix == '.md':
        print(f"⚠️ 警告: {file_path} 不是Markdown文件")

    # 读取文件
    with open(file_path, 'r', encoding='utf-8') as f:
        original = f.read()

    # 修复括号
    fixed = fix_parentheses(original)

    # 如果有变化，则写回文件
    if fixed != original:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(fixed)

        # 统计变化
        changes = sum(1 for a, b in zip(original, fixed) if a != b)
        print(f"✅ 已修复 {file_path}")
        print(f"   变更字符数: {changes}")

        # 显示部分示例
        if '(PET)' in original and '（PET）' in fixed:
            print("   示例: (PET) → （PET）")
    else:
        print(f"ℹ️ {file_path} 无需修复")


if __name__ == "__main__":
    main()
