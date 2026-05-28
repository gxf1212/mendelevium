#!/bin/bash
# 快速修复Markdown格式问题的shell脚本

if [ $# -eq 0 ]; then
    echo "用法: $0 <file_path>"
    exit 1
fi

file="$1"

if [ ! -f "$file" ]; then
    echo "❌ 文件不存在: $file"
    exit 1
fi

echo "正在修复: $file"

# 1. 删除列表项之间的空行（使用Python）
python3 -c "
import re
with open('$file', 'r', encoding='utf-8') as f:
    lines = f.readlines()

new_lines = []
i = 0
while i < len(lines):
    if (i < len(lines) - 2 and
        re.match(r'^- ', lines[i]) and
        lines[i + 1].strip() == '' and
        re.match(r'^- ', lines[i + 2])):
        new_lines.append(lines[i])
        new_lines.append(lines[i + 2])
        i += 3
    else:
        new_lines.append(lines[i])
        i += 1

with open('$file', 'w', encoding='utf-8') as f:
    f.writelines(new_lines)
"

# 2. 修复Unicode引号
sed -i 's/"/"/g' "$file"
sed -i 's/"/"/g' "$file"

# 3. 修复百分比加粗
sed -i 's/\*\*\([0-9]\+\)%\*\*/**\1**%/g' "$file"

# 4. 修复引用格式中的中文标点
sed -i '/^-\*\*引用格式\*\*:.*$/s/（/(/g' "$file"
sed -i '/^-\*\*引用格式\*\*:.*$/s/）/)/g' "$file"
sed -i '/^-\*\*引用格式\*\*:.*$/s/。/./g' "$file"
sed -i '/^-\*\*引用格式\*\*:.*$/s/：/:/g' "$file"

echo "✓ 修复完成: $file"
