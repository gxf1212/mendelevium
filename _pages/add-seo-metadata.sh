#!/bin/bash

# 批量为文章添加 SEO 元数据的脚本
# 使用方法: ./add-seo-metadata.sh <filename>

if [ -z "$1" ]; then
    echo "Usage: $0 <file-path>"
    echo "Example: $0 \"_pages/Some Category/some-article.md\""
    exit 1
fi

FILE_PATH="$1"

if [ ! -f "$FILE_PATH" ]; then
    echo "Error: File not found: $FILE_PATH"
    exit 1
fi

# 读取现有的 frontmatter
read -r -d '' CONTENT < "$FILE_PATH"

# 提取标题
TITLE=$(echo "$CONTENT" | grep '^title:' | head -1 | sed 's/^title: *//' | sed 's/^"//' | sed 's/"$//')

if [ -z "$TITLE" ]; then
    echo "Error: No title found in frontmatter"
    exit 1
fi

# 检查是否已有 description 或 image
if echo "$CONTENT" | grep -q "^description:"; then
    echo "Description already exists in $FILE_PATH"
else
    # 根据标题生成描述（这里可以手动修改）
    echo "Title found: $TITLE"
    echo "Please enter a description (50-160 characters):"
    read DESCRIPTION

    if [ -z "$DESCRIPTION" ]; then
        # 如果没有输入，使用默认描述
        DESCRIPTION="Detailed analysis of $TITLE in computational chemistry and molecular dynamics research."
    fi

    # 添加 description
    sed -i "/^date:/a description: \"$DESCRIPTION\"" "$FILE_PATH"
    echo "Added description to $FILE_PATH"
fi

if echo "$CONTENT" | grep -q "^image:"; then
    echo "Image already exists in $FILE_PATH"
else
    # 显示可用的图片选项
    echo "Available images:"
    echo "1. /assets/img/thumbnail/bricks.webp - 科学主题"
    echo "2. /assets/img/thumbnail/book.jpg - 教程类文章"
    echo "3. /assets/img/thumbnail/nightgardenflower.jpg - 生物/化学主题"
    echo "4. /assets/img/thumbnail/sample.png - 默认图片"
    echo "5. /assets/img/thumbnail/empty.jpg - 通用空白"

    echo "Select image number (1-5) or enter custom path:"
    read choice

    case $choice in
        1) IMAGE="/assets/img/thumbnail/bricks.webp" ;;
        2) IMAGE="/assets/img/thumbnail/book.jpg" ;;
        3) IMAGE="/assets/img/thumbnail/nightgardenflower.jpg" ;;
        4) IMAGE="/assets/img/thumbnail/sample.png" ;;
        5) IMAGE="/assets/img/thumbnail/empty.jpg" ;;
        *) IMAGE="$choice" ;;
    esac

    # 添加 image
    sed -i "/^tags:/i image: \"$IMAGE\"" "$FILE_PATH"
    echo "Added image to $FILE_PATH"
fi

echo "SEO metadata added successfully to $FILE_PATH"
echo "Updated file preview:"
head -15 "$FILE_PATH"