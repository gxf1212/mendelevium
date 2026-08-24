#!/usr/bin/env python3
import os
import re
from pathlib import Path

def add_last_modified_at(file_path):
    """Add last_modified_at field after date field if not present"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Check if file has date: field
        date_match = re.search(r'^date:\s*(.+)$', content, re.MULTILINE)
        if not date_match:
            return False  # No date field found

        # Check if last_modified_at already exists
        if re.search(r'^last_modified_at:', content, re.MULTILINE):
            return False  # Already has last_modified_at

        # Get the date value
        date_value = date_match.group(1).strip()

        # Insert last_modified_at after date line
        new_content = re.sub(
            r'(^date:\s*.+$)',
            r'\1\nlast_modified_at: ' + date_value,
            content,
            count=1,
            flags=re.MULTILINE
        )

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(new_content)

        return True

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False

def main():
    # Find all .md files
    md_files = list(Path('.').rglob('*.md'))

    # Exclude certain directories
    exclude_dirs = {'.git', '_site', '.jekyll-cache', 'vendor', '.history'}
    md_files = [f for f in md_files if not any(excl in str(f) for excl in exclude_dirs)]

    modified_count = 0
    for md_file in md_files:
        if add_last_modified_at(md_file):
            print(f"Updated: {md_file}")
            modified_count += 1

    print(f"\nTotal files modified: {modified_count}")

if __name__ == '__main__':
    main()
