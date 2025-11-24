#!/usr/bin/env python3
"""
Search PDF text extracted by pymupdf.
Usage: python3 search_pdf_text.py <pdf_file> <keyword> [context_lines]
"""

import sys
import fitz
from pathlib import Path

def extract_and_search_pdf(pdf_path, keyword, context_lines=3):
    """Extract PDF text and search for keyword"""

    # Open PDF and extract text
    doc = fitz.open(pdf_path)
    page_count = len(doc)

    full_text = ""
    page_mapping = {}  # Map line numbers to page numbers

    for page_num in range(page_count):
        page = doc[page_num]
        text = page.get_text()
        start_line = len(full_text.split('\n'))
        full_text += f"\n=== PAGE {page_num + 1} ===\n{text}"
        end_line = len(full_text.split('\n'))
        for line_num in range(start_line, end_line):
            page_mapping[line_num] = page_num + 1

    doc.close()

    # Search for keyword
    lines = full_text.split('\n')
    matches = []

    for i, line in enumerate(lines):
        if keyword.lower() in line.lower():
            start = max(0, i - context_lines)
            end = min(len(lines), i + context_lines + 1)
            context = '\n'.join(lines[start:end])
            page_num = page_mapping.get(i, "Unknown")
            matches.append((i, page_num, context))

    # Display results
    if matches:
        print(f"\n{'='*70}")
        print(f"Found {len(matches)} matches for '{keyword}' in {Path(pdf_path).name}")
        print(f"{'='*70}\n")
        for idx, (line_num, page_num, context) in enumerate(matches, 1):
            print(f"--- Match {idx} (Page {page_num}, line {line_num}) ---")
            print(context)
            print()
    else:
        print(f"No matches found for '{keyword}'")

    return matches

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 search_pdf_text.py <pdf_file> <keyword> [context_lines]")
        sys.exit(1)

    pdf_file = sys.argv[1]
    keyword = sys.argv[2]
    context_lines = int(sys.argv[3]) if len(sys.argv) > 3 else 3

    if not Path(pdf_file).exists():
        print(f"Error: PDF file not found: {pdf_file}")
        sys.exit(1)

    extract_and_search_pdf(pdf_file, keyword, context_lines)
