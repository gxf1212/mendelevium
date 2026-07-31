#!/usr/bin/env python3
"""Read and search PDF text using PyMuPDF (fitz)."""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Read and search PDF text with PyMuPDF")
    parser.add_argument("pdf", help="Path to PDF file")
    parser.add_argument("keyword", nargs="?", default="", help="Keyword to search for")
    parser.add_argument("-C", "--context", type=int, default=0,
                        help="Lines of context before/after each match")
    parser.add_argument("-l", "--lines", type=int, default=0,
                        help="Show first N lines of the PDF (no search)")
    parser.add_argument("-p", "--page", type=int,
                        help="Show full text of a specific page (1-indexed)")
    parser.add_argument("-P", "--page-range",
                        help="Page range to display, e.g. '3-5' or '8-' or '-10' (1-indexed)")
    parser.add_argument("-r", "--regex", action="store_true",
                        help="Treat keyword as regex pattern")
    parser.add_argument("-i", "--ignore-case", action="store_true",
                        help="Case-insensitive search")
    parser.add_argument("--pages", action="store_true",
                        help="Show page numbers before each match (use with keyword)")
    parser.add_argument("--page-num", action="store_true",
                        help="Show page number before each line (use with -p/-P/--lines)")
    args = parser.parse_args()

    if not args.keyword and not args.lines and not args.page and not args.page_range:
        parser.print_help()
        return

    import re
    import fitz

    try:
        doc = fitz.open(args.pdf)
    except Exception as e:
        print(f"Error: cannot open {args.pdf}: {e}", file=sys.stderr)
        return

    # ---- --lines: first N lines ----
    if args.lines:
        count = 0
        for i in range(len(doc)):
            for line in doc[i].get_text().split("\n"):
                if args.page_num:
                    print(f"[{i+1}] {line}")
                else:
                    print(line)
                count += 1
                if count >= args.lines:
                    doc.close()
                    return
        doc.close()
        return

    # ---- -p / -P: specific page(s) ----
    if args.page or args.page_range:
        pages = []
        if args.page:
            pages.append(args.page)
        if args.page_range:
            parts = args.page_range.split("-")
            start = max(1, int(parts[0]) if parts[0] else 1)
            end = len(doc) if parts[1] == "" else int(parts[1])
            for p in range(start, min(end + 1, len(doc) + 1)):
                pages.append(p)

        pages.sort()
        for p in pages:
            if p < 1 or p > len(doc):
                print(f"Error: page {p} out of range (1-{len(doc)})", file=sys.stderr)
                continue
            page = doc[p - 1]
            text = page.get_text()
            for line in text.split("\n"):
                if args.page_num:
                    print(f"[{p}] {line}")
                else:
                    print(line)
        doc.close()
        return

    # ---- keyword search ----
    pattern = args.keyword
    flags = re.IGNORECASE if args.ignore_case else 0
    if args.regex:
        compiled = re.compile(pattern, flags)
    else:
        compiled = re.compile(re.escape(pattern), flags)

    total = 0
    for i in range(len(doc)):
        lines = doc[i].get_text().split("\n")
        for j, line in enumerate(lines):
            if compiled.search(line):
                total += 1
                start = max(0, j - args.context)
                end = min(len(lines), j + args.context + 1)
                for k in range(start, end):
                    marker = ">>>" if k == j else "   "
                    if args.page_num:
                        print(f"[{i+1}] {marker} {lines[k]}")
                    else:
                        print(f"{marker} {lines[k]}")
    doc.close()
    if args.keyword:
        print(f"\n--- {total} match(es) ---", file=sys.stderr)


if __name__ == "__main__":
    main()
