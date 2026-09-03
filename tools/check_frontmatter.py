#!/usr/bin/env python3
"""Validate that _pages blog posts have complete frontmatter.

Modes:
  (default)  check staged .md files under _pages/  (for git pre-commit hook)
  --file PATH   validate a single file            (for testing / manual use)

Exits non-zero if any required field is missing or inconsistent, so a
pre-commit hook can block the commit until the frontmatter is complete.

Required fields (see .claude/skills/blog/modules/structure/frontmatter.md):
  title, date, last_modified_at, tags, description,
  image, thumbnail, author, lang
"""
import os
import re
import sys
import subprocess

REQUIRED = [
    "title", "date", "last_modified_at", "tags",
    "description", "image", "thumbnail", "author", "lang",
]

# Paths that should NOT carry frontmatter, per project rules.
EXCLUDE_DIRS = ("/archive/", "/diary/")
EXCLUDE_FILES = ("about.md", "index.md")


def is_excluded(path):
    low = path.lower().replace("\\", "/")
    for d in EXCLUDE_DIRS:
        if d in low:
            return True
    if os.path.basename(low) in EXCLUDE_FILES:
        return True
    return False


def parse_frontmatter(text):
    """Return dict of frontmatter key->raw value, or None if absent."""
    if not text.lstrip().startswith("---"):
        return None
    lines = text.split("\n")
    end = None
    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            end = i
            break
    if end is None:
        return None
    data = {}
    for line in lines[1:end]:
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*):\s*(.*)$", line)
        if m:
            key = m.group(1)
            val = m.group(2).strip()
            if len(val) >= 2 and val[0] in "\"'" and val[-1] == val[0]:
                val = val[1:-1]
            data[key] = val
    return data


def check_file(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            text = f.read()
    except OSError as e:
        return [f"无法读取文件: {e}"]
    data = parse_frontmatter(text)
    errors = []
    if data is None:
        errors.append("缺少 frontmatter（文件不以 --- 开头）")
        return errors
    for key in REQUIRED:
        if key not in data or data[key] == "":
            errors.append(f"缺少字段 `{key}`")
    img = data.get("image", "")
    thumb = data.get("thumbnail", "")
    if img and thumb and img != thumb:
        errors.append("image 与 thumbnail 不一致")
    for fld in ("image", "thumbnail"):
        if "empty.jpg" in data.get(fld, "").lower():
            errors.append(f"{fld} 使用了禁止的 empty.jpg")
    return errors


def main():
    if "--file" in sys.argv:
        target = sys.argv[sys.argv.index("--file") + 1]
        errs = check_file(target)
        if errs:
            print(f"[FAIL] {target}")
            for e in errs:
                print("  - " + e)
            sys.exit(1)
        print(f"[OK] {target}")
        sys.exit(0)

    # staged mode (for git pre-commit)
    try:
        root = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True).strip()
    except (subprocess.CalledProcessError, OSError):
        sys.exit(0)
    try:
        out = subprocess.check_output(
            ["git", "diff", "--cached", "--name-only",
             "--diff-filter=ACMR"], text=True).strip()
    except (subprocess.CalledProcessError, OSError):
        out = ""
    candidates = []
    for f in out.split("\n"):
        f = f.strip()
        if not f.lower().endswith(".md"):
            continue
        if "_pages" not in f.lower().replace("\\", "/"):
            continue
        if is_excluded(f):
            continue
        candidates.append(f)

    if not candidates:
        sys.exit(0)

    bad = {}
    for rel in candidates:
        errs = check_file(os.path.join(root, rel))
        if errs:
            bad[rel] = errs

    if bad:
        print("=" * 64)
        print("PRE-COMMIT 拦截：以下 _pages 文章 frontmatter 不完整，提交已阻止。")
        print("参考 .claude/skills/blog/modules/structure/frontmatter.md")
        print("=" * 64)
        for rel, errs in bad.items():
            print(f"\n{rel}")
            for e in errs:
                print("  - " + e)
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
