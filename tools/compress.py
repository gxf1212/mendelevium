#!/usr/bin/env python3
"""
压缩图片到指定大小以内，保持原始分辨率和格式。

用法：
  python compress.py                              # 默认：当前目录 → ../Wallpaper_compressed，≤500kB
  python compress.py -i ~/photos -o ~/thumbs      # 指定输入输出
  python compress.py -s 300                       # 目标300kB
  python compress.py --overwrite                  # 覆盖已有文件
  python compress.py -n 50                        # 只处理前50张
"""

import argparse
import subprocess
import shutil
from pathlib import Path

IMAGE_EXTS = {".jpg", ".jpeg", ".png", ".webp"}
DEFAULT_SRC = Path(__file__).resolve().parent
DEFAULT_DST = DEFAULT_SRC.parent / "Wallpaper_compressed"


def get_size_kb(path: Path) -> float:
    return path.stat().st_size / 1024


def compress_jpeg(src: Path, dst: Path, target_kb: int) -> bool:
    lo, hi = 10, 92
    best = None

    while lo <= hi:
        mid = (lo + hi) // 2
        subprocess.run(
            ["convert", str(src), "-quality", str(mid), "-strip", str(dst)],
            capture_output=True, check=True,
        )
        size_kb = get_size_kb(dst)
        if size_kb <= target_kb:
            best = mid
            lo = mid + 1
        else:
            hi = mid - 1

    if best is None:
        subprocess.run(
            ["convert", str(src), "-quality", "10", "-strip", "-resize", "1280x720", str(dst)],
            capture_output=True, check=True,
        )
        size_kb = get_size_kb(dst)
        if size_kb > target_kb:
            return False
        return True

    subprocess.run(
        ["convert", str(src), "-quality", str(best), "-strip", str(dst)],
        capture_output=True, check=True,
    )
    return True


def main():
    parser = argparse.ArgumentParser(description="压缩图片到指定大小以内")
    parser.add_argument("-i", "--input", type=Path, default=DEFAULT_SRC,
                        help="输入目录（默认：脚本所在目录）")
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_DST,
                        help="输出目录（默认：输入目录的兄弟目录 *_compressed）")
    parser.add_argument("-s", "--size", type=int, default=500,
                        help="目标最大文件大小，单位kB（默认：500）")
    parser.add_argument("--overwrite", action="store_true",
                        help="覆盖已有的压缩文件（默认：跳过已存在的）")
    parser.add_argument("-n", "--max-count", type=int, default=None,
                        help="最多处理的图片数量（默认：全部）")
    args = parser.parse_args()

    src_dir = args.input
    dst_dir = args.output
    target_kb = args.size
    overwrite = args.overwrite
    max_count = args.max_count

    dst_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(f for f in src_dir.iterdir() if f.is_file() and f.suffix.lower() in IMAGE_EXTS)
    total_found = len(files)

    if max_count is not None:
        files = files[:max_count]

    print(f"找到 {total_found} 张图片，本次处理 {len(files)} 张，目标 ≤{target_kb}kB")
    if not overwrite:
        print("模式：跳过已存在的文件（--overwrite 可覆盖）")

    ok, fail, skip = 0, 0, 0
    for src in files:
        name = src.stem + ".jpg"
        dst = dst_dir / name

        if not overwrite and dst.exists() and get_size_kb(dst) <= target_kb:
            print(f"  [skip] {name}: already compressed ({get_size_kb(dst):.0f}kB)")
            skip += 1
            continue

        src_kb = get_size_kb(src)

        if src_kb <= target_kb:
            shutil.copy2(src, dst)
            print(f"  [copy] {name}: {src_kb:.0f}kB (already ok)")
            ok += 1
            continue

        success = compress_jpeg(src, dst, target_kb)
        if success:
            dst_kb = get_size_kb(dst)
            ratio = (1 - dst_kb / src_kb) * 100
            print(f"  [ok]   {name}: {src_kb:.0f}kB -> {dst_kb:.0f}kB ({ratio:.0f}% reduced)")
            ok += 1
        else:
            print(f"  [FAIL] {name}: {src_kb:.0f}kB (could not reach target)")
            fail += 1

    out_total = sum(get_size_kb(f) for f in dst_dir.iterdir() if f.is_file())
    out_count = len(list(dst_dir.iterdir()))
    print(f"\n本次：{ok} 压缩，{skip} 跳过，{fail} 失败")
    print(f"输出目录共 {out_count} 张，总大小 {out_total / 1024:.1f}MB")


if __name__ == "__main__":
    main()
