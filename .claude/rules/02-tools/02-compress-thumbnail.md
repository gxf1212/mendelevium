# 图片压缩工具

## 工具位置

`tools/compress.py`

## 功能

将大尺寸图片压缩到指定文件大小以内，保持原始分辨率和格式。适用于生成博客缩略图、壁纸等需要控制文件大小的场景。

## 用法

```bash
# 基本用法：压缩当前目录所有图片到 ../Wallpaper_compressed，目标 ≤500kB
python3 tools/compress.py

# 指定输入输出路径
python3 tools/compress.py -i /mnt/f/Lenovo/Pictures/4K图片1080P/4K美女壁纸图片1080P -o assets/img/4K_1080P_compressed

# 指定目标大小（单位：kB）
python3 tools/compress.py -s 300

# 只处理前 N 张图片
python3 tools/compress.py -n 50

# 覆盖已存在的压缩文件
python3 tools/compress.py --overwrite
```

## 输入输出目录对照

| 输入目录 | 输出目录 | 用途 | 当前数量 |
|----------|----------|------|----------|
| `/mnt/f/Lenovo/Pictures/4K图片1080P/4K美女壁纸图片1080P` | `assets/img/4K_1080P_compressed` | 博客缩略图、封面图 | 250张 |
| `/mnt/f/Lenovo/Pictures/Wallpaper` | `assets/img/Wallpaper_compressed` | 壁纸类图片 | 216张 |
| `assets/img/thumbnail` | `assets/img/thumbnail` | 通用缩略图（bricks.webp等） | 5张 |
| `assets/img/thumbnail_mine` | `assets/img/thumbnail_mine` | 个人收藏缩略图 | 45张 |

## 参数

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| 输入目录 | `-i` / `--input` | 脚本所在目录 | 图片源目录 |
| 输出目录 | `-o` / `--output` | 输入目录的 `*_compressed` 子目录 | 压缩后图片存放位置 |
| 目标大小 | `-s` / `--size` | 500 | 最大文件大小，单位 kB |
| 覆盖 | `--overwrite` | false | 是否覆盖已有的压缩文件 |
| 最大数量 | `-n` / `--max-count` | 全部 | 最多处理的图片数量 |

## 工作模式

- **智能二分搜索**：自动寻找最佳 JPEG 质量参数，在满足文件大小限制的前提下保留最高画质
- **跳过已压缩**：默认只处理输出目录中不存在的文件（使用 `--overwrite` 强制覆盖）
- **自动复制小图**：源文件已小于目标大小的直接复制，不重复压缩

## 注意事项

- 依赖 `ImageMagick`（`convert` 命令）
- 输出统一为 `.jpg` 格式
- 压缩失败（无法达到目标大小）时自动降级到 1280×720 分辨率
- 处理大量图片时建议分批运行（如 `-n 50`）以避免内存压力