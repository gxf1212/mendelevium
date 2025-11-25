---
title: "Ubuntu Virtual Memory (Swap) Setup Tutorial: Enhance System Performance"
date: "2025-10-08"
tags: [ubuntu, swap, virtual-memory, system-administration, memory-management, server-management]
description: "详细Ubuntu系统虚拟内存(Swap)配置教程，涵盖交换空间检查、文件创建、启用挂载等完整步骤，为内存不足服务器提供性能优化方案"
thumbnail: "/assets/img/thumbnail_mine/wh-g78rle.jpg"
image: "/assets/img/thumbnail_mine/wh-g78rle.jpg"
---
# 在 Ubuntu 中增加虚拟内存（Swap）教程

在 Ubuntu 系统中增加虚拟内存（即交换空间，Swap）可以有效提升系统在内存不足时的性能。以下是详细的操作步骤：

## 一、检查当前交换空间

首先，您需要检查当前系统的交换空间情况。打开终端并运行以下命令：

```bash
sudo swapon --show
```

如果命令没有输出，说明当前系统没有启用交换空间。如果有输出，则会显示现有交换文件或分区的信息（例如 `/swapfile`）。

## 二、创建新的交换文件

### 方法一：使用 `fallocate` 命令（推荐）

运行以下命令创建一个新的交换文件：

```bash
sudo fallocate -l 4G /swapfile_new
```

- `-l 4G`：指定交换文件大小为 4GB。您可以根据需求调整大小，例如使用 `8G` 表示 8GB。
- `/swapfile_new`：新交换文件的路径。您可以自定义文件名，但需确保后续步骤中路径一致。

### 方法二：使用 `dd` 命令（若 `fallocate` 不可用）

如果 `fallocate` 命令不可用，可以使用 `dd` 命令创建交换文件：

```bash
sudo dd if=/dev/zero of=/swapfile_new bs=1G count=4
```

- `bs=1G`：每次写入 1GB 数据。
- `count=4`：写入 4 次，生成 4GB 文件。

## 三、设置交换文件的权限

为了安全起见，设置交换文件的权限，使其仅限 root 用户访问：

```bash
sudo chmod 600 /swapfile_new
```

## 四、格式化交换文件

将创建的文件标记为交换空间：

```bash
sudo mkswap /swapfile_new
```

## 五、启用交换文件

运行以下命令启用新创建的交换文件：

```bash
sudo swapon /swapfile_new
```

## 六、验证交换空间

检查新增的交换空间是否生效：

```bash
sudo swapon --show
```

您还可以查看内存使用情况以确认交换空间的变化：

```bash
free -h
```

## 七、配置开机自动挂载

为了使交换文件在系统重启后仍然有效，需要将其添加到 `/etc/fstab` 文件中：

1. 打开 `/etc/fstab` 文件进行编辑：

    ```bash
    sudo nano /etc/fstab
    ```

2. 在文件末尾添加以下内容：

    ```
    /swapfile_new none swap sw 0 0
    ```

3. 保存并退出编辑器（在 `nano` 中，按 `Ctrl+O` 保存，按 `Ctrl+X` 退出）。

## 注意事项

- **调整交换文件大小**：根据系统需求和使用场景调整交换文件的大小。一般建议交换文件大小为物理内存的 1-2 倍，但具体大小取决于您的应用场景。
- **权限管理**：确保交换文件的权限设置正确，避免非授权访问。
- **性能考量**：虽然增加交换空间可以缓解内存不足的问题，但过度依赖交换空间可能会降低系统性能，因为磁盘 I/O 速度远低于内存。

通过以上步骤，您可以成功增加 Ubuntu 系统的虚拟内存（Swap），从而提升系统的整体性能和稳定性。

希望这份教程对您有所帮助！如果您在操作过程中遇到任何问题，欢迎随时提问。





# Pandoc 生成 PDF 时字体问题解决方案教程

## 一、问题概述

在使用 Pandoc 将 Markdown 文件生成 PDF 时，如果指定使用 Times New Roman 字体，可能会遇到错误。这是因为 Times New Roman 是 Windows 系统的默认字体，在 Linux 或 macOS 上默认未安装。此外，对于中文支持，也需要确保系统中存在相应的中文字体。

## 二、检查字体是否安装

### 在 Linux 系统中

打开终端，运行以下命令查看系统中已安装的字体：

```bash
fc-list :lang=zh  # 查看中文字体
fc-list | grep "Times New Roman"  # 查找 Times New Roman 字体
```

如果没有输出，说明系统中未安装该字体。

### 在 macOS 系统中

使用 Font Book 应用程序检查字体是否安装。

### 在 Windows 系统中

打开“字体”文件夹（通常在 C:\Windows\Fonts），查找“Times New Roman”字体。

## 三、安装所需字体

### 安装 Times New Roman 字体

#### 对于 Ubuntu/Debian 系统：

运行以下命令安装 Microsoft 核心字体，其中包含 Times New Roman：

```bash
sudo apt-get update
sudo apt-get install ttf-mscorefonts-installer
```

在安装过程中，可能需要接受许可协议。安装完成后，运行以下命令刷新字体缓存：

```bash
sudo fc-cache -fv
```

#### 对于 CentOS/RHEL 系统：

使用以下命令安装字体：

```bash
sudo yum install curl curl-devel
sudo rpm -Uvh http://li.nux.ro/download/fedora/epel/5/i386/epel-release-5-4.noarch.rpm
sudo yum install ttf-mscorefonts-installer
```

#### 对于 macOS 系统：

从官方渠道下载并安装 Microsoft Office for Mac，它会附带安装 Times New Roman 字体。或者，您可以手动下载字体文件并安装。

### 安装中文支持字体

如果您需要在 PDF 中显示中文，还需要安装中文字体。例如，在 Ubuntu/Debian 系统上，可以安装 `texlive-lang-chinese` 包：

```bash
sudo apt install texlive-lang-chinese
```

该包包含中文支持的宏包（如 ctex），是 Debian 官方维护的包，具有良好的兼容性。

## 四、配置 Pandoc 使用正确字体

在 Pandoc 命令中指定字体时，确保使用的字体名称与系统中实际存在的字体名称完全匹配。例如：

```bash
pandoc input.md -o output.pdf --pdf-engine=xelatex --css style.css -V mainfont="Times New Roman" -V CJKmainfont="AR PL UMing CN"
```

- `mainfont`：指定西文字体。
- `CJKmainfont`：指定中文字体。

## 五、生成 PDF 的 Python 函数示例

以下是一个使用 Pandoc 生成 PDF 的 Python 函数示例，确保路径和字体名称正确：

```python
import subprocess
import logging
from pathlib import Path

log = logging.getLogger(__name__)

def generate_pdf_with_pandoc(md_path: Path, css_path: Path, output_pdf_path: Path) -> bool:
    """
    使用 Pandoc 和 XeLaTeX 生成 PDF 文件。
    
    参数:
        md_path: 输入的 Markdown 文件路径。
        css_path: CSS 文件路径（可选）。
        output_pdf_path: 输出的 PDF 文件路径。
    
    返回:
        PDF 生成成功返回 True，失败返回 False。
    """
    log.info(f"Attempting PDF generation with Pandoc for {md_path}.")
    pandoc_cmd = [
        'pandoc',
        str(md_path),
        '-o',
        str(output_pdf_path),
        '--pdf-engine=xelatex',
        '--css',
        str(css_path),
        '-V',
        'mainfont=Times New Roman',
        '-V',
        'CJKmainfont=AR PL UMing CN'
    ]
    result = subprocess.run(pandoc_cmd, capture_output=True, text=True, encoding='utf-8')
    if result.returncode != 0:
        log.error(f"Pandoc failed. Stderr: {result.stderr}")
        return False
    log.info(f"Successfully generated PDF with Pandoc at {output_pdf_path}")
    return True
```

## 六、验证和测试

1. **验证字体安装**：
   - 运行 `fc-list` 命令，检查是否列出了 Times New Roman 和中文字体。
   - 确保字体名称与 Pandoc 命令中指定的名称完全一致。

2. **测试 PDF 生成**：
   - 使用上述 Python 函数或直接运行 Pandoc 命令生成 PDF。
   - 打开生成的 PDF 文件，检查字体显示是否正确。

## 七、总结

通过以上步骤，您可以解决 Pandoc 在生成 PDF 时找不到指定字体的问题。确保系统中安装了所需的字体，并在 Pandoc 命令中正确指定字体名称。对于中文支持，安装 `texlive-lang-chinese` 包是一个推荐的解决方案。希望这份教程能帮助您顺利完成 PDF 生成任务。

如果您在操作过程中遇到任何问题或需要进一步的帮助，欢迎随时提问。