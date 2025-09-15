---
title: "【笔记整理|2025-07】系统管理与开发工具实用指南"
date: "2025-09-12"
tags: [linux, system-admin, development-tools, vscode, git, ssh, claude-code, latex, python, conda, docker]
---

# 【笔记整理|2025-07】系统管理与开发工具实用指南

基于2025年7月以来的技术讨论，整理系统管理和开发环境配置的实用技巧和解决方案。

## SSH密钥管理

### Linux/macOS SSH密钥操作
```bash
# 生成SSH密钥
ssh-keygen

# 查看公钥内容
cat ~/.ssh/id_rsa.pub
```

### Windows SSH密钥操作
```powershell
# 方法1：使用type命令
type %USERPROFILE%\.ssh\id_rsa.pub

# 方法2：使用PowerShell
Get-Content $env:USERPROFILE\.ssh\id_rsa.pub
```

## Git版本控制配置

### 基础配置
```bash
# 全局用户配置
git config --global user.name "username"
git config --global user.email "email@example.com"

# 查看当前配置
git config --list
```

### 输出重定向技巧
```bash
# 正确的stdout和stderr重定向方式
command > logfile 2>&1

# 简化语法
command &> logfile
```

## 开发环境配置

### LaTeX完整安装（Ubuntu）
```bash
# 基础LaTeX包
sudo apt install texlive-latex-extra

# XeLaTeX支持（支持中文排版）
sudo apt install texlive-xetex

# BibTeX支持
sudo apt install texlive-bibtex-extra

# 中文支持包
sudo apt install texlive-lang-chinese

# 推荐：完整安装
sudo apt install texlive-full
```

**中文文档错误解决**：
如果遇到`LaTeX Error: File 'ctexbook.cls' not found`错误，说明缺少CTEX包，需要安装中文支持包。

### Python环境管理

#### Conda缓存清理
```bash
# 查看pip缓存路径
pip cache dir

# 清理pip缓存
pip cache purge

# 查看conda环境路径
conda info --envs

# 清理Python缓存目录
find . -name "__pycache__" -type d -exec rm -rf {} +
```

#### 开发模式安装
```bash
# 以开发模式安装包
pip install -e .[dev]
```

## VSCode开发配置

### 基础搜索功能
- **快捷键**：`Ctrl+Shift+F`
- **菜单路径**：查找 -> 查找
- **支持功能**：正则表达式搜索、跨文件搜索

### 缓存管理（Windows）
```cmd
# 将VSCode缓存迁移到其他盘符
mklink /d C:\Users\Lenovo\.vscode F:\cache\.vscode
```

### 代码导航快捷键
- `Ctrl+B`：跳转到声明或用法
- `Ctrl+Click`：点击跳转到定义

## Claude Code AI工具配置

### 安装和配置流程
```bash
# 1. 安装Claude Code CLI
sudo npm install -g @anthropic-ai/claude-code
```

### NonoCode配置
```bash
function nonocode {
    export ANTHROPIC_BASE_URL="http://claude.nonocode.cn/api/"
    export ANTHROPIC_API_KEY="your_api_key"
    export ANTHROPIC_AUTH_TOKEN="your_token"
}
```

注册地址：[NonoCode注册](http://claude.nonocode.cn/register?ref=3S5UMR)：http://claude.nonocode.cn/register?ref=3S5UMR

## 系统管理工具

### 用户管理
```bash
# 查看用户详细信息
getent passwd username

# 删除用户及其相关文件
userdel -r username

# 查看系统版本
gnome-shell --version
```

### 文件传输和压缩

#### 保持目录结构的文件传输
```bash
# 使用tar+ssh保持目录结构上传
tar czf - mapping.png | ssh user@remote "cd /remote/path && tar xzf -"
```

#### 压缩和解压操作
```bash
# 创建zip压缩文件
zip archive.zip folder/ -r

# 解压tar到指定目录并去除顶级目录
tar -xf archive.tar -C target-dir --strip-components=1

# 查看压缩文件内容（不解压）
tar -tf archive.tar
```

### 文本处理技巧

#### sed批量处理
```bash
# 批量修复markdown转义字符
sed -i 's/\\\\/\\/g' *.md    # 修复反斜杠
sed -i 's/\\_/_/g' *.md      # 修复下划线
sed -i 's/\\\*/\*/g' *.md    # 修复星号
```

### Python缓存清理
```bash
# 递归删除所有__pycache__目录
find . -name "__pycache__" -type d -exec rm -rf {} +
```

## GNOME桌面环境

### 扩展程序修复
修复Dash to Panel扩展：[Dash to Panel](https://extensions.gnome.org/extension/1160/dash-to-panel/)：https://extensions.gnome.org/extension/1160/dash-to-panel/

## pKa计算工具推荐

- **PypKa**：[GitHub页面](https://github.com/mms-fcul/PypKa)：https://github.com/mms-fcul/PypKa（开源推荐）
- **DelphiPka**：[GitHub页面](https://github.com/delphi001/DelphiPka)：https://github.com/delphi001/DelphiPka（较老版本）
- **Rowan pKa Calculator**：[Rowan在线工具](https://rowansci.com/tools/pka)：https://rowansci.com/tools/pka（商业工具）

## 高级开发技巧

### 快速目录搜索
```bash
# 使用grep进行目录名匹配（恰好包含两个连字符的示例）
grep "^.*-.*-.*$" directory_list
```

### PyCharm使用注意事项
PyCharm作为IDE本身不能直接渲染localhost:8501等网页内容，需要在外部浏览器中打开Web服务。

### 文档服务本地部署
```bash
# 本地启动Sphinx文档服务
sphinx-build -b html source build
python -m http.server 8000 -d build
```

### 开发模式安装技巧
```bash
# Python包开发模式安装
pip install -e .[dev]
```

## 系统磁盘扩展

参考京东云文档：[扩展文件系统无分区](https://docs.jdcloud.com/cn/iavm/expand-filesystem-no-partition)：https://docs.jdcloud.com/cn/iavm/expand-filesystem-no-partition

---

*本文基于2025年7-9月技术讨论记录整理，涵盖系统管理和开发环境配置的实际操作经验*