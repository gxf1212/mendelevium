---
title: "【笔记整理|2023-09+2024年上半年】系统运维与故障排除实用指南"
date: "2024-06-30"
tags: [system-administration, troubleshooting, linux, remote-access, command-line, desktop-environment, operations]
description: "汇总Linux系统运维、远程连接、桌面环境配置及常见故障排除的实用技巧和解决方案，涵盖系统监控、显示器问题、键盘输入等多方面"
thumbnail: "/assets/img/thumbnail_mine/wh-1kdv6v.jpg"
image: "/assets/img/thumbnail_mine/wh-1kdv6v.jpg"
---

# 【笔记整理|2023-09+2024年上半年】系统运维与故障排除实用指南

本文汇总了Linux系统运维、远程连接、桌面环境配置以及常见故障排除的实用技巧和解决方案。

## 系统监控与性能诊断

### 系统兼容性问题识别

#### 软件兼容性检查
如果PyMOL和ChimeraX都有问题，通常是系统级别的问题，需要检查：
- 显卡驱动是否正常
- OpenGL支持是否完整
- 系统库文件是否缺失

#### 键盘输入问题
在某些终端环境下，VMD无法正常响应上下左右键，这通常与gnome terminal的设置有关。

### 显示器相关问题
每次关闭显示器后，dash to panel任务栏会消失，系统默认的会显示，这可能是扩展与电源管理的兼容性问题。

## 远程连接解决方案

### ToDesk使用体验
ToDesk在Linux环境下的特点：
- 无法在Pop!_OS中自动调整布局，但能记住布局设置
- Linux版本不支持复制粘贴功能
- 与Windows版本功能有差异

### AnyDesk配置管理

#### 安装问题解决
[Fedora AnyDesk安装问题](https://discussion.fedoraproject.org/t/cannot-install-anydesk/73854): https://discussion.fedoraproject.org/t/cannot-install-anydesk/73854

#### 自启动管理
[Ubuntu禁用AnyDesk自启动](https://devicetests.com/disable-anydesk-autostart-ubuntu): https://devicetests.com/disable-anydesk-autostart-ubuntu

建议直接禁用自启动功能，按需启动。

## 命令行工具技巧

### 跨平台命令对比

#### Windows PowerShell替代方案
在Windows系统中，没有与Linux系统中的tac命令完全相同的命令。可以使用PowerShell中的Get-Content命令和-Reverse参数来实现类似功能。

#### findstr命令使用
findstr命令类似于Unix系统中的grep，用于在文件中进行文本搜索：
```cmd
findstr "xxx" filename
```

### 文件批量处理

#### sed批量替换
批量文件名处理时，Linux命令更高效：
```bash
# 批量替换文件中的路径
sed -i 's/E:\\GitHub-repo\\notes\\research\\/https\:\/\/cdn.jsdelivr.net\/gh\/username\/notes\@master\/research\//g' *.md

# 批量替换assets路径
sed -i 's/assets\\/assets\//g' *.md
```

#### ZIP压缩操作
[Linux ZIP命令教程](https://www.runoob.com/linux/linux-comm-zip.html): https://www.runoob.com/linux/linux-comm-zip.html

## 桌面环境配置与故障排除

### GNOME扩展管理

#### 扩展兼容性问题
检查GNOME版本兼容性：
```bash
gnome-shell --version
```

某些扩展可能在特定版本的GNOME下存在兼容性问题。

#### Dash to Panel配置
[Dash to Panel扩展](https://extensions.gnome.org/extension/1160/dash-to-panel/): https://extensions.gnome.org/extension/1160/dash-to-panel/

配置注意事项：
- 检查GNOME Shell版本兼容性
- 避免与其他任务栏扩展冲突
- 注意电源管理对扩展的影响

### 工作区管理

#### 动态工作区设置
```bash
# 禁用动态工作区，使用固定数量
gsettings set org.gnome.mutter dynamic-workspaces false
```

建议设置1-4个固定工作区，而不是使用默认的Home设置。

#### 窗口管理优化
[Ubuntu单击任务栏图标最小化窗口](https://cn.linux-console.net/?p=17727): https://cn.linux-console.net/?p=17727

### 多显示器配置
工作区管理在多显示器环境下的注意事项：
- 不是在所有监视器上都显示工作区
- 可以设置主显示器和辅助显示器的不同行为

## Web服务故障排除

### 端口占用问题
```bash
# 检查端口占用情况
sudo apt-get update

# 释放被占用的端口
```

[端口释放指南](https://medium.com/@antonrosh/address-already-in-use-a-simple-guide-to-freeing-up-ports-fbc6a3822983): https://medium.com/@antonrosh/address-already-in-use-a-simple-guide-to-freeing-up-ports-fbc6a3822983

### WebView错误处理

**常见错误**：`Error loading webview: Error: Could not register service workers: TypeError: Failed`

[WebView错误解决方案](https://stackoverflow.com/questions/67698176/error-loading-webview-error-could-not-register-service-workers-typeerror-fai): https://stackoverflow.com/questions/67698176/error-loading-webview-error-could-not-register-service-workers-typeerror-fai

## 网络代理与连接问题

### 代理配置管理
```bash
# 手动设置代理
export http_proxy="http://127.0.0.1:7890"
```

#### CFW代理配置

**使用经验**：
- 现在CFW不影响conda，配置manual proxy即可
- 无法在重启后CFW缓慢启动前连接网络，但手动配置可以工作

### 网络连接故障排除
重启后网络连接问题的解决方案：
1. 检查网络服务状态
2. 验证代理配置
3. 测试DNS解析
4. 检查防火墙设置

## 开发工具集成

### Python外部程序调用
```python
import subprocess

# 调用外部程序的标准方法
def run_external_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout, result.stderr
```

#### 包管理集成
使用subprocess调用系统包管理器：
```python
# 调用antechamber等外部工具
def call_antechamber(input_file, output_file):
    cmd = f"antechamber -i {input_file} -o {output_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result
```

### JSON数据处理
```python
import json

# 加载JSON数据的标准方法
with open('data.json', 'r') as f:
    data = json.load(f)
```

## 系统文档与术语

### 技术术语翻译
- **de facto**：事实上的标准
- **Software Development Kit (SDK)**：软件开发工具包

### 编程概念
#### Arrow Functions
[JavaScript箭头函数](https://developer.mozilla.org/zh-CN/docs/Web/JavaScript/Reference/Functions/Arrow_functions): https://developer.mozilla.org/zh-CN/docs/Web/JavaScript/Reference/Functions/Arrow_functions

## 数据库与版本控制

### Git版本控制扩展
[基于Git版本控制的关系型数据库Dolt](https://jasonkayzk.github.io/2024/01/21/%E5%9F%BA%E4%BA%8EGit%E7%89%88%E6%9C%AC%E6%8E%A7%E5%88%B6%E7%9A%84%E5%85%B3%E7%B3%BB%E5%9E%8B%E6%95%B0%E6%8D%AE%E5%BA%93Dolt/): https://jasonkayzk.github.io/2024/01/21/%E5%9F%BA%E4%BA%8EGit%E7%89%88%E6%9C%AC%E6%8E%A7%E5%88%B6%E7%9A%84%E5%85%B3%E7%B3%BB%E5%9E%8B%E6%95%B0%E6%8D%AE%E5%BA%93Dolt/

这种新型数据库结合了版本控制的优势。

## LaTeX与文档处理

### LaTeX环境配置

#### 基础安装
```bash
# 安装LaTeX基础包
sudo apt install texlive-latex-extra

# 安装XeLaTeX
sudo apt install texlive-xetex

# 安装BibTeX支持
sudo apt install texlive-bibtex-extra
```

[Linux LaTeX安装指南](https://linuxconfig.org/how-to-install-latex-on-ubuntu-20-04-focal-fossa-linux): https://linuxconfig.org/how-to-install-latex-on-ubuntu-20-04-focal-fossa-linux

#### 中文支持
处理"LaTeX Error: File `ctexbook.cls' not found"错误：
这个错误表明缺少CTEX包，该包用于LaTeX中文文档的排版。需要安装相应的中文支持包。

### Markdown到PDF转换
[VSCode Markdown PDF插件](https://github.com/yzane/vscode-markdown-pdf?tab=readme-ov-file#usage): https://github.com/yzane/vscode-markdown-pdf?tab=readme-ov-file#usage

## Docker容器化

### Docker配置问题
[Linux Docker配置](https://blognas.hwb0307.com/linux/docker/654): https://blognas.hwb0307.com/linux/docker/654

容器化部署在开发环境中的重要性日益增加。

## API与Web开发

### GitHub相关服务
[GitHub Discussions快速入门](https://docs.github.com/zh/discussions/quickstart): https://docs.github.com/zh/discussions/quickstart

[GitHub Apps Giscus](https://github.com/apps/giscus): https://github.com/apps/giscus

### Vue.js开发
[Vue.js组合式函数](https://cn.vuejs.org/guide/reusability/composables): https://cn.vuejs.org/guide/reusability/composables

## 故障排除最佳实践

### 系统问题诊断流程
1. **问题重现**：确认问题的可重现性
2. **日志检查**：查看系统和应用程序日志
3. **资源监控**：检查CPU、内存、磁盘使用情况
4. **服务状态**：验证相关服务的运行状态
5. **配置验证**：检查关键配置文件
6. **权限确认**：验证文件和目录权限

### 网络问题排查
1. **连通性测试**：ping、traceroute
2. **端口检查**：netstat、ss命令
3. **DNS解析**：nslookup、dig命令
4. **防火墙状态**：iptables、ufw检查
5. **代理配置**：环境变量和应用配置

### 桌面环境问题解决
1. **重启服务**：重启显示管理器
2. **重置配置**：备份后重置用户配置
3. **扩展管理**：禁用可疑扩展
4. **兼容性检查**：验证软件版本兼容性

## 监控与维护

### 系统健康检查
定期进行系统健康检查：
- 磁盘空间使用情况
- 系统更新状态
- 服务运行状态
- 网络连接质量
- 安全更新应用

### 预防性维护
- 定期清理临时文件
- 更新系统软件包
- 检查硬件健康状态
- 备份重要配置文件
- 监控系统性能指标

---

*本文基于2023年9月至2024年上半年的系统运维实践整理，涵盖常见运维问题的诊断方法和解决方案*