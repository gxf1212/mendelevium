---
title: "【笔记整理|2023-09】Linux科研开发环境配置和管理指南"
date: "2023-12-12"
description: "Linux科研开发环境的完整配置和管理指南。涵盖Fedora、Ubuntu系统设置，KDE桌面配置，远程连接、代理设置、VS Code开发环境等核心内容。"
tags: [linux, fedora, ubuntu, kde, development-environment, remote-desktop, proxy, vscode, system-administration]
thumbnail: "/assets/img/thumbnail/sorry.png"
image: "/assets/img/thumbnail/sorry.png"
---

# 【笔记整理|2023-09】Linux科研开发环境配置和管理指南

本文总结了在Linux环境下进行科研开发的实用配置技巧、常见问题解决方案和工具推荐。

## 跨平台文件同步和远程控制

### KDE Connect：跨设备无缝协作

#### 功能特色
KDE Connect是一个强大的跨平台设备协作工具，支持Windows、Linux、macOS、iOS、Android之间的无缝连接：

```bash
# 安装KDE Connect
sudo apt install kdeconnect    # Ubuntu/Debian
sudo dnf install kdeconnect    # Fedora
```

#### 主要功能
KDE Connect虽然不是投屏软件，但功能非常丰富：

1. **文件传输**：
   - 电脑文件右键直接发送至手机
   - 手机图片视频可发送到电脑指定文件夹
   - 无需蓝牙，只要在同一局域网即可

2. **远程控制**：
   - 手机作为电脑遥控器
   - 音乐视频播放控制（音量、进度、暂停等）
   - PPT演示时手机可作为翻页器

3. **通知同步**：
   - 手机电话、短信通知同步到电脑
   - 在电脑上让手机发出声音找手机

4. **剪贴板共享**：
   - 跨设备剪贴板同步
   - 复制粘贴无缝衔接

5. **命令执行**：
   - 预设Linux命令，手机远程执行
   - 支持关机、锁屏、自定义脚本等

### 远程桌面解决方案对比

#### ToDesk
- **官网**：[ToDesk Linux版](https://www.todesk.com/linux.html)：https://www.todesk.com/linux.html
- **优点**：免费，跨平台支持好
- **缺点**：
  - Linux不支持复制粘贴功能
  - 任务栏显示问题（特别是全屏模式）
  - 输入法切换可能有问题

#### AnyDesk
- **安装和配置**：
  ```bash
  # 禁用自启动
  # 参考：[AnyDesk禁用自启动指南](https://devicetests.com/disable-anydesk-autostart-ubuntu)：https://devicetests.com/disable-anydesk-autostart-ubuntu
  
  # 会话管理
  # 参考：[AnyDesk会话管理](https://support.anydesk.com/knowledge/disconnecting-sessions)：https://support.anydesk.com/knowledge/disconnecting-sessions
  ```
- **Fedora安装问题**：[Fedora论坛讨论](https://discussion.fedoraproject.org/t/cannot-install-anydesk/73854)：https://discussion.fedoraproject.org/t/cannot-install-anydesk/73854

## 代理和网络配置

### 代理软件配置

#### electron-ssr配置
```bash
# 启动命令（解决沙盒问题）
/usr/bin/electron-ssr --no-sandbox

# Fedora38环境下使用
# electron-ssr在conda环境中不会报错
```

相关问题讨论：[Electron-SSR GitHub问题](https://github.com/shadowsocksrr/electron-ssr/issues/126)：https://github.com/shadowsocksrr/electron-ssr/issues/126

#### Clash for Windows Linux版
- 配置指南：[Linux Clash配置教程](https://bestoko.cc/p/linux-clash-for-windows/)：https://bestoko.cc/p/linux-clash-for-windows/

#### 其他代理工具
- go-proxy-bingai设置：[GitHub项目](https://github.com/adams549659584/go-proxy-bingai)：https://github.com/adams549659584/go-proxy-bingai

### 网络连接问题诊断

#### Fedora镜像源问题
```bash
# 常见错误：无法连接到Fedora镜像源
Failed to search for file: cannot update repo 'fedora': 
Cannot prepare internal mirrorlist: Curl error (7): 
Couldn't connect to server for https://mirrors.fedoraproject.org/metalink?repo=fedora-38&arch=x86_64
[Failed to connect to 127.0.0.1 port 12333 after 0 ms: Couldn't connect to server]
```

**解决方案**：
1. 检查代理设置是否正确
2. 尝试更换镜像源
3. 检查防火墙和网络配置

## 开发工具配置

### Visual Studio Code

#### 扩展开发
```bash
# VSCode扩展路径
/home/user/.vscode/extensions/md-highlighter-0.0.1

# 发布token配置
vscode token: your_token_here
```

#### 已知问题
- **虚拟桌面恢复**：VSCode或Firefox无法在Fedora KDE中将窗口恢复到正确的虚拟桌面，这是已知问题
- **调试配置**：缺少.vscode文件夹可能导致调试扩展无法识别

相关讨论：[VSCode VSCE GitHub问题](https://github.com/microsoft/vscode-vsce/issues/419)：https://github.com/microsoft/vscode-vsce/issues/419

### Linux原生应用

#### 微信支持
```
现在优麒麟下有Linux原生的微信，虽然功能简陋了一些，
但是有比没有强，基本的聊天需求是可以被满足的。
```

#### 文件权限管理
```bash
# 设置可执行权限
chmod +x /path/to/Multiwfn_3.8_dev_bin_Linux/Multiwfn

# 批量权限设置
find . -type f -exec chmod a+x {} \;
```

## 系统优化和故障排除

### 桌面环境配置

#### KDE Plasma优化
- **启动速度**：Plasma启动需要25-40秒（可能与NVIDIA显卡有关）
- **应用启动器**：左下角的"f"图标（Plasma application launchers）
- **窗口恢复**：重启后只有Firefox能够恢复窗口状态

#### 虚拟桌面管理
目前VSCode和Firefox在Fedora KDE中无法正确恢复虚拟桌面窗口位置。

### 编译工具链配置

#### Devtoolset（CentOS/RHEL）
```
Devtoolset是一个用于在Red Hat Enterprise Linux (RHEL)和CentOS系统上
安装和使用多个版本的编译器和开发工具的软件集合。
它提供了更新的编译器版本，以便开发人员可以使用最新的功能和优化。
```

#### 包管理问题
```bash
# Conda包损坏问题
InvalidArchiveError("Error with archive /home/user/anaconda3/pkgs/gxx_impl_linux-64-10.4.0-h7ee1905_16.tar.bz2. 
You probably need to delete and re-download or re-create this file.")

# 解决方案：清理并重新下载
conda clean -a
```

## 端口和服务管理

### 端口占用检查
```bash
# 查看端口12333的使用情况
sudo lsof -i :12333

# 识别占用端口的进程
lsof -i :port_number
```

### 系统服务管理
```bash
# 检查系统版本
gnome-shell --version

# 网络服务诊断
ping -c 4 mirrors.fedoraproject.org
```

## 文档和教程资源

### Bash编程
- [Bash序列表达式](https://linuxize.com/post/bash-sequence-expression/)：https://linuxize.com/post/bash-sequence-expression/
- [Python字典排序](https://www.golinuxcloud.com/python-sort-dictionary-by-key/)：https://www.golinuxcloud.com/python-sort-dictionary-by-key/

### 网络分析工具
```python
import networkx as nx
```
NetworkX文档：[NetworkX算法文档](https://networkx.org/documentation/stable/reference/algorithms/traversal.html)：https://networkx.org/documentation/stable/reference/algorithms/traversal.html

### Ubuntu系统资源
- [Amber22安装指南](http://archive.ambermd.org/202302/att-0090/Amber_22_and_Tools_22_install_Ubuntu_22.pdf)：http://archive.ambermd.org/202302/att-0090/Amber_22_and_Tools_22_install_Ubuntu_22.pdf

## 性能优化建议

### GPU计算支持
配置GPU支持的计算环境：
```
Quick package for Hartree-Fock and DFT electronic stucture calculations, with GPU support. 
Quick is integrated into sander for QM/MM simulations, and AmberTools23 contains significant performance improvements, 
a new geometry optimizer, and support for spin-unrestricted calculations.
```

### 跨平台兼容性
注意Linux和Windows之间文件格式的兼容性：
```
是因为在Linux里面读Windows的chk？
```
某些二进制文件在不同操作系统间可能存在兼容性问题。

## 故障排除检查清单

### 网络连接问题
1. **检查代理设置**
   - electron-ssr是否正常运行
   - 端口12333是否被占用
   - 防火墙设置是否正确

2. **包管理问题**
   - conda缓存是否损坏
   - 镜像源是否可访问
   - 网络连接是否稳定

### 桌面环境问题
1. **显示相关**
   - NVIDIA驱动是否正确安装
   - Plasma启动时间是否异常
   - 虚拟桌面功能是否正常

2. **应用兼容性**
   - VSCode扩展是否正确安装
   - .vscode配置文件夹是否存在
   - 权限设置是否正确

### 开发环境问题
1. **编译工具**
   - GCC版本是否兼容
   - 开发库是否完整安装
   - 环境变量是否正确设置

2. **Python环境**
   - conda环境是否激活
   - 包依赖是否满足
   - 路径配置是否正确

## 推荐的Linux发行版选择

### 科研用途推荐
1. **Ubuntu LTS**：稳定性好，社区支持强
2. **Fedora**：新技术支持好，适合开发
3. **优麒麟**：中文支持好，有原生微信

### 桌面环境选择
1. **KDE Plasma**：功能丰富，可定制性强
2. **GNOME**：简洁美观，资源占用相对较低

---

*本文基于2023年9-12月技术讨论记录整理，涵盖Linux环境下科研开发的实际经验和解决方案*