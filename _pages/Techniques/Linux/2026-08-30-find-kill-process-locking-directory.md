---
title: "用 Sysinternals Handle 查杀占用目录的进程"
date: "2026-08-30"
last_modified_at: "2026-08-30"
tags: [windows, sysinternals, process-handle, file-lock, troubleshooting, command-line]
description: "介绍用微软 Sysinternals 的 Handle 与 Process Explorer 查出占用目录的进程，并通过 taskkill 或关闭句柄解除占用，便于重命名或移动被锁定的文件夹"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/4K_1080P_compressed/081025rf9KH.jpg"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/4K_1080P_compressed/081025rf9KH.jpg"
author: Xufan Gao
lang: zh-CN
---

# 用 Sysinternals Handle 查杀占用目录的进程

最方便的是用微软官方 Sysinternals 的 `Handle` 或 `Process Explorer` 查哪个进程占用了这个目录。`Handle` 就是专门查“哪个进程打开了某个文件或目录”的。

如果你喜欢命令行，先以管理员身份打开 PowerShell/CMD，然后：

```bat
handle "E:\graduate_study\research\NanoMedicine\OP_unit"
```

它会输出类似：

```text
Code.exe pid: 12345 ...
explorer.exe pid: 6789 ...
python.exe pid: 24680 ...
```

然后直接杀进程：

```bat
taskkill /PID 12345 /F
```

如果是多个进程，就分别杀掉。之后再重命名文件夹即可。

如果你还没安装 `handle.exe`，可以从[微软官方 Sysinternals 下载](https://learn.microsoft.com/en-us/sysinternals/downloads/handle)并添加到环境变量。它要求以管理员权限运行。

如果查出来是 `explorer.exe`，通常不用杀死整个桌面，直接：

```bat
taskkill /f /im explorer.exe
start explorer.exe
```

或者任务管理器里右键“Windows 资源管理器”→“重新启动”。

还有一个更暴力但不太建议的办法：

```bat
handle "OP_unit"
```

得到具体 handle 编号和 PID 后，可以直接关闭那个句柄：

```bat
handle -c <HANDLE_ID> -p <PID>
```

比如：

```bat
handle -c 3F4 -p 12345
```

但微软明确警告，**强行关闭 handle 可能导致程序不稳定或数据损坏**，所以一般优先 `taskkill` 整个占用进程（参考 [Process Explorer 官方下载](https://learn.microsoft.com/en-us/sysinternals/downloads/process-explorer)）。
