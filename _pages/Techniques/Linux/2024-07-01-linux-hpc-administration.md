---
title: "【笔记整理|2024-07】Linux系统管理与HPC集群运维：从基础命令到SLURM作业调度"
date: "2024-07-01"
tags: [linux, hpc, slurm, system-administration, ssh, cloud-computing, technical-notes]
description: "详细整理Linux系统管理和HPC集群运维的关键知识，涵盖用户组管理、SSH配置、SLURM作业调度、性能监控等技术要点，为计算科学工作者提供实用运维指南"
thumbnail: "/assets/img/thumbnail_mine/wh-6o2m9q.jpg"
image: "/assets/img/thumbnail_mine/wh-6o2m9q.jpg"
---

# 【笔记整理|2024-07】Linux系统管理与HPC集群运维：从基础命令到SLURM作业调度

## 引言

Linux系统管理和HPC集群运维是计算科学研究的基石。无论是本地工作站还是大型计算集群，掌握Linux系统管理技能都是必不可少的。本文整理了从技术讨论中提取的Linux系统管理和HPC集群运维的关键知识和实用技巧，涵盖从基础命令到高级作业调度的各个方面。

## Linux基础命令与系统管理

### 系统信息查看

了解系统基本信息是系统管理的第一步。

**有趣的知识：** usr代表Unix System Resources，而不是user！

### 用户与组管理

Linux系统中的用户和组管理是多用户环境下的基础操作。

要在Linux系统中查看用户组，可以使用以下命令。usermod命令是一个用于修改用户属性的强大工具，其中包括将用户添加到现有用户组的功能。

**用户组管理的重要性：** 操作系统具有拥有完全权限的用户。然而，由于该用户不能与登录到系统的人员共享，因此他们临时与其他用户共享部分权限。

### SSH密钥管理

SSH密钥是远程管理和自动化任务的核心。

执行ssh-keygen命令生成密钥对。我们为每个人只存储一个SSH公钥。公钥可以与世界上的任何人共享（因此称为公钥）。只有您应该访问您的私钥。

### 虚拟内存管理

Linux系统的虚拟内存管理对于保证大规模计算任务的稳定运行至关重要。

在Linux中，当物理内存被耗尽时，会使用swap的虚拟内存（较慢）。当物理内存和虚拟内存都耗尽时就会出现程序跑不起来、启动这个进程会杀死另外一个进程的情况，以保证程序的良好运行。

### 包管理

不同的Linux发行版使用不同的包管理系统。面对如此多样的指令集结构，软件开发者想要为每一种架构都编译一份软件包十分困难。因此，在Linux生态中，源代码是最通用的软件分发形式。

**Zlib包安装问题处理：**

zlib的官网打不开，apt-get install zlib也找不到软件包，貌似不在软件源里。解决方法是打开Ubuntu Software Center，搜索zlib，找到zlib1g-dev这个包，安装成功。

使用APT安装Zlib：
```bash
sudo apt install zlib1g
# 如果需要开发文件（头文件和静态库）
sudo apt install zlib1g-dev
```

### 模块管理系统

在HPC环境中，模块管理系统是软件环境配置的关键。

```bash
module avail  # 显示可以使用的模块
```

## SLURM作业调度系统

### 作业提交与资源管理

SLURM是最常用的HPC作业调度系统之一，合理配置作业参数可以显著提高计算效率。

```bash
#SBATCH --exclude=node4,node5,node7,node8,node9
```

**节点选择策略：** --nodelist只能指定一个节点，但#SBATCH --exclude=node[1-16]这种范围表示法是可行的。

### 作业依赖与流程管理

复杂的计算流程通常需要作业之间的依赖关系管理：

**SLURM依赖作业提交指南：**
https://bioinformaticsworkbook.org/Appendix/HPC/SLURM/submitting-dependency-jobs-using-slurm.html#gsc.tab=0

### 作业状态监控

实时监控作业状态是集群管理的重要功能：

```bash
sacct --starttime=2024-06-29 --format=JobID%10,User%20,Partition,Submit,Start,Elapsed,AllocTRES%50 -X
```

### 作业控制

作业的暂停、恢复和取消是日常管理操作：

```bash
scontrol suspend jobid
```

### 用户账户管理

在SLURM集群中管理用户账户是系统管理员的职责：

```bash
sacctmgr add user User=${u} Account=urgent
```

## 云计算与远程服务

### AWS EC2使用

AWS EC2是常用的云计算平台，掌握基本操作非常重要，包括文件上传和下载等操作。

### 环境变量配置

合理配置环境变量可以简化日常操作：

```bash
export TZ='Asia/Shanghai'
```

## 文件系统与数据管理

### 文件压缩与解压

数据压缩和归档是数据管理的必备技能。

要清理pip的缓存，可以使用以下命令：

```bash
pip cache purge  # 清除所有缓存
pip cache remove <package_name>  # 清除特定包的缓存
pip cache dir  # 查看缓存路径
```

**参考：** [Zip文件操作指南](https://docs.hostdime.com/hd/command-line/how-to-tar-untar-and-zip-files)

### 文件搜索与过滤

高效的文件搜索和过滤可以大大提高工作效率。要查找恰好包含两个连字符的目录名，需要将grep模式"锚定"以匹配整行。排除特定文件可以使用-X选项。

### Git版本控制

Git是现代科研项目的标准版本控制工具。合理配置.gitignore规则可以避免提交不必要的文件。

## 编译与开发环境

### 编译系统理解

理解编译系统的工作原理有助于解决编译问题。

gcc的编译其实是四个过程的集合，分别是预处理（preprocessing）、编译（compilation）、汇编（assembly）、链接（linking），分别由cpp、cc1、as、ld这四个程序完成，gcc是它们的封装。

### C++编程技巧

掌握C++编程技巧可以提高开发效率。

在C++中，字符"*"是一个指针，包含变量的值。++i有时可以比i++更快，并且永远不会更慢。对于基本数据类型，编译器很可能会修复并优化掉任何不必要的复制。对于迭代器这更困难，对于用户定义类型可能完全不可能。

### Makefile编写

Makefile是自动化编译的重要工具，可以将多个C++源文件分别编译成不同的可执行文件。

### LaTeX排版系统

LaTeX是科学文档排版的标准工具。

可以使用apt命令安装LaTeX：

```bash
sudo apt install texlive-latex-extra
sudo apt install texlive-xetex  # XeLaTeX
sudo apt install texlive-bibtex-extra  # BibTeX支持
```

**中文字体支持问题：** 错误"LaTeX Error: File `ctexbook.cls' not found"表示缺少CTEX包，该包是LaTeX中用于排版中文文档的文档类文件。

**参考：** [LaTeX安装指南](https://linuxconfig.org/how-to-install-latex-on-ubuntu-20-04-focal-fossa-linux)

## 系统诊断与性能优化

### 系统监控工具

系统监控是保证服务稳定运行的关键。

**参考：** [VS Code缓存清理](https://wz.anoms.top/2025/02/23/vscode-cache-diskspace-clean/)

### 软件安装问题解决

解决软件安装过程中的常见问题，如"No rule to make target 'X'"通常表示文件缺失。

## 云原生与容器技术

### 虚拟化技术

虚拟化技术是现代云计算的基础。

Hypervisor（也称为虚拟机监视器或VMM）是创建和运行虚拟机（VM）的软件。

**虚拟化类型：**

- **Type 1 hypervisor：** 直接在主机硬件上运行以控制硬件并管理客户操作系统。例如VMware ESXi、Microsoft Hyper-V和Xen。

### Linux发行版选择

选择合适的Linux发行版对于特定应用场景很重要。

netinst版本是一个小型ISO镜像，仅包含启动安装所需的文件。DVD-1版本是一个大型ISO镜像，包含桌面环境、应用程序和其他软件。

## 总结与最佳实践

1. **基础命令**：掌握Linux基础命令是系统管理的基础，理解命令的内部工作原理有助于问题排查
2. **用户管理**：合理配置用户和组权限，确保系统的安全性和可管理性
3. **SSH密钥**：妥善管理SSH密钥，建立安全的远程访问机制
4. **虚拟内存**：合理配置swap空间，避免因内存不足导致的程序异常
5. **SLURM调度**：熟练掌握SLURM作业调度系统，优化计算资源使用
6. **版本控制**：建立良好的Git使用习惯，确保研究过程的可追溯性
7. **编译环境**：理解编译原理，能够独立解决编译和链接问题
8. **监控诊断**：建立系统监控体系，及时发现和解决潜在问题

通过这些系统管理和集群运维技能的掌握，可以为计算科学研究提供稳定、高效的计算环境支持。

## 参考资源

- [SLURM依赖作业提交指南](https://bioinformaticsworkbook.org/Appendix/HPC/SLURM/submitting-dependency-jobs-using-slurm.html#gsc.tab=0)
- [文件压缩操作指南](https://docs.hostdime.com/hd/command-line/how-to-tar-untar-and-zip-files)
- [Linux系统监控指南](https://wz.anoms.top/2025/02/23/vscode-cache-diskspace-clean/)
- [SLURM环境变量文档](https://nscc.mrzhenggang.com/faqs/slurm-built-in-environment-variables/)
- [LaTeX在Ubuntu上安装指南](https://linuxconfig.org/how-to-install-latex-on-ubuntu-20-04-focal-fossa-linux)