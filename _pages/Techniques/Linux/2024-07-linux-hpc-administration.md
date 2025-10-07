---
title: "【笔记整理|2024-07】Linux系统管理与HPC集群运维：从基础命令到SLURM作业调度"
date: "2024-07-01"
tags: [linux, hpc, slurm, system-administration, ssh, cloud-computing, technical-notes]
---

# 【笔记整理|2024-07】Linux系统管理与HPC集群运维：从基础命令到SLURM作业调度

## 引言

Linux系统管理和HPC集群运维是计算科学研究的基石。无论是本地工作站还是大型计算集群，掌握Linux系统管理技能都是必不可少的。本文整理了从技术讨论中提取的Linux系统管理和HPC集群运维的关键知识和实用技巧，涵盖从基础命令到高级作业调度的各个方面。

## Linux基础命令与系统管理

### 系统信息查看

了解系统基本信息是系统管理的第一步：

> Fun fact: usr stands for Unix System Resources, not user!

### 用户与组管理

Linux系统中的用户和组管理是多用户环境下的基础操作：

> 要在Linux系统中查看用户组，可以使用以下命令：

> usermod 命令是一个用于修改用户属性的强大工具，其中包括将用户添加到现有用户组的功能。下面是使用 usermod 命令将用户添加到现有用户组的示例操作：

**用户组管理的重要性：**
> Operating systems have a user with full privileges. However, since this user cannot be shared with the people logged into that system, they temporarily share some of their privileges with other users.

### SSH密钥管理

SSH密钥是远程管理和自动化任务的核心：

> Execute the ssh-keygen command

> We only store one SSH public key per person. The public key can be shared with anyone in the world (thus the name public). Only you should have access to your private key.

**SSH密钥格式示例：**
> ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABgQC8qLBAyfzgREFH4hMXwCkCTqu2/6p9JfFQmlW6FbD91iaG9EXnvuSxuZ7T+AlsmldFCAH9LaByvdXFZOkuwas0IcNGeb1HomSPLnS73JI8NBgIL/wbfYydGYesXXi9e13BL+/a3m35IajMBpra1K9tEvOOOovx4HRWzDhDWdqQQRKeJr+KHvVgJSr9MPU16kdJUDrlBlRT99F0hOF8DApNmoNI922wGXYkZDK171Qez0YAXIqtNZO4cy3kdjICsr6EbVWtzXAqVsusDxyDMvnpiEQ1cuVv5P5syqz7E56xVoxdCaXIMFa9LxR+VodkCrTzCx/ucqWPZkqSallJ7feyRxkFNW+OnH1qpEYsRSFjKC+ZW1zC556g/327GP7vW7K+6yZe4ReEt4OiW36Empb7/jML7X14nsKceTECd42J7qBANq1Pb/3Kqz3LQOjvZhgcgXpQ8L/MBiAZpQ1n8jbtl2muqfNKfMqFMn6x1EoQ2pcwmJuKjS2udmkl5PkycLM= lenovo@Xufan-Legion-y9000p

### 虚拟内存管理

Linux系统的虚拟内存管理对于保证大规模计算任务的稳定运行至关重要：

> Linux设置虚拟内存

> 在linux中，当mem物理内存被耗尽时，会使用swap的虚拟内存（较慢），当物理内存和虚拟内存都耗尽时就会出现程序跑不起来，启动这个进程会杀死另外一个进程的情况，已保证程序的良好运行

**创建swap空间示例：**
> Setting up swapspace version 1, size = 32 GiB (34359734272 bytes)

### 包管理

不同的Linux发行版使用不同的包管理系统：

> 面对如此多样的指令集结构，软件开发者想要为每一种架构都编译一份软件包十分困难。因此，在 Linux 生态中，源代码是最通用的软件分发形式。

**Zlib包安装问题处理：**
> zlib的官网打不开，apt-get install zlib也找不到软件包，貌似不在软件源里？解决方法是打开ubuntu software center，搜索zlib，找到zlib1g-dev这个包，安装成功。

> To install Zlib using APT, follow these steps: sudo apt install zlib1g If you require the development files for Zlib, which include header files and static libraries necessary for compiling programs that use Zlib, install the zlib1g-dev package as well: sudo apt install zlib1g-dev

### 模块管理系统

在HPC环境中，模块管理系统是软件环境配置的关键：

> module avail显示可以使用的模块

## SLURM作业调度系统

### 作业提交与资源管理

SLURM是最常用的HPC作业调度系统之一，合理配置作业参数可以显著提高计算效率：

> #SBATCH --exclude=node4,node5,node7,node8,node9

**节点选择策略：**
> we can only specify one for --nodelist, but #SBATCH --exclude=node[1-16] works

### 作业依赖与流程管理

复杂的计算流程通常需要作业之间的依赖关系管理：

**SLURM依赖作业提交指南：**
https://bioinformaticsworkbook.org/Appendix/HPC/SLURM/submitting-dependency-jobs-using-slurm.html#gsc.tab=0

### 作业状态监控

实时监控作业状态是集群管理的重要功能：

> sacct --starttime=2024-06-29 --format=JobID%10,User%20,Partition,Submit,Start,Elapsed,AllocTRES%50 -X

### 作业控制

作业的暂停、恢复和取消是日常管理操作：

> scontrol suspend jobid

### 用户账户管理

在SLURM集群中管理用户账户是系统管理员的职责：

> sacctmgr add user User=${u} Account=urgent

## 云计算与远程服务

### AWS EC2使用

AWS EC2是常用的云计算平台，掌握基本操作非常重要：

> AWS EC2 使用时的一些指令：上传文件、下载文件

### 环境变量配置

合理配置环境变量可以简化日常操作：

> export TZ='Asia/Shanghai'

## 文件系统与数据管理

### 文件压缩与解压

数据压缩和归档是数据管理的必备技能：

> 要清理pip的缓存，可以使用pip cache purge命令。这将清除pip缓存的所有内容，包括已下载但未安装的包和已安装但未使用的包的缓存。如果只想清除特定包的缓存，可以使用pip cache remove <package_name>命令，将package_name替换为要清除缓存的包名。

**查看pip缓存路径：**
> 要查看pip的缓存路径，可以使用pip cache dir命令。在命令行或终端中输入该命令，pip会显示其缓存的目录。

**Zip文件操作指南：**
https://docs.hostdime.com/hd/command-line/how-to-tar-untar-and-zip-files

### 文件搜索与过滤

高效的文件搜索和过滤可以大大提高工作效率：

> To find directory names with exactly two hyphens, you need to "anchor" your grep pattern to match the whole line.

**排除特定文件：**
> Exclude .batch and .extern by providing the -X option:

### Git版本控制

Git是现代科研项目的标准版本控制工具：

> 🔍 PDBFixer的工作机制：SEQRES与ATOM记录的对话

**Git忽略文件配置：**
> .gitignore 规则覆盖问题

## 编译与开发环境

### 编译系统理解

理解编译系统的工作原理有助于解决编译问题：

> gcc 的编译其实是四个过程的集合，分别是预处理（preprocessing）、编译（compilation）、汇编（assembly）、链接（linking）， 分别由 cpp、cc1、as、ld 这四个程序完成，gcc 是它们的封装。

### C++编程技巧

掌握C++编程技巧可以提高开发效率：

> In C++, the character "*" is a pointer which contains the value in a variable.

> To be accurate: ++i can sometimes be faster than i++ and is never slower. For fundamental data types, the compiler will very likely fix your mistake and optimise away any unneeded copying. For iterators this is more difficult and for user-defined types it may very well be impossible.

### Makefile编写

Makefile是自动化编译的重要工具：

> To create a Makefile that compiles two C++ source files (add.cpp and multi.cpp) separately into two different executables, you can follow the steps below. Here's an example of what your Makefile might look like:

### LaTeX排版系统

LaTeX是科学文档排版的标准工具：

> Regardless of your package choice you can install LaTeX by use of the apt command. The following linux command will install the LaTeX package: texlive-latex-extra. Replace the package name with the one you wish to install, open up terminal and enter:

**XeLaTeX安装：**
> The command to install XeLaTeX is:

> install the texlive-bibtex-extra package.

**中文字体支持问题：**
> The error message "LaTeX Error: File `ctexbook.cls' not found" in Ubuntu indicates that the ctexbook.cls document class file, which is part of the CTEX package for typesetting Chinese documents in LaTeX, is missing from your TeX Live installation.

**LaTeX安装指南：**
https://linuxconfig.org/how-to-install-latex-on-ubuntu-20-04-focal-fossa-linux

## 系统诊断与性能优化

### 系统监控工具

系统监控是保证服务稳定运行的关键：

> https://wz.anoms.top/2025/02/23/vscode-cache-diskspace-clean/

### 软件安装问题解决

解决软件安装过程中的常见问题：

> No rule to make target 'X' when X is simply missing.

## 云原生与容器技术

### 虚拟化技术

虚拟化技术是现代云计算的基础：

> A hypervisor, also known as a virtual machine monitor or VMM, is software that creates and runs virtual machines (VMs).

**虚拟化类型：**
> Type 1 hypervisor: Runs directly on the host's hardware to control the hardware and to manage guest operating systems. Examples include VMware ESXi, Microsoft Hyper-V, and Xen.

### Linux发行版选择

选择合适的Linux发行版对于特定应用场景很重要：

> The netinst version is a small ISO image that contains only the necessary files to start the installation. The DVD-1 version is a large ISO image that contains desktop environments, applications, and other software.

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