---
title: "几TB科研数据怎么安全搬家？"
date: "2026-11-21"
last_modified_at: "2026-11-21"
tags: [linux, rsync, data-migration, backup, system-administration, ntfs, technical-notes]
description: "几TB科研数据迁移备份的精简实录，含软链接展开、断点续传、百万小文件、NTFS损坏四个坑，Conda排除规则和一套可复用的rsync命令"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-g75zg3.jpg"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-g75zg3.jpg"
author: Xufan Gao
lang: zh-CN
published: false
---

# 几TB科研数据怎么安全搬家？

最近实验室调整服务器，要把几台旧机器的数据集中迁到新机器。原以为复制`/home`和`/data`、配置好再恢复就行，真做起来问题接连不断：软链接让数据重复复制，中断后要续传，几百万个小文件拖垮机械盘，最后一块刚买的4TB移动硬盘连NTFS文件系统都损坏了。

几天下来反复用`scp`、`rsync`、`du`、`iostat`、`smartctl`。本文只留这次迁移里最该记住的几个坑和处理命令，供以后复用。

## 坑一：scp把软链接展开，备份体积翻倍

最初用的是这条命令：

```bash
scp -r -P 2222 /home/* user@backup-host:/path/to/backup/home/
```

传输中断后检查，其中一台服务器的`home`备份竟占了约2TB，`data`备份另外占了约3TB。查一个用户目录就明白了：

```text
data -> /data/userA/
```

原本该保留成软链接的`data`，被复制成了真实目录。原因是OpenSSH的`scp -r`会跟随递归中遇到的软链接，把指向的内容一并复制。源数据若已在`/data`单独备份，这里等于又复制了一遍。

**备份Linux用户目录、尤其是有大量软链接时，应优先用rsync并保留链接本身**。`rsync -a`默认保留符号链接，只有`-L`（`--copy-links`）才主动展开。

```bash
rsync -aH --exclude='/*/data' /home/ /path/to/backup/home/
```

清理已被展开的重复副本前，务必先确认独立的`/data`备份完整，别在没核对时就批量删除`home/*/data`。

## 坑二：断点续传，参数选错等于白跑

几TB数据在SSH里传输，连接中断、终端关闭、网络波动都可能让任务中止。rsync的厉害之处是**可以反复运行同一条命令**：它比较两端，只传需要更新的文件。

```bash
rsync -aH --partial --info=progress2 -e 'ssh -p 2222' /data/ user@backup-host:/path/to/backup/server-c/data/
```

如果之前用`cp -r`或`scp -r`传过一部分，它们没有完整保留时间戳，后续rsync会把不少文件重新判为需要传输。补传阶段可以用`--size-only`按大小跳过：

```bash
rsync -rlt --size-only --partial --info=progress2 /home/ /path/to/backup/home/
```

这里`-r`是递归、`-l`保留软链接、`-t`保留修改时间。有一次只写了`-t`，输出停在`skipping directory .`，原因就是漏了`-r`，rsync根本没进目录。

几个参数的适用范围：

| 参数 | 行为 | 适用情况 |
|---|---|---|
| 默认比较 | 按大小和修改时间判断 | 正常续传 |
| `--size-only` | 大小一致就跳过 | 接续未保留时间戳的cp/scp复制 |
| `-c` | 读内容算校验和 | 核查数据是否一致 |
| `--ignore-existing` | 目标存在就跳过 | 明确不允许覆盖 |

但`--size-only`只是快速续传策略，**不能证明备份正确**。文件大小恰好一致、内容却损坏时也会被跳过。`--delete`会删除目标端多余内容，`--force`可能为替换软链接而删掉冲突的非空目录，这两个都不该当默认参数。

我们确实遇到过：先前被`scp`展开成真实目录的`data`，在rsync试图恢复软链接时报错`cannot delete non-empty directory`。稳妥做法是确认冲突目录只是重复备份后，先重命名暂存，同步完成再删。

## 坑三：三个rsync进程，速度却越来越慢

其中一台服务器通过USB连一块4TB的WD My Passport移动硬盘。前期速度还行，越往后越慢，百分比甚至从80%倒退回78%。

倒退不等于文件丢失。`ir-chk`表示rsync仍在增量扫描目录，发现的文件越多，待处理总数越大，整体进度估算跟着变。`top`里出现三个rsync也不奇怪，rsync内部用generator、sender、receiver三个角色协作，不代表跑了三份任务。

问题出在I/O表现上。当时的一条`iostat`结果：

| 指标 | 实际观察 |
|---|---|
| USB连接速率 | `480M` |
| 系统USB 3.x控制器 | 支持`5000M` |
| 磁盘`%util` | 长时间100% |
| 写入速度 | 部分采样仅0.12至0.48MB/s |
| `w_await` | 高达数百毫秒 |

`480M`是480Mbit/s，属于USB 2.0高速模式；USB 3.x SuperSpeed通常显示`5000M`或更高。这块盘支持USB 3.x，却实际跑在USB 2.0链路上。

但USB 2.0还解释不了每秒几百KB的写入。盘上还有几百万个小文件，机械盘频繁寻道、NTFS元数据更新和FUSE开销叠加在一起，才是主因。`du -d 1`同样受累，它要遍历下面的文件，遇到几十万文件的目录时会很久没有新输出。

**对海量小文件而言，文件数量往往比总容量更能影响备份耗时**。确认不需要的目录直接排除，比调rsync参数有效；必须完整保留的，可以先打包成tar归档，减少目标文件系统的元数据操作。

## 坑四：4TB移动硬盘的NTFS损坏了

先把结论放前面：**备份盘只要长期在Linux上用，就该一开始就格式化成ext4，别直接往NTFS分区写大量数据**。Linux靠FUSE或`ntfs3`驱动写NTFS，元数据操作开销大，权限和软链接也保不全。这块盘出厂就是NTFS，我们没提前重新格式化，后面才踩到这个坑。

复制继续时冒出`du: cannot access ... Input/output error`。原服务器对应目录正常，备份盘上的同一路径却无法访问。

先怀疑物理损坏，查SMART：

```bash
smartctl -H -A /dev/sdc
smartctl -l error /dev/sdc
```

结果反而很好：health `PASSED`，重映射扇区、待处理扇区、不可校正扇区全是0，没有错误日志，通电仅约11小时。

但SMART正常不代表文件系统没问题。重新连接时设备能识别，NTFS分区却挂不上：

```text
Failed to load $MFT:
Input/output error
```

`$MFT`是NTFS的主文件表，记录文件及其元数据。连它都加载不了，文件系统自然无法挂载。现有证据能确认NTFS元数据访问异常，但无法判定是文件系统软件问题、USB链路故障，还是底层存储曾出现未被SMART记录的错误。

如果盘上有唯一的重要数据，先停止写入、优先保护数据，再用Windows的`chkdsk`或专业工具处理。这块盘的备份还能从原服务器重取，我们重新格式化成ext4再备份。

**重新格式化只能重建文件系统，不能证明硬盘硬件一定可靠**。若设备再次出现I/O错误，必须重新检查硬盘、USB接口与线缆，不能继续把它当可信的唯一备份。

## 该跳过什么：Conda环境、编译软件和IDE服务端

多用户服务器里很多空间不是研究数据。Conda环境、编译安装的GROMACS或OpenMPI、Python虚拟环境、`node_modules`、各种IDE远程服务端文件，动辄几十上百GB，而且全是小文件，复制起来极慢。这些都能在新机器上重建，值得在备份前就排除。

排除规则写成独立文件`rsync-exclude.txt`，比在命令行堆十几个`--exclude`省事：

```text
miniconda*/
anaconda*/
miniforge*/
.conda/
.venv/
venv/
__pycache__/
node_modules/
.vscode-server/
```

调用时用`--exclude-from=rsync-exclude.txt`。注意已复制到备份端的同名目录不会被自动删除，要手动清理。

两条排除规则要小心。`gromacs*/`、`openmpi*/`可能匹配到用户自己修改过的源码甚至科研项目目录，不建议无条件排除；**Hugging Face和PyTorch的模型缓存尤其要小心**，模型默认在`~/.cache/huggingface/`和`.cache/torch/`，`.safetensors`、`.pt`这类权重也常藏在项目目录里。所以`.cache/`、`models/`、`checkpoints/`、`results/`都没有整体排除，是否备份要看实际内容。

## 直接能用的命令备忘

下面是整理后的安全流程，适用于源端Linux、目标端Linux文件系统、走SSH传输的情况。

开始前先确认设备：

```bash
lsblk -o NAME,SIZE,MODEL,FSTYPE,UUID,MOUNTPOINTS
df -hT
findmnt /data
```

要区分磁盘、分区和挂载点，别因为重启前后都叫`/dev/sdb1`就认定是同一块盘，核对型号、序列号和UUID更可靠。

首次备份分成`/home`和`/data`两条：

```bash
ID=server-c
DEST="user@backup-host:/path/to/backup/$ID"
rsync -aH --no-owner --no-group --partial --exclude='/*/data' --info=progress2 -e 'ssh -p 2222' /home/ "$DEST/home/"
rsync -aH --no-owner --no-group --partial --exclude-from=rsync-exclude.txt --info=progress2 -e 'ssh -p 2222' /data/ "$DEST/data/"
```

接收端是普通用户权限，所以没保留属主。若目的是完整恢复服务器，最好用有权限的备份账户加`-aHAX --numeric-ids`。排除规则用上一节整理的`rsync-exclude.txt`，内容按实际需要调整。

传输后用`-c`按内容校验、`-n`只检查不修改：

```bash
rsync -rltcni --exclude-from=rsync-exclude.txt -e 'ssh -p 2222' /data/ "$DEST/data/"
```

正常时同步范围内不该再有需要传输的差异。文件数可能达到几百万，完整校验会很慢。若源端在备份期间仍有作业写入，要等数据稳定后再做最终校验。长任务建议放进`tmux`，按`Ctrl+B`再按`D`退出而不中断。

## 这次留下的几条经验

正式复制前先检查用户目录结构，尤其是软链接和挂载点。对Linux服务器，`scp -r`不适合大规模备份，优先用rsync能省掉很多后续麻烦。

复制开始前就要确定哪些软件环境可以重建、哪些研究数据必须保留。提前规划排除规则，比传到一半再调整高效得多。

别只盯终端里的速度和百分比。文件系统类型、USB协商速率、磁盘等待时间和I/O是否还在推进，更能解释实际性能。

还有一条：**文件复制完成，不等于数据已经可靠备份**。目标盘可能损坏，软链接可能失效，权限可能丢失，部分内容可能在传输中被跳过。只有检查日志、核对内容、确认能恢复，才算完成备份。

> 数据能从备份中完整恢复，这次备份才算做完。少核对一个环节，之前花的时间都可能白费。
