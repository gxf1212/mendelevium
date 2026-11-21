---
title: "几TB数据、百万个文件和一块罢工的硬盘：Linux服务器迁移备份实录"
date: "2026-11-21"
last_modified_at: "2026-11-21"
tags: [linux, rsync, data-migration, backup, system-administration, ntfs, technical-notes]
description: "记录数TB科研数据从多台Linux服务器迁移的完整过程，梳理scp展开软链接、rsync断点续传与校验、百万小文件拖慢机械盘、NTFS文件系统损坏和硬盘重分配等问题的排查思路与可复用命令"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-g75zg3.jpg"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/thumbnail_mine/wh-g75zg3.jpg"
author: Xufan Gao
lang: zh-CN
---

# 几TB数据、百万个文件和一块罢工的硬盘：Linux服务器迁移备份实录

最近实验室调整服务器，准备将几台旧机器上的硬盘和数据集中迁移到新服务器。原以为把`/home`和`/data`复制到备份目录，等新机器配置好再恢复就行。实际操作起来，问题比预想的多得多：软链接导致数据重复复制，传输中断后需要续传，数百万个小文件让机械硬盘不堪重负，最后一块刚买不久的4TB移动硬盘甚至损坏了NTFS文件系统。

前后折腾了好几天，反复使用`scp`、`rsync`、`du`、`iostat`、`smartctl`等工具。趁着过程还记得清楚，把遇到的问题、原因和处理命令整理下来，以后迁移Linux服务器时也能直接参考。

## 一、最初的任务：把几台服务器的数据搬走

这次涉及server-a、server-b、server-c、server-d等多台实验室服务器。它们的硬件和使用情况不同，但用户数据大致遵循相同的目录结构：

```text
/home/
├── user1/
│   ├── .bashrc
│   ├── .cache/
│   ├── projects/
│   └── data -> /data/user1
└── user2/
    └── data -> /data/user2

/data/
├── user1/
└── user2/
```

`/home`一般保存用户配置、代码、环境和少量工作文件；`/data`则存放体积较大的科研数据。很多用户在家目录下创建了指向`/data`的软链接，便于访问。

最初的方案是分别备份`/home`和`/data`，通过SSH传到另一台服务器的共享存储，目录按照原服务器编号区分：

```text
/path/to/backup/
├── server-a/
│   ├── home/
│   └── data/
├── server-c/
│   ├── home/
│   └── data/
└── server-d/
    ├── home/
    └── data/
```

这样即使不同服务器上存在同名用户，备份也不会相互覆盖。

当时部分数据规模如下，均为迁移过程中的快照，并非最终备份容量：

| 服务器 | 源数据规模 | 当时的情况 |
|---|---|---|
| server-a | `/data`约3TB | `/home`复制过程中出现大量重复数据 |
| server-d | `/data`约1.8TB，`/home`约247GB | 一度只完成633GB和197GB |
| server-c | `/data`约2TB | 个别用户目录超过1TB，另有数百万个文件 |

麻烦从第一次`scp -r`开始。

## 二、第一个坑：scp把软链接展开了

最早用的是这条命令：

```bash
scp -r -P 2222 /home/* user@backup-host:/path/to/backup/home/
scp -r -P 2222 /data/* user@backup-host:/path/to/backup/data/
```

传输中断后检查备份，发现一个奇怪现象：210服务器的`home`备份竟然占了约2TB，而`data`备份另外占了约3TB。

检查其中一个用户：

```bash
ls -ld /home/userA/data
```

发现原目录实际上是：

```text
data -> /data/userA/
```

问题就出在这里。OpenSSH的`scp -r`会跟随递归遍历时遇到的软链接，将链接指向的实际内容复制过去。

于是本来只应该保留一个软链接的`home/userA/data`，在备份中变成了巨大的真实目录。源数据如果在`/data`中已经单独备份，相当于又复制了一遍。

**以后备份Linux用户目录，尤其是涉及大量软链接时，应当优先使用rsync，并保留链接本身**。

`rsync -a`会保留符号链接，而`-L`（`--copy-links`）才会主动展开它们。如果`/home`与`/data`分别备份，还可以直接排除用户目录中的`data`链接。

```bash
rsync -aH --exclude='/*/data' /home/ /path/to/backup/home/
```

这里的`/*/data`匹配传输根目录下各用户的`data`，不会误排除其他层级的同名目录。

至于已经被`scp`展开的重复副本，需要先确认独立的`/data`备份完整，再清理对应目录。尤其不要在没核对之前，直接批量删除所有`home/*/data`。

## 三、改用rsync：怎么跳过已有文件、继续传输？

服务器之间传输几TB数据，SSH连接中断、终端关闭、网络波动都可能造成任务中止。

`rsync`的优势是可以反复运行同一条命令：它会比较源端和目标端，只传输需要更新的文件。

例如，从server-c备份到另一台服务器：

```bash
rsync -aH --partial --info=progress2 -e 'ssh -p 2222' /data/ user@backup-host:/path/to/backup/server-c/data/
```

实际使用时，如果远端只有普通用户权限，可以添加`--no-owner --no-group`，避免尝试恢复原文件属主与组。

### 已经用cp或scp复制过一部分，怎么办？

这里涉及三个容易混淆的参数：

| 参数 | 行为 | 适用情况 |
|---|---|---|
| 默认比较 | 根据文件大小和修改时间判断 | 正常的rsync续传 |
| `--size-only` | 文件大小一致就跳过 | 接续先前未保留时间戳的cp/scp复制 |
| `-c` | 读取文件内容计算校验和 | 核查已复制数据是否一致 |
| `--ignore-existing` | 目标已存在就跳过 | 明确不允许覆盖现有文件 |

我们之前用过`cp -r`和`scp -r`，它们没有完整保留原文件时间戳，因此后续rsync可能把不少文件重新判断为需要传输。

所以在补传阶段采用了：

```bash
rsync -rlt --size-only --partial --info=progress2 /home/ /path/to/backup/home/
```

其中`-r`负责递归，`-l`保留软链接，`-t`保留修改时间。

有一次命令只写了`-t`，运行后提示：

```text
skipping directory .
0 100% (xfr#0, to-chk=0/0)
```

漏掉了递归选项`-r`，**rsync根本没有进入目录**。这个细节很容易被忽略。

需要注意，`--size-only`只是快速续传策略，不能用来证明备份内容完全正确。如果某个文件大小恰好一致、内容却已经损坏，它也会被跳过。`--partial`可以保留中断时的部分传输内容，但不意味着一定从上次中断的字节位置直接追加。

### rsync会不会覆盖原来的备份？

默认情况下：

- 目标端没有的文件，会新增。
- 路径相同、文件名相同且被判断需要更新的文件，会替换。
- 目标端多余的文件，在没有`--delete`时会保留。
- `--exclude`只是跳过匹配路径，并不会自动清理目标端已有内容。

因为每台服务器都放在独立的备份目录中，所以不同服务器之间不会发生路径冲突。

但有两个参数需要谨慎使用：`--delete`会删除目标端多余的内容；`--force`可能为了替换文件或软链接，删除目标端冲突的非空目录。

我们之前就遇到过：

```text
cannot delete non-empty directory: data
could not make way for new symlink: data
rsync error: code 23
```

这正是原先`scp`把软链接展开成真实目录后，rsync试图恢复软链接时产生的冲突。

更安全的处理是先确认冲突目录仅为重复备份，再将它重命名暂存，重新同步后才删除。不要把`--force`当作普通备份命令的默认参数。

## 四、怎么跳过Conda、编译软件等不必要的内容？

多用户科研服务器里，很多空间并不是研究数据。

常见的有Conda环境、编译安装的GROMACS或OpenMPI、Python虚拟环境、Node.js依赖，以及各种IDE远程服务端文件。

这些目录可能占用几十甚至上百GB，更重要的是包含大量小文件，复制起来非常慢。

我们后来整理了一套排除规则。相比把十几个`--exclude`都写在命令行里，更推荐保存成独立文件，例如`rsync-exclude.txt`：

```text
miniconda*/
anaconda*/
miniforge*/
mambaforge*/
.conda/
.venv/
venv/
.tox/
.pixi/
__pycache__/
.pytest_cache/
.mypy_cache/
.ruff_cache/
.ipynb_checkpoints/
node_modules/
.npm/
.vscode-server/
.cursor-server/
.qoder-server/
```

调用时只需要：

```bash
rsync -aH --partial --exclude-from=rsync-exclude.txt --info=progress2 /data/ /path/to/backup/data/
```

这些规则会匹配相应目录并跳过其内容。原来已经复制到备份端的同名目录不会被自动删除。

当时还用过`gromacs*/`和`openmpi*/`，但它们可能匹配到用户自己修改过的源码甚至科研项目目录，因此不建议无条件排除。

**Hugging Face和PyTorch的模型缓存尤其要小心**。

Hugging Face模型默认可能保存在`~/.cache/huggingface/`，PyTorch也可能使用`.cache/torch/`。项目目录里还可能有`.safetensors`、`.pt`、`.pth`、`.ckpt`等模型权重或训练checkpoint。

因此我们特意没有整体排除`.cache/`、`models/`、`checkpoints/`和`results/`。这些目录究竟是否需要备份，要根据里面的实际内容判断。

同样，`rsync -rlt`可以保留软链接，但不会跟随链接复制目标文件。如果链接指向的路径没有单独备份，最终恢复时仍然可能缺数据。

## 五、为什么三个rsync进程，复制速度却越来越慢？

这是这次备份中耗费时间最多的部分。

其中一台服务器最初尝试通过USB连接一块4TB的WD My Passport移动硬盘。前期速度还不错，复制推进后速率逐渐下降。

某次终端输出大致是：

```text
658,663,791,729  80%  12.52MB/s  13:56:04
(xfr#1340945, ir-chk=1000/1699982)
```

后来又变成：

```text
661,701,440,479  78%  3.22MB/s  14:51:20
(xfr#1342618, ir-chk=1007/1701881)
```

一个很奇怪的现象是，百分比居然从80%退回了78%。

原因在于`ir-chk`表示rsync仍在增量扫描目录，还没有得到完整文件列表。发现的文件越多，待处理文件总数越大，整体进度估算也跟着变。

所以百分比倒退并不意味着文件丢失。官方文档也说明，`--info=progress2`在增量递归扫描时显示的是不断更新的整体进度。

### 为什么top里面有三个rsync？

当时查看进程，发现同一任务出现了三个PID：

```text
30353 rsync
30354 rsync
30355 rsync
```

这是正常现象。rsync内部使用generator、sender和receiver三个角色协作完成文件比较与传输，相关进程并不代表运行了三份备份。

问题出在I/O表现上。

### 一次iostat揭示了性能瓶颈

当时运行：

```bash
iostat -dx 2 5 /dev/sdc
lsusb -t
ps -eo pid,comm,%cpu,stat,wchan:25 | grep -E 'rsync|ntfs'
```

结果是：

| 指标 | 实际观察 |
|---|---|
| USB连接速率 | `480M` |
| 系统USB 3.x控制器 | 支持`5000M` |
| 磁盘`%util` | 长时间100% |
| 写入速度 | 部分采样只有0.12–0.48MB/s |
| `w_await` | 高达数百毫秒 |
| NTFS进程 | `balance_dirty_pages`等待 |

原来，移动硬盘虽然支持USB 3.x，却实际运行在USB 2.0链路上。

`480M`表示480Mbit/s，是USB 2.0高速模式；正常的USB 3.x SuperSpeed连接通常应该显示`5000M`或更高。

不过，USB 2.0也不足以单独解释每秒只有几百KB的写入速度。移动硬盘还在处理数百万个小文件，机械硬盘频繁寻道、NTFS元数据更新和FUSE文件系统开销都会影响性能。

这也解释了为什么文件复制总量增长很慢，而`top`里进程仍然活着。

`du -d 1`也一样。它要遍历下面的文件才能算出目录容量，遇到包含几十万文件的目录时，可能很久都没有新的输出。

**对海量小文件而言，文件数量往往比总容量更能影响备份耗时**。

如果确认某些目录完全不需要备份，直接排除比调整rsync参数更有效。对于需要完整保留的海量小文件目录，也可以考虑先打包为tar归档，以减少目标文件系统上的元数据操作。

## 六、最麻烦的一次：4TB移动硬盘的NTFS损坏了

先把结论放在前面：**备份盘只要长期在Linux上使用，就应该一开始就格式化为ext4，别直接向NTFS分区写大量数据**。Linux通过FUSE或`ntfs3`驱动写入NTFS，元数据操作开销更大，权限和软链接也难以原样保留。这块移动硬盘出厂就是NTFS，我们没有提前重新格式化，才在后续备份中踩到这个问题。

复制继续，又遇到一个新错误：

```text
du: cannot access '.../include/sodium':
Input/output error
```

检查后确认，原服务器上的对应目录正常，但移动硬盘备份中的相同路径无法访问。

起初怀疑是硬盘物理损坏，于是检查SMART：

```bash
smartctl -H -A /dev/sdc
smartctl -l error /dev/sdc
```

结果反而很好：

```text
SMART overall-health self-assessment test result: PASSED
Reallocated_Sector_Ct: 0
Current_Pending_Sector: 0
Offline_Uncorrectable: 0
No Errors Logged
```

硬盘通电时间只有约11小时，温度30°C，没有发现明确的坏扇区记录。

但SMART正常并不能保证文件系统没有问题。后来重新连接硬盘时，设备仍能识别，NTFS分区却无法挂载。

```bash
lsblk -f
ntfs-3g.probe --readwrite /dev/sdc1
```

返回：

```text
ntfs_mapping_pairs_decompress() failed:
Input/output error

Failed to load $MFT:
Input/output error
```

内核日志也提示：

```text
ntfs3: sdc1: Mark volume as dirty due to NTFS errors
ntfs3: sdc1: Failed to load $MFT.
ntfs3: sdc1: It is recommended to use chkdsk.
```

即使尝试只读挂载，仍然失败。

这里的`$MFT`是NTFS的Master File Table（主文件表），用于记录文件系统中的文件及其元数据信息。连主文件表都无法正常加载，文件系统自然无法完成挂载。

现有证据已经确认NTFS元数据访问异常，但还不能断定究竟是文件系统软件问题、USB连接故障，还是底层存储曾出现未被SMART记录的错误。

如果磁盘上有唯一的重要数据，首先应该停止写入，优先保护数据，再使用Windows的`chkdsk`或专业恢复工具处理。

幸运的是，这块硬盘里的备份数据仍然可以从原服务器重新取得。我们重新格式化为ext4后再备份。

不过这次经历也说明，**重新格式化只能重建文件系统，不能证明硬盘硬件一定可靠**。如果设备再次出现I/O错误，就必须重新检查硬盘、USB接口与线缆，而不能继续把它当成可信的唯一备份。

## 七、迁移结束前，硬盘还得重新分配

数据备份之外，这次还有一项任务：把旧服务器拆下来的硬盘集中到其中一台，提高整体存储利用率。

管理员最初准备将两块2.4TB硬盘组成RAID1，再将两块600GB硬盘组成另一组RAID1。

RAID1会镜像保存相同数据，因此两块2.4TB硬盘组成RAID1后，实际可用容量仍然只有约2.4TB。

考虑到实验室的数据增长速度，我们后来决定将新增的两块2.4TB磁盘改为独立使用，换取约4.8TB总容量；600GB磁盘组仍然保留RAID1。

当时这台服务器的原始`/data`已经使用到约95%：

```text
Filesystem  Size  Used  Avail  Use%
/dev/sdb1   2.4T  2.2T   131G   95%
```

新增两块2.4TB磁盘一度挂载在`/data4`和`/data5`，600GB RAID1挂载在`/data3`。随后又讨论了重新调整挂载点名称，因此最终映射需要以配置完成后的`lsblk`和`findmnt`为准。

值得强调的是，新增加的两块2.4TB单盘也没有镜像冗余。用户原来的UID/GID可以随ext4文件系统保留下来，但跨机器使用前，需要确认新旧服务器的用户UID/GID对应关系一致。

未来如果继续多人共用数据盘，还可以配置ext4 quota，对不同用户设置容量限制。毕竟一块2.4TB的数据盘，在多个科研项目同时运行时，很容易再次被占满。

而且RAID1本身也不能替代独立备份：误删文件或发生逻辑损坏时，镜像可能同步保留错误状态。

## 八、以后直接复用的命令备忘

前面是实际经历，这里整理成较为安全、简洁的操作流程。以下命令适用于源端Linux服务器、目标端Linux文件系统，通过SSH传输的情况。

### 1. 开始前确认磁盘和目录

```bash
lsblk -o NAME,SIZE,MODEL,FSTYPE,UUID,MOUNTPOINTS
df -hT
findmnt /data
du -sh /home /data
```

要区分磁盘设备、分区和挂载点，不要因为重启前后都叫`/dev/sdb1`就直接认定它是同一块盘。更可靠的识别方式是检查型号、序列号和UUID。

### 2. 首次正式备份

在源服务器上运行，按实际服务器编号修改`ID`：

```bash
ID=server-c
DEST="user@backup-host:/path/to/backup/$ID"
```

备份`/home`，跳过各用户的`data`软链接：

```bash
rsync -aH --no-owner --no-group --partial --exclude='/*/data' --info=progress2 -e 'ssh -p 2222' /home/ "$DEST/home/"
```

备份`/data`，排除已确认不需要的软件环境：

```bash
rsync -aH --no-owner --no-group --partial --exclude-from=rsync-exclude.txt --info=progress2 -e 'ssh -p 2222' /data/ "$DEST/data/"
```

这里假设接收端是普通用户，因此没有保留原始所有权。若目的是日后完整恢复Linux服务器，最好由具备适当权限的备份账户接收，使用`-aHAX --numeric-ids`保留属主、组、ACL及扩展属性。

`rsync-exclude.txt`应当保存在源服务器，并根据实际需要调整；它不会自动包含本文列出的所有规则。

### 3. 中途断开后续传

通常直接重新运行原来的rsync命令即可。

如果先前是`scp`或`cp`复制，且目标端文件的时间戳与源端不同，可以在补传时考虑增加：

```bash
--size-only
```

这个选项只适合快速判断，不能取代最终数据校验。

如果准备运行很久，建议把任务放进`tmux`：

```bash
tmux new -s backup
```

运行rsync后，按`Ctrl+B`再按`D`即可退出会话而不终止任务。之后通过`tmux attach -t backup`重新进入。

### 4. 传输后验证完整性

对于需要严格检查的数据，使用`-c`按内容计算校验和，并通过`-n`预览差异：

```bash
rsync -rltcni --exclude-from=rsync-exclude.txt -e 'ssh -p 2222' /data/ "$DEST/data/"
```

`-n`表示只检查、不实际修改。正常情况下，同步范围内的文件应当不再出现需要传输的内容差异。由于文件数量可能达到数百万，完整校验也可能花费很长时间。

此外，还应该检查软链接、排除规则、权限和传输错误日志，特别是之前出现过`code 23`的目录。若源端在备份期间仍有作业持续写入，就应该在数据稳定后再做最终校验。

### 5. 磁盘变慢时检查什么？

```bash
iostat -dx 2 5
lsusb -t
ps -eo pid,stat,%cpu,wchan:25,args | grep -E 'rsync|ntfs'
journalctl -k -b | grep -Ei 'I/O error|usb|sdc|ntfs' | tail -30
```

`iostat`重点看实际吞吐、`w_await`和`%util`；`lsusb -t`可以确认USB是否工作在预期速率。

如果根本没有I/O进展，或者反复出现`Input/output error`，应当优先排查存储设备和文件系统，而不是继续增加并发任务。

### 6. 移动硬盘重新挂载

```bash
lsblk -f
blkid /dev/sdX1
findmnt -S /dev/sdX1
```

对于已有数据的ext4磁盘，如果只希望先检查内容，可以使用：

```bash
mkdir -p /mnt/diskcheck
mount -t ext4 -o ro,noload /dev/sdX1 /mnt/diskcheck
```

其中`sdX1`必须替换为确认过的实际分区，`ro,noload`可以避免读取检查时进行ext4日志回放写入。

检查完毕后再正常卸载，不要把未经确认的设备直接格式化或加入`/etc/fstab`。

## 九、这次迁移留下的几条经验

回头看，这次备份耗费大量时间，主要问题集中在几个环节。

正式复制前应该检查用户目录结构，尤其是软链接和挂载点。对于Linux服务器，`scp -r`并不是理想的大规模备份工具，默认优先选择rsync会省去很多后续麻烦。

复制开始前还要确定哪些软件环境可以重建，哪些研究数据必须保留。几百万个小文件可能让备份时间远超预期，提前规划排除规则，比传输到一半再调整效率高得多。

不能只看终端里的传输速度和百分比。文件系统类型、USB协商速率、磁盘等待时间，以及进程是否仍在正常进行I/O，往往更能解释实际性能。

还有一条：**文件复制完成，不等于数据已经可靠备份**。目标盘可能损坏，软链接可能失效，文件权限可能丢失，部分内容也可能在传输中被跳过。只有检查日志、核对内容并确认能够恢复，才算真正完成了备份。

> 数据能从备份中完整恢复，这次备份才算做完。少核对一个环节，之前花的时间都可能白费。

以后再做类似迁移，我会优先采用“盘点目录与文件系统—确定排除规则—rsync传输—校验完整性—再拆盘迁移”的顺序。可以少一些无谓的重复复制，也能避免在硬盘出问题时，才发现备份本身还有缺口。

这次最大的收获，大概就是终于把那些平时看似简单的Linux文件操作，在几TB真实科研数据上完整经历了一遍。

---

### 参考资料

1. [Rsync官方手册：参数、软链接、过滤规则与校验](https://rsync.samba.org/ftp/rsync/rsync.1)
2. [How Rsync Works：rsync的三个进程及工作机制](https://rsync.samba.org/how-rsync-works.html)
3. [OpenSSH scp手册：递归复制和软链接行为](https://man7.org/linux/man-pages/man1/scp.1.html)
4. [GNU Coreutils：df与du的统计差异](https://www.gnu.org/software/coreutils/faq/coreutils-faq.html)
5. [Linux Kernel：ext4文件系统及挂载选项](https://www.kernel.org/doc/html/latest/admin-guide/ext4.html)
6. [Microsoft Learn：chkdsk及NTFS文件系统检查](https://learn.microsoft.com/en-us/windows-server/administration/windows-commands/chkdsk)

*记录时间：2026年11月。文中的容量、传输进度及磁盘状态来自本次服务器迁移过程中的实际终端输出，用户名、主机地址和路径已做脱敏处理。*