---
title: "CentOS 7 升级到 Rocky Linux 8/9 完整指南"
date: "2025-10-08"
tags: [centos, rocky-linux, leapp, elevate, system-upgrade, linux-migration]
description: "从CentOS 7到Rocky Linux 8/9的系统迁移完整指南，详解ELevate和Leapp框架的使用方法，涵盖系统备份、软件源修复、内核升级和故障排查全流程"
thumbnail: "/assets/img/thumbnail_mine/wh-dp5x3l.jpg"
image: "/assets/img/thumbnail_mine/wh-dp5x3l.jpg"
---

# CentOS 7 升级到 Rocky Linux 8/9 完整指南

## 概述

随着 CentOS 7 于 2024 年 6 月 30 日正式停止维护，众多企业面临系统迁移的紧迫需求。Rocky Linux 作为 CentOS 的完美替代方案，为用户提供了稳定可靠的企业级解决方案。

本文将手把手教你使用 ELevate 项目和 Leapp 框架，实现从 CentOS 7 到 Rocky Linux 8/9 的无痛升级。

**项目主页链接**：
- ELevate项目主页：https://almalinux.org/elevate/
- Rocky Linux官方网站：https://rockylinux.org/
- CentOS官方网站：https://www.centos.org/

**重要警告**：
- 升级前务必备份所有重要数据并创建系统快照
- 生产环境推荐采用全新安装而非原地升级
- 原地升级存在一定风险，请先在测试环境验证

## 第一阶段：CentOS 7 升级到 Rocky Linux 8

### 准备工作

#### 1. 系统备份

创建系统快照（虚拟机环境）或备份重要配置和数据目录。

#### 2. 修复损坏的软件源（如需要）
```bash
cd /etc/yum.repos.d
mkdir bak
mv *.repo bak/

# 使用阿里云标准源
# 阿里云CentOS镜像：http://mirrors.aliyun.com/repo/
curl -o /etc/yum.repos.d/CentOS-Base.repo http://mirrors.aliyun.com/repo/Centos-7.repo

# 修复EPEL源配置（如已安装）
if [ -f "bak/epel.repo" ]; then
    cp bak/epel.repo /etc/yum.repos.d/
    sed -i 's/^metalink=/#metalink=/g' /etc/yum.repos.d/epel.repo
    sed -i 's/^#baseurl=/baseurl=/g' /etc/yum.repos.d/epel.repo
    sed -i 's|download.fedoraproject.org/pub|mirrors.aliyun.com|g' /etc/yum.repos.d/epel.repo
fi

yum clean all
yum makecache
yum install epel-release -y
yum update -y
```

#### 3. 检查系统状态
```bash
# 检查内核版本
rpm -qa | grep kernel
uname -r
cat /etc/redhat-release

# 清理旧内核（如需要）
# sudo yum remove kernel-3.10.0-1127.el7.x86_64 kernel-devel-3.10.0-1127.el7.x86_64
```

### 安装升级工具

#### 1. 安装 ELevate 和 Leapp
```bash
# 下载并安装 ELevate 仓库
# ELevate仓库地址：https://repo.almalinux.org/elevate/
curl -k -L -o /tmp/elevate-release-latest-el7.noarch.rpm https://repo.almalinux.org/elevate/elevate-release-latest-el7.noarch.rpm
yum localinstall -y /tmp/elevate-release-latest-el7.noarch.rpm

# 修复ELevate仓库SSL证书问题
cat > /etc/yum.repos.d/ELevate.repo << 'EOF'
# ELevate project repo for el7

[elevate]
name=ELevate
baseurl=https://repo.almalinux.org/elevate/el7/$basearch/
gpgcheck=1
enabled=1
priority=90
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-ELevate
sslverify=0

## Sources
[elevate-source]
name=ELevate - Source
baseurl=https://repo.almalinux.org/elevate/el7/SRPMS/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-ELevate
sslverify=0
EOF

# 清理缓存并安装升级工具
yum clean all
yum install -y leapp-upgrade leapp-data-rocky

# 移除可能冲突的包
yum remove javapackages-tools -y   # CUDA 相关包会被移除
```

### 执行升级

#### 1. 预升级检查
```bash
leapp preupgrade
```

**说明**：预升级检查会生成报告文件 `/var/log/leapp/leapp-report.txt`，包含所有潜在问题和解决方案。

#### 2. 执行升级
```bash
leapp upgrade
```

这需要一段时间，大概十几分钟，请保持网络稳定。即使远程运行也可以实时查看`/var/log/leapp/leapp-upgrade.log`来知道安装进度

![alt text](image.png)

> 这里用cc远程弄了一次

**注意**：Conda 环境可能在重启后出现冲突，但通常多次重启后会自动解决。

#### 3. 重启系统
```bash
reboot
```

![alt text](2b4a4bb327fc8cc5af9e44f9a9380a8e.jpeg)

![alt text](801587ec97f1b1e747f2018733f3c850.jpeg)

这里又得等好一会，执行剩下的升级。系统会自动重启两次，完成后 GRUB 菜单会显示 Rocky Linux 条目。

### 常见问题解决

#### 1. yum锁定问题
如遇到 "Another app is currently holding the yum lock" 错误：
```bash
# 强制杀掉相关进程
pkill -9 yum
pkill -9 PackageKit
pkill -9 packagekitd

# 删除锁文件
rm -f /var/run/yum.pid

# 彻底停用PackageKit
systemctl stop packagekit
systemctl disable packagekit
systemctl mask packagekit

# 禁用PackageKit插件
echo 'enabled=0' > /etc/yum/pluginconf.d/refresh-packagekit.conf
```

#### 2. EPEL仓库metalink错误
```bash
# 修复EPEL仓库配置
sed -i 's/^metalink=/#metalink=/g' /etc/yum.repos.d/epel.repo
sed -i 's/^#baseurl=/baseurl=/g' /etc/yum.repos.d/epel.repo
sed -i 's|download.fedoraproject.org/pub|mirrors.aliyun.com|g' /etc/yum.repos.d/epel.repo
```

#### 3. BTRFS 相关错误
如遇到 "btrfs has been removed from anolis8" 错误，这通常是正常的，因为大多数系统并未使用 btrfs 分区。

#### 4. 外部仓库包冲突
```bash
# 移除可能冲突的 EPEL 包
yum remove <package_name>
```

#### 5. Leapp升级网络连接失败
如果遇到 `Failed to synchronize cache for repo 'rocky8-*'` 错误：

```bash
# 创建leapp专用DNF配置（解决代理兼容性问题）
mkdir -p /etc/leapp/files
cat > /etc/leapp/files/dnf.conf << 'EOF'
[main]
gpgcheck=1
installonly_limit=3
clean_requirements_on_remove=True
best=True
skip_if_unavailable=False
proxy=socks5://127.0.0.1:1080
sslverify=0
timeout=300
retries=10
EOF

# 修复软件源配置（使用直接URL代替mirrorlist）
cp /etc/leapp/files/leapp_upgrade_repositories.repo{,.backup}
sed -i 's|^mirrorlist=.*|#&|; s|^#baseurl=.*|baseurl=https://download.rockylinux.org/pub/rocky/8/BaseOS/x86_64/os/|' /etc/leapp/files/leapp_upgrade_repositories.repo
sed -i '/^\[/a sslverify=0' /etc/leapp/files/leapp_upgrade_repositories.repo

# 清理并重试
rm -rf /var/lib/leapp/* /tmp/leapp_*
leapp upgrade --no-rhsm --target 8.10
```

#### 6. 升级中断恢复
```bash
# 尝试恢复升级（如果支持）
leapp upgrade --resume --no-rhsm

# 如果不支持恢复，清理后重新开始
rm -rf /var/lib/leapp/* /tmp/leapp_*
leapp upgrade --no-rhsm --target 8.10
```

### GRUB引导故障修复指南

#### 常见GRUB问题
升级过程中可能遇到GRUB引导失败，这是leapp升级的已知问题。典型症状包括：
- 重启后进入GRUB命令行模式
- 内核文件丢失或无法找到
- 系统无法正常启动

**参考资源**：
- Red Hat GRUB问题解决方案：https://access.redhat.com/solutions/7004146
- GRUB修复指南：https://phoenixnap.com/kb/grub-rescue
- CentOS GRUB救援命令：https://linuxhint.com/grub_rescue_commands_centos/

#### 方法1：使用Rocky Linux 8救援模式修复
```bash
# 1. 从Rocky Linux 8安装ISO启动
# 2. 选择 "Troubleshooting" -> "Rescue a Rocky Linux system"
# 3. 选择挂载文件系统选项 "1"

# 4. 进入chroot环境
chroot /mnt/sysimage

# 5. 重新安装GRUB和内核
grub2-install /dev/sda
grub2-mkconfig -o /boot/grub2/grub.cfg

# 6. 如果内核丢失，重新安装
dnf install kernel

# 7. 退出并重启
exit
reboot
```

#### 方法2：GRUB命令行紧急启动
如果在grub>提示符下，尝试以下命令：
```bash
# 加载LVM模块
grub> insmod lvm

# 查看可用分区
grub> ls

# 设置根分区（根据实际情况调整）
grub> set root=(lvm/centos-root)

# 手动加载内核（版本号需要根据实际情况调整）
grub> linux /boot/vmlinuz-4.18.0-553.el8_10.x86_64 root=/dev/mapper/centos-root ro

# 加载initrd
grub> initrd /boot/initramfs-4.18.0-553.el8_10.x86_64.img

# 启动系统
grub> boot
```

#### 方法3：预防措施
```bash
# 升级前检查磁盘空间（GRUB需要至少1024KB空间）
df -h /boot

# 备份当前GRUB配置
cp /boot/grub2/grub.cfg /boot/grub2/grub.cfg.backup

# 确保系统更新到最新
yum update -y
```

## 第二阶段：Rocky Linux 8 升级到 Rocky Linux 9

### 准备工作

参考链接：
- Phoenix NAP 升级指南：https://phoenixnap.com/kb/upgrade-rocky-linux-8-to-9
- Vultr 升级文档：https://docs.vultr.com/how-to-upgrade-from-rocky-linux-8-to-rocky-linux-9#upgrade-rocky-linux-8-to-rocky-linux-9
- ZJU 镜像站文档：https://mirrors.zju.edu.cn/docs/rocky/

#### 1. 安装 Rocky Linux 9 GPG 密钥
```bash
# 浙江大学Rocky Linux镜像：https://mirrors.zju.edu.cn/rocky/
wget https://mirrors.zju.edu.cn/rocky/9.5/BaseOS/x86_64/os/Packages/r/rocky-gpg-keys-9.5-1.2.el9.noarch.rpm
sudo rpm -ivh rocky-gpg-keys-9.5-1.2.el9.noarch.rpm
```

#### 2. 备份软件源配置
```bash
cp -r /etc/yum.repos.d/ /etc/yum.repos.d.bak8
```

#### 3. 更新软件源为浙大镜像

**方法 1：批量更新 EPEL 源**
```bash
for repo_file in /etc/yum.repos.d/*.repo; do
  sed -e 's!^metalink=!#metalink=!g' \
      -e 's!^#baseurl=!baseurl=!g' \
      -e 's!https://download\.example/pub/epel/!https://mirrors.zju.edu.cn/epel/!g' \
      -e 's!https://mirrors\.fedoraproject\.org/metalink!#https://mirrors.fedoraproject.org/metalink!g' \
      -i "$repo_file"
done
```

**方法 2：更新 Rocky 官方源**
```bash
sed -e 's|^mirrorlist=|#mirrorlist=|g' \
    -e 's|^baseurl=http://dl.rockylinux.org/$contentdir|baseurl=https://mirrors.zju.edu.cn/rocky|g' \
    -i.bak \
    /etc/yum.repos.d/Rocky-AppStream.repo \
    /etc/yum.repos.d/Rocky-BaseOS.repo \
    /etc/yum.repos.d/Rocky-Extras.repo \
    /etc/yum.repos.d/Rocky-PowerTools.repo
```

#### 4. 清理和准备升级环境
```bash
# 备份 Elevate.repo 文件
mv /etc/yum.repos.d/ELevate.repo /etc/yum.repos.d/Elevate.repo.bak

# 清除缓存并重建
sudo yum clean all
sudo yum makecache
sudo dnf upgrade --refresh

# 重新安装 xl2tpd（用于 ZJU 网络）
sudo yum remove xl2tpd
sudo yum install xl2tpd
```

### 执行升级到 Rocky 9

#### 1. 修改 DNF 配置
```bash
# 取消 exclude 行的注释
sudo sed -i 's/^exclude=/#exclude=/g' /etc/dnf/dnf.conf
```

#### 2. 移除旧版本包和依赖
```bash
# 移除 Leapp 相关包
sudo dnf remove leapp leapp-upgrade-el7toel8 python2-leapp

# 移除其他冲突包
sudo dnf -y remove rpmconf yum-utils epel-release
sudo rm -rf /usr/share/redhat-logos

# 移除 Python 2 相关包
sudo dnf remove --skip-broken --nobest python2 python2-libs python2-pip python2-setuptools \
  python2-requests python2-pytz python2-coverage python2-idna python2-backports python2-lxml \
  python2-backports-ssl_match_hostname python2-ipaddress pygobject2 python2-pysocks \
  python2-urllib3 python2-pyyaml python2-chardet python2-six python2-cairo

# 移除其他冲突组件
sudo dnf remove make-devel iptables-ebtables
```

#### 3. 安装升级工具
```bash
sudo dnf install dnf-plugin-system-upgrade
```

#### 4. 强制移除遗留包
```bash
sudo rpm -e --nodeps leapp leapp-upgrade-el7toel8 python2-leapp python2 python2-libs
```

#### 5. 导入 GPG 密钥
```bash
# Rocky Linux官方GPG密钥：https://dl.rockylinux.org/pub/rocky/
sudo rpm --import https://dl.rockylinux.org/pub/rocky/RPM-GPG-KEY-Rocky-9
```

#### 6. 执行系统升级
```bash
sudo dnf -y --releasever=9 --allowerasing --setopt=deltarpm=false distro-sync
sudo rpm --rebuilddb
sudo reboot
```

### 升级后验证和清理

#### 1. 验证系统版本
```bash
cat /etc/redhat-release
cat /etc/os-release
```

#### 2. 完成系统更新
```bash
sudo dnf update --allowerasing
```

#### 3. 重新安装 CUDA（如需要）
```bash
sudo dnf module install nvidia-driver:latest
dnf install nvidia-driver-cuda -y

# NVIDIA CUDA官方仓库：https://developer.download.nvidia.com/compute/cuda/repos/
sudo dnf config-manager --add-repo https://developer.download.nvidia.com/compute/cuda/repos/rhel9/x86_64/cuda-rhel9.repo
sudo dnf clean all
sudo dnf -y install cuda-toolkit-12-6
```

### EPEL 9 仓库配置

创建或更新 `/etc/yum.repos.d/epel.repo`：

```ini
[epel]
name=Extra Packages for Enterprise Linux $releasever - $basearch
baseurl=https://mirrors.zju.edu.cn/epel/$releasever/Everything/$basearch/
#mirrorlist=https://mirrors.fedoraproject.org/metalink?repo=epel-$releasever&arch=$basearch&infra=$infra&content=$contentdir
enabled=1
gpgcheck=1
countme=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-$releasever

[epel-debuginfo]
name=Extra Packages for Enterprise Linux $releasever - $basearch - Debug
baseurl=https://mirrors.zju.edu.cn/epel/$releasever/Everything/$basearch/debug/
#mirrorlist=https://mirrors.fedoraproject.org/metalink?repo=epel-debug-$releasever&arch=$basearch&infra=$infra&content=$contentdir
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-$releasever
gpgcheck=1

[epel-source]
name=Extra Packages for Enterprise Linux $releasever - $basearch - Source
baseurl=https://mirrors.zju.edu.cn/epel/$releasever/Everything/source/tree/
#mirrorlink=https://mirrors.fedoraproject.org/metalink?repo=epel-source-$releasever&arch=$basearch&infra=$infra&content=$contentdir
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-$releasever
gpgcheck=1
```

## 总结

本指南提供了从 CentOS 7 到 Rocky Linux 9 的完整升级路径：
1. **阶段一**：CentOS 7 → Rocky Linux 8（使用 Leapp）
2. **阶段二**：Rocky Linux 8 → Rocky Linux 9（使用 DNF 系统升级）

**最佳实践建议**：
- 生产环境建议采用全新安装而非原地升级
- 升级前充分测试并制定回滚计划
- 定期备份系统和数据
- 关注官方文档更新

## 相关教程和参考资料

### 官方文档
- Rocky Linux 迁移官方指南：https://docs.rockylinux.org/guides/migrate2rocky/
- ELevate 项目快速入门指南：https://wiki.almalinux.org/elevate/ELevate-quickstart-guide.html
- AlmaLinux ELevate 项目主页：https://almalinux.org/elevate/
- ELevate 项目 - CloudLinux：https://cloudlinux.com/elevate/
- Red Hat Leapp 升级文档：https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/8/html/upgrading_from_rhel_7_to_rhel_8/

### 详细教程
- 从 CentOS 7 迁移到 Rocky Linux 8 详细指南 - Linuxiac：https://linuxiac.com/migrating-from-centos-7-to-rocky-linux-8/
- CentOS 7 到 Rocky Linux 8 转换教程 - First2Host：https://first2host.co.uk/blog/migrate-centos-7-rocky-linux-8/
- CentOS 7.x 原地升级到 Rocky Linux 8.x - JetPatch：https://kc.jetpatch.com/hc/en-us/articles/28894194238989-In-place-upgrade-from-CentOS-7-x-to-Rocky-linux-8-x
- CentOS 7 迁移到 Rocky Linux 9 - phoenixNAP：https://phoenixnap.com/kb/migrate-centos-to-rocky-linux
- CentOS 7 迁移到 Rocky Linux 9 指南 - Medium：https://medium.com/@redswitches/how-to-migrate-centos-7-to-rocky-linux-9-bc00db9e4ee7

### Rocky 8 到 9 升级
- Phoenix NAP 升级指南：https://phoenixnap.com/kb/upgrade-rocky-linux-8-to-9
- Vultr 升级文档：https://docs.vultr.com/how-to-upgrade-from-rocky-linux-8-to-rocky-linux-9#upgrade-rocky-linux-8-to-rocky-linux-9
- Rocky Linux 8 到 9 升级 - Linuxiac：https://linuxiac.com/upgrade-rocky-linux-8-to-rocky-linux-9/
- Rocky Linux 8 到 9 升级 - Shapehost：https://shape.host/resources/how-to-upgrade-from-rocky-linux-8-to-rocky-linux-9

### 技术资源和工具
- GitHub - CentOS 7 升级到 8 脚本：https://gist.github.com/Trogvars/d93f8e370e9d01d4afc6e2a7e8c69ab2
- Linux Notes: ELevate - leapp 迁移工具：https://neilrieck.net/docs/linux_notes_leapp.html
- CentOS 到 Rocky Linux 迁移规划 - OpenLogic：https://www.openlogic.com/blog/planning-centos-rocky-linux-migration
- CIQ Ascender CentOS 7 到 Rocky 8 迁移：https://ciq.com/blog/ascender-migrates-host-from-centos-7-to-rocky-8/

### GRUB故障排查资源
- GRUB修复指南 - Phoenix NAP：https://phoenixnap.com/kb/grub-rescue
- GRUB救援模式修复 - HowToForge：https://www.howtoforge.com/tutorial/repair-linux-boot-with-grub-rescue/
- CentOS GRUB救援命令：https://linuxhint.com/grub_rescue_commands_centos/
- Red Hat GRUB问题解决方案：https://access.redhat.com/solutions/7004146

### 镜像源和下载
- 浙江大学 Rocky Linux 镜像站：https://mirrors.zju.edu.cn/docs/rocky/
- 阿里云 CentOS 镜像源：http://mirrors.aliyun.com/repo/
- 阿里云 EPEL 镜像：https://mirrors.aliyun.com/epel/
- Rocky Linux 官方下载：https://download.rockylinux.org/
- ELevate 项目仓库：https://repo.almalinux.org/elevate/

### 社区支持
- Rocky Linux 官方论坛：https://forums.rockylinux.org/
- CentOS 官方论坛：https://forums.centos.org/
- AlmaLinux 社区聊天室（~migration 频道）：https://chat.almalinux.org/
- Server Fault 社区：https://serverfault.com/
- Red Hat 客户门户：https://access.redhat.com/

### 官方网站和文档
- Rocky Linux 官方文档：https://docs.rockylinux.org/
- ELevate 项目：https://wiki.almalinux.org/elevate/
- 浙江大学镜像站：https://mirrors.zju.edu.cn/docs/rocky/
- CentOS 官方文档：https://docs.centos.org/
- NVIDIA CUDA 官方文档：https://docs.nvidia.com/cuda/
