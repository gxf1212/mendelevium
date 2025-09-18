
---
title: "CentOS 7升级Rocky Linux 8无网络环境解决方案"
date: "2025-09-17"
tags: [centos, rocky-linux, 无网络升级, ssh隧道, 代理配置, 系统升级, linux运维]
---

# CentOS 7升级Rocky Linux 8无网络环境解决方案

## 无网络环境升级解决方案

### 场景说明
在实际生产环境中，很多服务器出于安全考虑无法直接访问互联网。比如你浙，zjunet已经无法在老机子上安装了，Linux上无线网卡又不好使。本文将详细介绍如何在这种受限网络环境下，成功完成CentOS 7到Rocky Linux 8的平滑升级。

本文档参考了ELevate项目官方文档和社区最佳实践，ELevate项目主页：https://almalinux.org/elevate/

### 方法一：SSH动态代理隧道

#### 适用场景
- 有一台可以联网的跳板机
- 待升级机器可以SSH连接到跳板机
- 跳板机可以SSH连接到待升级机器

参考SSH隧道配置指南：https://www.ssh.com/academy/ssh/tunneling

#### 配置步骤

**1. 在待升级机器上建立SSH隧道**
```bash
# 在待升级机器上执行（后台运行）
ssh -D 1080 user@跳板机IP &
# 示例：ssh -D 1080 gxf1212@10.77.14.189 &
```

**2. 配置yum使用SOCKS5代理**
```bash
# 在yum.conf中添加代理配置
echo 'proxy=socks5://127.0.0.1:1080' >> /etc/yum.conf

# 验证代理是否工作
curl --proxy socks5://127.0.0.1:1080 -I http://www.baidu.com
```

**3. 解决常见代理问题**
```bash
# 如果遇到SSL证书过期问题，跳过SSL验证
# 对于需要下载的rpm包
curl --proxy socks5://127.0.0.1:1080 -k -L -o /tmp/package.rpm https://example.com/package.rpm

# 为ELevate仓库添加SSL跳过设置
echo 'sslverify=0' >> /etc/yum.repos.d/ELevate.repo
```

### 方法二：离线软件包准备

#### 适用场景
- 完全无网络环境
- 需要预先在联网机器上准备软件包

#### 准备软件包（在联网机器上执行）

**1. 下载ELevate相关包**
```bash
# 创建下载目录
mkdir -p /tmp/centos7-upgrade-packages

# 下载ELevate仓库包
# ELevate仓库地址：https://repo.almalinux.org/elevate/
curl -k -L -o /tmp/centos7-upgrade-packages/elevate-release-latest-el7.noarch.rpm \
    https://repo.almalinux.org/elevate/elevate-release-latest-el7.noarch.rpm

# 配置临时ELevate仓库
yum install -y /tmp/centos7-upgrade-packages/elevate-release-latest-el7.noarch.rpm

# 下载leapp相关包及其依赖
yumdownloader --resolve --destdir=/tmp/centos7-upgrade-packages \
    leapp-upgrade leapp-data-rocky
```

**2. 传输软件包到目标机器**
```bash
# 使用scp传输软件包目录
scp -r /tmp/centos7-upgrade-packages/ root@target-server:/tmp/

# 或使用rsync
rsync -avz /tmp/centos7-upgrade-packages/ root@target-server:/tmp/centos7-upgrade-packages/
```

#### 在目标机器上安装（无网络环境）

```bash
# 安装ELevate仓库
yum localinstall -y /tmp/centos7-upgrade-packages/elevate-release-latest-el7.noarch.rpm

# 安装所有下载的包
yum localinstall -y /tmp/centos7-upgrade-packages/*.rpm

# 继续正常的升级流程
leapp preupgrade
leapp upgrade
reboot
```

### 软件源配置文件模板

#### CentOS 7 基础源配置
创建 `/etc/yum.repos.d/CentOS-Base.repo`：

阿里云CentOS镜像源：https://mirrors.aliyun.com/centos/
```ini
[base]
name=CentOS-7 - Base - mirrors.aliyun.com
baseurl=http://mirrors.aliyun.com/centos/7/os/$basearch/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7

[updates]
name=CentOS-7 - Updates - mirrors.aliyun.com
baseurl=http://mirrors.aliyun.com/centos/7/updates/$basearch/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7

[extras]
name=CentOS-7 - Extras - mirrors.aliyun.com
baseurl=http://mirrors.aliyun.com/centos/7/extras/$basearch/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7

[centosplus]
name=CentOS-7 - Plus - mirrors.aliyun.com
baseurl=http://mirrors.aliyun.com/centos/7/centosplus/$basearch/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7
```

#### EPEL 7 源配置（修复版）
创建 `/etc/yum.repos.d/epel.repo`：

EPEL项目主页：https://fedoraproject.org/wiki/EPEL
阿里云EPEL镜像源：https://mirrors.aliyun.com/epel/
```ini
[epel]
name=Extra Packages for Enterprise Linux 7 - $basearch
baseurl=http://mirrors.aliyun.com/epel/7/$basearch
#mirrorlist=https://mirrors.fedoraproject.org/metalink?repo=epel-7&arch=$basearch
failovermethod=priority
enabled=1
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-7

[epel-debuginfo]
name=Extra Packages for Enterprise Linux 7 - $basearch - Debug
baseurl=http://mirrors.aliyun.com/epel/7/$basearch/debug
#mirrorlist=https://mirrors.fedoraproject.org/metalink?repo=epel-debug-7&arch=$basearch
failovermethod=priority
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-7
gpgcheck=1

[epel-source]
name=Extra Packages for Enterprise Linux 7 - $basearch - Source
baseurl=http://mirrors.aliyun.com/epel/7/SRPMS
#mirrorlist=https://mirrors.fedoraproject.org/metalink?repo=epel-source-7&arch=$basearch
failovermethod=priority
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-7
gpgcheck=1
```

#### ELevate 源配置（修复SSL问题版）
创建 `/etc/yum.repos.d/ELevate.repo`：

ELevate项目仓库：https://repo.almalinux.org/elevate/
```ini
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
```

### 代理配置最佳实践

#### yum.conf 代理配置
```ini
[main]
cachedir=/var/cache/yum/$basearch/$releasever
keepcache=0
debuglevel=2
logfile=/var/log/yum.log
exactarch=1
obsoletes=1
gpgcheck=1
plugins=1
installonly_limit=5
bugtracker_url=http://bugs.centos.org/set_project.php?project_id=23&ref=http://bugs.centos.org/bug_report_page.php?category=yum
distroverpkg=centos-release

# 代理配置（根据实际情况选择一种）
# HTTP代理
#proxy=http://proxy-server:port
#proxy_username=username
#proxy_password=password

# SOCKS5代理（推荐用于SSH隧道）
proxy=socks5://127.0.0.1:1080
```

### 关键修复步骤

#### leapp升级网络失败修复

如果遇到 `Failed to synchronize cache for repo 'rocky8-*'` 错误：

```bash
# 1. 创建leapp专用DNF配置
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

# 2. 更新系统DNF配置
cat > /etc/dnf/dnf.conf << 'EOF'
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

# 3. 修复Rocky Linux 8软件源配置（使用直接URL而非mirrorlist）
cp /etc/leapp/files/leapp_upgrade_repositories.repo /etc/leapp/files/leapp_upgrade_repositories.repo.backup
cat > /etc/leapp/files/leapp_upgrade_repositories.repo << 'EOF'
[rocky8-baseos]
name=Rocky Linux 8 - BaseOS
baseurl=https://download.rockylinux.org/pub/rocky/8/BaseOS/x86_64/os/
gpgcheck=1
enabled=1
gpgkey=file:///etc/leapp/repos.d/system_upgrade/common/files/rpm-gpg/8/RPM-GPG-KEY-Rocky-8
sslverify=0

[rocky8-appstream]
name=Rocky Linux 8 - AppStream
baseurl=https://download.rockylinux.org/pub/rocky/8/AppStream/x86_64/os/
gpgcheck=1
enabled=1
gpgkey=file:///etc/leapp/repos.d/system_upgrade/common/files/rpm-gpg/8/RPM-GPG-KEY-Rocky-8
sslverify=0

[rocky8-extras]
name=Rocky Linux 8 - Extras
baseurl=https://download.rockylinux.org/pub/rocky/8/extras/x86_64/os/
gpgcheck=1
enabled=1
gpgkey=file:///etc/leapp/repos.d/system_upgrade/common/files/rpm-gpg/8/RPM-GPG-KEY-Rocky-8
sslverify=0
EOF

# 4. 清理leapp状态并重试
rm -rf /var/lib/leapp/* /tmp/leapp_*
leapp upgrade --no-rhsm --target 8.10
```

#### 升级中断恢复

如果升级过程中断：

```bash
# 1. 检查leapp状态
ls -la /var/lib/leapp/

# 2. 尝试恢复升级（如果支持）
leapp upgrade --resume --no-rhsm

# 3. 如果resume不支持，清理后重新开始
rm -rf /var/lib/leapp/* /tmp/leapp_*
leapp upgrade --no-rhsm --target 8.10
```

### GRUB引导修复指南

#### 常见GRUB问题
如果升级后遇到GRUB引导问题，这是leapp升级的已知问题。参考解决方案：

**Red Hat GRUB修复文档**：https://access.redhat.com/solutions/7004146
**Rocky Linux救援模式指南**：https://docs.rockylinux.org/guides/
**GRUB修复社区文档**：https://phoenixnap.com/kb/grub-rescue

```bash
# 从Rocky Linux 8救援盘启动后执行：
chroot /mnt/sysimage
grub2-install /dev/sda
grub2-mkconfig -o /boot/grub2/grub.cfg
```

#### GRUB命令行紧急启动
```bash
# 在grub>提示符下执行：
insmod lvm
set root=(lvm/centos-root)
linux /boot/vmlinuz-4.18.0-553.el8_10.x86_64 root=/dev/mapper/centos-root ro
initrd /boot/initramfs-4.18.0-553.el8_10.x86_64.img
boot
```

**GRUB救援命令参考**：https://linuxhint.com/grub_rescue_commands_centos/

### 故障排查指南

#### 1. 验证网络连通性
```bash
# 测试基本网络连接
ping -c 3 8.8.8.8

# 测试域名解析
nslookup mirrors.aliyun.com

# 测试HTTP连接
curl -I http://mirrors.aliyun.com/

# 测试HTTPS连接
curl -I https://mirrors.aliyun.com/

# 测试代理连接
curl --proxy socks5://127.0.0.1:1080 -I http://www.baidu.com
```

#### 2. 检查SSH隧道状态
```bash
# 检查SSH隧道进程
ps aux | grep "ssh -D"

# 检查监听端口
netstat -tlnp | grep 1080

# 重新建立SSH隧道（如果断开）
ssh -D 1080 user@jumphost &
```

#### 3. 清理和重试
```bash
# 清理yum缓存
yum clean all
rm -rf /var/cache/yum/*

# 重新生成缓存
yum makecache

# 测试软件源
yum repolist
```

### 注意事项

1. **SSH隧道稳定性**：确保SSH隧道在整个升级过程中保持稳定，建议使用screen或tmux来管理长时间运行的任务。

2. **带宽考虑**：升级过程需要下载大量软件包，确保网络带宽足够。

3. **防火墙设置**：检查跳板机和目标机器的防火墙配置，确保必要的端口开放。

4. **备份重要性**：无网络环境下出现问题更难修复，务必在升级前做好完整备份。

5. **测试环境**：建议先在类似的测试环境中验证整个流程。

### 相关资源链接

**官方文档**：
- ELevate项目主页：https://almalinux.org/elevate/
- Rocky Linux官方文档：https://docs.rockylinux.org/
- Red Hat leapp工具文档：https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/8/html/upgrading_from_rhel_7_to_rhel_8/
- CentOS官方文档：https://docs.centos.org/

**镜像源**：
- 阿里云CentOS镜像：https://mirrors.aliyun.com/centos/
- 阿里云EPEL镜像：https://mirrors.aliyun.com/epel/
- Rocky Linux官方下载：https://download.rockylinux.org/
- ELevate项目仓库：https://repo.almalinux.org/elevate/

**技术支持**：
- Rocky Linux论坛：https://forums.rockylinux.org/
- CentOS论坛：https://forums.centos.org/
- Server Fault社区：https://serverfault.com/
- Red Hat客户门户：https://access.redhat.com/

**故障排查资源**：
- GRUB修复指南：https://phoenixnap.com/kb/grub-rescue
- SSH隧道配置：https://www.ssh.com/academy/ssh/tunneling
- GRUB救援命令：https://linuxhint.com/grub_rescue_commands_centos/
- Linux Foundation GRUB指南：https://www.linuxfoundation.org/blog/blog/classic-sysadmin-how-to-rescue-a-non-booting-grub-2-on-linux

---

*注：本文档整合了官方指南和社区最佳实践，包含了无网络环境的完整解决方案，适用于 2025 年的系统迁移需求。如遇到GRUB引导问题，建议优先考虑使用Rocky Linux 8救援盘进行修复，或在测试环境中先验证升级流程的完整性。*