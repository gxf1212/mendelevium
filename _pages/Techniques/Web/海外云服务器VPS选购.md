---
title: "【笔记整理|2026-07】VPS——海外便宜云服务器怎么选"
date: "2026-07-23"
last_modified_at: 2026-07-23
tags: [vps, cloud, server, oracle-cloud, google-cloud, linode, aws-lightsail, digitalocean, hetzner, web-techniques]
description: "海外云服务器VPS选购指南：Oracle免费方案、Google e2-micro、Akamai/Linode等主流方案对比，附配置推荐与安全设置"
image: "/assets/img/4K_1080P_compressed/074906xF9k4.jpg"
thumbnail: "/assets/img/4K_1080P_compressed/074906xF9k4.jpg"
author: Xufan Gao
lang: zh-CN
---

# 【笔记整理|2026-07】VPS——海外便宜云服务器怎么选

> **使用场景**：找一个海外云服务器，北美节点，能访问 Codex，配置尽量低（1 核 1 G Linux 档），长期使用。CPU 要求很低，本质上是**长期在线的 Linux 跳板/开发机**，用 SSH、VS Code Remote、Codex CLI，不跑重计算。但 **1 GB 内存只够纯终端；同时跑 VS Code Server、Codex、Node/Python 和语言服务器时容易 OOM**。建议至少 2 GB；坚持 1 GB 就配 2–4 GB swap。

> **注意**：北美（美国、加拿大）属于 OpenAI 支持地区，但“服务器位于支持地区”并不等于一定能用 Codex——还取决于你的 OpenAI 账号、套餐、验证状态，以及该云厂商 IP 是否被风控。详见 [OpenAI 支持国家列表](https://developers.openai.com/api/docs/supported-countries)。

## 结论先给

按学生党长期使用的优先级：

1. **Oracle Cloud Always Free A1：首选免费方案**
2. **Google Cloud e2-micro：最适合免费备用机**
3. **Akamai Cloud / Linode：最推荐的廉价付费方案，约 $5/月**
4. **DigitalOcean：约 $6/月，简单稳定**
5. **AWS Lightsail：1 GB 公网 IPv4 约 $7/月**
6. **Hetzner：性价比高，但美国节点价格通常不如欧洲节点便宜**

| 方案 | 大致配置 | 长期费用 | 北美节点 | 适合程度 | 主要问题 |
| --- | --- | --- | --- | --- | --- |
| Oracle Cloud Always Free A1 | ARM，免费额度等效 4 OCPU、24 GB RAM（2025 年后已升级） | **$0** | 有美国区域 | **最推荐** | 注册、容量和账户风控较严格 |
| Google Cloud Free Tier | e2-micro，共享核、约 1 GB RAM、30 GB 磁盘 | **$0** | 美国 3 个区域 | 推荐作备用 | 内存偏紧，免费出站流量仅 1 GB/月 |
| Akamai Cloud / Linode | Shared CPU，最低约 1 GB RAM | **$5/月起** | 美国、加拿大 | **最推荐付费** | 1 GB 跑 VS Code 偏紧 |
| DigitalOcean | Basic Droplet，1 vCPU、1 GB RAM | 通常约 **$6/月** | 美国、加拿大 | 简单省心 | 性价比一般 |
| AWS Lightsail | 2 vCPU、1 GB RAM、40 GB 磁盘、2 TB 流量 | **$7/月** | 多个美国区域 | 稳定易用 | 比同配置小厂贵 |
| Hetzner Cloud | 通常从 2 vCPU、4 GB RAM 起 | 欧洲低价档约 €5.99 起；美国需看控制台 | 弗吉尼亚、俄勒冈 | 需要更多内存时划算 | 美国配置/价格与欧洲不同 |

> **重要提示**：价格、资源和免费政策都会变化，实际下单时应以控制台为准。以下信息截至 2026 年 7 月。

## 第一推荐：Oracle Cloud Always Free

> Oracle 在 **2025 年中**对 Always Free 额度进行了升级，A1 Ampere 免费额度从原来的 2 OCPU / 12 GB RAM 提升至 **4 OCPU / 24 GB RAM**（对应 3,000 OCPU-hours + 18,000 GB-hours）。此外仍提供最多两台 AMD `VM.Standard.E2.1.Micro` 免费实例。详见 [Oracle Always Free 官方文档](https://docs.oracle.com/en-us/iaas/Content/FreeTier/freetier_topic-Always_Free_Resources.htm)。

对于 Codex 跳板用途，这个免费资源完全足够：

```text
VM.Standard.A1.Flex
1–2 OCPU（4 个额度中只用 1–2 个）
6–12 GB RAM
Ubuntu ARM64
北美区域
```

**优点**：

- 长期费用为零；
- 12–24 GB 内存足够 VS Code Remote、Codex、Node、Python；
- ARM64 上普通 Python、Node.js、Git、Docker 基本都能用；
- 适合长期挂 SSH、代理、博客构建和轻量开发。

**缺点**：

- 注册时经常要求真实银行卡；
- 免费 A1 实例有时显示容量不足（特别是热门区域）；
- Home Region 创建后不能更换；
- 免费账户如果长期极低负载，存在资源回收或账户审核风险；
- 某些仅提供 x86_64 二进制的软件在 ARM 上不能直接运行。

**建议**：选美国西部或东部区域，创建 **1 OCPU + 6 GB RAM**——这已经足够，给免费额度留有余量。

> Oracle 免费额度官方仍标为无限期可用的 Always Free 资源，而不是仅 30 天试用。30 天 / $300 是另一套 [Free Trial credit](https://www.oracle.com/cloud/free/faq/)，不要混淆。

## 第二推荐：Google Cloud Free Tier

Google Cloud 免费层目前包括：

- 每月一台非抢占式 `e2-micro`
- 可选美国俄勒冈 `us-west1`、爱荷华 `us-central1`、南卡罗来纳 `us-east1`
- 30 GB 标准持久磁盘
- 每月 1 GB **从北美区域传出**的免费出站流量；前往中国和澳大利亚不计入这项免费范围。详见 [Google Cloud Free Tier 文档](https://docs.cloud.google.com/free/docs/free-cloud-features)

它更适合：

```text
纯 SSH
Codex CLI
Git
简单脚本
轻量反向代理
```

不太适合：

```text
VS Code Remote + Python language server + Node language server
大型 Jekyll build
Docker 多容器
```

因为 e2-micro 是共享核且只有约 1 GB 内存。创建后立即加 2 GB swap：

```bash
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

Google 的免费额度足以让一台 e2-micro 整月持续运行。详见 [Sustained Use Discounts 文档](https://docs.cloud.google.com/compute/docs/sustained-use-discounts)。

**需要注意**：Google Cloud 必须绑定付款方式，免费额度外的磁盘、快照、静态公网 IPv4 或额外流量可能产生费用。建议设置：

```text
Billing budget：$1
Billing alert：50%、90%、100%
```

预算提醒通常不会自动停机，只负责通知。

## 最推荐的付费方案：Akamai Cloud / Linode

Akamai Cloud（原 Linode）的 Shared CPU 最低档官方标价为 **$5/月起**，资源采用固定月价，CPU、RAM、存储和传输通常打包计费。详见 [Akamai Shared CPU 页面](https://www.linode.com/products/shared/)。

**优点**：

- 正规大厂，长期价格比较稳定；
- 后台和文档比廉价年付 VPS 完善；
- 美国节点多；
- 不容易出现"促销到期后翻几倍"的情况；
- 支持快照、防火墙、反向 DNS；
- 删除机器后计费容易理解。

若最低档仍是 1 GB，建议先买最低档试用并加 swap：

```bash
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

如果 VS Code 经常卡死或 language server 被杀，再升到 2 GB。

## DigitalOcean：最简单，但不是最低价

DigitalOcean 的 Basic Droplet 是典型的学生/开发者入门 VPS。Droplet 支持按秒计费，存在月度价格上限；**关闭机器但不销毁仍继续收费**，因为资源仍被保留。详见 [DigitalOcean 定价文档](https://docs.digitalocean.com/products/droplets/details/pricing/)。

常见最低可用档：

```text
1 vCPU
1 GB RAM
约 25 GB SSD
约 $6/月
```

北美有纽约、旧金山、多伦多、亚特兰大等区域。详见 [DigitalOcean 区域文档](https://docs.digitalocean.com/products/droplets/details/availability/)。

**适合**：创建 Ubuntu 很快、SSH key 配置简单、网络和控制台稳定、文档清晰、VS Code Remote 兼容性通常没问题。

**缺点**：同样的钱，配置不如 Oracle Free 或一些欧洲 VPS。

## AWS Lightsail：稳定但略贵

AWS Lightsail 当前官方 Linux 公网 IPv4 套餐：

- 0.5 GB：$5/月
- 1 GB：$7/月，2 vCPU、40 GB 磁盘、2 TB 流量
- 2 GB：$12/月。详见 [AWS Lightsail 文档](https://docs.aws.amazon.com/lightsail/latest/userguide/amazon-lightsail-bundles.html)

不要为了省 $2 选 0.5 GB：

```text
0.5 GB：基本不适合 VS Code Server
1 GB：勉强可用，需要 swap
2 GB：比较舒服
```

**优势**：AWS 网络与账户体系成熟、静态 IP/快照/防火墙操作简单、北美区域充足、价格固定（账单比直接用 EC2 容易理解）。

但对学生党来说，$7/月并不算特别便宜。IPv6-only 版 1 GB 是 $5/月，但不建议作为第一台服务器——本地网络、学校网络和某些服务对纯 IPv6 的支持可能不完整。

## Hetzner 是否值得选

Hetzner 在美国有：

- Ashburn, Virginia
- Hillsboro, Oregon。详见 [Hetzner 区域文档](https://docs.hetzner.com/cloud/general/locations/)

欧洲区的成本优化型配置很强，例如 2 vCPU、4 GB、40 GB 的入门档。官方当前页面显示成本优化系列约 €5.99/月起，但这类资源容量有限，且具体型号是否能在美国区创建需要以控制台为准。详见 [Hetzner 成本优化页面](https://www.hetzner.com/cloud/cost-optimized)。

**Hetzner 适合**：需要 4 GB RAM、跑 VS Code + Docker + Jekyll、愿意支付约 €6–12/月。

但如果你只需要一个 1 GB Codex 跳板，它不一定是最便宜的北美选项。

## 不太建议的"超低价年付 VPS"

确实有一些 LowEndTalk、RackNerd、ColoCrossing 系商家会提供：

```text
1 GB RAM
15–25 GB SSD
年付 $10–25
```

折合可能只有每月 $1–2。

**但我不建议把它作为唯一长期环境**，原因是：

- CPU 严重超售；
- 磁盘 I/O 波动；
- IP 段信誉可能较差；
- 云厂商共享 IP 段可能被 GitHub / OpenAI / Google 风控；
- 年付退款困难；
- 商家或套餐可能消失；
- 备份、快照、救援系统不完善。

用来做临时 SSH 跳板可以；存放唯一代码、密钥和研究数据不合适。

## 最适合你的选择

### 完全不想花钱

先试：

```text
Oracle Cloud Always Free
美国区域
A1 Flex
1 OCPU + 6 GB RAM
Ubuntu 24.04 ARM64
```

注册或容量失败，再用：

```text
Google Cloud e2-micro
us-west1 / us-central1 / us-east1
Ubuntu 24.04
30 GB standard disk
2 GB swap
```

### 愿意每月付 5 美元

选：

```text
Akamai Cloud / Linode
美国西海岸或东海岸
最低 Shared CPU
Ubuntu 24.04
```

这是我认为**价格、稳定性、长期可用性最均衡**的方案。

### 需要 VS Code 使用舒服

不要买 1 GB，直接找：

```text
2 vCPU
2–4 GB RAM
30–50 GB SSD
```

优先级：

```text
Oracle Free A1
Hetzner 4 GB
Akamai 2 GB
Lightsail 2 GB
DigitalOcean 2 GB
```

## 创建后最小安全配置

无论选哪家，建议立即执行：

```bash
# 更新系统
sudo apt update && sudo apt upgrade -y

# 新建普通管理用户
sudo adduser dev
sudo usermod -aG sudo dev

# 安装常用工具
sudo apt install -y git curl wget tmux htop fail2ban ufw

# 只开放 SSH
sudo ufw allow OpenSSH
sudo ufw enable

# 启用 fail2ban
sudo systemctl enable --now fail2ban
```

确认 SSH key 登录可用以后，再关闭密码登录：

```bash
sudoedit /etc/ssh/sshd_config
```

设置：

```text
PasswordAuthentication no
PermitRootLogin no
```

然后：

```bash
sudo systemctl restart ssh
```

**我的最终推荐：先申请 Oracle A1；申请不到就买 Akamai / Linode $5 档。** Google e2-micro 适合作为免费的备用跳板，但长期运行 VS Code Remote 会明显局促。
