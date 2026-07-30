---
title: "VPS——海外便宜云服务器怎么选"
date: "2026-07-23"
last_modified_at: 2026-07-25
tags: [vps, cloud, server, oracle-cloud, google-cloud, linode, aws-lightsail, digitalocean, hetzner, web-techniques]
description: "海外云服务器VPS选购指南：Oracle免费方案、Google e2-micro、Akamai/Linode等主流方案对比，附配置推荐与安全设置"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/4K_1080P_compressed/074906xF9k4.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/4K_1080P_compressed/074906xF9k4.jpg"
author: Xufan Gao
lang: zh-CN
---

# VPS——海外便宜云服务器怎么选

> **使用场景**：找一个海外云服务器，北美节点，能访问 Codex，配置尽量低（1 核 1 G Linux 档），长期使用。CPU 要求很低，本质上是**长期在线的 Linux 跳板/开发机**，用 SSH、VS Code Remote、Codex CLI，不跑重计算。但 **1 GB 内存只够纯终端；同时跑 VS Code Server、Codex、Node/Python 和语言服务器时容易 OOM**。建议至少 2 GB；坚持 1 GB 就配 2–4 GB swap。

> **注意**：北美（美国、加拿大）属于 OpenAI 支持地区，但“服务器位于支持地区”并不等于一定能用 Codex——还取决于你的 OpenAI 账号、套餐、验证状态，以及该云厂商 IP 是否被风控。详见 [OpenAI 支持国家列表](https://developers.openai.com/api/docs/supported-countries)。

## 结论先给

按学生党长期使用的优先级：

1. **AWS 12 个月 Free Tier（t2.micro）：现有账号优先用完**
2. **Oracle Cloud Always Free A1：首选长期免费方案**
3. **Google Cloud e2-micro：最适合免费备用机**
4. **Akamai Cloud / Linode：最推荐的廉价付费方案，约 $5/月**
5. **DigitalOcean：约 $6/月，简单稳定**
6. **AWS Lightsail：1 GB 公网 IPv4 约 $7/月**
7. **Hetzner：性价比高，但美国节点价格通常不如欧洲节点便宜**

| 方案 | 大致配置 | 长期费用 | 北美节点 | 适合程度 | 主要问题 |
| --- | --- | --- | --- | --- | --- |
| **AWS Free Tier（t2.micro）** | t2.micro，1 vCPU、1 GB RAM、8 GB EBS、30 GB 磁盘、15 GB 出流量（12 个月免费） | **$0**（12 个月内） | 全球多个区域 | **现有账号优先用完** | 12 个月后按量计费；t2.micro 性能偏低；免费仅 t2.micro 实例类型 |
| Oracle Cloud Always Free A1 | ARM Ampere A1，2 OCPU、12 GB RAM（可配置为 1–2 台 VM）；外加 2 台 AMD 小实例 | **$0** | 有美国区域 | **最推荐** | 注册、容量和账户风控较严格 |
| Google Cloud Free Tier | e2-micro，2 vCPU、约 1 GB RAM、30 GB 磁盘 | **$0**（免费额度内） | 美国 3 个区域 | 推荐作备用 | 内存偏紧，免费出站流量仅 1 GB/月；超出免费额度约 $6.11/月 |
| Akamai Cloud / Linode | Shared CPU，最低约 1 GB RAM | **$5/月起** | 美国、加拿大 | **最推荐付费** | 1 GB 跑 VS Code 偏紧 |
| DigitalOcean | Basic Droplet，1 vCPU、1 GB RAM | 通常约 **$6/月** | 美国、加拿大 | 简单省心 | 性价比一般 |
| AWS Lightsail | 2 vCPU、1 GB RAM、40 GB 磁盘、2 TB 流量 | **$7/月** | 多个美国区域 | 稳定易用 | 比同配置小厂贵 |
| Hetzner Cloud | 通常从 2 vCPU、4 GB RAM 起 | 欧洲低价档约 €5.99 起；美国需看控制台 | 弗吉尼亚、俄勒冈 | 需要更多内存时划算 | 美国配置/价格与欧洲不同 |

## AWS 12 个月 Free Tier：现有账号先用完

AWS 免费套餐每月提供 **t2.micro 实例（1 vCPU、1 GB RAM、8 GB EBS）** 共 750 小时，外加 30 GB 标准 EBS 磁盘、15 GB 出流量，**免费 12 个月**。详见 [AWS Free Tier 文档](https://aws.amazon.com/cn/free/)。

t2.micro 性能偏低，跑 VS Code Remote + Codex 会明显局促，但当纯 SSH 跳板或 Codex CLI 跳板完全够用。如果你的 AWS 账号还在免费期内，**先把它用完**——这是零成本的现成资源，不用白不用。

### 查看 Free Tier 剩余时间

AWS 没有公开的 API 直接返回 Free Tier 到期日期，需要从控制台或通过 Cost Explorer 估算。

**方法一：控制台直接查看（最准确）**

登录 AWS Console，进入 **Billing and Cost Management** → **Free Tier** 页面，可以查看各项服务剩余免费额度：

- 找到 **Amazon Elastic Compute Cloud**，看 `BoxUsage:freetier.micro` 的当前使用量和总免费量
- **EC2** 的 750 小时/月是从账号激活起连续 12 个月；750 小时/月 ≈ 整月运行（30 天 × 24 小时 = 720 小时，留 30 小时余量）
- **EBS 磁盘**：30 GB-Mo（12 个月内累计）
- **出流量**：15 GB/月

**方法二：通过 AWS CLI 查询**

```bash
# 查询当前账号的 Free Tier 使用情况（需 cost-explorer 权限）
aws ce get-cost-and-usage \
  --time-period Start=$(date -u +%Y-%m-01),End=$(date -u +%F) \
  --granularity MONTHLY \
  --metrics UnblendedCost \
  --query 'ResultsByTime[0].Total.UnblendedCost' \
  --output table
```

返回结果中 `Estimated: True` 说明是估算值，实际金额可能有延迟。如果本月费用接近 0，说明 Free Tier 还在抵扣。

```bash
# 查看每日成本明细
aws ce get-cost-and-usage \
  --time-period Start=2026-07-25,End=2026-07-26 \
  --granularity DAILY \
  --metrics UnblendedCost \
  --group-by Type=DIMENSION,Key=SERVICE \
  --output table
```

返回结果中 `Amazon Elastic Compute Cloud - Compute` 为 0 说明 EC2 费用在免费额度内。

```bash
# 查看账号创建时间
aws account get-account-information \
  --region us-east-2 \
  --query 'AccountInformation[0].{Name:AccountName,Created:AccountCreatedDate,State:AccountState}' \
  --output table
```

**方法三：根据账号创建时间估算**

AWS 账号创建时间从 `AccountCreatedDate` 获取（UTC）。12 个月免费期从账号激活之日起计算。例如：

```text
AccountCreatedDate = 2025-07-07T16:40:48+00:00
→ 12 个月免费期截止约 2026-07-08（UTC）
```

注意：**账号创建日期 ≠ 实际激活日期**，实际 Free Tier 开始时间可能略晚 1-2 天。

### Free Tier 到期后会发生什么

12 个月免费期结束后，继续使用 t2.micro 实例将按标准费率计费（约 **$0.0116/小时 ≈ $8.50/月**）。**到期前建议的操作**：

- 在免费期结束前 1-2 个月主动关闭或销毁实例，避免意外计费
- 或提前迁移到 Oracle Always Free、Akamai / Linode 等免费或低价方案
- 关闭实例后**记得删除 EBS 快照和弹性 IP**，否则可能产生额外费用

### AWS Free Tier 版本差异

> **注意**：2025 年 7 月 15 日之后注册的 AWS 账号使用的是新版 Free Tier，免费额度和使用规则与旧版不同。旧版账号（2025-07-15 之前注册）仍适用 12 个月 t2.micro 免费规则；新版账号免费额度通过积分（credits）方式分配，规则略有不同。详见 [AWS Free Tier FAQ](https://aws.amazon.com/cn/free/faq/)。

旧版账号可以通过 AWS CLI 确认：

```bash
# 新版账号有 account-plan 数据；旧版账号返回 ResourceNotFoundException 是正常的
aws freetier get-account-plan-state \
  --region us-east-1 \
  --query '{Plan:accountPlanType, Status:accountPlanStatus, Expires:accountPlanExpirationDate}' \
  --output table
```

### 当前使用量示例

以下是某账号 2026 年 7 月 27 日的 Free Tier 使用情况截图：

> Amazon Virtual Private Cloud：750 小时，已用 621 小时（约 82.8%）
> Amazon Elastic Compute Cloud（t2.micro）：750 小时，已用 620 小时（约 82.6%）
> EBS 磁盘：30 GB-Mo，已用 7 GB-Mo
> 出流量：100 GB，已用 4 GB

该账号创建于 2025-07-07，Free Tier 理论上已于 2026-07-08 左右到期。当前账单仍显示 USD 0.00 可能是因为：

1. AWS 账单有延迟，Cost Explorer 中的数据可能标注 `Estimated: True`，实际费用需要 1-2 天才完全生成；
2. Free Tier 在部分结算周期内仍被计算为"有效"，直到月底才最终应用；
3. 账号实际激活时间可能晚于 `AccountCreatedDate` 1-2 天。

**最安全的做法**：假设免费期已到期，当前显示 0 美元不证明接下来不会产生费用。建议在免费期结束前主动关闭实例。

## 第一推荐：Oracle Cloud Always Free

> Oracle 官方文档显示，Always Free 的 Ampere A1 免费额度为每月 **1,500 OCPU-hours 和 9,000 GB-hours**，等效于持续运行 **2 OCPU + 12 GB RAM**，可配置为 1 台或 2 台 VM。此外还提供最多两台 AMD `VM.Standard.E2.1.Micro` 免费实例（每台 1/8 OCPU + 1 GB 内存）。总磁盘空间 200 GB，月出流量 10 TB。详见 [Oracle Always Free 官方文档](https://docs.oracle.com/en-us/iaas/Content/FreeTier/freetier_topic-Always_Free_Resources.htm)。

> Oracle 另外提供 **300 美元 Free Trial credit**，30 天内可以体验所有 Oracle 云服务，但这个和永久免费的 Always Free 是两套独立资源。详见 [Oracle Free FAQ](https://www.oracle.com/cloud/free/faq/)。

对于 Codex 跳板用途，这个免费资源完全足够。推荐配置（只占一半额度，留有余量）：

| 项目 | 配置 |
| --- | --- |
| 实例类型 | VM.Standard.A1.Flex |
| OCPU | 1（额度中只用 1 个） |
| RAM | 6 GB |
| 系统 | Ubuntu ARM64 |
| 区域 | 北美区域（Phoenix 等热门地区配额充足） |

**优点**：长期费用为零；12 GB 内存足够 VS Code Remote、Codex、Node、Python；ARM64 上普通 Python、Node.js、Git、Docker 基本都能用；适合长期挂 SSH、代理、博客构建和轻量开发。

**缺点**：注册时经常要求真实银行卡；免费 A1 实例有时显示容量不足（特别是热门区域）；Home Region 创建后不能更换；免费账户如果长期极低负载，存在资源回收或账户审核风险；某些仅提供 x86_64 二进制的软件在 ARM 上不能直接运行。

**建议**：选美国西部或东部区域，创建 **1 OCPU + 6 GB RAM**——这已经足够，给免费额度留有余量。

> Oracle 免费额度官方仍标为无限期可用的 Always Free 资源，而不是仅 30 天试用。30 天 / $300 是另一套 [Free Trial credit](https://www.oracle.com/cloud/free/faq/)，不要混淆。

**注册注意事项**：Oracle 免费实例需要保活，否则可能被释放（释放后对应区域无配额就尴尬了）。Oracle 也随时可能 ban 号（俗称"杀龟"），有时候新号活不过一个月，有时候几年的老号突然就没了。建议只当测试机，不要投放生产环境。注册时不需要代理，地址建议与信用卡账单地址一致，推荐手机 + 5G + Chrome 无痕模式 + 全中文地址。信用卡验证会扣几毛钱然后冲正，二次验证也类似，需要留几块钱。

## 第二推荐：Google Cloud Free Tier

Google Cloud 免费层每月提供一台非抢占式 `e2-micro`，可选美国俄勒冈 `us-west1`、爱荷华 `us-central1`、南卡罗来纳 `us-east1` 三个区域，附带 30 GB 标准持久磁盘，以及每月 1 GB **从北美传出**的免费出站流量（不计入中国和澳大利亚）。详见 [Google Cloud Free Tier 文档](https://docs.cloud.google.com/free/docs/free-cloud-features)。

**免费时长限制**：每月 744 小时（约等于整月不停机），受时长上限而非实例数量上限约束。所有区域合并计算，不超过 744 小时，整月持续运行一台 `e2-micro` 完全没问题。

**出站流量才是真正的限制**。1 GB/月 ≈ 每天 33 MB，跑一个代理或 SOCKS 跳板，一个用户刷个网页就爆了。所以 e2-micro 作为中转站或代理用不了——不是性能不够，是流量额度太小。这点对"用量较大"的用户尤其要避开。

它更适合：纯 SSH、Codex CLI、Git、简单脚本。

不太适合：VS Code Remote + Python language server + Node language server、大型 Jekyll build、Docker 多容器、任何需要大量出站流量的场景。

因为 e2-micro 只有约 1 GB 内存，CPU 为 2 vCPU，性能偏低。创建后立即加 2 GB swap：

```bash
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

**需要注意**：Google Cloud 必须绑定付款方式，免费额度外的磁盘、快照、静态公网 IPv4 或额外流量可能产生费用。建议设置：

```text
Billing budget：$1
Billing alert：50%、90%、100%
```

预算提醒通常不会自动停机，只负责通知。

> **超出免费额度后的费用**：e2-micro 按需定价约 **$6.11/月**（2 vCPU、1 GB RAM）。3 年 CUD 约 $3.30/月。所以免费额度外的磁盘、静态公网 IPv4 或额外流量才是真正产生费用的部分——e2-micro 本身在免费额度内是 $0，但一旦超出就要按这个价格计费。详见 [Google Compute Engine 定价](https://cloud.google.com/products/compute/pricing/general-purpose#e2-shared-core-machine-types)。

## 最推荐的付费方案：Akamai Cloud / Linode

Akamai Cloud（原 Linode）的 Shared CPU 最低档官方标价为 **$5/月起**，资源采用固定月价，CPU、RAM、存储和传输通常打包计费。详见 [Akamai Shared CPU 页面](https://www.linode.com/products/shared/)。

**优点**：正规大厂，长期价格比较稳定；后台和文档比廉价年付 VPS 完善；美国节点多；不容易出现"促销到期后翻几倍"的情况；支持快照、防火墙、反向 DNS；删除机器后计费容易理解。

若最低档仍是 1 GB，建议先买最低档试用并加 swap：

```bash
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

按内存档位的实际体验：

| 内存档位 | 适合的场景 | 主要问题 |
| --- | --- | --- |
| 1 GB | 纯 SSH、Codex CLI、Git、简单脚本 | 跑 VS Code Server 偏紧，需要 swap |
| 2 GB | VS Code Remote + Python/Node language server、Jekyll build | 性价比最高的档位 |
| 4 GB+ | Docker 多容器、跑多个服务、长时编译 | 价格明显上升 |

## DigitalOcean：最简单，但不是最低价

DigitalOcean 的 Basic Droplet 是典型的学生/开发者入门 VPS。Droplet 支持按秒计费，存在月度价格上限；**关闭机器但不销毁仍继续收费**，因为资源仍被保留。详见 [DigitalOcean 定价文档](https://docs.digitalocean.com/products/droplets/details/pricing/)。

常见最低档是 **1 vCPU、1 GB RAM、约 25 GB SSD、约 $6/月**。

北美有纽约、旧金山、多伦多、亚特兰大等区域。详见 [DigitalOcean 区域文档](https://docs.digitalocean.com/products/droplets/details/availability/)。

**适合**：创建 Ubuntu 很快、SSH key 配置简单、网络和控制台稳定、文档清晰、VS Code Remote 兼容性通常没问题。

**缺点**：同样的钱，配置不如 Oracle Free 或一些欧洲 VPS。

## AWS Lightsail：稳定但略贵

AWS Lightsail 当前官方 Linux 公网 IPv4 套餐：

| 套餐 | 月费 | 配置 |
| --- | --- | --- |
| 0.5 GB | $5/月 | 0.5 GB RAM、40 GB 磁盘、2 TB 流量 |
| 1 GB | $7/月 | 2 vCPU、40 GB 磁盘、2 TB 流量 |
| 2 GB | $12/月 | 2 vCPU、80 GB 磁盘、5 TB 流量 |

详见 [AWS Lightsail 文档](https://docs.aws.amazon.com/lightsail/latest/userguide/amazon-lightsail-bundles.html)。

**优势**：AWS 网络与账户体系成熟、静态 IP/快照/防火墙操作简单、北美区域充足、价格固定（账单比直接用 EC2 容易理解）。

但对学生党来说，\$7/月并不算特别便宜。IPv6-only 版 1 GB 是 $5/月，但不建议作为第一台服务器——本地网络、学校网络和某些服务对纯 IPv6 的支持可能不完整。

## Hetzner 是否值得选

Hetzner 在美国有 Ashburn, Virginia 和 Hillsboro, Oregon 两个节点。详见 [Hetzner 区域文档](https://docs.hetzner.com/cloud/general/locations/)。

欧洲区的成本优化型配置很强，例如 2 vCPU、4 GB、40 GB 的入门档。官方当前页面显示成本优化系列约 €5.99/月起，但这类资源容量有限，且具体型号是否能在美国区创建需要以控制台为准。详见 [Hetzner 成本优化页面](https://www.hetzner.com/cloud/cost-optimized)。

**Hetzner 适合**：需要 4 GB RAM、跑 VS Code + Docker + Jekyll、愿意支付约 €6–12/月。

但如果你只需要一个 1 GB Codex 跳板，它不一定是最便宜的北美选项。

## $4/月以内的低价 VPS 候选

下面这份清单来自 ChatGPT 联网查询的整理，截到 **2026 年 7 月**。低价套餐库存、续费价格和地区经常变化，下单页仍需再确认一次。

| 服务商 | 折算月费 | 付款方式 | 主要配置 | IPv4 | 主要限制 |
| --- | --- | --- | --- | --- | --- |
| **GreenCloud Budget KVM** | **$2.08/月** | $25/年 | 2 vCPU、4 GB RAM、35 GB SSD/NVMe、4 TB 流量 | 有 | 年付；无退款；部分地区缺货 |
| **CloudCone SSD VPS 3** | **$1.79/月** | $21.59/年 | 3 vCPU、3 GB RAM、41 GB SSD、3 TB 流量 | 有 | 年付；目前页面显示 Missouri |
| **RackNerd Specials** | **$3.00/月** | $35.99/年 | 2 vCPU、2 GB RAM、35 GB SSD、5 TB 流量 | 有 | 年付；共享 CPU，地区取决于库存 |
| **IONOS 入门 VPS** | **$2/月** | 3 年合同 | 1 vCPU、1 GB RAM、10 GB NVMe | 有 | 需要较长合同，资源较小 |
| **ServerHost** | **$3/月起** | 月付 | 官网称 VPS 从 $3/月起，流量不计量 | 需在订单页确认 | 官网公开搜索结果未完整列出入门配置 |
| **Vultr Cloud Compute** | **$3.50/月左右** | 按小时/月付 | 入门共享 CPU 实例 | $3.50 档通常有；$2.50 档仅 IPv6 | 低配内存较小；$2.50 套餐限 IPv6 |
| **AWS Lightsail** | **$3.50/月** | 按小时，月费封顶 | 2 vCPU、512 MB RAM、20 GB SSD、1 TB 流量 | **无，仅 IPv6** | 512 MB 很小；有 IPv4 的版本超过 $4 |
| **DigitalOcean** | **$4/月起** | 按小时/月付 | 基础 Droplet | 通常有 | 卡在预算上限，入门资源较少 |

### GreenCloud：性价比最高，但要挑有货地区

GreenCloud 的 Budget KVM 页面目前显示 **4 GB RAM、2 核、35 GB SSD/NVMe、4 TB 月流量、1 个 IPv4**，年付 $25 折合约 $2.08/月。详见 [GreenCloud 购物车](https://greencloudvps.com/billing/cart.php)。

当前 Jacksonville、Kansas City、Buffalo、Coventry、Amsterdam、Frankfurt 等部分 4 GB 或 8 GB 套餐仍显示有库存；洛杉矶、圣何塞和一些亚洲节点的低价档经常缺货。该系列明确写着**不退款**。

对部署 **Sub2API + PostgreSQL + Redis** 的场景来说，这一档配置比较合适——4 GB RAM 够 Nginx、Sub2API、PostgreSQL、Redis、Docker、少量日志和监控同时运行。但低价共享 CPU 的持续性能可能波动，不适合重型数据库或高并发 API。

### CloudCone：纸面性价比高，先验证库存和续费

CloudCone 的 [SSD VPS 3](https://cloudcone.com/vps/) 当前列出 **3 vCPU、3 GB RAM、41 GB SSD、3 TB 月流量、1 IPv4**，年付 $21.59 折合约 $1.79/月，地区 Missouri。下单前需要确认：

- 是否真的可以直接购买，而不是旧套餐页面；
- 次年是否仍按相同价格续费；
- IPv4 是否包含在最终订单中；
- 是否有备份、快照额外费用。

如果续费同价，它比 VoyraCloud $4.50/月便宜很多。

### RackNerd：配置和稳定性之间较均衡

[RackNerd Specials](https://www.racknerd.com/specials/) 页面当前有 **2 GB RAM、2 vCPU、35 GB SSD、5 TB 月流量、1 IPv4**，年付 $35.99 折合约 $3/月。另有 1 GB 套餐年付 $21.99。

RackNerd 还会不定期推出更低的活动套餐，例如 2026 年促销曾出现约 $10–19/年的 1–2.5 GB 套餐，但这类活动库存和入口不稳定，不能当作长期固定价格。详见 [LowEndBox 报道](https://lowendbox.com/blog/racknerd-ranks-90-on-2026-inc-regionals-pacific-list-3-years-in-a-row-check-out-their-latest-vps-deals/)。

相对而言，官网 Specials 的 $35.99/年更适合作为长期预算依据。

### IONOS：正规大厂，但 $2 套餐要签三年

[IONOS 官网](https://www.ionos.com/servers/cheap-vps) 目前宣传的 $2/月 VPS 是 **1 vCPU、1 GB RAM、10 GB NVMe，需要 3 年期限**。

**优点**：大型服务商，基础设施和账号体系比 LowEnd 商家规范，独立 IPv4 适合轻量反代、VPN、监控和小型静态服务。

**缺点**：三年合同，10 GB 磁盘很小，1 GB 内存不适合同时跑 PostgreSQL、Redis 和多个 Docker 容器，促销套餐与续费套餐需要区分。

IONOS 的另一类 [Cloud VPS 页面](https://www.ionos.com/servers/cloud-vps) 存在 "$2/月，前三个月" 的促销，但之后恢复到 $5/月，因此不能算长期不超过 $4。

### ServerHost：真正的 $3 月付候选

[ServerHost 官网](https://serverhost.com/) 表示 VPS 从 **$3/月**起，标注所有套餐提供不计量流量。

这类月付套餐的优势是不需要一次付一年、不满意时迁移成本低、适合短期测试。

但官网搜索结果没有完整展示 $3 套餐的 RAM、磁盘和 IPv4 条件。下单前必须核对是否包含独立 IPv4、实际端口速率、CPU 是否严重限频、可用地区、是否有安装费、是否允许 Docker、代理或 API 服务。

### Vultr：$3.50 可能有 IPv4，使用灵活

[Vultr FAQ](https://www.vultr.com/resources/faq/) 说明 **$2.50 套餐是 IPv6-only sandbox**，每个账户最多创建两个；$3.50 套餐与其类似但定位更接近正常入门实例。

Vultr 优点是按小时计费、地区多、API 成熟、快照/重装/防火墙方便。但这一档通常只有约 512 MB RAM，跑 Sub2API 全套服务不现实。更适合 Nginx 反代、WireGuard、frp、简单探针、极轻量网页这类场景。

### AWS Lightsail：$3.50 但仅 IPv6

[AWS Lightsail $3.50/月](https://aws.amazon.com/lightsail/pricing/) 的 Linux 套餐是 **512 MB RAM、2 vCPU、20 GB SSD、1 TB 流量，仅 IPv6**。

带公网 IPv4 的 Lightsail Linux 套餐目前从约 $5/月开始，超过 $4 预算。如果客户端、域名解析、反向代理或上游 API 仍依赖 IPv4，不建议购买这档。可用 [Lightsail 定价计算器](https://cloudburn.io/tools/amazon-lightsail-pricing-calculator) 复核。

### DigitalOcean：正好 $4，但性能价格比一般

[DigitalOcean 官网](https://www.digitalocean.com/solutions/vps-hosting) 仍标注 VPS/Droplet 从 **$4/月**起。优势是平台成熟、文档多、按小时计费；但 $4 预算下资源明显少于 GreenCloud、CloudCone 或 RackNerd。买它更多是为**平台可靠性、API、网络、文档、快照和生态**付费，而不是 RAM/价格比。

### VoyraCloud 是否偏贵

[VoyraCloud Cloud VPS](https://www.voyracloud.com/pricing) 当前最低是 **$4.50/月**，超预算 $0.50，纸面配置和同价位几家比没有明显优势。页面还列出 Residential IP VPS 从 $9/月、Windows VPS 从 $12/月、更换 IPv4 需要 $5/次。

如果需要香港、日本、新加坡或特定网络质量，$4.50 未必绝对不合理；如果只是部署美国节点的普通 API 服务，则没有明显购买理由。

### 按 Sub2API 用途的推荐

前面已经部署了 PostgreSQL、Redis 和 Sub2API，对内存和磁盘有最低要求：

| 要求 | 说明 |
| --- | --- |
| RAM | 2 GB 起步，推荐 4 GB |
| 磁盘 | 20 GB 起步，推荐 35 GB 以上 |
| IPv4 | 必须有 |
| 虚拟化 | KVM |
| 系统 | Ubuntu 22.04/24.04 |

按这个门槛筛选的优先级：

1. **GreenCloud 4 GB / $25 年付**：配置最匹配，确认地区、库存和不退款条款
2. **CloudCone 3 GB / $21.59 年付**：最便宜，先确认订单和续费价
3. **RackNerd 2 GB / $35.99 年付**：更保守，配置够用，但 PostgreSQL 内存要调小
4. **ServerHost $3 月付**：不想年付时值得测试
5. **IONOS $2**：适合长期轻量服务，但 1 GB 和三年合同限制较大

不建议用于 Sub2API 全套部署：**AWS Lightsail $3.50** 仅 IPv6 且只有 512 MB；**Vultr $2.50** 仅 IPv6；**DigitalOcean $4** 资源太少；**VoyraCloud $4.50** 超预算且纸面性价比一般。

购买年付 LowEnd VPS 前，最好先查清楚**续费是否同价、是否允许退款、IPv4 是否包含、服务器地区以及商家是否允许代理/API 转发业务**。这些因素比宣传的 vCPU 数量更重要。

## 最适合你的选择

[Oracle 免费套餐申请](https://www.oracle.com/cn/cloud/free/)（https://www.oracle.com/cn/cloud/free/）。

Oracle 注册成功率取决于信用卡信息、网络环境（不要用代理）、地址与信用卡账单地址一致等。推荐手机 + 5G + Chrome 无痕模式 + 全中文地址。

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
