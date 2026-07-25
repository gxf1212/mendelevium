---
title: "【笔记整理|2026-07】Sub2API：让 Claude Code、Codex CLI 都用一个中转 API"
date: "2026-07-25"
last_modified_at: 2026-07-26
tags: [sub2api, Claude-Code, Codex-CLI, API-gateway, Anthropic, OpenAI, technical-notes]
description: "Sub2API 一站式开源中转服务，34k stars。用一台 Linode 服务器搭好中转 API，Claude Code 和 Codex CLI 都能直接连，拼车共享还能摊薄订阅成本。附完整配置教程"
image: "/assets/img/thumbnail/bricks.webp"
thumbnail: "/assets/img/thumbnail/bricks.webp"
author: Xufan Gao
lang: zh-CN
---

# 【笔记整理|2026-07】Sub2API：让 Claude Code、Codex CLI 都用一个中转 API

白嫖了一台 Linode VPS，自带命令行真的巨难用，直到直接关掉防火墙了，ssh才能用了（）。

闲着想让它再干点活——把 Claude 和 OpenAI 的订阅接出来，统一中转给 Claude Code 和 Codex CLI 用，顺带拼车共享摊薄成本。

Claude Relay Service已经迁移，不好一键安装了。Sub2API 就是继承者。一个项目，34k stars，更新到 2026 年 7 月 25 日，issue 响应比较快。

## 这篇教程覆盖什么

下面按 Sub2API 官方文档和新手指引走一遍：先在 Linux 服务器上部署 Sub2API，启动后会引导到 Web 界面配置上游账号和 API Key；过程中需要先安装 PostgreSQL 和 Redis 两个依赖；账号、渠道、Key 这些配置都按新手指引填就行。下面把容易卡住的几步单独拎出来。

## Sub2API 是什么

Sub2API（[GitHub 仓库](https://github.com/Wei-Shaw/sub2api)）是一个一站式开源中转服务，把 Claude、OpenAI、Gemini、Grok 的订阅统一接入，支持拼车共享，原生工具能无缝使用。

对 Claude Code 和 Codex CLI 用户来说，它解决了一个真实存在的痛点：Claude Code 用 Anthropic Messages API，Codex CLI 用 OpenAI Responses API，路径在 `/v1/responses`，两边要分别配 `ANTHROPIC_BASE_URL` 和 `OPENAI_BASE_URL`，维护起来麻烦。Sub2API 把这两套 API 都封装到同一个 `http://xxx:8080` 后面，客户端只需要改一个基础地址。

## 安装

服务器装好 Linux（Ubuntu 24.04 比较稳妥），把端口开放或者临时关掉防火墙（SSH 不受影响）。

一行脚本搞定安装，然后启动和开机自启：

```bash
curl -sSL https://raw.githubusercontent.com/Wei-Shaw/sub2api/main/deploy/install.sh | sudo bash
sudo systemctl enable --now sub2api
```

验证一下端口是否正常响应：

```bash
curl http://你的IP:8080/
```

正常会返回 Sub2API 的欢迎页面，说明服务已经起来了。

## 数据库配置

Sub2API 启动后会自动引导进入 Web 配置页面，过程中需要 PostgreSQL 和 Redis 作为依赖。如果服务器上还没装，先一次性装好：

```bash
sudo apt update
sudo apt install -y postgresql postgresql-contrib redis
sudo systemctl enable --now postgresql redis-server
```

给默认 `postgres` 用户设置密码，避免明文写入命令历史：

```bash
sudo -u postgres psql
```

进入 psql 后修改密码，再退出：

```sql
\password postgres
\q
```

> 如果 Sub2API 和 PostgreSQL 都在同一台服务器上，保持 PostgreSQL 只监听本机回环即可，不要为了应用连接直接开放公网 5432 端口。

## API 端点

Sub2API 当前兼容以下等价入口，客户端任选一个就行：

| 工具 | 实际请求路径 | 说明 |
| --- | --- | --- |
| Codex CLI | `POST http://xxx:8080/v1/responses` | OpenAI Responses API，Codex 的公开入口 |
| Codex CLI（兼容） | `POST http://xxx:8080/responses` | Sub2API 兼容入口 |
| Codex CLI（兼容） | `POST http://xxx:8080/backend-api/codex/responses` | Sub2API 兼容入口 |
| Claude Code | `POST http://xxx:8080/v1/messages` | Anthropic Messages API，Claude Code 自动追加 |

Claude Code 只需要设基础地址：

```bash
export ANTHROPIC_BASE_URL="http://你的服务器IP:8080"
```

Codex CLI 同理：

```bash
export OPENAI_BASE_URL="http://你的服务器IP:8080"
```

> **注意**：Codex CLI 调用的是 OpenAI Responses API，路径是 `/v1/responses`，写成 `/openai/` 是错的。Claude Code 会自动追加 `/v1/messages`，基础地址不需要带路径后缀。

## Web 界面与账号配置

安装完成后访问 `http://服务器IP:8080`，看到的是 Sub2API 的 Web 界面。首页左侧是导航栏（账号、渠道、Key、订单），右侧是功能区域。配置好上游账号和 API Key 后，可以在"我的账户"里查看 API Key 列表，在"渠道管理"里查看各模型的连接状态。

Sub2API 是 Vue 3 + Go 后端的 Web 应用，登录方式目前只有注册邮箱验证。首页显示各模型（Claude、OpenAI、Gemini、Grok）的可用状态和延迟。

账号、渠道、Key 这些步骤直接跟着新手指引填就行。具体配置（包括上游账号如何接入、各模型分组怎么分配）可以参考这篇教程：[akokoi 的 Sub2API 配置教程](https://x.com/akokoi1/article/2041464436876345524)（https://x.com/akokoi1/article/2041464436876345524）。

## 小结

- **Sub2API**（[GitHub](https://github.com/Wei-Shaw/sub2api)）一站式中转 Claude、OpenAI、Gemini、Grok 订阅，34k stars，issue 响应比较快
- **部署**走官方安装脚本，一行命令带上开机自启，端口 8080 验证是否启动
- **依赖**PostgreSQL 和 Redis 都需要装，postgres 用户设密码，数据库保持本机回环访问
- **Claude Code** 设 `ANTHROPIC_BASE_URL` 指向 `http://xxx:8080`，自动追加 `/v1/messages`
- **Codex CLI** 设 `OPENAI_BASE_URL` 指向 `http://xxx:8080`，自动追加 `/v1/responses`
- **账号配置**按新手指引填即可，详细步骤参考 akokoi 那篇教程

