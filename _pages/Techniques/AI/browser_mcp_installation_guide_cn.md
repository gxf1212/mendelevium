---
title: "让 Claude Code 控制浏览器：Playwright MCP 完全配置指南"
date: "2025-10-17"
tags: ['tutorial', 'ai-tools', 'mcp', 'playwright', 'browser-automation', 'claude']
description: "想让 AI 直接帮你操作浏览器吗？**Model Context Protocol (MCP)** 让这一切成为现实。通过 MCP 服务器，Claude Code 可以像人类一样浏览网页、填写表单、截图、抓取数据，甚至生成自动化测试代码。"
image: "/assets/img/thumbnail/bricks.webp"
thumbnail: "/assets/img/La-Mancha.jpg"
author: Xufan Gao
lang: zh-CN
---

# 让 Claude Code 控制浏览器：Playwright MCP 完全配置指南

## 引言

想让 AI 直接帮你操作浏览器吗？**Model Context Protocol (MCP)** 让这一切成为现实。通过 MCP 服务器，Claude Code 可以像人类一样浏览网页、填写表单、截图、抓取数据，甚至生成自动化测试代码。

**Playwright MCP** 是微软官方推出的浏览器自动化 MCP 服务器，它采用**基于可访问性树的创新方法**，无需视觉模型即可让 LLM 理解网页结构。这意味着更快的响应速度、更低的资源消耗，以及更精准的页面交互。

本文将手把手教你如何在 Claude Code 中配置 Playwright MCP，让 AI 成为你的浏览器自动化助手。

## 什么是 MCP？

**Model Context Protocol (MCP)** 是 Anthropic 推出的开放协议，用于连接 AI 应用与外部数据源和工具。通过 MCP，LLM 可以：

- 访问文件系统、数据库、API
- 操作浏览器、执行代码
- 与 GitHub、Slack 等第三方服务集成

MCP 的设计理念是**标准化 AI 与工具的连接方式**，就像 USB 协议统一了设备连接标准一样。开发者只需实现一次 MCP 服务器，就能在所有支持 MCP 的 AI 应用中使用。

Playwright MCP 是 MCP 生态中最受欢迎的浏览器自动化工具之一，由微软官方维护，已被数千个项目使用。

## 实际应用场景

安装 Playwright MCP 后，你可以让 Claude Code 帮你：

**Web 开发调试**
- "访问我的本地开发服务器 localhost:3000 并截图"
- "检查页面控制台是否有错误信息"
- "点击登录按钮，填写测试账号并提交表单"

**数据抓取**
- "访问这个产品页面，提取所有商品标题和价格"
- "抓取这个表格的数据并整理成 CSV 格式"

**自动化测试**
- "生成这个登录流程的 Playwright 测试代码"
- "验证这个页面在不同屏幕尺寸下的布局"

**内容监控**
- "每天检查这个网站的首页内容变化"
- "监控竞品的价格更新"

## MCP 服务器对比

Claude Code 支持两种主流浏览器自动化 MCP 服务器：

- **Playwright MCP**（推荐）：微软官方出品，支持多浏览器（Chrome/Firefox/WebKit），无需图形界面，性能优异
- **Chrome DevTools MCP**：基于 Chrome DevTools Protocol，适合 Chrome 专用调试场景

安装后，只需在对话中提及浏览器操作（如"访问这个网址并截图"），Claude Code 会自动调用相应的 MCP 工具完成任务。

## 完整安装步骤（Ubuntu/Debian）

### 方案一：Playwright MCP（推荐）

```bash
# 1. 添加到 Claude Code（无头模式）
claude mcp add -s user playwright -- npx @playwright/mcp@latest --headless

# 2. 安装 Playwright 浏览器
npx playwright install chromium

# 3. 安装系统依赖
sudo apt-get update
sudo apt-get install -y \
    libnss3 libnspr4 libdbus-1-3 \
    libatk1.0-0 libatk-bridge2.0-0 \
    libcups2 libdrm2 libxkbcommon0 \
    libxcomposite1 libxdamage1 libxfixes3 \
    libxrandr2 libgbm1 libpango-1.0-0 \
    libcairo2 libasound2

# 4. 验证安装
npx playwright --version

# 5. 完成！现在可以在 Claude Code 中使用浏览器功能
```

**优点：**
- 无需图形界面（X Server）
- 支持多浏览器（Chrome、Firefox、WebKit）
- 系统依赖少
- 开箱即用

### 方案二：Chrome DevTools MCP（备选）

```bash
# 1. 添加到 Claude Code
claude mcp add chrome-devtools npx chrome-devtools-mcp@latest

# 2. 安装 Chrome 浏览器
wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
sudo apt install ./google-chrome-stable_current_amd64.deb

# 3. 安装 Puppeteer 系统依赖
sudo apt-get update
sudo apt-get install -y \
    ca-certificates fonts-liberation \
    libappindicator3-1 libasound2 \
    libatk-bridge2.0-0 libatk1.0-0 \
    libcairo2 libcups2 libdbus-1-3 \
    libgbm1 libglib2.0-0 libgtk-3-0 \
    libnspr4 libnss3 libpango-1.0-0 \
    libx11-6 libxcomposite1 libxdamage1 \
    libxext6 libxfixes3 libxrandr2 \
    libxrender1 libxss1 libxtst6 \
    xdg-utils wget

# 4. 如果无图形界面，安装 xvfb（虚拟显示）
sudo apt-get install -y xvfb

# 5. 验证安装
google-chrome --version

# 6. 完成！
```

安装完成后，你就可以开始让 AI 帮你自动化浏览器操作了！

**注意：**
- 需要更多系统依赖
- 在无图形界面的服务器上需要 xvfb
- 仅支持 Chrome/Chromium

## 使用方法

安装完成后，直接在对话中提及浏览器操作即可，例如：

```
你：请访问 http://localhost:8504 并截图
Claude：好的，我来访问这个地址... [自动调用 mcp__playwright__browser_navigate]

你：查看页面上的错误信息
Claude：我来检查控制台日志... [自动调用 mcp__playwright__browser_console_messages]

你：点击"Performance Analysis"标签
Claude：我来点击这个标签... [自动调用 mcp__playwright__browser_click]
```

Claude Code 会自动选择合适的 MCP 工具执行操作。

## 常见问题

### 1. Playwright 找不到浏览器

```bash
# 重新安装浏览器
npx playwright install --force chromium

# 或指定浏览器路径
export PLAYWRIGHT_BROWSERS_PATH=/path/to/browsers
npx playwright install
```

### 2. Chrome DevTools 报错："Missing X server"

这是因为服务器没有图形界面。解决方案：

```bash
# 方案 A：安装 xvfb（虚拟显示）
sudo apt-get install -y xvfb

# 方案 B：使用 Playwright（推荐）
# Playwright 默认无头模式，无需图形界面
claude mcp add -s user playwright -- npx @playwright/mcp@latest --headless
```

### 3. 权限错误

```bash
# 使用 sudo 安装系统依赖
sudo npx playwright install-deps

# 或修改 npm 全局目录权限
mkdir ~/.npm-global
npm config set prefix '~/.npm-global'
echo 'export PATH=~/.npm-global/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

### 4. 检查 MCP 是否安装成功

```bash
# 查看已安装的 MCP 服务器
claude mcp list

# 测试 Playwright
npx playwright --version

# 测试 Chrome
google-chrome --version
```

## 方案对比

| 特性 | Playwright MCP | Chrome DevTools MCP |
|------|---------------|---------------------|
| **安装难度** | 非常简单 | 中等 |
| **无头模式** | 默认支持 | 需要配置 |
| **多浏览器** | Chrome, Firefox, WebKit | 仅 Chrome |
| **系统依赖** | 少 | 多 |
| **需要 X Server** | 不需要 | 需要（或 xvfb） |
| **性能** | 快 | 中等 |
| **推荐场景** | 通用自动化、测试 | Chrome 专用调试 |

**推荐**：优先使用 Playwright MCP，特别是在无图形界面的服务器上。

## 参考资源

- Playwright MCP GitHub 仓库：https://github.com/microsoft/playwright-mcp
- Playwright 官方文档：https://playwright.dev
- Model Context Protocol 规范：https://modelcontextprotocol.io
- Claude Code MCP 文档：https://docs.claude.com/en/docs/claude-code
- Chrome DevTools Protocol：https://chromedevtools.github.io/devtools-protocol
- MCP Hub（发现更多 MCP 服务器）：https://mcphub.com

另外再推荐个小工具


cc相关的工具太多了，肯定学不完，随缘了。
