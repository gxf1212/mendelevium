---
title: "【笔记整理|2026-07】Claude Code 用第三方模型爆 API Error 400：上下文超限怎么办"
date: "2026-07-24"
last_modified_at: 2026-07-24
tags: [Claude-Code, API-error, context-limit, third-party-model, compact, technical-notes]
description: "Claude Code用SenseNova等第三方模型时，上下文超过模型上限，自动compact没触发，报API Error 400。分析原因，给出手动compact、环境变量配置和模型切换的完整方案"
image: "/assets/img/thumbnail/bricks.webp"
thumbnail: "/assets/img/thumbnail/bricks.webp"
author: Xufan Gao
lang: zh-CN
---

# 【笔记整理|2026-07】Claude Code 用第三方模型爆 API Error 400，上下文超限怎么办

前几天用 SenseNova 的 `sensenova-6.7-flash-lite` 跑一个长任务，突然弹出一个错误：

```text
API Error: 400
{
  "type": "error",
  "error": {
    "type": "invalid_request_error",
    "message": "the input prompt token len 303876 + max_new_tokens 32000 > 262144"
  }
}
```

**303,876 个输入 token 加上 32,000 个输出 token，合计 335,876，超过了模型 262,144 的上限**。第三方服务在请求进入模型前直接拒绝，Claude Code 根本来不及做任何处理。

这个错误的核心在于**Claude Code 对第三方模型的上下文能力判断，与第三方 API 实际限制不一致**。

## 为什么自动 compact 没触发

Claude Code 的自动 compact 是客户端行为：它需要先知道当前模型的真实上下文窗口，然后在接近限制时主动总结历史消息。

问题在于，使用 `ANTHROPIC_BASE_URL` 指向第三方 Anthropic-compatible API 时，Claude Code 未必能正确获取模型能力。官方仓库已有相关问题：第三方 Base URL 可能无法进入官方模型能力检测流程，随后回退到默认值或错误的上下文配置，导致 compact 触发阈值与实际模型不一致。

具体就是两种相反的情况：

- Claude Code 以为模型支持更长上下文，迟迟不 compact，但第三方模型实际只有 262,144 token；
- Claude Code 以为模型只有较短上下文，过早 compact。

我遇到的显然更接近第一种：**Claude Code 允许上下文增长到 303,876 token，但后端模型实际只能接收 262,144 token**。

这与官方 Claude API 的行为差异有关。Anthropic 较新的官方模型 API 在部分情况下允许 `input + max_tokens` 超过窗口，只在实际生成撞到窗口时停止；较早的模型或第三方兼容接口可能直接返回 400。所以同一套 Claude Code 在官方 Claude 模型上正常，换第三方模型后报错。

## 报错时先救场

当前会话里先手动执行 `/compact` 把上下文压缩一下：

```text
/compact
```

如果 `/compact` 本身也因为上下文太长而失败，就执行 `/clear` 清空当前会话，或者退出当前会话重新开启一个。

Claude Code 官方的建议是：继续当前任务用 `/compact`，切换任务用 `/clear`。自动 compact 应该在接近上下文上限时运行，官方也明确保留了这两个手动命令作为兜底。

在长任务中，不要等到报错再 compact。对于 262,144 token 的模型，可以在 180,000–210,000 token 左右主动执行一次：

```text
/context
/compact
```

## 更根本的解决方法

### 方法一：把上下文上限设成第三方模型的真实值

如果当前使用的 Claude Code 版本支持，可以在启动前设置环境变量：

```bash
export CLAUDE_CODE_MAX_CONTEXT_TOKENS="262144"
```

永久设置在 `.bashrc` 或 `.zshrc` 里加：

```bash
export CLAUDE_CODE_MAX_CONTEXT_TOKENS="262144"
```

然后关闭并重新打开终端及 Claude Code。

不同 Claude Code 版本对该变量的支持可能有差异。设置后通过以下命令确认 Claude Code 显示的上下文上限：

```text
/context
/status
```

若仍显示错误值，说明当前版本没有读取该配置，或者第三方代理覆盖了模型能力配置。

### 方法二：降低第三方代理传入的 max_new_tokens

当前每次请求预留了 32,000 个输出 token：

```text
输入 token + 32000 <= 262144
允许的最大输入 = 230144
```

所以输入达到约 230,000 token 后就已经危险。如果第三方网关允许配置，可以把最大输出降到 8,192：

```text
262144 - 8192 = 253952
```

但不要卡着极限使用。工具定义、系统提示、MCP 返回值等也占上下文，更稳妥的目标是让输入低于约 200,000 token。

这个参数通常不在 Claude Code 的普通设置里修改，而是在第三方代理、网关或其模型映射配置里调整。

### 方法三：确认模型名映射和上下文设置

第三方服务常把一个 Claude 模型名映射到其他模型，比如：

```bash
export ANTHROPIC_MODEL="sensenova-6.7-flash-lite"
```

Claude Code 可能按照"Claude Sonnet"的能力估计上下文，但服务端实际运行的是一个 256,000 token 模型。需要检查第三方服务的实际模型名称、上下文窗口、最大输出 token、Anthropic-compatible 接口如何映射 `max_tokens`、是否向 Claude Code 暴露模型能力、是否有网关侧自动截断或压缩功能。

不要只看 Claude Code 中的显示名称。

### 方法四：更新 Claude Code

先检查版本，然后更新：

```bash
claude --version
npm update -g @anthropic-ai/claude-code
```

第三方 provider、模型发现和上下文窗口处理一直在调整。官方文档也说明，通过 LLM gateway 时，Claude Code 未必能验证扩展上下文能力，需要显式选择或配置模型。[Claude Code LLM Gateway 文档](https://code.claude.com/docs/en/llm-gateway)：https://code.claude.com/docs/en/llm-gateway

## 上下文为什么会快速膨胀到 30 万 token

Claude Code 每一轮都会重新携带整个对话历史、读取过的文件及结果、`CLAUDE.md`、工具定义、shell 输出、MCP 返回、diff 和测试日志、以及当前输入。长调试会话中读大量文件、生成大量 diff 后，每轮都会继续带着这些内容。

常见的上下文膨胀来源包括：大体积日志、完整构建输出、`package-lock.json` 等锁文件、大型 JSON、重复读取整个项目、MCP 返回过多内容、过长的 `CLAUDE.md`、连续处理多个无关任务。

建议养成几个习惯：用 `/context` 随时查看各部分占用；完成一个独立任务后执行 `/clear`；长任务中阶段性执行 `/compact`；引用文件时让 Claude 按路径按需读取，不要直接把整个大文件粘贴进聊天。官方特别提醒，使用 `@文件` 可能把整个文件及相关 `CLAUDE.md` 注入上下文，仅写文件路径更利于按需读取。

## 实际配置：SenseNova + MiniMax 双模型切换

我用 SenseNova 和 MiniMax 两个第三方模型，各配了一个 shell 函数。SenseNova 需要额外的上下文限制环境变量，切换到 MiniMax 时要用 `unset` 删掉，避免环境变量污染。

SenseNova 函数（262,144 token 上限，compact 窗口设 220,000）：

```bash
function sensenova {
    export ANTHROPIC_BASE_URL="https://token.sensenova.cn"
    export ANTHROPIC_API_KEY=""
    export ANTHROPIC_AUTH_TOKEN="sk-xxx"
    export ANTHROPIC_MODEL="sensenova-6.7-flash-lite"

    # 让 Claude Code 按 220K 上下文容量计算自动 compact
    export CLAUDE_CODE_AUTO_COMPACT_WINDOW="220000"

    # 第三方模型检测异常时的兼容性覆盖
    export CLAUDE_CODE_MAX_CONTEXT_TOKENS="262144"

    claude
}
```

切回 MiniMax 时用 `unset` 删掉覆盖变量：

```bash
function minimaxcc {
    export ANTHROPIC_API_KEY=""
    export ANTHROPIC_AUTH_TOKEN="sk-cp-xxx"
    export ANTHROPIC_BASE_URL="https://api.minimaxi.com/anthropic"

    unset ANTHROPIC_MODEL
    unset CLAUDE_CODE_AUTO_COMPACT_WINDOW
    unset CLAUDE_CODE_MAX_CONTEXT_TOKENS

    claude
}
```

这里最好同时 `unset ANTHROPIC_MODEL`。否则先执行 `sensenova`，再执行 `minimaxcc` 时，shell 中仍然可能保留 `ANTHROPIC_MODEL=“sensenova-6.7-flash-lite”`，导致 MiniMax 网关收到错误的模型名称。

## 检查变量是否切换正确

执行 `sensenova` 后，检查环境变量：

```bash
env | grep -E 'ANTHROPIC|CLAUDE_CODE'
```

应该能看到 SenseNova 的四个变量都设好了。退出后执行 `minimaxcc`，再次检查，这三个覆盖变量应该消失。

进入 Claude Code 后还可以运行 `/context`、`/status`、`/model`，确认当前模型、上下文占用和模型映射。

## 总结

- **报错本质**：Claude Code、第三方兼容层和实际模型三者对上下文窗口的认知不一致。自动 compact 只有在 Claude Code 正确判断临界点时才会及时触发；一旦后端先按自己的限制拒绝请求，Claude Code 就没有机会完成这一轮调用。
- **报错时**：先 `/compact`，失败就 `/clear` 或开新会话。长任务中 180,000–210,000 token 时主动 compact。
- **配置时**：设 `CLAUDE_CODE_AUTO_COMPACT_WINDOW=220000`（最稳妥），同时 `CLAUDE_CODE_MAX_CONTEXT_TOKENS=262144` 做兼容；切模型用 `unset` 清理环境变量。
- **日常习惯**：`/context` 随时看占用，`/clear` 切任务，`@文件` 改用路径按需读取。

[Claude Code 环境变量文档](https://code.claude.com/docs/en/env-vars)：https://code.claude.com/docs/en/env-vars
