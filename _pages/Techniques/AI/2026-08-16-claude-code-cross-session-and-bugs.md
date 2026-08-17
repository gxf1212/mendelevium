---
title: "Claude Code 跨 session 通信、版本更新与常见报错"
date: "2026-08-16"
last_modified_at: 2026-08-16
tags: [claude-code, cross-session, peer-messaging, compact, debugging, ai-tools]
description: "Claude Code 跨 session 通信机制说明，以及版本更新失败、peer messaging 不可用、thinking block 解析错误、频繁 compact 等问题的排查记录"
image: "/assets/img/thumbnail_mine/wh-1kdv6v.jpg"
thumbnail: "/assets/img/thumbnail_mine/wh-1kdv6v.jpg"
author: Xufan Gao
lang: zh-CN
---

# Claude Code 跨 session 通信、版本更新与常见报错

## 1. Cross-session messaging 机制

整理自 [Claude Code 官方文档](https://code.claude.com/docs/en/cross-session-messaging)：https://code.claude.com/docs/en/cross-session-messaging，结合近期使用中遇到的三个实际问题。

Claude Code v2.1.224 起支持跨 session 通信，macOS 和 Linux 上只要版本满足就自动启用，不需要额外配置。核心机制如下：

- 一个 session 里的 Claude 可以给另一个 session 发送文本消息，但不能传递对话历史或文件。要转移整个上下文，应该用 resume 恢复 session
- 依赖两个工具：**ListAgents** 发现可达的 agent，**SendMessage** 按名称发送消息。同一个 SendMessage 工具还可以给同一 session 内的 subagent 和 agent team 成员发消息
- 消息在一个 Claude 和另一个 Claude 之间传递，不经过用户中转

**典型使用场景**：

- 一个 session 发现了影响另一个 session 工作的变更，自动通知对方
- 多个 session 在不同 worktree 处理同一仓库，一个提交后通知其他 session
- 让长时间运行的任务（如迁移、测试）向当前正在看的 session 回报进展
- 跨机器或跨网页端给其他 session 发消息

**消息投递状态**：接收方 session 会对每条消息做入站检查，结果有三种。

- **Delivered**：消息正常送达
- **Held**：消息被暂存，等你批准或后续设置变更后才会送达
- **Refused**：消息被直接丢弃

消息一旦送达，会计入用量，接收方的 Claude 可以同样方式回复发送方（跨机器的单向场景除外）。权限边界是 per-session 的，Claude 被禁止要求另一个 session 执行它自己 session 中被拒绝的操作。

**session 命名规则**：用 `/rename` 或 `--name` 设置的名称就是它接收消息的地址。如果不设置，Claude Code 会从工作目录名称自动生成，比如 `my-app` 目录下自动命名为 `my-app-3f`。用 `/list-agents` 命令可以看到当前可达的所有 session。

> cross-session messaging 适用于你自己启动和操控的独立 session。要继续另一个对话，用 resume；要协调 Claude 自己派生的 team，用 agent teams；要从一个地方监控多个 session，用 agent view；要从手机操控，用 Remote Control；要推送外部事件进 session，用 channels。

## 2. 版本更新踩坑

参考：[微信公众号文章](https://mp.weixin.qq.com/s/u8TREHQCWP73DNqFPsK-FA)：https://mp.weixin.qq.com/s/u8TREHQCWP73DNqFPsK-FA

打开终端跑 `claude update`，直接报错说无法从 npm registry 获取最新版本号。第一反应是代理问题，但用 `npm view @anthropic-ai/claude-code version` 单独探测同一个包，三秒就返回了 `2.1.226`——registry 明明是通的。

**根因**：`claude update` 自己那套检测最新版本号的逻辑坏了，跟网络无关。处理方式：

```bash
npm install -g @anthropic-ai/claude-code@latest
claude --version
```

前后不到 5 分钟，从 2.1.222 跳到了 2.1.226。

**注意**：更新后 `/home/gxf1212/.local/bin/claude` 可能还是旧版本的路径，直接删掉。确认当前用的是正确的：

```bash
which claude
# /home/gxf1212/.nvm/versions/node/v22.18.0/bin/claude
claude --version
# 2.1.233 (Claude Code)
```

## 3. 问题一：peer messaging 不可用

### 现象

尝试让一个 session 向另一个 session 发消息，但 SendMessage 工具只看见了当前 team 里的 agent，看不见其他 peer session。各种方式都搞不定。

### 根因排查

这台机器上 Claude Code 没有成功启用 peer messaging：

- `~/.claude/sessions/*.json` 里没有注册 `messagingSocketPath`
- `tengu_fleetview_peers` 是关的
- TIOCSTI 键盘注入也禁用了

### 结论

**Claude Code 2.1.233 本身有 cross-session peer messaging 能力**，官方 SDK 文档已经出现了 `cross-session peer` 相关定义。参考 [Claude Platform Docs](https://docs.anthropic.com/en/docs/claude-code/sdk/sdk-typescript)：https://docs.anthropic.com/en/docs/claude-code/sdk/sdk-typescript

当前更像是 **session 初始化时没有成功注册 peer messaging**，不是版本不支持。参考 [GitHub issue #42737](https://github.com/anthropics/claude-code/issues/42737)：https://github.com/anthropics/claude-code/issues/42737，SendMessage 确实存在过“未暴露、不可用、静默失败”的 bug。

### 排查命令

```bash
ls -lah ~/.claude/sessions/
ls -lah /tmp/cc-socks/
pgrep -af claude
```

检查各 session 是否有 `messagingSocketPath`：

```bash
jq . ~/.claude/sessions/*.json
```

如果某个 session 没有 `messagingSocketPath` 或对应 socket，**重启那个 session** 最值得先试。

进一步确认 discovery 是否正常：

```bash
claude agents --json
```

如果这里能看到那些 session，但各 session JSON 仍然没有 `messagingSocketPath`，说明 discovery 正常但 peer inbox 没 bind，SendMessage 自然发不出去。

## 4. 问题二：`i.thinking.length` 报错

### 现象

个别 session 执行任何操作都报同一个错：

```
Error: undefined is not an object (evaluating 'i.thinking.length')
Error during compaction: undefined is not an object (evaluating 'i.thinking.length')
```

普通 Bash 和 `/compact` 都触发。

### 根因

**Claude Code 自己的 thinking 块解析 / TUI bug**。Claude Code 在读取当前 session 历史时遇到了一个 `thinking` 块，但其 `thinking` 字段是 `undefined/null`，代码却直接做了 `.length`。参考 [GitHub issue #65240](https://github.com/anthropics/claude-code/issues/65240)：https://github.com/anthropics/claude-code/issues/65240，已有高度相似的报告：`thinking` 为空导致 TUI 报 `thinking.trim`、thinking summary 为空、compaction 被 thinking 块弄坏。

前面刚跑过 background agent / peer 相关操作，很可能是某个 agent 返回写进 session 时生成了异常 thinking block，后续所有需要遍历 history 的功能一起炸。

### 解决方案

```bash
claude --resume
```

选这个 session 看是否恢复，不行就新开 session。参考 [Reddit 讨论](https://www.reddit.com/r/ClaudeAI/comments/1pfd5sr/fix_claude_code_error_error_during_compaction/)：https://www.reddit.com/r/ClaudeAI/comments/1pfd5sr/fix_claude_code_error_error_during_compaction/，已有案例里重启或重新 resume 能绕过部分 compaction 异常。

想定位具体是哪条坏了：

```bash
grep -n '"type":"thinking"' ~/.claude/projects/*/*.jsonl | tail -20
```

找缺失 `thinking` 字段的块：

```bash
grep '"type":"thinking"' SESSION.jsonl | jq -c '
  select(
    .message.content[]? |
    .type == "thinking" and
    (.thinking == null)
  )
'
```

**结论**：不是 Markdown server 的问题，先救 Claude session。

## 5. 问题三：频繁 compact

### 现象

新版本 Claude Code 出现频繁 compact，“还没说两句就 compact 了”。

### 根因

Claude Code **v2.1.223** 开始改了 auto-compact 逻辑：

- 对**无法识别的 model ID**，不再允许它使用声明的大 context，而是按 Claude Code 内部“推测的 context window”强制 compact
- `CLAUDE_CODE_DISABLE_1M_CONTEXT` 现在会把所有原生 1M 模型都通过 auto-compaction 限制在 200K
- 社区已有大量报告“1M 模型却在约 76K / 200K / 400K 就 compact”

参考 [Claude Code Releases](https://github.com/anthropics/claude-code/releases)：https://github.com/anthropics/claude-code/releases 和 [GitHub issue #34332](https://github.com/anthropics/claude-code/issues/34332)：https://github.com/anthropics/claude-code/issues/34332

### 解决方案

如果在使用第三方模型：

```bash
ANTHROPIC_BASE_URL=...
ANTHROPIC_MODEL="sensenova-..."
# 或 glm / minimax / deepseek 等
```

Claude Code 可能把这个模型视为 **unknown model ID**，新版本会更保守地限制它。

先查当前环境变量：

```bash
env | grep -E 'CLAUDE_CODE.*CONTEXT|COMPACT|ANTHROPIC_MODEL'
```

试官方给的逃生开关：

```bash
export CLAUDE_CODE_DISABLE_UNKNOWN_MODEL_WINDOW_ENFORCEMENT=1
```

再重新启动 `claude`。

### 实际环境

当前用的 sensenova 配置：

```bash
function sensenova {
    export ANTHROPIC_BASE_URL="https://token.sensenova.cn"
    export ANTHROPIC_API_KEY=""
    export ANTHROPIC_MODEL="sensenova-6.8-flash-lite"
    export ANTHROPIC_DEFAULT_SONNET_MODEL="sensenova-6.8-flash-lite"
    export ANTHROPIC_DEFAULT_HAIKU_MODEL="sensenova-6.8-flash-lite"
    export ANTHROPIC_DEFAULT_OPUS_MODEL="sensenova-6.8-flash-lite"
    export CLAUDE_CODE_AUTO_COMPACT_WINDOW="200000"
    export CLAUDE_CODE_MAX_CONTEXT_TOKENS="262144"
}
```

设置了 `CLAUDE_CODE_AUTO_COMPACT_WINDOW=200000` 和 `CLAUDE_CODE_MAX_CONTEXT_TOKENS=262144`，但仍然频繁 compact。最可能的原因：**2.1.223+ 把 sensenova 的 model ID 判成 unknown，然后按内部较小窗口频繁 compact**，覆盖掉环境变量也没用。

## 经验总结

- `claude update` 坏了时不要怀疑网络，直接用 `npm install -g @anthropic-ai/claude-code@latest`。更新后检查 `which claude` 指向的路径，旧的 `.local/bin/claude` 可能还在
- peer messaging 不可用时，先查 `~/.claude/sessions/*.json` 里有没有 `messagingSocketPath`。没有的话重启对应 session，不是版本太旧
- `i.thinking.length` 报错是 thinking 块污染了 session 历史，`claude --resume` 或开新 session 即可
- 频繁 compact 大概率是第三方 model ID 被当成 unknown，设 `CLAUDE_CODE_DISABLE_UNKNOWN_MODEL_WINDOW_ENFORCEMENT=1`