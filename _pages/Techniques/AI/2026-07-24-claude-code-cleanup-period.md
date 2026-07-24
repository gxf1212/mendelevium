---
title: "Claude Code 默认30天清理会话？已经删完了，只剩 prompt了。。"
date: "2026-07-24"
last_modified_at: "2026-07-24"
tags: [Claude-Code, settings, session-cleanup, AI-tools, data-recovery]
description: "Claude Code默认30天自动清理会话，超过期限的项目对话只剩prompt。设置cleanupPeriodDays延长保留期，再配一个本地备份脚本"
image: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/Wallpaper_compressed/wallhaven-9mogqx.jpg"
thumbnail: "https://raw.githubusercontent.com/gxf1212/mendelevium/main/assets/img/Wallpaper_compressed/wallhaven-9mogqx.jpg"
author: Xufan Gao
lang: zh-CN
---

# Claude Code 默认30天清理会话？已经删完了，只剩 prompt了。。

今天回看几个之前跑的科研项目，发现 `~/.claude/projects/` 下对应的 session 都没了，甚至文件夹都没有了。。**对话jsonl文件全部/部分清空，只剩 prompt 还在总体jsonl文件**。

查了一下，这是 Claude Code 的官方行为。

## 30天默认清理的逻辑

Claude Code 每次启动都会清理一次 `~/.claude/projects/` 下早于设定天数的会话文件，孤立的 worktree 也按同一个阈值清理。

[官方文档](https://code.claude.com/docs/en/settings)（https://code.claude.com/docs/en/settings）：

> **cleanupPeriodDays** — 默认30天，最小值1。Claude Code 启动时删除早于该期限的会话文件和其他应用数据。设置成 0 会触发验证错误。

v2.1.203 之前有个坑：`settings.json` 解析失败时，清理会在30天默认值下运行，**配置里设了更长的保留期也会被绕过**。30天以内的文件不受影响。

## 如何把清理周期延长

打开 `~/.claude/settings.json`，加一行 `cleanupPeriodDays`：

```json
{
  "cleanupPeriodDays": 9999999,
  "permissions": { "allow": [] }
}
```

9999999 天大概是两万七千年。设置成 `0` 会直接报验证错误，所以用一个大数代替。

## 备份也要做

配置只能避免以后再丢，已经删掉的 session 暂时没有官方恢复方法。建议再加一个本地备份脚本：

```bash
#!/bin/bash
# backup_claude_projects.sh
set -euo pipefail

SOURCE="$HOME/.claude/projects"
BACKUP="$HOME/data2/backup/claude_projects"
LOG="$BACKUP/backup.log"
LOCK="/tmp/backup_claude.lock"

mkdir -p "$BACKUP"

# 防并发
if [ -f "$LOCK" ]; then
  echo "[$(date)] SKIP: another backup is running"
  exit 0
fi
echo $$ > "$LOCK"
trap 'rm -f "$LOCK"' EXIT

# rsync 增量同步，不删源端文件
rsync -av --backup --backup-dir=".bak/$(date +%Y%m%d)" \
  "$SOURCE/" "$BACKUP/" >> "$LOG" 2>&1

# 清理180天前的备份
find "$BACKUP/.bak" -maxdepth 1 -type d -mtime +180 -exec rm -rf {} + 2>/dev/null

echo "[$(date)] Backup done" >> "$LOG"
```

每天 cron 跑一次，增量同步到 `data2/backup/claude_projects/`，旧备份超过 180 天自动清理。`rsync --backup` 把覆盖的内容移到 `.bak/日期`，既节省空间，回滚也方便。

## 总结一下

- **30天清理是 Claude Code 的官方行为**，v2.1.203 之后逻辑更稳定
- **修改 settings.json 加 `cleanupPeriodDays`** 可以延长保留期；`0` 不行，用大数代替
- **本地备份要做**：已经删掉的 session 没有官方恢复途径，增量备份加 rsync `.bak` 是简单可行的方案

如果你的 session 也莫名消失，大概率就是这个机制在跑。建议先改 settings，再把备份脚本挂上 cron。
