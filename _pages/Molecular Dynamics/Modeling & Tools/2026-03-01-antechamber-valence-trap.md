---
title: "antechamber 的一个隐蔽坑：羧基键级被改写后的 valence 报错"
date: "2026-03-01"
tags: [molecular-dynamics, ambertools, antechamber, gaff, mol2, forcefield]
description: "记录一次 antechamber 报 Weird atomic valence 的排查过程，根因是中间文件丢失双键"
image: "/assets/img/thumbnail/nightgardenflower.jpg"
thumbnail: "/assets/img/thumbnail/nightgardenflower.jpg"
author: Xufan Gao
lang: zh-CN
---

# antechamber 的一个隐蔽坑：羧基键级被改写后的 valence 报错

下面是一段完整、可复现的排查故事。场景很常见：羧酸盐配体在自动化流程中报错，但单独跑 antechamber 又能过。

## 症状与第一眼判断

报错信息通常长这样：

```
Fatal Error!
Weird atomic valence (3) for atom (ID: 1, Name: C1).
Possible open valence.
Warning: This molecule has no hydrogens nor halogens.
```

第一反应往往是“结构不合理”或“键级没写对”。但这个案例里，**原始 mol2 的键级完全正确**。

## 复现路径

直接在命令行运行下列命令可以通过：

```
antechamber -i ligand.mol2 -fi mol2 -o ligand.prep -fo prepi -at gaff -nc -2
```

而在自动化流程里，通常会采用两步式处理：

```
antechamber -i ligand.mol2 -fi mol2 -o ligand_gaff.mol2 -fo mol2 -c gas -s 2 -at gaff -nc -2
```

```
antechamber -i ligand_gaff.mol2 -fi mol2 -o ligand.prep -fo prepi -at gaff -nc -2
```

报错发生在第二步。

## 关键证据：中间文件改写了双键

对比原始 mol2 与中间 mol2 的键级后发现，羧基双键被改写成了单键。对于 sp2 碳而言，这会让连接数降为 3，acdoctor 以连接数而非键级和判定 valence，于是直接终止。

这一点解释了两个看似矛盾的现象：
- 原始 mol2 能通过
- 中间 mol2 会触发 “Weird atomic valence (3)”

## 另一个会干扰判断的细节

如果在排查过程中手动加了 H 或更改质子化态，**务必同步更新 mol2 的部分电荷**。否则 `-nc` 与总电荷不一致，会把排查方向彻底带偏。这个问题和 valence 报错是两条独立链路，需要分别确认。

## 为什么文档会建议 `-s 2`

antechamber 会调用一系列子程序并生成多个中间文件，文档说明这些中间文件通常是全大写命名。遇到问题时，推荐用 `-s 2` 输出详细日志，逐步定位是哪一步把键级改写了。

在本例中，acdoctor 在预检查阶段就失败，**还没进入重新判断键级的流程**。这也是为什么调整 `-j` 并没有效果。

## 稳定修复方式

最稳妥的修复是跳过 acdoctor 诊断：

```
antechamber -i ligand_gaff.mol2 -fi mol2 -o ligand.prep -fo prepi -at gaff -nc -2 -dr no
```

`-dr no` 只是不做诊断，不改变实际参数化逻辑。对结构正常的分子来说，acdoctor 原本就全部通过，跳过与否结果一致。

## 一句话结论

**不是结构错，而是中间 mol2 丢了双键，acdoctor 又在最前面把流程截断了**。先看中间文件，再考虑化学结构。

## 避坑清单

- 先单独运行 antechamber，确认原始 mol2 是否能过
- 核对 mol2 的部分电荷总和与 `-nc` 是否一致
- 用 `-s 2` 输出详细日志，检查中间文件是否保留键级
- 若中间 mol2 丢双键，可用 `-dr no` 跳过 acdoctor 诊断
