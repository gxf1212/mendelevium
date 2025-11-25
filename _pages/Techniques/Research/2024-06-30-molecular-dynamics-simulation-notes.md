---
title: "【笔记整理|2024年上半年】分子动力学模拟实用技巧与经验总结"
date: "2025-10-08"
tags: [molecular-dynamics, gromacs, namd, plumed, umbrella-sampling, pmf, martini, simulation]
description: "分子动力学模拟的实用技巧与经验，涵盖GROMACS和NAMD轨迹处理、PLUMED增强采样、伞形采样PMF计算、Martini粗粒化力场应用等关键技术"
thumbnail: "/assets/img/thumbnail_mine/wh-6o2m9q.jpg"
image: "/assets/img/thumbnail_mine/wh-6o2m9q.jpg"
---

# 【笔记整理|2024年上半年】分子动力学模拟实用技巧与经验总结

## MD模拟技巧

### 轨迹分析与处理

#### Amber轨迹重启时间设置问题
```bash
ncdump -v time [path to your rst7 file]
```
当重启模拟时，初始时间从重启文件中读取。可以用上述命令检查重启文件中的时间设置。

#### Amber轨迹文件合并
使用cpptraj工具合并多个.nc轨迹文件：
```bash
cpptraj -p topology.prmtop
trajin file1.nc
trajin file2.nc
trajout combined.nc
```
cpptraj是AmberTools套件中处理轨迹文件的多功能程序，可以处理包括合并在内的各种操作。

### 温度耦合组优化设置

在GROMACS中，温度耦合组（tc-grps）的设置需要根据体系各组分的动力学特性进行合理分组，以平衡温度控制的精度和计算效率。针对脂双层膜-水-溶质体系，建议：
- 脂质分子单独成组
- 水分子单独成组
- 蛋白质/小分子溶质单独成组

#### 动态负载平衡设置
```bash
-dlb auto  # 默认开启
-dlb yes   # 显式指定
```
在粒子分布不均或相互作用强度不同的情况下动态调整域大小。

**注意**：在GPU常驻模式（使用-update gpu）时，动态负载平衡会被关闭。

## 伞形采样与PMF计算

### 拉动参数优化

#### 拉动力常数建议
拉动力常数建议设置在1000-5000之间比较合适，需要根据具体体系进行调试。

#### 收敛性检查
```bash
gmx wham -b 50000  # 只包含最后50ns
gmx wham -b 75000  # 只包含最后25ns
```
检查收敛性时，可以只包含每个模拟的最后50ns或25ns数据，通过-b选项控制。

#### PMF解读注意事项
PMF表面上最多计数的区域不一定对应能量最小值。这是因为PMF模拟施加了偏置势来采样特定区域，在能量计算时会去除这个偏置。如果用"无偏"模拟估算自由能，最小值才对应最大采样区域。

### 伞形采样窗口设置

#### 结合位点附近窗口密度
对于蛋白质-配体结合体系，可能需要在结合位点附近设置更多的窗口，而不是单纯延长每个窗口的模拟时间。

#### 长距离拉动设置
Direction-periodic选项应该只用于需要拉动超过半个盒子长度距离的情况。这种情况很少见，拉动大型聚合物可能是一个有效的使用场景。建议拉动距离略小于完整盒子尺寸，以避免周期性映像间的相互作用。

## Martini粗粒化力场

### Martini 3.0 参数和设置

#### Colvars使用
[Colvars](https://colvars.github.io/): https://colvars.github.io - 集合变量库，可用于增强采样和自由能计算。

#### Martini 3.0甾醇参数
[Martini 3.0甾醇参数](https://github.com/Martini-Force-Field-Initiative/M3-Sterol-Parameters/blob/main/martini_v3.0_sterols_v1.0.itp): https://github.com/Martini-Force-Field-Initiative/M3-Sterol-Parameters/blob/main/martini_v3.0_sterols_v1.0.itp

#### Martini 3.0脂质参数
[Martini脂质参数库](https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters): https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters

#### 镁离子表示
镁离子用一个TQ3p珠子表示，带电荷+1。

### 几何结合规则设置

#### vdWGeometricSigma参数
```
vdwGeometricSigma yes
```
在Martini力场中使用几何结合规则计算范德华相互作用参数。

## NAMD高级应用

### 多拷贝/副本交换设置

#### 多拷贝副本交换脚本接口
NAMD提供专门的脚本接口用于多拷贝/副本交换模拟设置。

#### 命令行参数传递
```bash
namd3 --outputenergies 100 --run 100
```
可以通过--keyword value参数对直接在命令行指定配置参数。

### 配置文件路径管理

工作目录自动切换
执行时NAMD会自动切换到包含配置文件的目录，使配置文件中的所有文件路径都相对于配置文件目录。可以指定多个配置文件，但所有文件路径都相对于第一个调用"run"命令的配置文件，或如果没有调用"run"则相对于最后一个配置文件。

## 轨迹可视化技巧

### ChimeraX使用技巧

正交投影设置
```
camera ortho
```
在ChimeraX中设置正交投影视图，便于科学可视化。

晶胞显示
```
unitcell outline
```
显示周期性边界条件的晶胞轮廓。

调整显示尺寸
参考[ChimeraX尺寸命令文档](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/size.html): https://www.cgl.ucsf.edu/chimerax/docs/user/commands/size.html

### PyMOL轨迹制作

PyMOL轨迹电影制作
参考[PyMOL电影制作教程](https://pymol.org/tutorials/moviemaking/): https://pymol.org/tutorials/moviemaking/

PyMOL正交投影设置
[PyMOL正交投影文档](https://pymolwiki.org/index.php/Orthoscopic): https://pymolwiki.org/index.php/Orthoscopic

## GROMACS选择语法

### 距离计算和选择

距离计算命令
```bash
gmx distance -s md_smd.tpr -f md_smd.xtc -n index.ndx -oav dist.xvg
```
计算指定原子组间的距离变化。

gmx select工具
```bash
gmx select  # 基本动态选择数据输出
gmx help selections  # 详细选择语法帮助
```
gmx select可以输出动态选择的基本数据，用于简单分析或与其他程序组合进行更复杂的计算。

## 编译与安装问题

### 库文件依赖解决

glibc库链接问题
```bash
ln -s /usr/lib64/libz.so.1 /path/to/glibc/lib
ln -s /usr/lib64/libstdc++.so.6 /path/to/glibc/lib
ln -s /usr/lib64/libgcc_s.so.1 /path/to/glibc/lib
```
编译安装新版本glibc时，需要手动链接系统中的其他必要库文件。

### CUDA兼容性

CUDA 12.2支持
CUDA版本12.2已被检测到，需要相应修改cmake/CudaConfig.cmake配置文件以确保兼容性。

## 相关资源

### GROMACS社区

[GROMACS论坛](https://gromacs.bioexcel.eu/): https://gromacs.bioexcel.eu - GROMACS官方技术支持论坛

[GROMACS PMF讨论](https://gromacs.bioexcel.eu/t/how-can-i-get-smooth-pmf-from-umbrella-sampling/3629): https://gromacs.bioexcel.eu/t/how-can-i-get-smooth-pmf-from-umbrella-sampling/3629

[伞形采样直方图问题](https://gromacs.bioexcel.eu/t/problem-with-umbrella-histograms/9216): https://gromacs.bioexcel.eu/t/problem-with-umbrella-histograms/9216

### 技术博客

[GROMACS分子间相互作用计算](https://jerkwin.github.io/2019/09/06/%E4%BD%BF%E7%94%A8GROMACS%E8%AE%A1%E7%AE%97%E5%88%86%E5%AD%90%E9%97%B4%E7%9B%B8%E4%BA%92%E4%BD%9C%E7%94%A8/): https://jerkwin.github.io/2019/09/06/%E4%BD%BF%E7%94%A8GROMACS%E8%AE%A1%E7%AE%97%E5%88%86%E5%AD%90%E9%97%B4%E7%9B%B8%E4%BA%92%E4%BD%9C%E7%94%A8/

## 小结

分子动力学模拟涉及众多技术细节，从参数设置到结果分析都需要丰富的经验积累。合理的温度耦合、动态负载平衡、以及针对性的采样策略是获得可靠结果的关键。同时，可视化工具的熟练使用能够帮助更好地理解模拟结果和发现问题。