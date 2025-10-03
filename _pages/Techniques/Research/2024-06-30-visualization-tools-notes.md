---
title: "【笔记整理|2024年上半年】科学可视化工具实用技巧集锦"
date: "2024-06-30"
tags: [visualization, vmd, pymol, chimerax, scientific-visualization, molecular-graphics, trajectory-analysis]
---

# 【笔记整理|2024年上半年】科学可视化工具实用技巧集锦

## VMD使用技巧

### 基本设置与渲染

**渲染模式优化**：
VMD默认使用称作Normal的Rendermode，但此时有些材质的显示效果很差，甚至Transparent材质根本没法正确显示出透明效果。因此通过以下命令将默认的Rendermode设为效果好得多的GLSL：
```tcl
display rendermode GLSL
```

### VMD脚本与命令

**TCL脚本中执行bash命令**：
可以在TCL脚本中直接执行bash命令：
```tcl
exec grep 'ATOM' ${i}.pdb >> ${outputFile}
```

**动画控制**：
```tcl
animate goto 296
```

**播放MVD文件**：
```tcl
play view.mvd
```

### VMD路径与集成

**Windows上的VMD路径**：
```bash
/mnt/c/Program\ Files/VMD/vmd.exe
```

在WSL中使用Windows版VMD：
```bash
alias vmd='vmd.exe'
```

**VMD插件路径**：
```
/lib/vmd/plugins/LINUXAMD64/bin/catdcd5.2
```

### VMD坐标变换

#### transabout命令详解

**语法和参数**：
```tcl
# 绕指定轴和向量旋转的变换矩阵
transabout v amount [deg|rad|pi]
```

**参数说明**：
- `v`：旋转轴向量，格式为 `{x y z}`，如 `{0 0 1}` 表示绕Z轴旋转
- `amount`：旋转角度的数值
- `deg|rad|pi`：角度单位，分别表示度、弧度或π的倍数

**实际应用示例**：
```tcl
# 绕Z轴旋转90度
set rot_matrix [transabout {0 0 1} 90 deg]

# 绕任意向量{1 1 1}旋转π/4弧度
set rot_matrix [transabout {1 1 1} 0.25 pi]

# 应用变换到原子选择
set sel [atomselect top "protein"]
$sel move $rot_matrix
```

**变换原理**：生成绕通过原点沿给定向量的轴逆时针旋转指定角度的4x4齐次变换矩阵，可以与其他变换（平移、缩放）组合使用。

[VMD变换命令文档](https://www.ks.uiuc.edu/Research/vmd/current/ug/node194.html): https://www.ks.uiuc.edu/vmd/current/ug/node194.html

#### 嵌套列表处理问题详解

**问题背景**：VMD中获取原子坐标时经常遇到嵌套列表格式问题，这是VMD Tcl脚本编程中的常见陷阱。

**问题表现**：
```tcl
# 错误的坐标格式（嵌套列表）
set coords [$atm get {x y z}]
# 结果: {{10.5 20.3 30.7}} - 注意双重大括号！

# 期望的格式（简单列表）
# 结果: {10.5 20.3 30.7} - 单层大括号
```

**为什么会出现嵌套列表**：
- VMD的`get`命令返回的是列表的列表
- 每个原子的坐标作为一个子列表存储
- 即使只有一个原子，也会返回包含一个元素的列表

**解决方案**：
```tcl
# 方法1：使用lindex提取第一个元素
set coord1 [lindex [$atm get {x y z}] 0]

# 方法2：处理多个原子的坐标
set sel [atomselect top "protein"]
set coords [$sel get {x y z}]
foreach coord $coords {
    set x [lindex $coord 0]
    set y [lindex $coord 1]
    set z [lindex $coord 2]
    # 处理单个原子坐标
}

# 方法3：计算两点间距离的完整示例
set sel1 [atomselect top "resid 1 and name CA"]
set sel2 [atomselect top "resid 10 and name CA"]
set coord1 [lindex [$sel1 get {x y z}] 0]
set coord2 [lindex [$sel2 get {x y z}] 0]
set distance [vecdist $coord1 $coord2]
```

[VMD用户邮件列表参考](https://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/2584.html): https://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/2584.html

#### 高级坐标变换技巧

**组合变换**：
```tcl
# 先平移再旋转
set trans_matrix [transoffset {5 0 0}]  # 沿X轴平移5埃
set rot_matrix [transabout {0 0 1} 45 deg]  # 绕Z轴旋转45度
set combined_matrix [transmult $rot_matrix $trans_matrix]
$sel move $combined_matrix
```

**分子对齐**：
```tcl
# 将分子质心移到原点，然后旋转
set sel [atomselect top "backbone"]
set center [measure center $sel]
set trans_to_origin [transoffset [vecscale -1 $center]]
$sel move $trans_to_origin
$sel move $rot_matrix
```

## PyMOL操作指南

### 基本操作

**菜单操作**：
启动VMD后按"Push Menus"

**蛋白质轨迹对齐**：
在PyMOL中，使用intra_fit命令将蛋白质轨迹对齐到第一帧：
```pymol
intra_fit
```

### PyMOL设置优化

**正交投影设置**：
```pymol
set orthoscopic, on
```

[PyMOL正交投影文档](https://pymolwiki.org/index.php/Orthoscopic): https://pymolwiki.org/index.php/Orthoscopic

### PyMOL轨迹制作

**电影制作教程**：
[PyMOL电影制作指南](https://pymol.org/tutorials/moviemaking/): https://pymol.org/tutorials/moviemaking/

PyMOL提供了完整的轨迹电影制作功能，适合制作高质量的分子动画。

## ChimeraX高级功能

### 视图设置

**正交视图**：
```
camera ortho
```

**相机设置文档**：
[ChimeraX相机命令](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/camera.html): https://www.cgl.ucsf.edu/chimerax/docs/user/commands/camera.html

### 晶胞显示

**显示晶胞轮廓**：
```
unitcell outline
```

这对于显示周期性边界条件下的分子动力学模拟结果特别有用。

### 尺寸控制

**对象尺寸调整**：
[ChimeraX尺寸命令文档](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/size.html): https://www.cgl.ucsf.edu/chimerax/docs/user/commands/size.html

### PBC盒子显示

在Chimera中显示蛋白质-配体系统周围的PBC盒子/单元晶胞，这对于MD模拟结果的可视化很重要。可以用于录制MD模拟后的影片。

### 螺旋圆柱显示

ChimeraX提供了螺旋圆柱显示功能，可以更好地展示蛋白质的二级结构。

[ChimeraX螺旋圆柱命令文档](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/spiral.html): https://www.cgl.ucsf.edu/chimerax/docs/user/commands/spiral.html

### 系统兼容性检查

**WSL中的显示问题**：
WSL中的VMD，display功能无法正常显示任何内容，建议使用原生Linux版本或Windows版本。

如果PyMOL和ChimeraX都有问题，那就是系统级别的问题。需要检查：
- 显卡驱动是否正常
- OpenGL支持是否完整
- 系统库文件是否缺失

## 分子结构文件处理

### 坐标文件转换

**从坐标和拓扑文件生成PDB**：
```bash
ambpdb -p topology-file < coordinates-file > filename.pdb
```

**Amber文件转换示例**：
```bash
ambpdb -p cram.prmtop -c min_qmmm.rst > min_qmmm.rst.pdb
```

### 分子重心平移

将mol2文件的质心平移到(0,0,0)是常见的分子预处理操作，可以通过坐标计算和平移实现。

## 轨迹分析与可视化

### 文件上传与目录结构保持

上传所有mapping.png文件并保持父目录结构时，简单的scp不太理想。更强大简洁的方法是结合tar和ssh：

```bash
tar -czf - mapping.png | ssh user@remote 'cd /target/dir && tar -xzf -'
```

### 主成分轴长度计算

在MDAnalysis中计算蛋白质三个主成分轴的长度：
```python
import MDAnalysis as mda

# 计算惯性张量和主成分轴
# 然后计算每个轴的长度
```

这对于分析蛋白质形状变化很有用。

## 数据可视化选择

### 图表库选择

在现代web开发中，推荐使用：
- **Tailwind CSS**：用于布局和样式设计
- **Chart.js**：用于标准图表
- **Plotly.js**：用于复杂图表，确保使用Canvas/WebGL渲染

所有图表和图示都应该避免使用SVG和Mermaid JS，转而使用HTML/CSS、Unicode字符或Canvas来实现。

### 分子网格显示

**mols2grid使用**：
```python
import mols2grid

# 显示和滚动浏览聚类样本
mols2grid.display(molecules)
```

这对于大量分子的筛选和比较非常有用。

## 小结

科学可视化工具的选择和配置对研究效率有重要影响。VMD适合复杂的轨迹分析和脚本化操作，PyMOL在分子图形制作方面表现出色，ChimeraX则提供了现代化的用户界面和强大的渲染能力。正确配置这些工具，结合合适的数据处理流程，能够显著提升科学研究中的可视化质量和效率。同时，了解跨平台兼容性问题和性能优化技巧，有助于构建稳定高效的可视化工作环境。