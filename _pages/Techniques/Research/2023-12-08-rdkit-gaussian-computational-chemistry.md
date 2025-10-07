---
title: "【笔记整理|2023-09+2024年上半年】RDKit和Gaussian计算化学工具使用经验"
date: "2023-12-08"
tags: [rdkit, gaussian, computational-chemistry, quantum-chemistry, molecular-descriptors, antechamber, resp-charges, deepchem, molecular-alignment, mcs]
description: "总结RDKit化学信息学处理与Gaussian量子化学计算的实用技巧，涵盖分子导入、片段连接、电荷计算、最大共同子结构分析等核心功能，为计算化学工作者提供工具使用经验"
thumbnail: "/assets/img/thumbnail_mine/wh-dpe6lm.jpg"
image: "/assets/img/thumbnail_mine/wh-dpe6lm.jpg"
---

# 【笔记整理|2023-09+2024年上半年】RDKit和Gaussian计算化学工具使用经验

本文总结了在使用RDKit进行化学信息学处理和Gaussian进行量子化学计算时的实用技巧、常见问题和解决方案。

## RDKit分子处理

### 基础分子操作

#### 分子导入和基本处理
```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

# 读取分子
mol = Chem.MolFromMol2File('molecule.mol2')
mol = Chem.AddHs(mol)  # 添加氢原子
```

#### 分子片段连接
RDKit提供了强大的分子片段连接功能：
```python
from rdkit.Chem import rdmolops

def connect_mols(mol1, mol2, atom1, atom2):
    # 连接两个分子片段的函数
    # atom1和atom2是连接点的原子索引
    pass

# 参考资源：[RDKit片段连接指南](https://iwatobipen.wordpress.com/2020/10/16/easy-way-to-connect-fragments-rdkit-tips-memo/)：https://iwatobipen.wordpress.com/2020/10/16/easy-way-to-connect-fragments-rdkit-tips-memo/
```

#### 分子片段处理
```python
# 获取分子片段
from rdkit.Chem.rdmolops import GetMolFrags

# 处理虚原子标记片段
# 在RDKit中，虚原子可以用来标记这是一个片段
```

#### 分子组合
```python
from rdkit.Chem import CombineMols

# 组合多个分子
combined_mol = CombineMols(mol1, mol2)
```

### 分子可视化和绘制

#### 网格图像生成
```python
from rdkit.Chem import Draw

# 生成分子网格图像
```

**注意**：目前rdkit.Chem.Draw.MolsToGridImage函数没有直接设置图例字体大小的选项。

#### 高级绘制选项
```python
# 分子绘制选项设置
from rdkit.Chem.Draw import MolDrawing, rdMolDraw2D

# 分子绘制选项
# 参考: [RDKit绘制选项文档](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.MolDrawing.html#rdkit.Chem.Draw.MolDrawing.DrawingOptions): https://www.rdkit.org/docs/source/rdkit.Chem.Draw.MolDrawing.html#rdkit.Chem.Draw.MolDrawing.DrawingOptions

# 分子2D绘制选项
# 参考: [RDKit 2D绘制选项](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions): https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
```

#### 多分子高亮显示

[RDKit高亮显示博客](https://greglandrum.github.io/rdkit-blog/posts/2021-08-07-rgd-and-highlighting.html): https://greglandrum.github.io/rdkit-blog/posts/2021-08-07-rgd-and-highlighting.html

**注意**：DrawMolsToGridImage()不支持多重高亮显示功能。

### 文件格式和兼容性

#### mol2文件处理
处理mol2文件时的常见问题：

**价态错误处理**：
- 如果遇到："Explicit valence for atom # 8 N, 4, is greater than permitted"
- 这通常是因为氮原子的价态设置不正确

#### 分子坐标处理
```python
# 将分子质心移动到原点(0,0,0)
def translate_mol_to_origin(mol):
    # 计算质心并进行平移变换
    pass
```

### 高级分子处理

#### 分子体积计算
```python
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.AllChem import ComputeMolVolume

# 计算分子体积
volume = ComputeMolVolume(mol)
```
[RDKit分子体积计算文档](https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html#rdkit.Chem.AllChem.ComputeMolVolume): https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html#rdkit.Chem.AllChem.ComputeMolVolume

#### 分子对齐与匹配
```python
from rdkit.Chem import rdMolAlign

# 分子对齐：提供原子映射，使用反向GetSubstructureMatch
match = mol.GetSubstructMatches(cmn_core)
```
[RDKit分子对齐文档](https://www.rdkit.org/new_docs/source/rdkit.Chem.rdMolAlign.html): https://www.rdkit.org/new_docs/source/rdkit.Chem.rdMolAlign.html

#### 最大公共子结构（MCS）
```python
# MCS计算
from rdkit.Chem import rdFMCS

# 计算最大公共子结构
mcs = rdFMCS.FindMCS([mol1, mol2])
```
[RDKit MCS文档](https://rdkit.org/docs/source/rdkit.Chem.MCS.html): https://rdkit.org/docs/source/rdkit.Chem.MCS.html

**3D MCS应用**：
[RDKit博客3D MCS文章](https://greglandrum.github.io/rdkit-blog/posts/2022-06-23-3d-mcs.html): https://greglandrum.github.io/rdkit-blog/posts/2022-06-23-3d-mcs.html

## Gaussian计算

### 环境配置和权限问题

#### 权限问题解决
Gaussian对文件权限要求非常严格：
```bash
# 运行时如果提示"files in the gaussian directory are world accessible. this must be fixed"
find . -type f -exec chmod a+x {} \;
# 或者使用
chmod 750 -R *
```

**原因**：Gaussian如果发现其可执行文件对所有用户都可访问时就会拒绝运行，这是Gaussian的一个固执特点。

### 输入文件生成

#### 从mol2文件生成Gaussian输入
```bash
# 常见需求：从mol2文件生成包含连接信息的Gaussian输入文件
# 可以使用antechamber进行转换
antechamber -i input.mol2 -fi mol2 -o output.gjf -fo gcrt
```

#### 连接信息处理

**注意**：antechamber/G16猜测连接列表时，键序不一定正确，但需要保证合理性。

### 量子化学计算类型

#### RESP电荷计算
RESP (Restrained Electrostatic Potential) 电荷是分子动力学中常用的原子电荷：

```bash
# 使用antechamber计算RESP电荷
antechamber -fi gout -fo ac -i pet.log -o pet.ac -c resp -pf y

# 分离运行RESP计算
run resp separately....
```

#### AM1-BCC电荷方法

AM1-BCC stands for Austin Model 1 with Bond Charge Correction. 它是计算原子电荷的半经验方法。AM1方法是一种半经验量子化学方法，使用拟合到实验数据的参数集。BCC方法是对AM1电荷的修正，提高其准确性。

#### 电荷约束设置
在antechamber或Multiwfn中手动指定电荷约束：

**示例**：残基末端的电荷为0

**参考**：Multiwfn手册 4.7.7.4 Example 4: 天冬氨酸残基的原子电荷评估，包含等价和电荷约束的示例。

### 高级计算设置

#### 连接信息和拓扑

**问题**：Gaussian默认不提供连接信息，是否可能获得MD模拟的准确键、角度？

这是一个常见问题，通常需要：
1. 使用其他工具（如antechamber）推断连接
2. 手动指定键连接信息
3. 使用分子编辑器预处理

## 文件格式处理

### mol2格式详解

#### TRIPOS格式理解

TRIPOS mol2格式示例：
```
@<TRIPOS>MOLECULE
lig
45 47 0 0 0
SMALL
GASTEIGER
```

#### 常见格式问题
- Gview导出时坐标格式的一致性
- 不同软件之间mol2格式兼容性
- 原子类型和电荷信息的处理

### antechamber工具深度应用

#### 基本用法
```bash
# 从Gaussian输出文件生成mol2
antechamber -i bay.log -fi gout -o bay.mol2 -fo mol2

# 支持的文件格式
# .mc文件支持：antechamber accept .mc file?
```

#### Python集成
```python
# 在Python中调用antechamber
import subprocess

def run_antechamber(input_file, output_file, input_format, output_format):
    cmd = f"antechamber -i {input_file} -fi {input_format} -o {output_file} -fo {output_format}"
    subprocess.run(cmd, shell=True)
```

## 力场参数优化

### CGenFF参数优化器

#### 自动优化功能
CGenFF Parameter Optimizer提供自动优化可旋转二面角的功能：

1. **用户指定**待优化的二面角
2. **QM数据生成**：协调生成量子力学目标数据
3. **参数拟合**：使用LSFitPar最小二乘拟合程序
4. **多重度优化**：
   - 初始多重度由CGenFF程序分配
   - 自动尝试多重度1, (1,2), (1,2,3), (1,2,3,6)
   - 如果RMSE改善超过阈值（默认10%），选择更好的参数

#### QM计算集成
- 首先生成Psi4 QM任务
- 收集QM二面角扫描数据
- 拟合力场参数到这些目标数据

## 实用工具和脚本

### Multiwfn应用
```bash
# Multiwfn可执行文件权限设置
chmod +x /path/to/Multiwfn_3.8_dev_bin_Linux/Multiwfn
```

### ACPYPE工具
结合AmberTools + ACPYPE + Gaussian创建小分子GAFF力场的拓扑文件：
- 参考：[ACPYPE GAFF力场创建指南](https://jerkwin.github.io/2015/12/08/使用AmberTools+ACPYPE+Gaussian创建小分子GAFF力场的拓扑文件/)：https://jerkwin.github.io/2015/12/08/使用AmberTools+ACPYPE+Gaussian创建小分子GAFF力场的拓扑文件/

### 在线工具和资源

#### RESP电荷计算工具
- **R.E.D.** (RESP ESP charge Derive)：在线RESP电荷计算程序
- 虽然界面设计较旧，但功能齐全
- 更新状态：Last update of the R.E.D. Home Page: June 16th, 2017

#### 文档和教程
- [RESP电荷计算指南](https://jamesmccarty.github.io/research-wiki/RESP)：https://jamesmccarty.github.io/research-wiki/RESP
- [RDKit讨论区](https://sourceforge.net/p/rdkit/mailman/)：https://sourceforge.net/p/rdkit/mailman/
- [mol2格式说明](http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf)：http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf

## 常见错误和解决方案

### RDKit相关错误

#### 价态问题
```
reading mol2: Explicit valence for atom # 8 N, 4, is greater than permitted
```
**解决方案**：
- 检查mol2文件中氮原子的键连接
- 确认原子类型设置正确
- 必要时手动调整分子结构

#### 导入问题
- 确保mol2文件格式正确
- 检查原子坐标和连接表的一致性
- 注意不同软件生成的mol2文件格式差异

### Gaussian相关错误

#### 权限错误
最常见的Gaussian错误之一，严格按照权限设置要求执行：
```bash
chmod 750 -R gaussian_directory/
```

#### 连接猜测问题
- Gaussian的连接猜测算法有时不准确
- 建议使用其他工具预处理分子结构
- 或手动指定连接信息

## 工作流程建议

### 典型的小分子参数化流程
1. **结构优化**：Gaussian几何优化
2. **电荷计算**：RESP或AM1-BCC电荷
3. **参数生成**：antechamber生成力场参数
4. **验证检查**：RDKit验证分子结构合理性
5. **MD准备**：转换为MD程序所需格式

### 质量控制检查点
- 分子几何的合理性
- 电荷分布的物理意义
- 力场参数的完整性
- 与实验数据的一致性

## 深度学习与化学信息学

### DeepChem应用

#### 基础使用
```python
import deepchem as dc

# DeepChem是用于药物发现和化学信息学的深度学习库
# 提供分子特征化、模型训练和预测功能
```

DeepChem是专门为药物发现和化学信息学设计的深度学习库，集成了多种分子表示方法、模型架构和评估指标。

### 分子可视化扩展工具

#### Mols2grid网格显示
```python
import mols2grid

# 显示和滚动浏览聚类样本
mols2grid.display(molecules)
```

mols2grid提供了交互式的分子网格显示功能，特别适合大量分子的筛选和比较工作。

## 集成化学信息学工作流

### 现代化学信息学技术栈

1. **RDKit**: 核心分子处理和计算
2. **DeepChem**: 深度学习模型开发
3. **Gaussian**: 量子化学计算
4. **Mols2grid**: 交互式分子可视化
5. **antechamber**: 力场参数生成

### 推荐的集成工作流程

1. **分子预处理**: RDKit标准化和验证
2. **特征提取**: 结合传统描述符和深度学习特征
3. **量子计算**: Gaussian优化和性质计算
4. **模型开发**: DeepChem构建预测模型
5. **结果可视化**: mols2grid交互式展示

---

*本文基于2023年9-12月和2024年上半年技术讨论记录整理，涵盖计算化学和化学信息学工具使用中的实际问题和解决方案*