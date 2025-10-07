---
title: "【笔记整理|2024-07】Python开发环境构建与性能优化：从编码规范到科学计算"
date: "2024-07-01"
tags: [python, performance-optimization, development-environment, scientific-computing, data-analysis, technical-notes]
---

# 【笔记整理|2024-07】Python开发环境构建与性能优化：从编码规范到科学计算

## 引言

Python作为科学计算和数据科学的主要编程语言，其开发环境的配置和性能优化直接影响研究效率。本文整理了从技术讨论中提取的Python开发环境构建、性能优化和科学计算的实用技巧，涵盖从编码规范到高级性能优化的各个方面。

## Python编码最佳实践

### 属性访问与动态操作

Python提供了灵活的属性访问和动态操作机制：

> 在Python中，如果你想要根据传递的变量动态地设置对象的属性值，可以使用setattr函数。

> In Python, the __getattr__ method is called when you try to access an attribute that does not exist, but it's not a standard way to access attributes. Instead, you typically access attributes using the dot notation (e.g., object.attribute).

### 迭代器优化

在Python中，前置和后置增量操作的性能差异值得注意：

> To be accurate: ++i can sometimes be faster than i++ and is never slower. For fundamental data types, the compiler will very likely fix your mistake and optimise away any unneeded copying. For iterators this is more difficult and for user-defined types it may very well be impossible.

### 排序算法与数据结构

Python内置的排序算法和数据结构特性：

> 在Python 中，内置的 sorted() 函数使用的是双轴快排算法（timsort）来对序列进行排序。 这种算法的时间复杂度在最坏情况下是O(n * log n)，平均情况下是O(n * log n)

> n * log n + a * log n ≈ n * log n < n * a

> log2 or ln?

### 哈希表与集合操作

理解Python中集合和字典的内部实现有助于性能优化：

> 是的，集合和字典在Python中都是通过哈希表实现的。对于集合和字典的元素或键的查找，时间复杂度通常是O(1)，这是因为哈希表使得元素的位置可以快速定位。

## Python性能分析

### 代码性能分析工具

Python提供了多种性能分析工具来识别性能瓶颈：

> 在Python中，你可以使用cProfile模块来分析每个函数的执行时间123。以下是一个示例：

**性能分析输出示例：**
> update_results    2    96579    66

### 实际性能对比

实际测试显示不同运行环境下的性能差异：

> pycharm profile says 71s,simply debug 56s, cmd just 31s

> other causes, fit: ~8s; concat: 6s

> process_dict 11.6s, including the two?

## 数据处理与优化策略

### Pandas数据处理

Pandas是Python数据分析的核心库，掌握其高级功能非常重要：

> df = df_input.copy(deep=True)  # Use pandas' built-in copy method

### 大规模数据优化

处理大规模数据时，性能优化尤为重要：

> # Optimal pipeline for huge data: fast_histogram + memory mapping

> fast_histogram doesn't require parallel processing as it's already optimized internally

### 字符串处理

字符串处理在数据分析中经常是性能瓶颈：

> transform the code into a clean, efficient, and maintainable analysis framework.

## 科学计算环境配置

### 包管理工具

合理的包管理策略可以避免依赖冲突：

> pip install -e .[dev]

### Conda环境管理

Conda是科学计算环境管理的首选工具：

> conda install conda-forge::libmamba

> 要查看pip的缓存路径，可以使用pip cache dir命令。在命令行或终端中输入该命令，pip会显示其缓存的目录。

### 环境共享

在多用户环境中共享conda环境可以提高效率：

> 看起来你想将用户 xucx 的 boltz2 Conda 环境共享给其他用户，让大家都能方便地通过 conda activate boltz2 来使用。最直接且对原用户影响较小的方式是创建符号链接。

## Python科学计算生态

### 科学计算库

Python拥有丰富的科学计算库生态系统：

> import deepchem as dc

### 数据可视化

数据可视化是科学计算的重要组成部分：

> In Matplotlib, the axes can be easily hidden by calling the set_visible() method on the axes object and setting it to False. This can be done either by using the axes object itself or by looping through the list of axes in a figure.

### 色彩映射与数据表达

合适的色彩映射可以增强数据的可读性：

> In the context of seaborn.diverging_palette(), h_neg and h_pos refer to the anchor hues that define the endpoints of the color spectrum for the diverging palette. These hues are specified in the HUSL (Hue, Saturation, Lightness) color space, where hue is an angle on the color wheel ranging from 0 to 360 degrees.

### 高级可视化技术

高级可视化技术可以更好地展示复杂数据：

> https://medium.com/@alexbelengeanu/getting-started-with-raincloud-plots-in-python-2ea5c2d01c11

## 开发工具与环境配置

### 代码编辑器配置

合适的代码编辑器配置可以提高开发效率：

> 1. 打开VSCode，并在左侧的文件资源管理器中选择你要检索字符串的项目文件夹。 2. 使用快捷键Ctrl+Shift+F，或者点击顶部菜单栏中的"查找" -> "查找"来打开查找面板。 3. 在查找面板的文本输入框中输入你要搜索的字符串。 你可以使用普通的文本字符串进行搜索，也可以使用正则表达式进行更高级的搜索。

> PyCharm 本身是一个代码编辑器（IDE），而不是一个网页浏览器。所以它不能像 Chrome 或 Edge 那樣直接"打开"并渲染 localhost:8501 的页面内容。

### 前端开发与后端集成

Python在现代Web开发中也有广泛应用：

> 我将使用Tailwind CSS进行布局和样式设计，并采用Chart.js（用于标准图表）和Plotly.js（如果需要更复杂的图表，并确保使用Canvas/WebGL渲染）来创建可视化内容。所有图表和图示都将严格遵守无SVG和无Mermaid JS的要求，转而使用HTML/CSS、Unicode字符或Canvas来实现。

> I designed a frontend to manage the analysis and figures. here's the overview. understand it

## Python包管理与发布

### 包缓存管理

合理管理包缓存可以节省磁盘空间并提高安装速度：

> 要清理pip的缓存，可以使用pip cache purge命令。这将清除pip缓存的所有内容，包括已下载但未安装的包和已安装但未使用的包的缓存。如果只想清除特定包的缓存，可以使用pip cache remove <package_name>命令，将package_name替换为要清除缓存的包名。

### Git与代码版本控制

版本控制是现代软件开发的标准实践：

> git config advice.addIgnoredFile false

> git config --global user.name "gxf1212"

## 文档生成与部署

### Sphinx文档系统

Sphinx是Python项目文档生成的标准工具：

> How do I serve `sphinx` documentation locally?

> 用claude code写文案可能会有点过于浪费了

### 静态网站生成

现代文档部署通常使用静态网站生成器：

> 📚 Complete Workflow: Public Documentation with Private

## 高级编程技巧

### 正则表达式应用

正则表达式是文本处理的强大工具：

> 要查找目录名中恰好包含两个连字符的目录，需要将grep模式"锚定"以匹配整行。

### 代码重构与优化

代码重构是提高代码质量的重要手段：

> transform the code into a clean, efficient, and maintainable analysis framework.

### 函数设计与最佳实践

良好的函数设计是高质量代码的基础：

> The most straightforward and conventional method is to prefix each line of the desired comment block with the hash symbol (#).

## Python与AI集成

### AI辅助开发

AI工具正在改变Python开发的方式：

> Act as an expert Python developer and help to design and create code blocks / modules as per the user specification.

> I asked ChatGPT about this, it says:

### Claude Code集成

Claude Code为Python开发提供了AI辅助：

> https://www.yuque.com/beihu-iq2oo/zlyf06/vlg45fk72pu9gmtk?singleDoc#%20%E3%80%8AClaude%20Code%EF%BC%9A%E8%AE%A1%E8%B4%B9%E4%B8%8E%E8%AE%A2%E9%98%85%E3%80%8B

> Claude Code：计费与订阅

> AICodemirror，必须curl -fsSL https://download.aicodemirror.com/env_deploy/env-deploy.sh | bash -s -- "sk-ant-api03-JQBd6V2vGYfPrl20II1Y3mGvRoK52kP7BJKUPSh4jCSoou4Jxw7ctQ3lVFJQ36tTO10cypFIIU8MYgbQ_78E3g"之后才能用

> What the Script Does: After setting the environment variables, the script finds your API key, takes the last 20 characters of it, and uses the jq command to add this snippet to a list inside the ~/.claude.json file. Specifically, it adds it to the customApiKeyResponses.approved array.

> must do this after sudo npm install -g @anthropic-ai/claude-code

### 环境配置脚本

自动化环境配置脚本可以简化开发环境搭建：

> (cat ~/.claude.json 2>/dev/null || echo 'null') | jq --arg key "${ANTHROPIC_API_KEY: -20}" '(. // {}) | .customApiKeyResponses.approved |= (.[], $key) | unique)' > ~/.claude.json.tmp && mv ~/.claude.json.tmp ~/.claude.json

## 实用编程技巧

### 文件操作技巧

高效的文件操作是数据处理的基础：

> Working with Zip Files

> zip s.zip software-copyright/ -r

### 系统命令集成

Python与系统命令的集成可以扩展功能：

> 03:14:40  |base|gxf1212@gxf-pop-os file-transfer → gnome-shell --version

> to fix https://extensions.gnome.org/extension/1160/dash-to-panel/

### 条件判断与逻辑

良好的条件判断逻辑可以提高代码的健壮性：

> for what it's worth

## 总结与最佳实践

1. **编码规范**：遵循Python编码规范，使用合适的属性访问方式和动态操作
2. **性能优化**：熟练使用性能分析工具，理解Python内部数据结构的实现原理
3. **环境管理**：合理使用conda和pip管理Python环境，解决依赖冲突
4. **科学计算**：掌握Python科学计算生态，包括数据处理、可视化和分析工具
5. **开发工具**：配置合适的开发环境，使用现代化的编辑器和工具链
6. **版本控制**：建立良好的Git使用习惯，确保代码的可追溯性
7. **文档生成**：使用Sphinx等工具生成高质量的项目文档
8. **AI集成**：合理利用AI工具提高开发效率，但不过度依赖

通过这些Python开发技巧的掌握，可以显著提高科学计算和数据处理的效率和质量。

## 参考资源

- [雨云图Python教程](https://medium.com/@alexbelengeanu/getting-started-with-raincloud-plots-in-python-2ea5c2d01c11)
- [Claude Code使用指南](https://www.yuque.com/beihu-iq2oo/zlyf06/vlg45fk72pu9gmtk?singleDoc#%20%E3%80%8AClaude%20Code%EF%BC%9A%E8%AE%A1%E8%B4%B9%E4%B8%8E%E8%AE%A2%E9%98%85%E3%80%8B)
- [文件压缩操作指南](https://docs.hostdime.com/hd/command-line/how-to-tar-untar-and-zip-files)
- [GNOME扩展修复](https://extensions.gnome.org/extension/1160/dash-to-panel/)
- [VS Code搜索功能文档](https://docs.github.com/zh/discussions/quickstart)