

https://playbooks.com/mcp/ndchikin/reference-mcp 
如何安装到全局？
● reference-mcp - get_citeas_data (MCP)(resource: " https://doi.org/10.1038/s41598-018-20204-8 ")
  ⎿  @journal-article{ITEM1, title={Measurements of AMPs in stratum corneum of atopic dermatitis and healthy skin–tape stripping technique},       
     journal={Scientific Reports},
     volume={8},
     … +5 lines (ctrl+o to expand)

● 🎉 成功了！ reference-mcp现在可以用了！工具名称是 mcp__reference-mcp__get_citeas_data（注意是连字符 - 而不是下划线 _）。
https://blog.csdn.net/gitblog_00094/article/details/154972394

  🔬 Context7的优势

  - ✅ 快速获取最新文档：比搜索引擎更精准
  - ✅ 代码示例丰富：直接可用的代码
  - ✅ API文档完整：参数、返回值、用法说明
  - ✅ 多库支持：不仅是Biopython，还有Uni-Mol等分子建模工具
https://github.com/anthropics/life-sciences
# Add the marketplace
/plugin marketplace add  https://github.com/anthropics/life-sciences.git 

# Install MCP servers
/plugin install pubmed@life-sciences
/plugin install biorender@life-sciences
/plugin install synapse@life-sciences
/plugin install wiley-scholar-gateway@life-sciences
/plugin install benchling-mcp@life-sciences
/plugin install 10x-genomics@life-sciences

# Install skills
/plugin install single-cell-rna-qc@life-sciences
  PubMed工具可用功能

  根据测试，PubMed MCP服务器提供的工具包括：
  - ✅ search_articles - 搜索文献（支持关键词、作者、期刊、日期过滤）
  - ✅ get_article_metadata - 获取详细文章信息（标题、摘要、作者、DOI等）
  - 📚 还有14个其他工具用于更高级的PubMed操作
 1. 已安装的组件

  - ✓ uv/uvx (v0.9.18) - Python 包管理器
    - 安装位置: /home/gxf1212/.local/bin/uvx
    - 已添加到 PATH (永久生效)
  - ✓ reference_mcp - CiteAssist MCP 服务器
    - 配置类型: stdio
    - 配置文件: ~/.claude.json (用户级别)

  2. 当前已配置的 MCP 服务器

  [
    "context7",       // 实时代码文档
    "filesystem",     // 文件系统访问
    "github",         // GitHub 集成
    "playwright",     // 浏览器自动化
    "reference_mcp"   // BibTeX 引用获取 ✨ 新增
  ]

  3. 服务器配置详情

  {
    "type": "stdio",
    "command": "/home/gxf1212/.local/bin/uvx",
    "args": ["reference-mcp"],
    "env": {}
  }

  可用工具

  🔍 get_citeas_data

  从 CiteAs 获取 BibTeX 格式的引用数据
  - 输入: DOI、URL 或关键词
  - 输出: BibTeX 格式的引用

  示例用法:
  请使用 get_citeas_data 获取这篇文章的引用信息：
  DOI: 10.1038/s41586-021-03819-2

  📚 get_scholar_data

  从 Google Scholar 获取 BibTeX 格式的引用数据
  - 输入: 搜索查询 + 结果数量(可选)
  - 输出: BibTeX 格式的引用列表
 示例用法:
  请使用 get_scholar_data 搜索关于 "molecular dynamics GROMACS" 的文献，
  返回前 5 条结果的 BibTeX 引用
能不能在用户级别安装并测试： https://github.com/blazickjp/arxiv-mcp-server 

● arxiv - search_papers (MCP)(query: "coarse-grained molecular dynamics polymer self-assembly", max_results: 5, categories: ["cond-mat.soft","physics.chem-ph"])
  ⎿  {                                                                                                                                                                                                                       
       "total_results": 5,
       "papers": [
     … +101 lines (ctrl+o to expand)

● 🎉 arXiv MCP 服务器测试成功！

  我成功搜索到了 5 篇与粗粒化分子动力学和聚合物自组装相关的论文：

  1. 1312.4416v1 (2013) - "Relation between the grafting density of liquid crystal macromolecule and the symmetry of self-assembled bulk phase: coarse-grained molecular dynamics study"
    - 作者: J. M. Ilnytskyi
    - 关于液晶大分子自组装的粗粒化 MD 研究
  2. 2308.12620v1 (2023) - "Modeling for heterogeneous oxidative aging of polymers using coarse-grained molecular dynamics"
    - 作者: Takato Ishida 等
    - 使用粗粒化 MD 模拟聚合物氧化老化过程
  3. cond-mat/0311358v1 (2003) - "Automatic Coarse Graining of Polymers"
    - 作者: Roland Faller
    - 聚合物的自动粗粒化方法

  可用的 MCP 工具：
  - ✅ search_papers - 搜索论文（已测试，工作正常）
  - ✅ download_paper - 下载论文
  - ✅ list_papers - 列出已下载论文
  - ✅ read_paper - 阅读论文内容
● ✅ Serena MCP 服务器安装成功！

  安装状态

  serena: ✓ Connected

  配置详情

  - 作用范围: User scope（所有项目可用）
  - 完整命令:
  uvx --from git+ https://github.com/oraios/serena  serena start-mcp-server \
    --context ide-assistant \
    --project-from-cwd
  关键修复

  添加了 --project-from-cwd 参数，使Serena能够：
  - 自动从当前工作目录检测项目
  - 搜索 .serena/project.yml 或 .git
  - 回退到当前工作目录（CWD）
  - 专为CLI代理（如Claude Code）设计

  Serena功能特性

  根据官方文档，Serena提供：

  1. 语义代码操作
    - find_symbol - 查找符号定义
    - find_referencing_symbols - 查找引用
    - insert_after_symbol - 智能插入代码
  2. IDE集成优化
    - ide-assistant context专为IDE内编码辅助优化
    - 支持VSCode、Cursor、Cline等
  3. 复杂代码库导航
    - 避免读取整个文件
    - 避免grep式搜索
    - 使用语义理解进行精确操作