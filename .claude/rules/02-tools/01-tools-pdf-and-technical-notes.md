# 可复用工具

## 可复用工具（tools/）

> **重要提示**：**不要使用 `tools/fix_markdown.py` 脚本**。该脚本的自动修复可能破坏格式、误改内容。格式修复应该通过Read和Edit工具手动完成，确保准确性。

- **PMC/PubMed全文下载（best-effort）**：`tools/fetch_pmc_assets.py`
  - 用途：给定PMCID或PMC链接，自动收集全文PDF与补充材料链接并尝试下载，输出一个可追溯的`index.json`。
  - 重要限制（2026-01现状）：PMC网页端的PDF/补充材料下载可能触发**Proof-of-Work反爬页面**（“Preparing to download ... / cloudpmc-viewer-pow”）。脚本会自动改走官方`oa.fcgi`返回的FTP直链（或下载`oa_package/*.tar.gz`并解包提取PDF）来绕开网页端拦截；若文章不在OA范围内，则仍可能需要浏览器手动下载。
  - 用法示例：
    - `python3 tools/fetch_pmc_assets.py PMC12519464 --out downloads/pmc`

# PDF搜索工具

## 推文

现在有了一个方便的工具来快速搜索PDF内容。我已经为你准备好了：

使用方式：`python3 tools/search_pdf_text.py <pdf_file> <keyword> [context_lines]`

  例子：
```python
# 搜索"MCS"，显示前后3行上下文
python3 tools/search_pdf_text.py "_pages/Free Energy/fep-the-end-of-parameter-tuning.pdf" "MCS" 3
```
必要的时候还是要直接读PDF全文

**格式清理工具**：`tools/remove_extra_blank_lines.py` —— 清理多余空行，`tools/fix_all_punctuation.py` —— 标点符号修复，`tools/convert_quotes.py` —— 引号修复（支持`CHECK_QUOTES=1`环境变量检查引号配对问题）

# 技术记录

### 技术记录

- 这种笔记整理，除了[GROMACS论坛](https://gromacs.bioexcel.eu/)，还需要把链接显式地写出来，否则微信识别不了，比如：[GROMACS论坛](https://gromacs.bioexcel.eu/)：https://gromacs.bioexcel.eu。
- 没有代码、都是link的，不要放代码框。。正常的文字不要放quote 
- pip cache purge这种放在代码框或``包裹
- 不要泄露隐私内容。
- 对于这些文档，title: 和第一个# 要加点【笔记整理|2025-07】之类的标签，文件名不用加。这个时间是笔记实际的大致时间，不是创建文档的今天。
- 创建的文档不要保留奇怪的时间戳，[02:56] 这种
- 不要把参考资源：这种文字类的放在```里面，没必要代码框的文字说明就不要代码框，如“目前rdkit.Chem.Draw.MolsToGridImage函数没有直接设置图例字体大小的选项”就是个经验，文字就行
- 检查最新的一些md文件，把frontmatter中的date: "2025-11-27"这种东西改成最后修改日期，文件名也需要相应地改。多检查一些，不限于git没提交的；附录的日期得跟着主文档吧
- 反正每次我给你提啥问题都是希望你按照我的要求来修改文档，都尽量看文献原文
