# 翻译类文章规范

翻译他人英文博客 / 文章（如 Pat Walters、其他课题组博文）与原创「文献精读」是两类不同的工作流。本文定义翻译稿的处理规范；精读稿仍按 `07-article-structure/*` 与 `04-article-formatting/*` 执行。

---

## 1. 文首注明出处

正文（frontmatter 之后）第一行写一句翻译说明，采用引用块：

> 本文翻译自 [作者名] 的博客 [原文标题](原文URL)，原文发表于 YYYY-MM-DD。以下为直译，版权归原作者所有。

只写这一句，不要展开、不要加「译者按」长文。URL 用 `[]()` 形式，禁止尖括号。

## 2. 不要加精读稿才有的章节

翻译稿**只忠实呈现原文结构**，不要自行添加以下章节（这些属于原创精读稿）：

- `## 摘要` / `### 核心结论` / `### 创新点`
- `## 背景` / `### 关键科学问题`
- `## 关键结论与批判性总结` / `## 参考链接`（原文若在文末给了裸链接，按原文位置保留即可，不必提成独立「参考链接」标题）
- 任何原文没有的表格、小结、批判性评论

如果原文本身有小结/表格，照翻；原文没有的，绝不补。

## 3. 直译，保留第一人称

- 原文是第一人称（I / we），译文中用「我 / 我们」，**不要改成第三人称「作者说 / 作者指出 / 作者认为」式叙述**。
- 不增译、不意译发挥、不补背景知识。遇到原文口语/俏皮表达，尽量贴近原味译出，不做学术化润色。
- 章节标题可译为贴合原文内容的中文，但应对应原文存在的标题，不要新造。

## 4. 专业名词

- 专有技术名词保留英文，首次出现用「中文（English）」或「English（中文）」括注一次，其后可只用英文。
- 模型名、数据集名、方法名、工具名一律保留英文原写法（如 MEGA-CL、Monroe、Mol-JEPA、CheMeleon、ChemProp、TabPFN、TabICL、ADMET、ECFP、PM6、PCBA、ChEMBL、TDC、InChIKey、SMILES、GRIT、JEPA、herdr、Tukey HSD、MAE）。
- 与博客标签纯中文的约定不冲突：正文名词保留英文，tags 字段仍用纯中文小写连字符。

## 5. frontmatter 与图片

- 九字段照常补全：`author` 仍填博客主 `Xufan Gao`（出处靠文首翻译说明体现，不在 author 字段改原作者）；`image`/`thumbnail` 用本仓库图床（GitHub raw URL），不用原博客图。
- 内文配图：**优先把图下载到文章同级文件夹**（如 `emlmf/`、`agents_benchmarking_figs/`），用相对路径引用（如 `agents_benchmarking_figs/fig1.png`），与博客「图片放 md 同级子文件夹、不散放、不依赖远程」的惯例一致。下载前按 `07-article-structure/01-title-and-paper-info.md` 的「写入前必须访问验证链接可达」核验原图链接。**禁止写「（此处原图）」「图X略」之类占位文字**。
- 确实要直引原博客远程图时，须先核验链接可达，并在交付说明里注明「远程图未入库」；但这只是兜底，默认仍下载本地下放仓库。
- 链接一律裸链接，禁止尖括号 `<https://...>`。翻译时**保留原文的内联超链接**：原文在论文标题、数据集/工具名、arXiv编号、webinar等处的链接，逐一以 `[文字](URL)` 形式补回对应位置（如 Pat Walters 博文里的《Practically significant method comparison protocols for machine learning in small molecule drug discovery》论文、herdr、TabPFN、各模型 arXiv 编号等）；URL 用 `[]()` 形式、禁止尖括号，且与相邻中文之间不空格。

## 6. 校验

翻译稿写完后同样跑 `tools/check_frontmatter.py --file <path>` 与 `tools/check_blog_quality.py <path>`；若因特殊原因用了远程图而被标为「图不存在」，在交付说明里注明原因即可，但默认仍应把图下载到本地、放进取仓库后再引用相对路径。
