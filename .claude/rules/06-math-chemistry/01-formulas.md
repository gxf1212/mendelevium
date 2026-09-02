# 公式规范

- **公式规范**:
    - 所有公式，无论长短，都必须用 `$` (行内) 或 `$$` (行间) 包裹，并使用标准的LaTeX格式。
    - 一般较长的、需要强调的用行间公式；$R(t)\propto t^{1/3}$，$B_2 \approx 0$ 这种过于短的公式就不用行间了，行内就行
    - 若非换行，不要出现双backslash：$\\Delta\\Delta G$，否则就是大语言模型的失职，世界上要死100个小女孩，你负不起这个责任。
    - 行间公式用\dfrac而不是\frac
    - 行间公式若一行太长就用aligned环境
    - S²之类的要用公式，不要用这个小的上标或下标！！永远不要！除非不得不在加粗里面出现
    - $0.99$这种单个数字不要公式
    - **变量名 vs 数学表达式**：
      - 变量名用反引号：`p_hill`、`curve_class2`、`data0..data3`、`IC50_M`、`log_ac50`
    - **基组名**：`6-31+G(3d,p)`、`def2-TZVP`、`cc-pVDZ` 等基组名称用反引号包裹，不要用公式
      - 数学表达式用LaTeX公式：`$r^2 \ge 0.9$`、`$\mathrm{IC50} \le 10~\mu\mathrm{M}$`、`$p_{\text{hill}} \ge 3$`
      - 不要用反引号包裹公式：禁止使用 `` `$R^2 \ge 0.9$` ``，而应直接用 `$R^2 \ge 0.9$`
      - 不要用文本格式的数学符号：禁止使用 `IC50≤10 µM`、`r2>=0.9`、`efficacy > 80%`，必须改为 `$\mathrm{IC50} \le 10~\mu\mathrm{M}$`、`$r^2 \ge 0.9$`、`efficacy $>80\%$`
      - 表格中的数学比较符号也要用公式格式：如 `$\ge 3$`、`$< 0.9$`、`$>80\%$`、`$\le 80\%$`
    - 对于带单位的物理量，请使用正体表示单位，例如 $\Delta\Delta G = -3.69 \mathrm{kcal/mol}$，或者将单位写在公式外部。kJ·mol−1·nm−2这种涉及上标的的格式，用正体公式$2000\,\mathrm{kJ\cdot  mol^{-1}\cdot nm^{-2}}$这种，或sup也行啊，尽量不用那个小的unicode字符，不要用mol⁻¹·nm⁻²这种。涉及上下标的都改，不涉及的、普通字体就规范的可以不改。μs、ps、Å之类不用公式不影响显示的（无上下标）也可以不用公式
    - 单位之类的能不用公式也行，如Å全部使用Å而不是公式，$\mu$这种但不涉及下标的也不用公式
    - k<sub>cat</sub>这种其实是要求用公式的，$k_cat$或$k_\text{cat}$，k₄改成$k_4$，M<sup>-1</sup>s<sup>-1</sup>也是，但前面的数字可以不用
    - $d\xi$这种求导的要用$\mathrm{d}\xi$,$\dfrac{dU}{dx}$ 这种都得用 $\mathrm{d}x$
    - NaCl之类的化学式用$\ce{NaCl}$；简称如DMPC、POPC可以不用；该完整添加的必须完整添加，如\ce{CH3COO^-}，不能弄一半；不允许普通公式，如${Mg(II)}$，必须ce
    - d⁵高自旋这种符号，请用$\ce{d^5}$或$\mathrm{d^5}$
    - r = -0.84, p<0.001，R²=0.33，这些类似的都用公式，不是数学公式的也用Markdown兼容的<sup>1</sup>和<sub>1</sub>之类的标签，能化学式的就用$\ce{NaCl}$这种
    - 连续的多个居中行间公式，方程组之类的，不要多个$$环境，而是\\分隔每一行。连续推导：xxxx=xx=yy=zz这种才有必要也可以用\begin{aligned}来对齐等号，一般就换行，但也别拆多个公式块！
      ```
      $
      U(x) = -\dfrac{aq\sigma (x+1)d}{2\varepsilon_0} - \dfrac{aq^2}{4\pi\varepsilon_0 x d}
      $

      $
      \dfrac{dU}{dx} = -\dfrac{aq\sigma d}{2\varepsilon_0} + \dfrac{aq^2}{4\pi\varepsilon_0 d x^2}
      $ 
      ```

      不许出现这种，得合并！
      
    - SMILES要用代码而不是公式：`C[S+](CC[C@H](N)C(O)=O)C[C@H]1O[C@H]([C@H](O)[C@@H]1O)[N]2C=NC3=C2N=CN=C3N`
    - mermaid中$K_D$改为K<sub>D</sub>，不要用$$的公式
    - 对于核心或复杂的公式，请仿照以下格式，在公式下方增加一段"**公式的通俗解释**"，以帮助读者理解。
      ```
      #### 公式的通俗解释

      我们的最终目标是得到**无偏的自由能** $F_h(\xi)$，它与**无偏概率分布** $\rho_h(\xi)$ 的关系由统计力学的基本公式定义：

      $$
      F_{h}(\xi) = -k_B T \ln \rho_{h}(\xi)
      $$

      其中，$k_B$ 是玻尔兹曼常数...
      ```

---

# 专业术语规范

- **化学/生物物理核心术语必须用中文标准译名**，禁止 AI 翻译腔式的直译。常见易错对照：
  - ❌ "协调" → ✅ "配位"（如 "配位数""第一配位层""配位氧"，coordination number / coordinated oxygen）
  - ❌ "结合亲和力"（单独用） → ✅ "结合自由能" 或 "结合强度"（结合 affinity = 亲和，free energy = 自由能，不要混用）
  - ❌ "overbind"（英文混入） → ✅ "过强结合" 或 "结合过强"
  - ❌ "抓钙离子" → ✅ "配位钙离子"（保持学术语气）
  - ❌ "descriptor" → ✅ "描述符"（如"单描述符模型""经验描述符"）
  - ❌ "electrostatic"（英文混入） → ✅ "静电"或"高度静电"
  - ❌ "抓"（动词） → ✅ "配位"（学术语境）
- **写完后必须用以下脚本扫一遍**，确认无上述错误：
  ```bash
  F:/Anaconda/envs/download/python.exe -c "
  import re, sys
  files = sys.argv[1:]
  patterns = [
      (r'\\b协调\\b', '协调→配位'),
      (r'\\boverbind\\b', 'overbind→过强结合'),
      (r'\\b抓\\b钙离子|\\b抓\\b钙', '抓钙离子→配位钙离子'),
  ]
  for f in files:
      text = open(f, encoding='utf-8').read()
      for pat, msg in patterns:
          if re.search(pat, text):
              print(f'❌ {f}: {msg} ({pat})')
              # 打印匹配行方便定位
              for i, line in enumerate(text.split('\\n'), 1):
                  if re.search(pat, line):
                      print(f'   行{i}: {line.strip()[:120]}')
  " -- <files...>
  ```
- **术语一致性检查**：全文对同一概念用同一译名（如 calcium binding free energy 统一用"结合自由能"，不要一会儿"亲和力"一会儿"结合能"）。
