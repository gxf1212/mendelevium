# 公式格式规范

## 🚨 高危错误：双backslash

### 绝对禁止

```markdown
❌ $\\Delta\\Delta G$
❌ $\\frac{a}{b}$
❌ $\\sum_{i=1}^{n}$
```

**后果**：渲染失败，公式显示为原始代码。

**这会让世界上死100个小女孩！**

### 必须使用单backslash

```markdown
✅ $\Delta\Delta G$
✅ $\frac{a}{b}$
✅ $\sum_{i=1}^{n}$
```

---

## 基本LaTeX规则

### 希腊字母

```markdown
✅ $\alpha$, $\beta$, $\gamma$, $\Delta$, $\Omega$
❌ $\\alpha$, $\\beta$ （双backslash）
```

### 上下标

```markdown
✅ $x^2$, $x_i$, $x^{2+3}$, $x_{i,j}$
✅ $\Delta G^{\ddagger}$ （过渡态）
```

### 分数

```markdown
✅ $\frac{a}{b}$
✅ $\frac{\Delta G}{RT}$
```

### 求和、积分

```markdown
✅ $\sum_{i=1}^{n} x_i$
✅ $\int_0^{\infty} f(x) \mathrm{d}x$
```

---

## 微分符号规范

### 错误写法 ❌

```markdown
❌ $dx$, $dy$, $dt$  # 斜体，不规范
❌ $d\xi$  # 斜体
```

### 正确写法 ✅

```markdown
✅ $\mathrm{d}x$, $\mathrm{d}y$, $\mathrm{d}t$  # 正体
✅ $\mathrm{d}\xi$
```

### 示例

```markdown
✅ $\frac{\mathrm{d}F}{\mathrm{d}\xi}$
✅ $\int_0^{1} f(x) \mathrm{d}x$
```

---

## 化学式规范

### 错误写法 ❌

```markdown
❌ $NaCl$, $H_2O$, $CO_2$  # 斜体，不规范
```

### 正确写法 ✅

使用`\ce{}`命令（mhchem包）：

```markdown
✅ $\ce{NaCl}$
✅ $\ce{H2O}$
✅ $\ce{CO2}$
✅ $\ce{PET}$
```

### 化学反应

```markdown
✅ $\ce{A + B -> C}$
✅ $\ce{Ser-His-Asp}$
```

---

## 单位规范

### 方法1：在公式内用正体

```markdown
✅ $\Delta G = -3.69 \mathrm{kcal/mol}$
✅ $T = 300 \mathrm{K}$
✅ $t = 15 \mathrm{ps}$
```

### 方法2：单位写在公式外

```markdown
✅ $\Delta G = -3.69$ kcal/mol
✅ $T = 300$ K
✅ $t = 15$ ps
```

### 常见单位

- 能量：`\mathrm{kcal/mol}`, `\mathrm{kJ/mol}`
- 温度：`\mathrm{K}`
- 时间：`\mathrm{ps}`, `\mathrm{ns}`, `\mathrm{fs}`
- 距离：`\mathrm{Å}`, `\mathrm{nm}`
- 角度：`\mathrm{deg}`, `^\circ`

---

## 行内公式 vs 行间公式

### 行内公式

```markdown
文中提到 $\Delta G = -3.69$ kcal/mol。
```

**使用场景**：
- 简短的公式
- 不超过一行
- 不是文章的核心公式

### 行间公式

```markdown
$$
\Delta G^{\ddagger} = -RT \ln\frac{k_{\mathrm{cat}} h}{k_{\mathrm{B}} T}
$$
```

**使用场景**：
- 重要的公式
- 较长的公式
- 需要突出展示的公式

---

## 多行公式

### 使用aligned环境

```markdown
$$
\begin{aligned}
\Delta G^{\ddagger} &= -RT \ln\frac{k_{\mathrm{cat}} h}{k_{\mathrm{B}} T} \\
&= -0.603 \times 303 \ln\frac{0.08 \times 6.626 \times 10^{-34}}{1.381 \times 10^{-23} \times 303} \\
&= 18.6 \mathrm{kcal/mol}
\end{aligned}
$$
```

**注意**：
- 使用`&`对齐等号
- 使用`\\`换行（行间公式中可以，行内不行）

---

## 公式解释模板

对于核心或复杂的公式，添加"公式的通俗解释"：

```markdown
#### 公式的通俗解释

我们的最终目标是得到**无偏的自由能** $F_h(\xi)$，它与**无偏概率分布** $\rho_h(\xi)$ 的关系由统计力学的基本公式定义：

$$
F_{h}(\xi) = -k_B T \ln \rho_{h}(\xi)
$$

其中：
- $k_B$ 是玻尔兹曼常数（$1.381 \times 10^{-23}$ J/K）
- $T$ 是温度（单位：K）
- $\rho_h(\xi)$ 是归一化的概率密度
- $\xi$ 是反应坐标（单位：Å）
```

---

## 特殊符号

### 常用符号

```markdown
✅ $\pm$ （加减号）
✅ $\times$ （乘号，不要用×）
✅ $\approx$ （约等于）
✅ $\leq$, $\geq$ （小于等于、大于等于）
✅ $\neq$ （不等于）
```

### 希腊字母表

| 小写 | 大写 | LaTeX |
|------|------|-------|
| α | Α | `\alpha` |
| β | Β | `\beta` |
| γ | Γ | `\gamma` |
| δ | Δ | `\delta`, `\Delta` |
| ξ | Ξ | `\xi`, `\Xi` |
| ρ | Ρ | `\rho` |
| σ | Σ | `\sigma`, `\Sigma` |
| ω | Ω | `\omega`, `\Omega` |

---

## 检查方法

### 自动检查

使用`tools/check_blog_quality.py`会检测：
- 双backslash（高危错误）
- 可能的微分符号问题

### 手动检查

在编辑器中搜索：
1. `\\\\` - 查找双backslash（在编辑器中显示为\\）
2. `\$[^\$]*\bd[A-Za-z]` - 查找可能的微分符号

---

## 常见错误汇总

| 错误写法 | 正确写法 | 错误原因 |
|---------|---------|---------|
| `$\\Delta$` | `$\Delta$` | 双backslash |
| `$dx$` | `$\mathrm{d}x$` | 微分符号斜体 |
| `$NaCl$` | `$\ce{NaCl}$` | 化学式斜体 |
| `$3.69 kcal/mol$` | `$3.69 \mathrm{kcal/mol}$` 或 `$3.69$ kcal/mol` | 单位斜体 |
| `$H_2O$` | `$\ce{H2O}$` | 化学式斜体 |

---

## 完整示例

```markdown
本研究使用**伞形采样方法**计算了PETase催化反应的自由能曲线。酰化步骤的活化自由能为：

$$
\Delta G^{\ddagger}_{\mathrm{acylation}} = 20.0 \pm 0.5 \mathrm{kcal/mol}
$$

这个值与实验测得的表观活化能（$18.0 \mathrm{kcal/mol}$）非常接近，差异在DFT方法的预期误差范围内（$\pm 3 \mathrm{kcal/mol}$）。

#### 公式的通俗解释

活化自由能 $\Delta G^{\ddagger}$ 描述了反应物转变为过渡态所需的能量壁垒。该值越低，反应越容易发生。本研究计算得到的 $20.0 \mathrm{kcal/mol}$ 表明，在室温（$300 \mathrm{K}$）下，PETase的催化效率为：

$$
k_{\mathrm{cat}} = \frac{k_B T}{h} \exp\left(-\frac{\Delta G^{\ddagger}}{RT}\right) \approx 0.08 \mathrm{s}^{-1}
$$

其中 $k_B$ 是玻尔兹曼常数，$h$ 是普朗克常数，$R$ 是气体常数。
```
