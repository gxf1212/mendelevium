# David Tong Fluid Mechanics Note -- 1

##  1. 前提与符号约定

- 空间变量：$\mathbf{x}$，时间：$t$。

- 速度场：$\mathbf{u}(\mathbf{x},t)$。

- 密度：$\rho(\mathbf{x},t)$。

- 压强：$p(\mathbf{x},t)$。

- 物质导数（随流体团块的导数），场$\phi(x(t), t)$对时间的导数：
  $$
  \frac{d}{dt}\phi(x(t),t)=\partial_t\phi+\dot{\mathbf x}\cdot \nabla\phi
  $$


  因此将物质导数的算符写为：$$\displaystyle \frac{D}{Dt}=\frac{\partial}{\partial t}+\mathbf{u}\cdot\nabla.$$

- 本章节主要讨论**无粘性 (inviscid)** 流体；在大部分段落假设**不可压** (incompressible)，i.e.  $\rho= const$ ，且 $\nabla\cdot\mathbf{u}=0$。

------

## 2. 质量守恒（连续方程）

### 2.1 积分形式与局域微分形式

固定控制体 $V$ 内总质量：
$$
M(t)=\int_V \rho(\mathbf{x},t),dV
$$
对时间求导并用通量概念得：
$$
\frac{dM}{dt}=\int_V \frac{\partial\rho}{\partial t},dV = -\int_V \nabla\cdot(\rho\mathbf{u}),dV= -\int_S \rho\mathbf{u}\cdot d\mathbf{S}
$$
因为对任意 $V$ 成立，得到微分形式：
$$
\boxed{\displaystyle \frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf{u})=0}
$$
改写为物质导数形式：
$$
\boxed{\frac{D}{Dt}\rho - \rho\nabla\cdot\mathbf u=0}
$$

### 2.2 不可压缩极限

若 $\rho$ 恒定或随团块不变，同时因为质量守恒，物质导数 $D\rho/Dt=0$，由上式可以得：
$$
\boxed{\nabla\cdot\mathbf{u}=0}
$$

------

## 3. 动量守恒与 Euler 方程

### 3.1 流体中的牛顿第二定律写法（局域形式）

对单位体积的动量守恒（外力来源主要为压强梯度）有：
$$
\rho\frac{D\mathbf{u}}{Dt}= -\nabla p + \mathbf{f}_{\text{body}}
$$
其中 $\mathbf{f}_{\text{body}}$ 表示体积力（例如重力 $\rho\mathbf{g}$）。在本节主要先讨论无体积力或将其单独处理的情形，所得无粘、无体积力的 Euler 方程为：
$$
\boxed{\frac{D}{Dt} \mathbf u = \displaystyle \frac{\partial\mathbf{u}}{\partial t}+ (\mathbf{u}\cdot\nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p}
$$

### 3.2 非线性项的恒等变形（引入旋度）

使用向量恒等式：
$$
(\mathbf{u}\cdot\nabla)\mathbf{u} = \nabla\left(\frac12|\mathbf{u}|^2\right) - \mathbf{u}\times(\nabla\times\mathbf{u})
$$
令涡量 $\boldsymbol{\omega}=\nabla\times\mathbf{u}$，代入Eular方程得到：
$$
\displaystyle \frac{\partial\mathbf{u}}{\partial t} - \mathbf{u}\times\boldsymbol{\omega} = -\nabla\left(\frac{p}{\rho} + \frac12|\mathbf{u}|^2\right)
$$
定义伯努利函数：
$$
\displaystyle H\equiv \frac{p}{\rho}+\frac12|\mathbf{u}|^2
$$
于是Eular方程可进一步简写为：
$$
\boxed{\displaystyle \frac{\partial\mathbf{u}}{\partial t} - \mathbf{u}\times\boldsymbol{\omega} = -\nabla H}
$$
**物理注解**：项 $-\mathbf{u}\times\boldsymbol{\omega}$ 改变速度方向但不改变速度大小（它与 $\mathbf{u}$ 垂直），故不对动能做功。

------

## 4. 伯努利定理与能量守恒（定常与无旋两种情形）

### 4.1 能量方程（把 Euler 与速度点乘）

对Eular方程的一般形式两边点乘 $\mathbf{u}$：
$$
\mathbf{u}\cdot\frac{\partial\mathbf{u}}{\partial t}+\mathbf u\cdot (\mathbf u \times \boldsymbol \omega)=\cdot\frac{\partial(\frac 12 |\mathbf{u}|^2)}{\partial t}+\mathbf u\cdot (\mathbf u \times \boldsymbol \omega)=-\mathbf u \cdot \nabla \mathcal H=-\mathbf u \cdot \nabla(\frac{p}{\rho}+\frac12|\mathbf{u}|^2)
$$
 $\mathbf{u}\cdot(\mathbf{u}\times\boldsymbol{\omega})=0$，上式可写为物质导数形式：
$$
\boxed{\displaystyle \frac{D}{Dt}\left(\frac12|\mathbf{u}|^2\right) = -\frac{1}{\rho}\mathbf{u}\cdot\nabla p}
$$
再使用恒等式 $\mathbf{u}\cdot\nabla p = \nabla\cdot(p\mathbf{u}) - p(\nabla\cdot\mathbf{u})$，在不可压 ($\nabla\cdot\mathbf{u}=0$) 情形：
$$
\displaystyle \frac{D}{Dt}\left(\frac12|\mathbf{u}|^2\right) = -\frac{1}{\rho}\nabla\cdot(p\mathbf{u})
$$
对任意体积分并利用散度定理可得动能守恒的积分形式（压力通量项表征压力所做的功）：
$$
\frac{d}{dt}\int_V \frac12\rho|\mathbf{u}|^2\ dV = -\int_S p\ \mathbf{u}\cdot d\mathbf{S}
$$

### 4.2 伯努利公式（定常流）

若为定常流 $\partial\mathbf{u}/\partial t=0$，则Eular公式变为
$$
\mathbf{u}\times\boldsymbol{\omega} = \nabla H
$$
两边同时点成速度场$\mathbf u$，沿流线取导：$\mathbf{u}\cdot\nabla H = 0$，因此：
$$
\boxed{H\text{ 在每条流线上为常数。}}
$$
若进一步假定 **无旋**（$\boldsymbol{\omega}=0$），则 $\nabla H=0$，得到
$$
\boxed{H=\text{常数（全域）}}
$$
即在无旋定常流中，伯努利函数在整个流场内恒为同一常数。

------

## 5. 可压缩流修正

若流体可压，$\rho$ 不再恒定，之前的 $p/\rho$ 项不再能表示单位质量的压强能量。更一般地，在恰当的热力学约束（例如不可热交换或等熵流）下，应使用**比焓** $h$：
$$
\displaystyle dh=\frac{dp}{\rho}
$$
从而广义的伯努利函数为：
$$
\boxed{H = h + \frac12|\mathbf{u}|^2}
$$
对于理想气体，$h=c_p T$，故可写为 $H=c_pT+\frac12|\mathbf{u}|^2$（在等熵下）。

------

## 6. Euler 方程的守恒形式与应力张量

将质量守恒与动量守恒结合:
$$
\nabla\cdot\mathbf{u}=0\\
\frac{D}{Dt} \mathbf u = \displaystyle \frac{\partial\mathbf{u}}{\partial t}+ (\mathbf{u}\cdot\nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p
$$
动量守恒的第二项可以改写为：
$$
(\mathbf{u}\cdot\nabla)\mathbf{u}=\nabla\cdot (\mathbf u \mathbf u)-\mathbf u \nabla \cdot \mathbf u
$$
带入质量守恒的表达式：
$$
(\mathbf{u}\cdot\nabla)\mathbf{u}=\nabla\cdot (\mathbf u \mathbf u)
$$
我们可以把动量方程写成严格的守恒律形式：
$$
\displaystyle \frac{\partial}{\partial t}(\rho\mathbf{u}) + \nabla\cdot(\rho\mathbf{u}\mathbf{u}+p\mathbf{I}) = 0
$$
此处动量通量张量 $T=\rho\mathbf{u}\mathbf{u}+p\mathbf{I}$，其中 $p\mathbf{I}$ 即压力贡献；该结构与场论中的应力-能量张量相似:
$$
T_{ij}=\rho u_i u_j+p_{ij}\delta_{ij}
$$
Euler方程的张量形式可以改写为：
$$
\rho \partial_t u_i=\partial_jT_{ij}
$$
上式应用Einstein约定求和，$j$为哑指标

------

## 7. 浮力（Archimedes 原理）的推导

对被流体完全浸没的固体体积 $V$，受到的压力力为：

$$\mathbf{F}=-\oint_{\partial V} p,d\mathbf{S} = -\int_V \nabla p,dV. $$

在静水（或静力平衡）条件下 $\nabla p=\rho\mathbf{g}$，因此$$\mathbf{F} = -\int_V \rho\mathbf{g},dV = -\rho_{\text{fluid}} V \mathbf{g}$$，即浮力等于被排开的流体的重量。
