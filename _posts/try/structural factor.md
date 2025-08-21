# 结构因子简介

## 1 什么是结构因子 *S(q)*  

结构因子 $$S(\mathbf q)$$ 描述波矢 $$\mathbf q$$（或对应空间尺度 $$2\pi/|\mathbf q|$$）上密度涨落的强弱，是连接模拟与小角 X-射线/中子散射 (SAXS/SANS) 实验的重要桥梁。对于各向同性体系，$$S(\mathbf q)$$ 仅依赖模长 $$q=|\mathbf q|$$，通常记为 $$S(q)$$。

---

## 2 从密度涨落到结构因子  

### 2.1 瞬时密度的傅里叶分量  

体系在时刻 $$t$$ 的单粒子密度  

$$
\rho(\mathbf r,t)=\sum_{i=1}^{N}\delta\!\bigl(\mathbf r-\mathbf r_i(t)\bigr)
$$

其在波矢 $$\mathbf q$$ 上的傅里叶分量定义为  

$$
\rho_{\mathbf q}(t)=\int_{\text{box}}
e^{-i\mathbf q\cdot\mathbf r}\,
\rho(\mathbf r,t)\,\mathrm d\mathbf r
\tag{1}
$$

### 2.2 结构因子  

对准静态体系在时间窗口 $$\Delta t$$ 内取平均，有  

$$
S(\mathbf q,t)=
\frac{\bigl\langle\rho_{-\mathbf q}(t)\,\rho_{\mathbf q}(t)\bigr\rangle_{\Delta t}}{N}
\tag{2}
$$

在各向同性体系中，对同模长壳平均得到标量函数 $$S(q,t)$$。

---

## 3 高分子体系的计算细节  

1. **粒子选择**  
   * 原子级模拟：常选重原子或质心；  
   * 粗粒化链：选代表性 *bead*（如主链 Cα）。  

2. **时间平均** 
   在 *quench* 模拟中，常取前 2000 帧作为统计时间段；可按 1/100 – 1/200 的轨迹长度滑动窗口取样以平滑噪声。  

3. **径向平均** 
   在离散 $$N_xN_yN_z$$ 网格上映射密度并计算 (1) 式后，将 $$\mathbf q$$-空间划壳并平均获得 $$S(q,t)$$。  

4. **数值实现** 
   使用 FFT (FFTW、cuFFT 等) 加速傅里叶变换，并去除 $$\mathbf q=0$$ 项以免整体漂移贡献。

---

## 4 结构因子的一阶矩与特征长度  

令  

$$
\langle q\rangle(t)=
\frac{\displaystyle\int_{0}^{\infty} q\,S(q,t)\,\mathrm d q}
     {\displaystyle\int_{0}^{\infty}   S(q,t)\,\mathrm d q}
\tag{3}
$$

其倒数（取 $2\pi$ 因子方便与实空间长度对应）给出瞬时特征长度  

$$
\ell(t)=\frac{2\pi}{\langle q\rangle(t)}
\tag{4}
$$

---

## 5 粗化动力学的幂律行为  

在相分离或网络重排过程中，$$\ell(t)$$ 往往遵循幂律 $$\ell(t)\propto t^{\alpha}$$，不同指数 $$\alpha$$ 反映不同主导机制：

| 指数 $$\alpha$$ | 物理机制                                                     | 典型体系/条件                                                |
| :-------------: | :----------------------------------------------------------- | :----------------------------------------------------------- |
|  $$\tfrac13$$   | **扩散受限 (Lifshitz–Slyozov–Wagner, LSW)**<br>界面张力驱动，质量通过体扩散在液滴间传输；<br>特征时间尺度 $$t_D\sim R^3/D$$ | 高黏度溶剂或无水动力学耦合的聚合物/胶体液滴相分离            |
|  $$\tfrac12$$   | **粘性-水动力耦合粗化**<br>内摩擦较小、流体可流动；<br>Navier–Stokes 粘性项主导界面演化 | 聚电解质共沉淀网络 (Yuan & Tanaka, 2025)、高分子-溶剂体系在中等雷诺数下 |
|      $$1$$      | **惯性-水动力主导**<br>惯性项不可忽略，界面张力直接驱动物料流 | 低黏度二元流体、相分离后期“涡旋合并”阶段                     |
|      $$0$$      | **自限化**<br>当玻色化或交联网限制进一步粗化时，长度饱和     | 化学交联凝胶、永久网络                                       |

### $$\alpha=\tfrac13$$ 幂律  

* 液滴半径 $$R(t)$$ 的演化由溶解-沉积 (i.e. 奥斯瓦尔德熟化过程) 控制：  
  $$
    \frac{\mathrm dR}{\mathrm dt} \propto \frac{1}{R^2},
    \qquad\Rightarrow\qquad R(t)\propto t^{1/3}.
  $$
* 结构因子主峰 $$q_\mathrm p(t)$$ 近似与逆特征长度成正比 ($$\ell\approx 2\pi/q_\mathrm p$$)，因此 $$q_\mathrm p(t)\propto t^{-1/3}$$，对应公式 (3) 的一阶矩 $$\langle q\rangle(t)$$ 同样呈 $$t^{-1/3}$$ 衰减。  
* 适用条件：  
  - 动量耗散快，水动力互作被抑制；  
  - 典型如**长链聚合物/胶体在良溶剂中的液滴-连续相分离初期**

---

## 6 研究案例  

**Yuan & Tanaka (2025)** 研究了带异号电荷的柔性聚电解质在溶剂中形成 **coacervate 网络**：  

* 体系首先产生贯通网络而非孤立液滴；  
* 采用公式 (4) 得到 $$\ell(t)\propto t^{1/2}$$，指示水动力耦合粗化；  
* 提升电荷不对称性导致 $$S(q)$$ 主峰向低 $$q$$ 移动且粗化速度减慢；  
* SAXS 实验观测到一致趋势，为模拟提供验证。

---

## 参考文献  

[1] **Yuan J.; Tanaka H.** *Network-forming phase separation of oppositely charged polyelectrolytes forming coacervates in a solvent*. **Nat. Commun.** **2025**, *16*, 1517.  

[2] Hansen, J.-P.; McDonald, I. R. *Theory of Simple Liquids*, 4th ed.; Academic Press, 2013.  