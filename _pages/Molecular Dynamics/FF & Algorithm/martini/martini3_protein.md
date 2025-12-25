  🔴 您的蛋白质散架的主要原因

  1. 最关键问题：缺少结构约束！

  您的protein_only.itp文件没有elastic network或Go-Martini约束！

  - ❌ 只有6个二硫键约束
  - ❌ 没有elastic network（ElNeDyn）
  - ❌ 没有Go-Martini约束
  - ❌ 没有位置限制（POSRES）

  根据网络搜索结果：

  "Martini proteins are often simulated in combination with an elastic network, such as Elnedyn, to maintain the overall structure" (https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/)

  "The defined network of restraints (elastic network) aims to maintain the structural and dynamical properties of a protein" (https://www.nature.com/articles/s41467-024-50647-9)

  2. MDP文件问题

  ;;;define = -DPOSRES ; position restrain  ← 被注释掉了！
  constraints = h-bonds                      ← 只约束氢键，对CG模型无效

  - ❌ 位置限制被禁用
  - ❌ constraints = h-bonds 对粗粒化模型无意义（CG模型没有氢原子）
  - ✅ epsilon_r = 78.4 是正确的（Martini 3推荐15，但78.4用于explicit水模型可接受）

  3. Martinize2命令问题

  martinize2 -f ... -dssp -cys auto -scfix -maxwarn 50

  缺少关键参数：
  - ❌ 没有 -elastic 或 -go 参数来生成elastic network
  - ❌ 没有 -maxwarn 应该避免（可能掩盖警告）

  ✅ 解决方案

  方案1：重新生成带Elastic Network的ITP（推荐）

  martinize2 -f protein_only.pdb \
    -ff martini3001 \
    -x protein_cg.pdb \
    -o protein.top \
    -elastic \                    # 添加elastic network
    -ef 500 \                     # elastic force constant (kJ/mol/nm²)
    -el 0.5 \                     # elastic lower cutoff (nm)
    -eu 0.9 \                     # elastic upper cutoff (nm)
    -from amber \
    -dssp \
    -cys auto

  方案2：使用Go-Martini（更灵活）

  martinize2 -f protein_only.pdb \
    -ff martini3001 \
    -x protein_cg.pdb \
    -o protein.top \
    -go \                         # 使用Go-Martini模型
    -goepsilon 9.414 \           # Go势能强度
    -from amber \
    -dssp \
    -cys auto


  Based on the latest research, here's a comprehensive comparison:

   Feature-by-Feature Comparison

  | Feature                      | Elastic Network (EN)        | Gō-Martini                     |
  |------------------------------|-----------------------------|--------------------------------|
  | Stability                    | ⭐⭐⭐⭐⭐ Most robust      | ⭐⭐⭐⭐ Good                  |
  | Flexibility                  | ⭐⭐ Limited                | ⭐⭐⭐⭐ High                  |
  | Conformational Changes       | ❌ Restricted               | ✅ Allowed                     |
  | Setup Complexity             | ✅ Simple                   | ⚠️ Requires tuning             |
  | High-throughput              | ✅ Excellent                | ⚠️ Limited                     |
  | Protein-protein interactions | ⚠️ May be inaccurate        | ✅ More realistic              |
  | Ligand binding studies       | ❌ Fixed tertiary structure | ✅ Captures structural changes |
  | Multi-domain proteins        | ✅ Works well               | ⚠️ Single monomers only        |

   Technical Details

  Elastic Network (ElNeDyn):
  - Mechanism: Harmonic bonds between backbone beads within cutoff distance
  - Recommended parameters (Martini 3):
    - Force constant: 700 kJ/mol/nm² (don't go below this!)
    - Cutoff: 0.8-0.9 nm
  - Key limitation: "The use of the elastic network restricts the use of the Martini force field for the study of large conformational changes" (https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/)

  Gō-Martini:
  - Mechanism: Contact map-based Lennard-Jones potentials (can break/reform)
  - Recommended parameters:
    - Potential depth: 10-12 kJ/mol
  - Advantage: "Unlike elastic networks which restrict the study of large conformational changes (e.g., folding), the GōMartini approach removes this limitation" (https://www.biorxiv.org/content/10.1101/2025.03.17.643608v2.full)

   Which Should You Choose?

  Use Elastic Network when:
  - ✅ You have a stable, folded protein with known structure
  - ✅ You want simple setup and robust performance
  - ✅ Running high-throughput screening
  - ✅ Studying protein-membrane interactions with fixed protein structure
  - ✅ Need maximum stability (e.g., membrane protein simulations)

  Use Gō-Martini when:
  - ✅ Studying conformational changes or flexibility
  - ✅ Investigating ligand-induced structural changes
  - ✅ Researching folding/unfolding processes
  - ✅ Need balance between stability and flexibility
  - ✅ Want more realistic protein-protein interactions
  - ⚠️ Only working with single monomeric units

  ⚠️ Important Warnings from Literature

  "Extensive tests have shown the combination of Martini with elastic networks may contribute to inaccurate protein-protein interactions" (https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/)

  "Systems with elastic networks show distributions shifted toward larger cluster sizes, and elastic networks can increase protein stickiness" (https://pmc.ncbi.nlm.nih.gov/articles/PMC6785803/)


• - 重新运行 martinize2 … -elastic -ef 700 -el 0.5 -eu 0.9 -eunit chain 后，
    protein_only.itp 的内容从“只包含 include”变成了完整的 MARTINI 3 定义，并新增了一个
    名为 Elastic Network 的 [ bonds ] 区段。该区段列出了所有需要保持原子距离的 bead
    对，每条线通常形如 i  j  6  0.700  0.500  0.900：
      - i / j：被弹性弹簧连接的两个 beads（例如 backbone BB 之间或跨链接触对）。
      - funct=6：GROMACS 中的 harmonic bond，但带有限制距离（ElNeDyn 形式）。
      - 0.700：力常数 k = 700 kJ mol⁻¹ nm⁻²（-ef 700）。
      - 0.500 / 0.900：下/上截止（-el 0.5, -eu 0.9，单位 nm）。当真实距离 r < 0.5 nm 时
        弹簧被“夹断”，0.5–0.9 nm 区间呈 Hooke 行为，>0.9 nm 则不再施力。这能保持折叠但
        又允许一定柔性。
  - 其它段落与原始 martinize2 输出相同，但现在 全部都在同一个文件中：

      - 若存在 [virtual_sites]、[exclusions] 等也在同一文件里，用于虚拟点、排除列表。
      - 文件头部注释里保留了 martinize2 命令行，方便追溯参数。
  
  1. 新增"Rubber band"部分：
  - 1489个弹性网络键
  - 每个键的格式：原子1 原子2 1 距离(nm) 力常数(kJ/mol/nm²)
  - 例如：1 7 1 0.60982 700.0
    - 连接原子1和原子7
    - 平衡距离：0.60982 nm
    - 力常数：700 kJ/mol/nm² (用户设置的-ef 700)
  - 
  - 总结：加入弹性网络后，protein_only.itp 多了上百条额外弹簧，这些弹簧在 0.5–0.9 nm 内
    提供 700 kJ mol⁻¹ nm⁻² 的恢复力，等价于 MARTINI 官方教程中的 ElNeDyn 设定，可防止
    KLK5 在 CG 模型下“散架”。同时我们把完整参数都放回 protein_only.itp，后续脚本只需
    #include "protein_only.itp" 就能获得所有定义，无需再额外 include PROA_0.itp。
  2. 调和势能函数（Harmonic potential）

  弹性网络键使用调和势：

  V(r) = ½ k (r - r₀)²

  其中：
  - k = 700 kJ/mol/nm² （力常数，-ef 700参数）
  - r₀ = 平衡距离（从初始结构测量）
  - r = 当前模拟中的实际距离

  3. 实际作用

  模拟过程中：
  - 如果距离 r ≈ r₀：势能很小，几乎没有约束力
  - 如果距离 r > r₀：产生向内的拉力（像拉伸的橡皮筋）
  - 如果距离 r < r₀：产生向外的推力（像压缩的弹簧）

  力的大小：F = -k(r - r₀)

  4. 为什么这样设计？

  - 保持初始结构的整体形状：通过记住初始结构中所有接近原子对的距离
  - 允许局部灵活性：力常数700比骨架共价键(4000)弱得多，允许一定程度的构象变化
  - 防止蛋白质展开：如果蛋白质开始展开，大量的弹性网络键会同时被拉伸，产生很大的回复力

  5. 与二硫键对比

  # 二硫键（共价键，非常强）
  12  13 1 0.24 None  # CYS-CYS, 距离固定在0.24 nm

  # 弹性网络键（弱约束）
  1   7 1 0.60982 700.0  # 距离可变，但偏离0.60982会有能量惩罚

  所以平衡距离确实就是从初始PDB结构测量的快照值，模拟过程中这个值不会改变，但实际距离会围绕这个平衡值波动！


• 弹性网络（Elastic Network / ElNeDyn）

  - 原理：在 martini 粗粒化蛋白中，真实的氢键与二级结构势能被削弱，因而需要在所有同链
    backbone (BB) 珠之间按空间 cutoff 加“弹簧”。martinize2 的 -elastic 选项会扫描初始结
    构，凡是彼此相距 <0.8–0.9 nm、且序列差≥2 的 BB 对都加一个谐振子，默认劲度常数
    700 kJ·mol⁻¹·nm⁻²，力只在 0.5–0.9 nm 内生效，超出则失效，从而锁住折叠但保留局部柔
    性。(vermouth-martinize.readthedocs.io (https://vermouth-martinize.readthedocs.io/
    en/stable/tutorials/elastic_networks.html?utm_source=openai))
  - 经典 ElNeDyn（2009 JCTC）把这种结构弹簧与 Martini 2.1 非键相互作用耦合，实现大型蛋
    白、病毒衣壳等微秒级模拟，并证明 EN+物理力场的组合能重现原子级的集体振动、蛋白结合
    等性质。(cgmartini.nl (https://cgmartini.nl/docs/downloads/tools/proteins-and-
    bilayers.html?utm_source=openai))
  - 局限：弹簧不可断，限制大尺度构象改变，且对多域/IDR 需要手动移除部分弹簧或降低劲度，
    否则会“绷死”蛋白。GoMartini 等后续方案正是为了解决此缺陷。(nature.com (https://
    www.nature.com/articles/s41467-025-58719-0?utm_source=openai))

  GōMartini / GōMartini 3

  - 结构—物理混合：先用 martinize2 生成蛋白，再用 contact map（OV+rCSU）筛出原生接触。
    新版 GōMartini 3（2025 Nature Communications）把每个 backbone 复制出虚拟站点，只在
    这些虚拟点之间施加 Lennard‑Jones 势，势谷深度 ε≈9.4 kJ·mol⁻¹，可恢复特定 native
    contacts 而不需要通用弹簧。(nature.com (https://www.nature.com/articles/s41467-025-
    58719-0?utm_source=openai))
  - 虚拟站点优势：这些 LJ 势走非键邻域，因此享受截断和并行效率，解决了旧版把 pair
    potential 当“键”处理、难以并行的瓶颈；同时还能附加“水偏置”，即调节虚拟站点与
    Martini 水 bead 的 ε，用来矫正 IDP 过度塌缩、跨膜肽易被水拉出的现象。(nature.com
    (https://www.nature.com/articles/s41467-025-58719-0?utm_source=openai))
  - 应用：通过 tuning BB–water 偏置（例如 α 螺旋 ε = −1.0 kJ·mol⁻¹、IDP ε =
    +0.5 kJ·mol⁻¹），可以在不改动主力场的前提下，稳住 WALP 跨膜构型、让凝聚相含水量接近
    实验，也能用“只加水偏置、不加 Gō”的模式，兼容 EN 或其他约束。(nature.com (https://
    www.nature.com/articles/s41467-025-58719-0?utm_source=openai))
  - 适用场景：需要可逆 unfolding/大幅构象变化（AFM 拉伸、signal transduction）时，
    GōMartini 比 EN 更灵活；但它依赖参考结构，若想模拟未知折叠或多稳态，需要改造成多盆
    Gō 或动态更新 contact map。(nature.com (https://www.nature.com/articles/s41467-025-
    58719-0?utm_source=openai))

  OLIVES（Hydrogen-Bond Native Contact Network）

  - 研究动机：Martini 3 缺乏显式氢键能，导致蛋白在无偏置时迅速散架；传统 EN/Gō 虽能稳住
    结构，但会引入大量弹簧或对所有接触“一视同仁”。OLIVES（2024 JCTC）提出自动识别 CG 蛋
    白中具有氢键潜力的原生接触，并只对这些对施加 Go-like LJ 势，势深来自 ab initio 氢键
    能（~2–5 kcal/mol），因此既补上缺失的氢键，又能大幅减少偏置项。
    (pubmed.ncbi.nlm.nih.gov (https://pubmed.ncbi.nlm.nih.gov/39235392/?
    utm_source=openai))
  - 算法：先用 martinize2 得到 CG 结构，再用 OLIVES 脚本扫描所有 donor/acceptor，判断是
    否形成 H-bond（结合距离、角度和可及性），把通过筛选的接触写回 itp；脚本还可输出独立
    的 .itp，供需要与 EN 或 GoMartini 混用时灵活 include。(github.com (https://
    github.com/Martini-Force-Field-Initiative/OLIVES?utm_source=openai))
  - 性能：OLIVES 用真实 H-bond 能量恢复了折叠稳定性，且因偏置数量更少，GPU 上比传统
    GoMartini 快约 30%；对多聚体测试也显示无需额外的四级结构偏置。
    (pubmed.ncbi.nlm.nih.gov (https://pubmed.ncbi.nlm.nih.gov/39235392/?
    utm_source=openai))
  - 何时选 OLIVES：当体系需要“只补氢键、不限制大尺度运动”（例如需要蛋白可张开/闭合）
    时，OLIVES 是介于纯 EN 与全接触 Go 的折中方案；若蛋白包含大量非氢键主导的疏水接口，
    仍可与 EN 或 Go 共用。

  对比与实务建议

  1. 维持稳定折叠：纯 EN 成本最低、设置最简单，适合刚体蛋白或只看局部相互作用。若蛋白需
     呼吸式运动，可改用 GoMartini 或 OLIVES，甚至把 EN 只保留在核心、在柔性环/多域处禁
     用。(nature.com (https://www.nature.com/articles/s41467-025-58719-0?
     utm_source=openai))
  2. 可逆构象变化：GoMartini 通过 LJ 原生接触而非硬弹簧，可允许接触断裂再形成；若只关注
     氢键驱动的结构，则 OLIVES 更轻量。(nature.com (https://www.nature.com/articles/
     s41467-025-58719-0?utm_source=openai))
  3. 速度 / 并行：虚拟站点版 GoMartini 与 OLIVES 都兼容 GROMACS 域分解并可 GPU 加速；EN
     因只添加谐振子，代价也低，但一旦大分子需要跨域弹簧，仍可能影响并行效率。
     (nature.com (https://www.nature.com/articles/s41467-025-58719-0?
     utm_source=openai))
  4. 实现流程：
      - EN：martinize2 ... -elastic -ef 700 -el 0.5 -eu 0.9。
      - GoMartini：martinize 时加 --govs-include/--govs-moltype，再用 create_goVirt 生
        成虚拟站点与接触；如需水偏置，追加 --bias_* 或 --idp_*。(github.com (https://
        github.com/Martini-Force-Field-Initiative/GoMartini?utm_source=openai))
      - OLIVES：常规 martinize 后运行 OLIVES_v2.0_M3.0.0.py -c protein_cg.pdb -i
        protein.itp，脚本会插入 H-bond 接触项，可与 EN/Go 叠加。(github.com (https://
        github.com/Martini-Force-Field-Initiative/OLIVES?utm_source=openai))

  依据体系目标（保持整体 vs. 允许大尺度运动）、对物性修正的需求（氢键补偿、水偏置）以及
  可接受的建模工作量，在 EN、GoMartini 及 OLIVES 之间做组合，可让 Martini 3 的蛋白模拟
  既稳定又不过度“捆绑”。
