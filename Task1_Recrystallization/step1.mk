## U-Mo 合金中裂变诱导气泡演化与再结晶相场模型摘要
该模型旨在模拟辐照 U-Mo 合金中由裂变引起的气泡形成、生长以及再结晶过程。模型结合了位错密度演化、多相场方法和气体扩散/反应动力学。

### 1. 核心变量
- 晶粒序参数 $\eta_i(r,t)$ ($i=1, \dots, p$) : 描述多晶结构，每个 $\eta_i$ 代表一个具有特定取向的晶粒。
- 气泡序参数 $\eta_j(r,t)$ ($j=p+1, \dots, q$) : 描述气泡相。
- 裂变气体 (Xe) 浓度 $c_g(r,t)$ : 描述氙原子在材料中的浓度分布。
- 平均位错密度 $\rho_d(t)$ : 描述材料中位错的累积程度，影响再结晶驱动力。
### 2. 控制方程
模型主要由以下耦合方程组描述：

- Allen-Cahn 方程 (序参数演化) : 描述晶粒和气泡相的演化。
  
  ```math
  \frac{\partial \eta_i}{\partial t} = -L_\eta \frac{\delta F}{\delta \eta_i}, \quad i = 1, 2, \ldots, q
   ```
  ```
  
  其中 $L_\eta$ 是界面动力学系数，$F$ 是体系的总自由能。
- Cahn-Hilliard 方程 (浓度演化) : 描述裂变气体 Xe 的扩散、产生和溶解。
  
  ```math
  \frac{\partial c_g}{\partial t} = \nabla \cdot \left( M(\eta) \nabla \frac{\delta F}{\delta c_g} \right) + \dot{G} - \dot{R}
   ```
  ```
  
  其中：
  
  - $M(\eta)$ 是与序参数相关的气体原子迁移率 (晶界处通常更高)。
  - $\frac{\delta F}{\delta c_g}$ 是化学势。
  - $\dot{G}$ 是裂变气体产生速率 (与裂变率 $\dot{f}$ 相关)。
  - $\dot{R}$ 是裂变诱导的气体溶解速率 (主要发生在气泡内部)。
- 位错密度演化 (简化模型) : 假设温度较低，忽略热退火，位错密度随裂变密度 $f_d = \dot{f} \cdot t$ 指数增长。
  
  ```math
  \rho_d(t) = \rho_d(0) \cdot \exp[\lambda \dot{f} \cdot t]
   ```
  ```
  
  其中 $\rho_d(0)$ 是初始位错密度，$\lambda$ 是位错累积速率。该 $\rho_d(t)$ 主要用于计算储能。
### 3. 总自由能泛函
体系的总自由能 $F$ 定义为：

```math
F(c_g, \eta_i) = \int \left[ f_{bulk}(c_g, \eta_i) + \sum_{i=1}^q \frac{\kappa_\eta}{2} |\nabla \eta_i|^2 + \frac{\kappa_c}{2} |\nabla c_g|^2 + f_{stored}(\eta_i) \right] dV
 ```
```

包含以下几项：

- $f_{bulk}$ : 体自由能密度，描述了各相的热力学稳定性以及它们之间的相互作用（晶粒-晶粒，晶粒-气泡，气体-界面）。其具体形式为：
  ```math
  \begin{aligned}
  f_{bulk}(c_{g},\eta_{1},\ldots,\eta_{q}) &= A\left[\sum_{i=1}^{q}\left(-\frac{1}{2}\eta_{i}^{2}+\frac{1}{4}\eta_{i}^{4}\right)+C_{p}\sum_{i=1}^{p}\sum_{j\neq i}^{p}\eta_{i}^{2}\eta_{j}^{2}\right] \\
  &+ C_q\sum_{i=p+1}^q\sum_{j=1}^p\eta_i^2\eta_j^2+\frac14+f_g^m\theta_m+f_g^b\theta_b
  \end{aligned}
   ```
  ```
   其中 $A, C_p, C_q$ 是常数，$f_g^m, f_g^b$ 是气体在基体和气泡中的化学自由能，$\theta_m, \theta_b$ 是基体和气泡的体积分数。
- $\frac{\kappa_\eta}{2} |\nabla \eta_i|^2$ : 序参数的梯度能项，与界面能相关。
- $\frac{\kappa_c}{2} |\nabla c_g|^2$ : 浓度的梯度能项。
- $f_{stored}$ : 位错储存能密度，是再结晶的主要驱动力之一。
  ```math
  f_{stored}(\eta_i) = \frac{1}{2} \eta_i^2 G b_v^2 \rho_d, \quad (i = 1, 2, \ldots, p)
   ```
  ```
   其中 $G$ 是剪切模量，$b_v$ 是柏氏矢量，$\rho_d$ 是当前位错密度 (假设各晶粒内相同)。该项仅对晶粒序参数 ($\eta_i, i=1..p$) 有贡献。
### 4. 关键物理过程建模
- 气体迁移率 $M(\eta)$ : 考虑了晶界扩散加速效应，形式为 $M(\eta)=M_b h(\eta)+M_g[1-h(\eta)]$，其中 $M_b, M_g$ 分别是基体和晶界的迁移率，$h(\eta)$ 是插值函数。
- 气体产生 $\dot{G}$ : 与裂变率 $\dot{f}$ 成正比，$\dot{G} = \varpi \cdot \text{Ran}$, $\varpi = \dot{f} \Omega$。
- 气体溶解 $\dot{R}$ : 发生在气泡内，$\dot{R} = \eta_j^2 \Lambda c_g$ ($j=p+1..q$)，$\Lambda$ 是溶解速率。
- 气泡形核 (可选) : 模型还描述了晶界处气泡形核的概率模型，基于经典形核理论，形核率 $J(t)$ 和形核概率 $P(t)$ 与局部过饱和度 $\Delta c$ 相关。