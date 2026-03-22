# Belousov–Zhabotinsky (BZ) Reaction Modeling Guide

This note gives a practical path from nonlinear chemistry to a C++ reaction–diffusion simulation.

## 1) Nonlinear kinetics: Oregonator as a reduced mechanism

A widely used reduced BZ model is the **Oregonator**, which keeps the core feedback loop:

- autocatalytic production of activator species,
- delayed inhibition,
- recovery process.

A common nondimensional 3-variable form is:

\[
\begin{aligned}
\frac{\partial u}{\partial t} &= \frac{1}{\epsilon}\left[u(1-u)-f v\frac{u-q}{u+q}\right] + D_u\nabla^2 u,\\
\frac{\partial v}{\partial t} &= u-v + D_v\nabla^2 v,\\
\frac{\partial w}{\partial t} &= \phi(u-w) + D_w\nabla^2 w.
\end{aligned}
\]

Interpretation (coarse-grained):

- \(u\): activator (roughly related to bromous acid / fast oxidation branch),
- \(v\): inhibitor (bromide-controlled suppression branch),
- \(w\): recovery/oxidized catalyst-like slow variable.

Key nonlinearities:

- \(u(1-u)\): growth + saturation,
- \(v\frac{u-q}{u+q}\): inhibition with threshold-like behavior,
- small \(\epsilon\): separation of timescales (stiff fast chemistry).

## 2) ODE model (well-mixed reactor)

For a perfectly stirred reactor (no spatial gradients), set diffusion terms to zero:

\[
\dot{u}=\frac{1}{\epsilon}\left[u(1-u)-f v\frac{u-q}{u+q}\right],\quad
\dot{v}=u-v,\quad
\dot{w}=\phi(u-w).
\]

This can produce:

- fixed points,
- limit-cycle oscillations,
- excitability (large transient response to perturbation).

These correspond to temporal color oscillation observed in a stirred BZ reaction.

## 3) Reaction–diffusion and propagating color fronts

In unstirred medium, local kinetics couples through diffusion. If one region oxidizes first,
nearby regions receive diffusive activator/inhibitor flux and can be triggered after a delay.
This creates traveling excitation fronts.

- Activator-dominated leading edge triggers neighbors,
- inhibitor/recovery behind the front enforces refractory behavior,
- refractory behavior prevents immediate re-excitation, so waves have finite width and speed.

## 4) Oxidation state and visible color

Color indicator depends on catalyst redox state:

- **Ferroin system**:
  - reduced ferroin \(\mathrm{Fe(phen)_3^{2+}}\): red,
  - oxidized ferriin \(\mathrm{Fe(phen)_3^{3+}}\): blue.
- **Cerium system**:
  - \(\mathrm{Ce^{3+}}\): colorless/pale,
  - \(\mathrm{Ce^{4+}}\): yellow.

In simulation, map one variable (often the slow/redox variable, e.g. \(w\)) to a color ramp.

## 5) ODE \(\rightarrow\) PDE extension

Given ODE \(\dot{y}=F(y)\), PDE extension is:

\[
\frac{\partial y}{\partial t}=F(y)+D\nabla^2 y.
\]

For multiple species, use one diffusion coefficient per field. Often inhibitor diffusion is smaller than activator in gels; in liquids they can be closer.

## 6) Why target and spiral waves appear

- **Target waves**: local pacemaker point periodically re-triggers concentric waves.
- **Spiral waves**: form when a wavefront breaks (obstacle/noise/heterogeneity), and the broken tip rotates around refractory medium.

Excitable media + refractory period + curvature-dependent wave speed is enough for robust spiral dynamics.

## 7) Computational implementation (2D grid)

Use uniform 2D grid with fields \(u,v,w\) over \((i,j)\).

1. Compute Laplacian with 5-point stencil:
   \[
   \nabla^2 u_{i,j}\approx\frac{u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1}-4u_{i,j}}{\Delta x^2}.
   \]
2. Compute reaction terms from local \(u,v,w\).
3. Update with explicit Euler or RK4.
4. Apply boundary conditions (Neumann/no-flux is common).

### Numerical options

- **Explicit Euler**: simplest, small \(\Delta t\) needed for stability.
- **Heun / RK4**: better accuracy for reaction terms.
- **Operator splitting**: reaction step + diffusion step.
- **Semi-implicit diffusion**: larger stable time step for stiff diffusion.

## 8) Practical parameter and visualization suggestions

A useful nondimensional starting point:

- \(\epsilon=0.04\)
- \(f=1.4\)
- \(q=0.002\)
- \(\phi=0.08\)
- \(D_u=1.0, D_v=0.6, D_w=0.2\)

Grid/time:

- \(N_x=N_y=200\), \(\Delta x=1.0\)
- \(\Delta t\in[10^{-3},10^{-2}]\) for explicit method
- simulate at least 5k–50k steps for wave development.

Initialization:

- set near resting state everywhere,
- seed a local patch with elevated \(u\) to trigger wave,
- add small random noise (\(10^{-4}\) to \(10^{-3}\)) to break symmetry.

Visualization:

- map \(w\) (or oxidized fraction proxy) to color:
  - ferroin style: red \(\rightarrow\) blue,
  - cerium style: pale \(\rightarrow\) yellow.
- save frames as PPM/PNG every 10–50 steps and stitch into a video.

---

This model is intentionally reduced, but it captures the core nonlinear mathematics: oscillation,
excitability, refractory recovery, and self-organized wave patterns.

## 9) 文献中常见的 BZ 波纹演化特征（用于调参参考）

基于 Oregonator / 可激发介质文献的典型现象：

- **Target waves（靶心波）**：由局部起搏点（pacemaker）周期触发，形成同心圆外扩。
- **Inward waves（内向波）**：边界或外环触发后，前沿可向中心收缩传播。
- **Spiral waves（螺旋波）**：波前断裂（噪声、障碍、参数突变）后，波头旋转形成螺旋。
- **Refractory period（不应期）**：每一圈波后留下恢复区，导致后续波有固定间距。

想得到更“像实验”的视频，可优先使用：

- 较短起搏周期（例如 140–200），
- `view=mix`（更能看到波前运动），
- 自动对比度（提升弱波前可见性），
- 适度噪声与更长演化时间（>10k steps）。

这些设置对应可激发反应扩散系统的常见图样学行为：

- Field, Körös, Noyes（FKN）机理与其简化（Oregonator）
- Tyson/Fife 等对 BZ 可激发波与反应扩散图样的数值研究
- Winfree 关于螺旋波与可激发介质动力学

## 10) 螺旋波的数学描述（和实验图样对应）

在 BZ 可激发介质中，螺旋波通常由“断裂的激发前沿 + 不应期”形成。常见数学框架：

1. **反应扩散方程（Oregonator / Barkley 类）**
   \[
   \partial_t \mathbf{y}=\mathbf{F}(\mathbf{y})+D\nabla^2\mathbf{y}
   \]
   其中 \(\mathbf{y}\) 可取 \((u,v,w)\) 或 \((u,v)\)。

2. **波前几何近似（曲率关系）**
   \[
   c_n \approx c_0 - D\kappa
   \]
   - \(c_n\)：法向传播速度，
   - \(\kappa\)：前沿曲率。

   螺旋尖端附近曲率大、速度低；远处曲率小、速度接近平面波速度，形成稳定旋转波臂。

3. **不应期 + 重新激发条件**
   每次激发后抑制变量尚未恢复时，局部不能再次被激发；因此波臂之间有间距，形成条纹/涡旋结构。

这和你提供的培养皿图像一致：淡色波前代表高激发区域，后方是恢复区，前沿连续外扩并卷曲成螺旋臂。

## 11) 经典论文（数学建模与波纹演化）

下面是 BZ 反应数学描述和螺旋波研究中最常引用的文献方向（可直接作为阅读入口）：

1. Field, R. J.; Körös, E.; Noyes, R. M. (1972). *Oscillations in chemical systems. II. Thorough analysis of temporal oscillation in the bromate-cerium-malonic acid system*. JACS.（FKN 机理基础）
2. Field, R. J.; Noyes, R. M. (1974). *Oscillations in chemical systems. IV. Limit cycle behavior in a model of a real chemical reaction*. J. Chem. Phys.（Oregonator 经典来源）
3. Tyson, J. J.; Fife, P. C. (1980). *Target patterns in a realistic model of the Belousov–Zhabotinskii reaction*. J. Chem. Phys.（靶心波与反应扩散模型）
4. Winfree, A. T. (1972). *Spiral waves of chemical activity*. Science.（螺旋波可视化与动力学早期经典）
5. Tyson, J. J. (1985). *Complex dynamics in the Belousov–Zhabotinskii reaction*. Lecture Notes / reviews.（复杂动力学综述）
6. Barkley, D. (1991). *A model for fast computer simulation of waves in excitable media*. Physica D.（可激发介质简化模型，螺旋波常用）

> 建议按顺序阅读：FKN/Oregonator \(\rightarrow\) 反应扩散靶心波 \(\rightarrow\) 螺旋波几何与可激发介质通用理论。

## 12) 让波“更慢、层数更多”的数学调参思路

如果你希望图样从“点状快速喷发”变成“缓慢、多层并行波前”，在反应扩散模型中可按下面方式调整：

1. **降低有效传播速度**
   - 适当减小扩散系数（尤其激发变量扩散 \(D_u\)），
   - 减小时间步长并增加总步数，让视频时间尺度更“慢”。

2. **增强多层波前结构**
   - 采用周期起搏（pacemaker），但每次起搏不是单点，而是“多层环形扰动”，
   - 对应数学上可看作在源项 \(S(x,y,t)\) 中叠加多个不同半径的环形脉冲。

3. **保持不应期**
   - 不应期由抑制变量恢复速度控制，恢复太快会让波前挤在一起失真；
   - 通过 \(\phi, f, q\) 的组合调节恢复窗口，可得到更稳定层状波纹。

在数值上，这些都可写成：
\[
\partial_t u = F_u(u,v,w) + D_u\nabla^2u + S_u(x,y,t),
\]
其中 \(S_u\) 可以是中心多层环形起搏项，用于构造多圈同步外扩前沿。

## 13) 让图样不“机械重复”的论文一致改法（本程序已采用）

实验中的 BZ 波纹往往不是严格同心重复，而是会出现慢变花瓣/螺旋臂。常见建模做法是把
“介质异质性 + 缓慢调制 + 局部起搏源”加入到 Oregonator 反应扩散框架：

\[
\partial_t u=\frac{1}{\epsilon(x,y)}\Big[u(1-u)-f v\frac{u-q}{u+q}\Big]+D_u\nabla^2u+S(x,y,t),
\]
\[
\partial_t v=u-v+D_v\nabla^2v,\qquad
\partial_t w=\phi(x,y,t)\,(u-w)+D_w\nabla^2w.
\]

其中：

- \(\epsilon(x,y)=\epsilon_0(1+\eta(x,y))\)：空间异质性（实验里对应局部浓度/温度/厚度差异）；
- \(\phi(x,y,t)=\phi_0\,[1+\alpha \eta(x,y)+\beta\sin(\omega t+\theta(x,y))]\)：慢时间调制；
- \(S(x,y,t)\)：移动或环绕的局部起搏 + 偶发波前断裂（更容易得到持续螺旋臂）。
- 若希望“从一个点向外衍射并逐渐复杂化”，可令 \(S\) 只有单中心起搏，再通过弱异质性与间歇断裂触发次级结构。
- 若目标是“大小不一圆形 + 不均匀线条”（接近培养皿图样），可将次级中心写成
  \[
  S(x,y,t)=S_0(x,y,t)+\sum_{k=1}^{N_c} a_k\,G_k(x,y)\,\Pi_k(t),
  \]
  其中 \(G_k\) 是局部高斯环源，\(\Pi_k\) 是不同周期/相位的脉冲。这样会自然形成不同半径圆环和互相切割的波前。

这样能避免“每一圈都完全相同”的单调感，更接近 BZ 视频里“波前演化、分叉、再组织”的观感。

### 建议重点阅读（与上式直接相关）

1. Field, Körös, Noyes (1972), JACS：FKN 机理基础。  
2. Field, Noyes (1974), J. Chem. Phys.：Oregonator 经典模型。  
3. Tyson, Fife (1980), J. Chem. Phys.：BZ 靶心波/反应扩散图样。  
4. Winfree (1972), Science：化学螺旋波经典实验图像。  
5. Zhabotinsky et al. / photosensitive Oregonator 系列：将外场与时空调制显式写入动力学，解释波前慢变与螺旋演化。
