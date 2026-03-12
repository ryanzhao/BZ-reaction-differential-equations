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
