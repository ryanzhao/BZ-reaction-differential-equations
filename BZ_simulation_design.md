# BZ Reaction Simulation Design (Oregonator + Reaction–Diffusion)

This note gives a practical mathematical and computational blueprint for a C++ simulation of the Belousov–Zhabotinsky (BZ) reaction.

## 1) Nonlinear chemical kinetics: Oregonator model

A widely used reduced model is the **Oregonator** (Field–Noyes). In non-dimensional form, a common 3-variable kinetic system is:

\[
\begin{aligned}
\frac{dX}{dt} &= \frac{1}{\epsilon}\left(qY - XY + X(1-X)\right),\\
\frac{dY}{dt} &= -qY - XY + fZ,\\
\frac{dZ}{dt} &= \delta\,(X-Z).
\end{aligned}
\]

Where:

- \(X\): activator-like intermediate (often related to bromous acid dynamics).
- \(Y\): inhibitor-like intermediate (related to bromide influence).
- \(Z\): oxidized catalyst fraction (or proxy variable tied to catalyst redox state).
- \(\epsilon, q, f, \delta\): kinetic parameters controlling timescale separation and nonlinearity.

This system is nonlinear through terms like \(XY\), \(X(1-X)\), and coupling between variables.

## 2) ODE time evolution of intermediates

At a single well-mixed point, the Oregonator ODEs are enough. Typical behavior:

- Fixed points at some parameter values.
- Stable oscillations (limit cycles) in the oscillatory regime.
- Excitable dynamics if perturbed from rest.

Interpretation:

- Activator \(X\) rises quickly when autocatalysis dominates.
- Inhibitor \(Y\) suppresses activation.
- Catalyst state \(Z\) follows more slowly, introducing lag and phase shift.

## 3) Reaction–diffusion for pattern propagation

To model visible waves in space, add diffusion terms:

\[
\begin{aligned}
\frac{\partial X}{\partial t} &= \frac{1}{\epsilon}(qY - XY + X(1-X)) + D_X\nabla^2X,\\
\frac{\partial Y}{\partial t} &= (-qY - XY + fZ) + D_Y\nabla^2Y,\\
\frac{\partial Z}{\partial t} &= \delta(X-Z) + D_Z\nabla^2Z.
\end{aligned}
\]

- \(D_X, D_Y, D_Z\): diffusion coefficients.
- Usually, not all species diffuse equally; catalyst diffusion may be small in some setups.

## 4) Catalyst oxidation states and color

Color changes come from catalyst redox switching:

- **Ferroin system**:
  - Reduced ferroin (Fe(II) complex): reddish.
  - Oxidized ferriin (Fe(III) complex): blue.
- **Cerium system**:
  - Ce(III): colorless (or weakly colored).
  - Ce(IV): yellow.

In simulation, map a variable (often \(Z\), oxidized fraction) to color intensity. For example:

- low \(Z\): red (ferroin reduced),
- high \(Z\): blue (oxidized).

## 5) ODE to PDE extension

Procedure:

1. Start with validated ODE kinetics (single cell).
2. Place one ODE system at each grid point \((i,j)\).
3. Couple neighboring points via discrete Laplacian:

\[
\nabla^2U_{i,j} \approx \frac{U_{i+1,j}+U_{i-1,j}+U_{i,j+1}+U_{i,j-1}-4U_{i,j}}{\Delta x^2}.
\]

4. Choose boundary conditions (Neumann/no-flux is common for closed dishes).

## 6) Why diffusion makes target and spiral waves

In an excitable medium:

- Local reaction provides threshold-triggered activation and refractory recovery.
- Diffusion spreads activation to neighbors.
- A localized trigger produces expanding **target waves**.
- Broken wavefronts can curl and rotate, forming **spirals**.

Spiral stability depends on kinetics, diffusion ratio, and perturbations/heterogeneities.

## 7) Computational architecture (2D grid)

Recommended minimal structure:

- Grid: `Nx x Ny` arrays for `X, Y, Z`.
- Time stepping:
  - Explicit Euler or RK4 for reaction terms.
  - Finite difference for diffusion.
  - For simplicity, explicit Euler for full PDE can work with small `dt`.
- Boundary handling:
  - No-flux via mirrored indices (Neumann).
- Initialization:
  - Near steady state + random noise.
  - Add local perturbation (central disk) to trigger waves.

## 8) Practical parameters, numerics, visualization

### Example dimensionless starting values

Try these as exploratory values (tune as needed):

- `epsilon = 0.05`
- `q = 0.002`
- `f = 1.2`
- `delta = 0.02`
- `Dx = 1.0`
- `Dy = 0.6`
- `Dz = 0.0` (or small)

These are not universal; BZ regimes are sensitive.

### Discretization advice

- Use 5-point Laplacian on a square grid.
- Start with `dx = 1.0`.
- Stability for explicit diffusion in 2D roughly requires:

\[
\Delta t \lesssim \frac{\Delta x^2}{4D_{\max}}.
\]

- Also satisfy reaction timescale constraints (often need smaller `dt`).

### Visualization

- Map `Z` or `X` to grayscale or RGB.
- For ferroin-like display, interpolate between red and blue:
  - `color = (1-Z)*red + Z*blue`.
- Save frames as PPM/PNG and stitch into video (`ffmpeg`).
- Use a colorbar and fixed value range for consistent interpretation.

## Suggested development roadmap

1. Validate ODE oscillations for a single cell.
2. Extend to 2D PDE with diffusion off (`D=0`) and confirm each cell behaves independently.
3. Turn on diffusion and trigger a pulse to observe wave propagation.
4. Create spiral by cutting a wavefront or asymmetric perturbation.
5. Add parameter sweep scripts for `(epsilon, q, f, Dx, Dy)`.

