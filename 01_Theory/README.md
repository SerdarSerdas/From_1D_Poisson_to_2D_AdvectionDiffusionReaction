# Mathematical Framework and Theory

> **Quick navigation:** For practical results, see [main README](../README.md). This document provides the complete theoretical foundation, including weak formulations, error estimates, singular perturbation analysis, and comprehensive literature review.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Mathematical Formulation](#2-mathematical-formulation)
3. [Standard Galerkin FEM](#3-standard-galerkin-fem-sfem)
4. [Least-Squares FEM](#4-least-squares-fem-lsfem)
5. [Convergence Theory](#5-convergence-theory)
6. [Singular Perturbation Analysis](#6-singular-perturbation-analysis)
7. [Adaptive Mesh Refinement](#7-adaptive-mesh-refinement)
8. [References](#8-references)

---

## 1. Introduction

The finite element method (FEM) represents a fundamental numerical technique for the approximate solution of partial differential equations (PDEs) governing diverse physical phenomena, including fluid dynamics, heat transfer, transport in porous media, and fluid-structure interaction. The classical Galerkin finite element method discretizes the weak (variational) formulation of PDEs through projection onto finite-dimensional approximation spaces, typically constructed from piecewise polynomial basis functions. 

While widely adopted and theoretically well-established, the Galerkin approach exhibits numerical instabilities in certain problem classes, particularly:

- **Convection-dominated transport phenomena** — Spurious oscillations arise when advection dominates diffusion (high Péclet numbers)
- **Saddle-point systems** — The incompressible Navier-Stokes equations require careful treatment of velocity-pressure coupling
- **Mixed formulations** — Must satisfy the Ladyzhenskaya-Babuška-Brezzi (LBB) condition through compatible finite element pairs

These challenges have motivated extensive research into stabilized and alternative formulation strategies.

### 1.1 Advantages of Least-Squares FEM

The Least-Squares Finite Element Method (LSFEM) offers several compelling advantages over the standard Galerkin approach:

1. **Equal-order interpolation spaces** — May be employed for all unknown fields, eliminating the need for inf-sup stable element pairs
2. **Symmetric positive definite matrices** — The discrete system matrices are always sparse, symmetric, and positive definite, regardless of the underlying PDE structure
3. **Robust iterative solvers** — The resulting algebraic systems admit efficient solution via conjugate gradient methods and preconditioners
4. **Natural error estimators** — A residual-based error estimator emerges directly from the least-squares functional, enabling straightforward adaptive mesh refinement (AMR) strategies
5. **Flux superconvergence** — For certain problem classes, the flux variable achieves higher-order convergence rates

### 1.2 LSFEM Methodology

The standard L² LSFEM framework follows a systematic procedure:

1. **Reformulate** the governing equations as a first-order system
2. **Define** least-squares functionals in the L²-norm
3. **Discretize** using C⁰-continuous finite element spaces

This approach avoids the need for stabilization terms and mixed formulations while maintaining optimal convergence properties.

### 1.3 Literature Context

Building upon these computational advantages, LSFEM has been successfully applied to a broad spectrum of partial differential equations [1,2,3,4]. A particularly active research area concerns general transport equations arising in heat transfer and fluid mechanics. Various first-order formulations and stabilization techniques have been proposed and rigorously analyzed [2,5,6,7,8,9,10,11,12].

While substantial theoretical developments and error estimates exist for this problem class, few studies systematically compare LSFEM variants against standard Galerkin methods across the full hierarchy from Poisson to advection-diffusion-reaction problems, particularly with adaptive refinement and scaling strategies. **This repository addresses this gap through comprehensive 1D/2D benchmarks.**

### 1.4 Problem Hierarchy

This investigation adopts a systematic, hierarchical progression toward the general advection-diffusion-reaction equation:

$$
-\nabla \cdot (\nu \nabla u + \mathbf{a} u) + c u = f
$$

The progression is:

1. **Poisson equation** ($\mathbf{a} = \mathbf{0}$, $c=0$) — Pure diffusion, well-conditioned
2. **Advection-diffusion** ($c=0$) — Introduces directional transport, potential instabilities
3. **Full ADR** — Complete coupling of all transport mechanisms

This hierarchical approach isolates the effects of advection and reaction, enabling precise characterization of each method's strengths and limitations.

### 1.5 Scope of Investigation

The study conducts a comprehensive comparison between the standard Galerkin finite element method (SFEM) and adaptive least-squares finite element method (LSFEM) applied to this sequence of steady-state linear transport problems. Within the LSFEM framework, residual-based error indicators guide adaptive mesh refinement strategies, thereby concentrating computational resolution in regions characterized by sharp gradients or elevated approximation errors.

**Three dimensions of analysis:**

1. **Discretization effects** — P₁ vs P₂ polynomial approximations on regular and perturbed meshes
2. **Solver performance** — Unpreconditioned CG, Jacobi-preconditioned CG, and algebraic multigrid (AMG)-preconditioned CG
3. **Mesh adaptivity** — Uniform refinement vs residual-based adaptive refinement (Dörfler marking)

**Solver configuration:**

The AMG preconditioner employs PyAMG smoothed aggregation [13,14] with:
- Symmetric Gauss-Seidel smoothers (4 pre/post-smoothing steps)
- W-cycle iteration
- Target tolerance: 10⁻¹⁰ relative residual

This enables quantitative assessment of both discretization and preconditioning efficacy across problem size, mesh adaptivity, element order, and extreme scaling regimes (singular perturbations with ε → 0).

---

## 2. Mathematical Formulation

### 2.1 Strong Form

The general steady advection-diffusion-reaction (ADR) equation in domain $\mathcal{B} \subset \mathbb{R}^d$ ($d=1,2,3$) is:

$$
\begin{aligned}
-\nabla \cdot (\nu \nabla u + \mathbf{a} u) + c u &= f &&\text{in } \mathcal{B}, \\
u &= g_D &&\text{on } \partial\mathcal{B}_D, \\
\mathbf{n} \cdot (\nu \nabla u + \mathbf{a} u) &= g_N &&\text{on } \partial\mathcal{B}_N,
\end{aligned}
$$

where:
- $u: \mathcal{B} \to \mathbb{R}$ is the scalar unknown (temperature, concentration, etc.)
- $\nu > 0$ is the diffusion coefficient (constant or piecewise constant)
- $\mathbf{a}: \mathcal{B} \to \mathbb{R}^d$ is the advection velocity field
- $c \geq 0$ is the reaction coefficient
- $f: \mathcal{B} \to \mathbb{R}$ is the source term
- $\mathbf{n}$ is the outward unit normal on the boundary
- $\partial\mathcal{B} = \partial\mathcal{B}_D \cup \partial\mathcal{B}_N$ (disjoint)

### 2.2 Special Cases

**Poisson equation** ($\nu = 1$, $\mathbf{a} = \mathbf{0}$, $c = 0$):

$$
-\Delta u = f \quad \text{in } \mathcal{B}
$$

**Diffusion-Reaction** ($\nu = 1$, $\mathbf{a} = \mathbf{0}$):

$$
-\Delta u + c u = f
$$

**Advection-Diffusion** ($c = 0$):

$$
-\nabla \cdot (\nu \nabla u + \mathbf{a} u) = f
$$


### 2.3 Dimensionless Parameters

**Péclet number** (ratio of advection to diffusion):

$$
\text{Pe} = \frac{|\mathbf{a}| h}{2\nu}
$$

where $h$ is the characteristic mesh size. Large Pe (Pe ≫ 1) indicates convection-dominated flow requiring stabilization or refinement.

**Damköhler number** (ratio of reaction to advection):

$$
\text{Da} = \frac{c L}{|\mathbf{a}|}
$$

where $L$ is a characteristic length scale.

---

## 3. Standard Galerkin FEM (SFEM)

### 3.1 Weak Formulation

Multiply the strong form by a test function $v \in H^1_0(\mathcal{B})$ and integrate by parts:

$$
\int_\mathcal{B} (\nu \nabla u + \mathbf{a} u) \cdot \nabla v \, dx + \int_\mathcal{B} c u v \, dx = \int_\mathcal{B} f v \, dx + \int_{\partial\mathcal{B}_N} g_N v \, ds
$$

**Variational problem:** Find $u \in H^1_{g_D}(\mathcal{B})$ such that

$$
a(u, v) = \ell(v) \quad \forall v \in H^1_0(\mathcal{B})
$$

where:
- $a(u, v) = \int_\mathcal{B} [(\nu \nabla u + \mathbf{a} u) \cdot \nabla v + c u v]  dx$
- $\ell(v) = \int_\mathcal{B} f v \, dx + \int_{\partial\mathcal{B}_N} g_N v  ds$
- $H^1_{g_D} = \{w \in H^1 : w|_{\partial\mathcal{B}_D} = g_D\}$

### 3.2 Discrete Formulation

Introduce a finite-dimensional subspace $V_h \subset H^1_0(\mathcal{B})$ spanned by basis functions $\{\phi_i\}_{i=1}^N$:

$$
u_h(x) = \sum_{i=1}^N u_i \phi_i(x)
$$

**Galerkin orthogonality:** Find $u_h \in V_h$ such that

$$
a(u_h, v_h) = \ell(v_h) \quad \forall v_h \in V_h
$$

Substituting the discrete representation yields the linear system:

$$
\mathbf{K} \mathbf{u} = \mathbf{f}
$$

where:
- $K_{ij} = a(\phi_j, \phi_i)$
- $f_i = \ell(\phi_i)$
- $\mathbf{u} = [u_1, \ldots, u_N]^T$

### 3.3 Matrix Structure

**For the Poisson equation** ($\mathbf{a} = \mathbf{0}$, $c = 0$):

$$
K_{ij} = \int_\mathcal{B} \nu \nabla \phi_j \cdot \nabla \phi_i \, dx
$$

This produces a symmetric positive definite (SPD) stiffness matrix.

**For advection-diffusion** ($\mathbf{a} \neq \mathbf{0}$):

$$
K_{ij} = \int_\mathcal{B} [(\nu \nabla \phi_j + \mathbf{a} \phi_j) \cdot \nabla \phi_i] \, dx
$$

The advection term introduces asymmetry. The matrix may lose positive definiteness for high Péclet numbers.

### 3.4 Flux Recovery

In SFEM, the flux (gradient) is obtained through post-processing:

$$
\mathbf{q}_h = -\nu \nabla u_h - \mathbf{a} u_h
$$

This derivative is typically discontinuous across element boundaries and achieves only O(h^p) accuracy for P_p elements.

---

## 4. Least-Squares FEM (LSFEM)

### 4.1 First-Order System Reformulation

Introduce the flux as an independent variable:

$$
\mathbf{q} = -\nu \nabla u - \mathbf{a} u
$$

The second-order PDE becomes a first-order system:

$$
\begin{aligned}
\mathbf{q} + \nu \nabla u + \mathbf{a} u &= \mathbf{0}, \\
-\nabla \cdot \mathbf{q} + c u &= f.
\end{aligned}
$$

### 4.2 Least-Squares Functional

Define the L² least-squares functional:

$$
\mathcal{J}(u, \mathbf{q}) = \frac{1}{2} \|\mathbf{q} + \nu \nabla u + \mathbf{a} u\|_0^2 + \frac{1}{2} \|-\nabla \cdot \mathbf{q} + c u - f\|_0^2
$$

where $\|\cdot\|_0$ denotes the L² norm:

$$
\|w\|_0^2 = \int_\mathcal{B} w^2 \, dx
$$

**Minimization principle:** Find $(u, \mathbf{q}) \in H^1_0 \times H(\text{div})$ such that

$$
\mathcal{J}(u, \mathbf{q}) = \min_{(v, \mathbf{r}) \in H^1_0 \times H(\text{div})} \mathcal{J}(v, \mathbf{r})
$$

### 4.3 Discrete Formulation

Introduce finite element spaces:
- $V_h \subset H^1_0$ for the scalar field $u$
- $\mathbf{Q}_h \subset H(\text{div})$ for the flux field $\mathbf{q}$

For equal-order approximation (e.g., P₁-P₁ or P₂-P₂):

$$
u_h = \sum_{i=1}^N u_i \phi_i, \quad \mathbf{q}_h = \sum_{i=1}^N q_i \boldsymbol{\psi}_i
$$

The discrete minimization problem is:

$$
\min_{u_h \in V_h, \mathbf{q}_h \in \mathbf{Q}_h} \mathcal{J}(u_h, \mathbf{q}_h)
$$

Taking variations with respect to $u_h$ and $\mathbf{q}_h$ yields the Euler-Lagrange equations:

$$
\begin{aligned}
(\mathbf{q}_h + \nu \nabla u_h + \mathbf{a} u_h, \nu \nabla v_h + \mathbf{a} v_h) + (-\nabla \cdot \mathbf{q}_h + c u_h - f, c v_h) &= 0, \\
(\mathbf{q}_h + \nu \nabla u_h + \mathbf{a} u_h, \mathbf{r}_h) + (-\nabla \cdot \mathbf{q}_h + c u_h - f, -\nabla \cdot \mathbf{r}_h) &= 0.
\end{aligned}
$$

for all $v_h \in V_h$ and $\mathbf{r}_h \in \mathbf{Q}_h$.

### 4.4 Matrix Form (1D Poisson Example)

For the 1D Poisson equation ($\nu = 1$, $\mathbf{a} = 0$, $c = 0$), the first-order system is:

$$
q = u', \quad -q' = f
$$

The least-squares functional simplifies to:

$$
\mathcal{J} = \frac{1}{2} \|q - u'\|_0^2 + \frac{1}{2} \|q' + f\|_0^2
$$

The discrete system becomes:

**Matrix form:** [M, -K^T; -K, M] [q; u] = [Mf; 0]

where:
- $M_{ij} = (\phi_i, \phi_j)$ is the mass matrix
- $K_{ij} = (\phi_i', \phi_j')$ is the stiffness matrix

**Key properties:**
- The system matrix is **symmetric positive definite** (always!)
- Block structure enables specialized solvers
- Both $u$ and $q$ are primary variables (not post-processed)

### 4.5 Advantages Over SFEM

1. **No inf-sup condition** — Equal-order spaces $V_h = \mathbf{Q}_h$ are permissible
2. **Always SPD** — Regardless of problem parameters (Pe, Da)
3. **Flux superconvergence** — $\mathbf{q}_h$ achieves O(h^{p+1}) accuracy
4. **Natural error estimator** — $\mathcal{J}(u_h, \mathbf{q}_h)$ provides element-wise error indicators

### 4.6 Disadvantages

1. **Higher dimension** — 2× DOFs compared to SFEM in 1D, up to (d+1)× in d dimensions
2. **Larger system** — More memory and potentially more iterations (though AMG mitigates this)
3. **Implementation complexity** — Requires careful treatment of flux boundary conditions

---

## 5. Convergence Theory

### 5.1 Standard FEM Error Estimates

For smooth solutions and quasi-uniform meshes, the standard Galerkin FEM satisfies:

**Solution convergence:**

$$
\|u - u_h\|_0 \leq C h^{p+1} |u|_{p+1}, \quad \|u - u_h\|_1 \leq C h^p |u|_{p+1}
$$

where:
- $p$ is the polynomial degree (P_p elements)
- $|\cdot|_k$ denotes the H^k seminorm
- $C$ is a constant independent of $h$

**Flux convergence:**

$$
\|\nabla u - \nabla u_h\|_0 \leq C h^p |u|_{p+1}
$$

Note: The flux achieves only O(h^p) convergence—one order lower than the solution.

### 5.2 LSFEM Error Estimates

For the least-squares formulation with equal-order approximation:

**Solution convergence:**

$$
\|u - u_h\|_0 \leq C h^{p+1} (|u|_{p+1} + |\mathbf{q}|_{p+1})
$$

**Flux superconvergence:**

$$
\|\mathbf{q} - \mathbf{q}_h\|_0 \leq C h^{p+1} (|u|_{p+1} + |\mathbf{q}|_{p+1})
$$

**Key observation:** Both $u$ and $\mathbf{q}$ achieve O(h^{p+1}) convergence! This is superconvergence for the flux compared to SFEM.

### 5.3 Conditions for Optimal Rates

Optimal convergence rates require:

1. **Regularity:** $u \in H^{p+1}(\mathcal{B})$
2. **Mesh quality:** Quasi-uniform triangulation with bounded aspect ratios
3. **Asymptotic regime:** $h/\lambda \ll 1$ where $\lambda$ is the characteristic solution length scale

For singularly perturbed problems (ε ≪ 1), the third condition is critical and often violated on uniform meshes.

---

## 6. Singular Perturbation Analysis

### 6.1 Boundary Layer Problems

Consider the 1D advection-diffusion equation:

$$
-\epsilon u'' - u' = 0, \quad u(0) = 0, \quad u(1) = 1
$$

where $\epsilon \ll 1$ is the diffusion coefficient.

**Analytical solution:**

$$
u(x) = \frac{e^{x/\epsilon} - 1}{e^{1/\epsilon} - 1}
$$

**Boundary layer structure:**
- **Outer region** ($x \ll 1$): $u \approx 0$
- **Boundary layer** (near $x=1$): Width $\lambda \sim \epsilon$, exponential transition

**Resolution requirement:**

$$
h < \lambda \sim \epsilon \quad \text{(uniform mesh)}
$$

For $\epsilon = 10^{-7}$, this requires $h < 10^{-7}$ → over 10 million elements!

### 6.2 Interior Layer Problems

For reaction-diffusion with localized source:

$$
-\epsilon u'' + u = f(x), \quad f(x) = \delta(x - x_0)
$$

The solution exhibits an interior layer near $x_0$ with width:

$$
\lambda \sim \sqrt{\epsilon}
$$

**Key difference:** Interior layers are wider than boundary layers ($\sqrt{\epsilon}$ vs $\epsilon$), but still require $h \ll \sqrt{\epsilon}$ for resolution.

### 6.3 Pre-Asymptotic vs Asymptotic Regimes

**Pre-asymptotic regime:** $h/\lambda \not\ll 1$
- Convergence rates are unreliable
- May observe negative or chaotic apparent rates
- Pollution effects dominate

**Asymptotic regime:** $h/\lambda \ll 1$
- Theoretical convergence rates apply
- Predictable error behavior

**Transition point:** Typically when $h/\lambda < 0.1$ (rule of thumb)

### 6.4 Uniform vs Adaptive Refinement

**Uniform refinement:**
- Requires $O(\epsilon^{-1/d})$ elements in $d$ dimensions
- Exponential cost as $\epsilon \to 0$
- Impractical for $\epsilon < 10^{-5}$

**Adaptive refinement:**
- Concentrates elements in layers
- Achieves exponential convergence in adapted regime
- Requires $O(|\log \epsilon|^d)$ elements (much better!)

---

## 7. Adaptive Mesh Refinement

### 7.1 Error Estimation

For LSFEM, the least-squares functional provides a natural error estimator:

**Global estimator:**

$$
\eta^2 = \mathcal{J}(u_h, \mathbf{q}_h) = \sum_{K \in \mathcal{T}_h} \eta_K^2
$$

**Element-wise indicator:**

$$
\eta_K^2 = \|\mathbf{q}_h + \nu \nabla u_h + \mathbf{a} u_h\|_{0,K}^2 + \|-\nabla \cdot \mathbf{q}_h + c u_h - f\|_{0,K}^2
$$

### 7.2 Dörfler Marking Strategy

**Bulk criterion [15]:** Given parameter $0 < \theta \leq 1$ (typically $\theta = 0.5$):

1. Sort elements by decreasing $\eta_K$
2. Mark minimal set $\mathcal{M} \subset \mathcal{T}_h$ such that:
3. $$
\sum_{K \in \mathcal{M}} \eta_K^2 \geq \theta \sum_{K \in \mathcal{T}_h} \eta_K^2
$$
4. Refine all marked elements

**Properties:**
- Guarantees geometric convergence $\eta^{(n)} \leq C \rho^n$ for some $\rho < 1$
- Balances refinement and computational cost
- Concentrates elements where error is largest

### 7.3 Refinement Strategies

**Bisection (1D):**
- Split marked elements into two children
- Maintain conformity automatically

**Red-green refinement (2D triangles):**
- Red: Split triangle into 4 congruent children
- Green: Split to maintain conformity (temporary)

**Quadrilateral refinement:**
- Isotropic: Split into 4 congruent children
- Anisotropic: Split in one direction (for directional layers)

### 7.4 AMR Algorithm

```
1. Initialize: Coarse mesh T_h^(0), tolerance tol
2. For n = 0, 1, 2, ...
   a. Solve: Compute (u_h, q_h) on T_h^(n)
   b. Estimate: Compute element indicators {η_K}
   c. Check: If ||η|| < tol, STOP
   d. Mark: Dörfler strategy with θ = 0.5
   e. Refine: Create T_h^(n+1) by bisection/red-green
   f. Transfer: Interpolate solution to new mesh
3. Return: (u_h, q_h) on final mesh
```

### 7.5 Advantages of AMR for Singular Perturbations

**Efficiency gains:**
- Uniform: $N = O(\epsilon^{-1})$ elements (1D)
- Adaptive: $N = O(|\log \epsilon|)$ elements (1D)

**Example:** For $\epsilon = 10^{-6}$:
- Uniform: ~1,000,000 elements
- Adaptive: ~50 elements (20,000× reduction!)

**Preservation of superconvergence:**
- LSFEM flux maintains O(h²) convergence under AMR
- Verified numerically in this repository

---

## 8. References

### Core LSFEM Theory

[1] Bochev, P. B., & Gunzburger, M. D. (1998). *Finite element methods of least-squares type*. SIAM Review, 40(4), 789-837. DOI: [10.1137/S0036144597321156](https://doi.org/10.1137/S0036144597321156)

[2] Bochev, P. B., & Gunzburger, M. D. (2009). *Least-Squares Finite Element Methods*. Springer. DOI: [10.1007/b13382](https://doi.org/10.1007/b13382)

[3] Jiang, B. N. (1998). *The Least-Squares Finite Element Method: Theory and Applications in Computational Fluid Dynamics and Electromagnetics*. Springer. DOI: [10.1007/978-3-662-03740-9](https://doi.org/10.1007/978-3-662-03740-9)

[4] Bochev, P., & Gunzburger, M. (2006). Least-squares methods for the velocity-pressure-stress formulation of the Stokes equations. *Computer Methods in Applied Mechanics and Engineering*, 195(7-8), 828-845. DOI: [10.1016/j.cma.2005.02.018](https://doi.org/10.1016/j.cma.2005.02.018)

### Transport Equations

[5] Pehlivanov, A. I., & Carey, G. F. (1994). Error estimates for least-squares mixed finite elements. *RAIRO-Modélisation mathématique et analyse numérique*, 28(5), 499-516. DOI: [10.1051/m2an/1994280504991](https://doi.org/10.1051/m2an/1994280504991)

[6] Pehlivanov, A. I., Carey, G. F., & Lazarov, R. D. (1994). Least-squares mixed finite elements for second-order elliptic problems. *SIAM Journal on Numerical Analysis*, 31(5), 1368-1377. DOI: [10.1137/0731071](https://doi.org/10.1137/0731071)

[7] Yang, S. Y. (2000). Analysis of optimal least-squares finite element methods for the Stokes equations on rectangular elements. *Journal of Computational Mathematics*, 18(3), 277-288.

[8] Carey, G. F., Pehlivanov, A. I., & Vassilevski, P. S. (1995). Least-squares mixed finite element methods for non-selfadjoint elliptic problems: I. Error estimates. *Numerische Mathematik*, 72(4), 501-522. DOI: [10.1007/s002110050180](https://doi.org/10.1007/s002110050180)

[9] Fialko, S. Y., Manko, V. P., & McConnell, P. W. (1998). First-order system least squares for the convection-diffusion problem. *International Journal for Numerical Methods in Engineering*, 43(2), 307-321.

[10] Hsieh, P. W., & Yang, S. Y. (2009). On efficient least-squares finite element methods for convection-dominated problems. *Computer Methods in Applied Mechanics and Engineering*, 199(1-4), 183-196. DOI: [10.1016/j.cma.2009.09.029](https://doi.org/10.1016/j.cma.2009.09.029)

[11] Lazarov, R., Tobiska, L., & Vassilevski, P. S. (1997). Streamline diffusion least-squares mixed finite element methods for convection-diffusion problems. *East-West Journal of Numerical Mathematics*, 5(4), 249-264.

[12] Cai, Z., Manteuffel, T. A., & McCormick, S. F. (1997). First-order system least squares for second-order partial differential equations: Part II. *SIAM Journal on Numerical Analysis*, 34(2), 425-454. DOI: [10.1137/S0036142994266066](https://doi.org/10.1137/S0036142994266066)

### Algebraic Multigrid

[13] Bell, W. N., Olson, L. N., & Schroder, J. B. (2022). PyAMG: Algebraic Multigrid Solvers in Python v5.0. Release 5.0. DOI: [10.5281/zenodo.596446](https://doi.org/10.5281/zenodo.596446)

[14] Virtanen, P., Gommers, R., Oliphant, T. E., et al. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261-272. DOI: [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)

### Adaptive Mesh Refinement

[15] Dörfler, W. (1996). A convergent adaptive algorithm for Poisson's equation. *SIAM Journal on Numerical Analysis*, 33(3), 1106-1124. DOI: [10.1137/0733054](https://doi.org/10.1137/0733054)

### Singular Perturbations

[16] Roos, H.-G., Stynes, M., & Tobiska, L. (2008). *Robust Numerical Methods for Singularly Perturbed Differential Equations* (2nd ed.). Springer Series in Computational Mathematics, Vol. 24. DOI: [10.1007/978-3-540-34467-4](https://doi.org/10.1007/978-3-540-34467-4)

[17] Shishkin, G. I. (1988). Grid approximation of singularly perturbed elliptic and parabolic equations. Second doctoral thesis, Keldysh Institute, Moscow. (In Russian)

[18] Linß, T. (2010). *Layer-Adapted Meshes for Reaction-Convection-Diffusion Problems*. Springer Lecture Notes in Mathematics, Vol. 1985. DOI: [10.1007/978-3-642-05134-0](https://doi.org/10.1007/978-3-642-05134-0)

### Convergence Theory

[19] Brenner, S. C., & Scott, L. R. (2008). *The Mathematical Theory of Finite Element Methods* (3rd ed.). Springer Texts in Applied Mathematics, Vol. 15. DOI: [10.1007/978-0-387-75934-0](https://doi.org/10.1007/978-0-387-75934-0)

[20] Ciarlet, P. G. (2002). *The Finite Element Method for Elliptic Problems*. SIAM Classics in Applied Mathematics, Vol. 40. DOI: [10.1137/1.9780898719208](https://doi.org/10.1137/1.9780898719208)

---

## Additional Resources

### Software Implementations

- **FEniCS Project**: [https://fenicsproject.org/](https://fenicsproject.org/)
- **deal.II**: [https://www.dealii.org/](https://www.dealii.org/)
- **PyAMG**: [https://github.com/pyamg/pyamg](https://github.com/pyamg/pyamg)

### Related Work

- **This repository**: [From 1D Poisson to 2D ADR](https://github.com/SerdarSerdas/From_1D_Poisson_to_2D_AdvectionDiffusionReaction)
- **Numerical results**: See [main README](../README.md) for comprehensive benchmarks
- **1D Poisson**: Detailed convergence studies in [1D_Poisson/](../1D_Poisson/)

---

**Document version:** 1.0  
**Last updated:** January 2025  
**Maintainer:** Serdar Serdas
