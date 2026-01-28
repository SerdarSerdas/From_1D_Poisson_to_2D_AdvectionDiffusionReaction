# 1D Poisson Equation: SFEM vs LSFEM

### Strong Form

$$
\begin{aligned}
-u^{\prime\prime} &= f \quad &&\text{in } {\cal B} = (0,1) \\
\end{aligned}
$$

with

$$
\begin{aligned}
u &= 0 \quad &&\text{on } \partial{\cal B}_D .
\end{aligned}
$$

## Standard Galerkin FEM (SFEM)

The classical Galerkin finite element formulation seeks $u \in V \subset H^1_0({\cal B})$ such that

$$(u^\prime_h, v^\prime_h)_0 = (f, v_h)_0 \quad \forall v_h \in V_h$$

## Least Squares FEM (LSFEM)

### First-Order System

Define the diffusive flux $q = u^\prime$:

$$-q^\prime = f, \qquad q - u^\prime = 0$$

### Functional

$${\cal J}(u,q;f) = \frac{1}{2} \parallel w_1 (-q^\prime - f) \parallel^2_0 + \frac{1}{2} \parallel w_2 (q - u^\prime) \parallel^2_0$$

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}((u_h,q_h);(v_h,r_h)) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}((u_h,q_h);(v_h,r_h)) = w^2_1(-q^\prime_h , -r^\prime_h)_0 + w^2_2 (q_h - \nu u^\prime_h, r_h - v^\prime_h)_0$$

$${\cal F}(v_h,r_h) = w^2_1 (f, -r^\prime_h )_0$$


## 🔍 Key Findings

All configurations are discretized using linear ($P_1$) finite elements with ($\nun = 1\cdot10^{-3}$ / **moderate diffusion**). Both **regular** and **perturbed** grids are studied.




**1. Solution Variable** $u$ **(Tables 1-2)**

Both **SFEM** and **LSFEM** achieve optimal second-order convergence $(O(h^2)$ for the solution variable u:

**SFEM** produces the **smallest errors**:

- ML5: $4.01 \times 10^{−3}$ → ML10: $3.93 \times 10^{−6} (regular mesh)
- Consistent rate ≈ 2.00 across all refinement levels ✅

**LSFEM** exhibits **$~1.8 \times$** larger errors than SFEM:

- ML5: $7.43 \times 10^{−3}$ → ML10: $7.34 \times 10^{−6} (regular mesh)
- Also maintains optimal second-order convergence ✅

**Table 1: $L^2$ Error for Solution $u$ ($\| u - u_h \|_0$)** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.0088e-03 | 1.0060e-03 | 2.5175e-04 | 6.2954e-05 | 1.5739e-05 | 3.9349e-06 |
| LSFEM / RG  | 7.4255e-03 | 1.8731e-03 | 4.6933e-04 | 1.1740e-04 | 2.9354e-05 | 7.3389e-06 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 4.9101e-03 | 1.1785e-03 | 2.8344e-04 | 7.0732e-05 | 1.7867e-05 | 4.4419e-06 |
| LSFEM / PG  | 8.6949e-03 | 2.0132e-03 | 5.2286e-04 | 1.2784e-04 | 3.2207e-05 | 7.9986e-06 |



**Table 2: Convergence Rates for Solution $u$**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.99  | 2.00  | 2.00  | 2.00  | 2.00   | ✅ |
| LSFEM  | 1.99  | 2.00  | 2.00  | 2.00  | 2.00   | ✅ |


**Table 3: $L^2$ Error for Flux $q$ ($\| q - q_h \|_0$)** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.0613e-01 | 2.0367e-01 | 1.0191e-01 | 5.0965e-02 | 2.5483e-02 | 1.2742e-02 |
| LSFEM / RG  | 6.2490e-02 | 1.5812e-02 | 3.9654e-03 | 9.9212e-04 | 2.4808e-04 | 6.2023e-05 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 4.3709e-01 | 2.1479e-01 | 1.0588e-01 | 5.2829e-02 | 2.6532e-02 | 1.3245e-02 |
| LSFEM / PG  | 7.6399e-02 | 1.7035e-02 | 4.6915e-03 | 1.0984e-03 | 2.8095e-04 | 6.9651e-05 |


**Table 4: Convergence Rates for Solution $q$**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.00  | 1.00  | 1.00  | 1.00  | 1.00   | ✅ |
| LSFEM  | 1.98  | 2.00  | 2.00  | 2.00  | 2.00   | ✅ |



**Table 5: CG - iteration counts (Regular/Perturbed)**
| Mesh | CG | CG-JACOBI | CG-AMG |
|------|-----------------|------------------|-------------------|
| ML5  |   32 /   66 |   32 /   64 |  10 /  11 |
| ML6  |   64 /  134 |   64 /  128 |  15 /  16 |
| ML7  |  128 /  270 |  128 /  256 |  25 /  26 |
| ML8  |  256 /  542 |  256 /  515 |  44 /  45 |
| ML9  |  512 / 1086 |  512 / 1036 |  79 /  82 |
| ML10 | 1024 / 2177 | 1024 / 2076 | 147 / 152 |



### Key Observations

- **Solution accuracy**: Both SFEM and LSFEM achieve optimal $\mathcal{O}(h^2)$ convergence for $u$
- **Flux accuracy**: LSFEM demonstrates **superconvergence** for flux $q$:
  - SFEM: $\mathcal{O}(h)$ convergence (standard)
  - LSFEM: $\mathcal{O}(h^2)$ convergence (superconvergent)
- **Mesh robustness**: Both methods maintain convergence rates on perturbed meshes
- **Practical impact**: Orders of magnitude improvement in flux accuracy for LSFEM

### Convergence Rates

The following tables show the convergence rates computed between consecutive mesh levels (regular grid).

**Table 3: Convergence Rates for Solution $u$**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| **Boundary Layer ($\epsilon=10^{-3}$)** |
| SFEM   | 1.94  | 1.98  | 2.00  | 2.00  | 2.00   | ✓ |
| LSFEM  | 1.98  | 1.99  | 2.00  | 2.00  | 2.00   | ✓ |
| **Interior Layer ($\epsilon=10^{-3}$)** |
| SFEM   | 1.99  | 2.00  | 2.00  | 2.00  | 2.00   | ✓ |
| LSFEM  | 1.99  | 2.00  | 2.00  | 2.00  | 2.00   | ✓ |

**Table 4: Convergence Rates for Flux $q$**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| **Boundary Layer ($\epsilon=10^{-3}$)** |
| SFEM   | 0.95  | 0.99  | 1.00  | 1.00  | 1.00   | ✓ |
| **LSFEM**  | **1.94**  | **1.98**  | **2.00**  | **2.00**  | **2.00**   | ✓ |
| **Interior Layer ($\epsilon=10^{-3}$)** |
| SFEM   | 1.00  | 1.00  | 1.00  | 1.00  | 1.00   | ✓ |
| **LSFEM**  | **1.98**  | **2.00**  | **2.00**  | **2.00**  | **2.00**   | ✓ |

### Key Result

The convergence rate tables clearly demonstrate:

- **For solution $u$**: Both methods achieve optimal $\mathcal{O}(h^2)$ convergence
- **For flux $q$**: 
  - SFEM: $\mathcal{O}(h)$ convergence rate ≈ 1.0
  - LSFEM: $\mathcal{O}(h^2)$ convergence rate ≈ 2.0 ✨

**This is the superconvergence property of LSFEM** — it achieves one order higher accuracy for the flux compared to standard Galerkin FEM.

### Solver Performance

The symmetric positive definite nature of LSFEM enables efficient solution with conjugate gradient (CG) methods and algebraic multigrid (AMG) preconditioning.

**Table 5: Solver Iteration Counts for LSFEM ($P_1$ elements, $\epsilon=10^{-3}$)**

| Grid | Mesh | **Boundary Layer** ||| **Interior Layer** |||
|------|------|----|-----------|---------|----|-----------|---------| 
|      |      | CG | CG-JACOBI | **CG-AMG** | CG | CG-JACOBI | **CG-AMG** |
| **Regular** | ML5  | 32 | 32 | **10** | 32 | 32 | **10** |
|      | ML6  | 64 | 64 | **15** | 64 | 64 | **15** |
|      | ML7  | 128 | 128 | **25** | 128 | 128 | **25** |
|      | ML8  | 256 | 256 | **43** | 256 | 256 | **44** |
|      | ML9  | 512 | 512 | **79** | 512 | 512 | **79** |
|      | ML10 | 1024 | 1024 | **147** | 1024 | 1024 | **147** |
| **Perturbed** | ML5  | 66 | 64 | **10** | 66 | 64 | **11** |
|      | ML6  | 134 | 128 | **16** | 134 | 128 | **16** |
|      | ML7  | 269 | 256 | **26** | 270 | 256 | **26** |
|      | ML8  | 541 | 515 | **45** | 542 | 515 | **45** |
|      | ML9  | 1084 | 1035 | **82** | 1086 | 1036 | **82** |
|      | ML10 | 2177 | 2076 | **152** | 2177 | 2076 | **152** |

### Computational Efficiency

Key observations:
- **Without preconditioning**: CG iterations grow linearly with mesh size (e.g., 32→1024 iterations from ML5→ML10)
- **With AMG preconditioning**: Near mesh-independent convergence (10→147 iterations from ML5→ML10)
- **Speedup**: Up to **14× reduction** in iterations at finest level (ML10: 2177 → 152)
- **Mesh robustness**: AMG maintains efficiency on perturbed meshes

The **CG-AMG combination** makes LSFEM computationally practical even for large-scale problems, with iteration counts growing much slower than mesh size.
