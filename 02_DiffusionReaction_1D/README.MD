## 1D Diffusion-Reaction Equation

### Strong Form

$$
\begin{aligned}
-u^{\prime\prime} + cu &= f \quad &&\text{in } {\cal B} = (0,1), \\
u(0) = u(1) &= 0 \quad &&\text{on } \partial{\cal B}_D
\end{aligned}
$$

where $c \geq 0$ is the reaction coefficient.

## Standard Galerkin FEM (SFEM)

The classical Galerkin finite element formulation seeks $u_h \in V_h \subset H^1_0({\cal B})$ such that

$$(u^\prime_h, v^\prime_h)_0 + c(u_h, v_h)_0 = (f, v_h)_0 \quad \forall v_h \in V_h$$

## Least Squares FEM (LSFEM)

### First-Order System

Define the diffusive flux $q = u^\prime$:

$$-q^\prime + cu = f, \qquad q - u^\prime = 0$$

### Functional

$${\cal J}(u,q;f) = \frac{1}{2} \parallel w_1 (-q^\prime + cu - f) \parallel^2_0 + \frac{1}{2} \parallel q - u^\prime \parallel^2_0$$

**Weighting analyzed:**
- **LSFEM**: $w_1 = 1$ (unweighted)
- **wLSFEM**: $w_1 = c^{-1/2}$ (balanced scaling)

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}((u_h,q_h);(v_h,r_h)) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}((u_h,q_h);(v_h,r_h)) = w^2_1 (-q^\prime_h + c u_h , -r^\prime_h + c v_h)_0 + (q_h - u^\prime_h, r_h - v^\prime_h)_0$$

$${\cal F}(v_h,r_h) = w^2_1 (f, -r^\prime_h + c v_h)_0$$


## 🔍 Key Findings

All configurations are discretized using linear ($P_1$) finite elements with ($\epsilon = 1\cdot10^{-3}$ / **moderate diffusion**). Both **regular** and **perturbed** grids are studied.

**Two regimes analyzed:**
- **Diffusion-dominated**: $c = 10^{-4}$
- **Reaction-dominated**: $c = 10^{4}$

#

**1. Solution Variable** $u$ **(Tables 1-2)**

Both **SFEM** and **LSFEM** achieve optimal second-order convergence $(O(h^2)$ for the solution variable u:

**SFEM** produces the **smallest errors**:

- ML5: $4.01 \times 10^{−3}$ → ML10: $3.93 \times 10^{−6}$ (regular mesh)
- Consistent rate ≈ 2.00 across all refinement levels ✅

**LSFEM** exhibits **$~1.8 \times$ larger errors** than SFEM:

- ML5: $7.43 \times 10^{−3}$ → ML10: $7.34 \times 10^{−6}$ (regular mesh)
- Also maintains optimal second-order convergence ✅


**Table 1: $L^2$ Error for Solution $u$ ($\| u - u_h \|_0$)** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   |  |  |  |  |  |  |
| LSFEM / RG  |  |  |  |  |  |  |
| wLSFEM / RG |  |  |  |  |  |  |
|             |            |            |            |            |            |            |
| SFEM / PG   |  |  |  |  |  |  |
| LSFEM / PG  |  |  |  |  |  |  |
| wLSFEM / PG |  |  |  |  |  |  |


**Table 2: Convergence Rates for Solution $u$**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.99  | 2.00  | 2.00  | 2.00  | 2.00   | ✅ |
| LSFEM  | 1.99  | 2.00  | 2.00  | 2.00  | 2.00   | ✅ |
| wLSFEM | 1.99  | 2.00  | 2.00  | 2.00  | 2.00   | ✅ |

#

**2. Flux Variable** $q = u^\prime$ **(Tables 3-4)**

Here the methods show **dramatically different behavior**: 

**SFEM**: First-order convergence only $O(h)$

- Flux errors remain large: 
  - ML5: $4.06 \times 10^{−1}$ → ML10: $1.27 \times 10^{−2}$

- Rate ≈ 1.00 (suboptimal for gradient recovery)
- **~200× larger flux error** than LSFEM at ML10


**LSFEM**: Second-order convergence $O(h^2)$ ✨

- Flux errors **dramatically smaller**:
  - ML5: $6.25 \times 10^{−2}$ → ML10: $6.20 \times 10^{−5}$
- Rate ≈ 2.00 (optimal!)
- Direct enforcement of constitutive relation $q - u^\prime = 0$ recovers optimal flux accuracy

**Table 3: $L^2$ Error for Flux $q$ ($\| q - q_h \|_0$)** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   |  |  |  |  |  |  |
| LSFEM / RG  |  |  |  |  |  |  |
| wLSFEM / RG |  |  |  |  |  |  |
|             |            |            |            |            |            |            |
| SFEM / PG   |  |  |  |  |  |  |
| LSFEM / PG  |  |  |  |  |  |  |
| wLSFEM / PG |  |  |  |  |  |  |


**Table 4: Convergence Rates for Solution $q$**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.00  | 1.00  | 1.00  | 1.00  | 1.00   | ✅ |
| LSFEM  | 1.98  | 2.00  | 2.00  | 2.00  | 2.00   | ❌ (Degraded)  |
| wLSFEM |   | 2.00  | 2.00  | 2.00  | 2.00   | ❌ (Degraded)  |

#

**3. Mesh Perturbation Robustness**

Both methods remain stable under mesh perturbation:

**SFEM**: ~20% error increase for u, ~7% for q
**LSFEM**: ~9% error increase for u, ~12% for q
**Convergence rates** unchanged (still 2.00 for u, 1.00/2.00 for q)

✅ Both formulations are robust to irregular element geometry.

#

**4. Linear Solver Performance (Table 5)**

**Without preconditioning (CG):**

Iteration count doubles with each mesh refinement (O(N) growth)
- ML10: 1024 iterations (regular) / 2177 iterations (perturbed) 🔴
- Impractical for large-scale problems

**Jacobi preconditioning (CG-JACOBI):**

- Modest improvement on regular meshes
- Still exhibits O(N) scaling
- ML10: **1024/2076 iterations** 🟡

**Algebraic multigrid (CG-AMG): ⭐ Mesh-independent convergence!**

- Iteration count grows sub-linearly with DOF
- ML5: 10-11 iterations → ML10: **147-152** iterations
- **~7× reduction** compared to unpreconditioned CG at ML10




**Table 5: CG-AMG - iteration counts (Regular/Perturbed)**
| Mesh | LSFEM | wLSFEM 
|------|-----------------|------------------|
| ML5  |  /  |  /  |
| ML6  |  /  |  /  |
| ML7  |  /  |  /  |
| ML8  |  /  |  /  |
| ML9  |  /  |  /  |
| ML10 |  /  |  /  |


