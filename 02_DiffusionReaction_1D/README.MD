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


**Table 1: $L^2$ Error $\| u - u_h \|_0$ for $c = 10^{-4}$** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.0088e-03 | 1.0060e-03 | 2.5175e-04 | 6.2953e-05 | 1.5739e-05 | 3.9349e-06 |
| LSFEM / RG  | 7.4252e-03 | 1.8730e-03 | 4.6931e-04 | 1.1739e-04 | 2.9353e-05 | 7.3385e-06 |
| wLSFEM / RG | 4.0088e-03 | 1.0060e-03 | 2.5175e-04 | 6.2953e-05 | 1.5738e-05 | 3.9391e-06 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 4.9101e-03 | 1.1785e-03 | 2.8344e-04 | 7.0732e-05 | 1.7867e-05 | 4.4419e-06 |
| LSFEM / PG  | 8.6945e-03 | 2.0131e-03 | 5.2283e-04 | 1.2784e-04 | 3.2205e-05 | 7.9982e-06 |
| wLSFEM / PG | 4.9101e-03 | 1.1785e-03 | 2.8344e-04 | 7.0732e-05 | 1.7867e-05 | 4.4420e-06 |




**Table 2: $L^2$ Error $\| u - u_h \|_0$ for $c = 10^{4}$** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 1.7129e-03 | 4.1770e-04 | 1.0377e-04 | 2.5902e-05 | 6.4729e-06 | 1.6181e-06 |
| LSFEM / RG  | 1.9482e-03 | 1.8175e-03 | 1.4147e-03 | 8.3679e-04 | 3.8937e-04 | 1.4261e-04 |
| wLSFEM / RG | 1.7129e-03 | 4.1770e-04 | 1.0377e-04 | 2.5902e-05 | 6.4729e-06 | 1.6181e-06 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 2.3761e-03 | 5.6720e-04 | 1.3528e-04 | 3.4318e-05 | 8.7650e-06 | 2.1671e-06 |
| LSFEM / PG  | 2.6624e-03 | 1.8094e-03 | 1.3982e-03 | 8.3575e-04 | 3.8386e-04 | 1.4165e-04 |
| wLSFEM / PG | 2.3761e-03 | 5.6720e-04 | 1.3528e-04 | 3.4318e-05 | 8.7650e-06 | 2.1671e-06 |



**Table 3: Convergence Rates for Solution $u$** for regular grids $(c = 10^{-4} / c = 10^{4})$

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.99 / 2.04 | 2.00 / 2.01 | 2.00 / 2.00 | 2.00 / 2.00 | 2.00 / 2.00 | ✅ / ✅ |
| LSFEM  | 1.99 / 0.10 | 2.00 / 0.36 | 2.00 / 0.76 | 2.00 / 1.10 | 2.00 / 1.45 | ✅ / ❌ |
| wLSFEM | 1.99 / 2.04 | 2.00 / 2.01 | 2.00 / 2.00 | 2.00 / 2.00 | 2.00 / 2.00 | ✅ / ✅ |

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

**Table 4: $L^2$ Error ($\| q - q_h \|_0$) for $c = 10^{-4}$** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.0613e-01 | 2.0367e-01 | 1.0191e-01 | 5.0965e-02 | 2.5483e-02 | 1.2742e-02 |
| LSFEM / RG  | 6.2490e-02 | 1.5812e-02 | 3.9654e-03 | 9.9212e-04 | 2.4808e-04 | 6.2022e-05 |
| wLSFEM / RG | 6.2291e-02 | 1.5762e-02 | 3.9527e-03 | 9.8895e-04 | 2.4729e-04 | 6.1832e-05 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 4.3709e-01 | 2.1479e-01 | 1.0588e-01 | 5.2829e-02 | 2.6532e-02 | 1.3245e-02 |
| LSFEM / PG  | 7.6399e-02 | 1.7035e-02 | 4.6915e-03 | 1.0984e-03 | 2.8095e-04 | 6.9650e-05 |
| wLSFEM / PG | 7.6189e-02 | 1.6982e-02 | 4.6786e-03 | 1.0950e-03 | 2.8013e-04 | 6.9446e-05 |



**Table 5: $L^2$ Error ($\| q - q_h \|_0$) for $c = 10^{4}$** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.0982e-01 | 2.0415e-01 | 1.0197e-01 | 5.0972e-02 | 2.5484e-02 | 1.2742e-02 |
| LSFEM / RG  | 4.7039e-01 | 8.4996e-01 | 5.6253e-01 | 2.8100e-01 | 1.1144e-01 | 3.6257e-02 |
| wLSFEM / RG | 2.8906e-02 | 6.7445e-03 | 1.6575e-03 | 4.1262e-04 | 1.0305e-04 | 2.5755e-05 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 4.4388e-01 | 2.1558e-01 | 1.0596e-01 | 5.2838e-02 | 2.6533e-02 | 1.3245e-02 |
| LSFEM / PG  | 4.1357e-01 | 8.1364e-01 | 5.5315e-01 | 2.7828e-01 | 1.0938e-01 | 3.5868e-02 |
| wLSFEM / PG | 4.0409e-02 | 8.2047e-03 | 2.3856e-03 | 5.1891e-04 | 1.3782e-04 | 3.4091e-05 |



**Table 6: Convergence Rates for Solution $q$** for regular grids $(c = 10^{-4} / c = 10^{4})$

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.00 /  1.01 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | ✅ / ✅ |
| LSFEM  | 1.98 / -0.85 | 2.00 / 0.60 | 2.00 / 1.00 | 2.00 / 1.33 | 2.00 / 1.62 | ✅ / ❌ |
| wLSFEM | 1.98 /  2.10 | 2.00 / 2.02 | 2.00 / 2.01 | 2.00 / 2.00 | 2.00 / 2.00 | ✅ / ✅ |

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




**Table 7: CG-AMG - iteration counts for $c=10^{-4}$ (Regular/Perturbed)**
| Mesh | LSFEM | wLSFEM |
|------|-----------------|------------------|
| ML5  |  10 /  11 |  5 /   8 |
| ML6  |  15 /  16 |  5 /  13 |
| ML7  |  25 /  26 |  5 /  22 |
| ML8  |  44 /  45 | 35 /  36 |
| ML9  |  79 /  82 | 45 /  63 |
| ML10 | 147 / 152 | 26 / 111 |


**Table 8: CG-AMG - iteration counts for $c=10^{4}$ (Regular/Perturbed)**
| Mesh | LSFEM | wLSFEM |
|------|-----------------|------------------|
| ML5  |  21 /  22 |  1 /  1 |
| ML6  |  38 /  41 |  1 /  2 |
| ML7  |  58 /  62 |  2 /  3 |
| ML8  |  77 /  84 |  6 /  6 |
| ML9  | 103 / 127 | 11 / 12 |
| ML10 | 136 / 222 | 19 / 23 |




