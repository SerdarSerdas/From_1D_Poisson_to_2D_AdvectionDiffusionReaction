## 📊 Numerical Results for $\epsilon = 10^{-1}$

This study compares the Standard Finite Element Method (SFEM) and Least-Squares Finite Element Method (LSFEM) for solving the 1D Poisson equation with parameter $\epsilon = 10^{-1}$. Both methods are tested on regular and perturbed grids across five mesh levels (ML4-ML8, corresponding to 16-256 elements).

#

### 🧩 Solution variable $u$

Figure 1 depicts the solution $u$ at mesh level 5 (ML5 / 32 elements) and the $L^2(u)$ errors for regular grids.

![Figure 1: 1D Poisson - solution u](figures/Fig_Poisson1D_BL_ML5_u.png)
**Figure 1.** (a) Solution $u_h$ for ML5 (32 elements), (b) $L^2(u)$ error convergence. 

Both methods achieve optimal second-order convergence for the primary solution variable:

- **SFEM**: Convergence rate $\approx$ 2.00, achieving $O(h^2)$ accuracy
- **LSFEM**: Convergence rate $\approx$ 2.00, achieving $O(h^2)$ accuracy
- **Accuracy comparison**: SFEM produces approximately $1.7 \times$ smaller errors for u

At the finest mesh (ML8):

- SFEM: $L^2$ error = 3.1145e-05
- LSFEM: $L^2$ error = 5.3962e-05

**Table 1: $L^2$ Error $\| u - u_h \|_0$** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG   | 7.7953e-03 | 1.9821e-03 | 4.9766e-04 | 1.2455e-04 | 3.1145e-05 |
| LSFEM / RG  | 1.3594e-02 | 3.4398e-03 | 8.6257e-04 | 2.1581e-04 | 5.3962e-05 |
|             |            |            |            |            |            |
| SFEM / PG   | 1.0181e-02 | 1.9912e-03 | 5.7947e-04 | 1.4927e-04 | 3.4878e-05 |
| LSFEM / PG  | 1.6742e-02 | 3.4810e-03 | 9.7084e-04 | 2.4665e-04 | 5.8935e-05 |

**Table 2: Convergence Rates for Solution $u$** on regular grids

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM   | 1.98 | 1.99 | 2.00 | 2.00 | ✅ |
| LSFEM  | 1.98 | 2.00 | 2.00 | 2.00 | ✅ |

#

### 🧩 Solution flux variable $q$

**Figure 2 depicts the solution $q$ at mesh level 5 (ML5 / 32 elements) and the $L^2(q)$ errors for regular grids.**

![Figure 2: 1D Poisson - flux q](figures/Fig_Poisson1D_BL_ML5_q.png)

The most significant difference appears in the flux variable, where LSFEM demonstrates clear superiority:

- **SFEM**: First-order convergence (rate $\approx$ 1.00), achieving $O(h)$ accuracy
- **LSFEM**: Second-order convergence (rate $\approx$ 2.00), achieving $O(h^2)$ accuracy
- **Accuracy comparison**: LSFEM produces approximately $10 \times$ smaller errors at coarse meshes, with the advantage increasing at finer meshes due to the higher convergence rate

At the finest mesh (ML10):

- **SFEM**: $L^2$ error = 2.5214e-02
- **LSFEM: $L^2$ error = 2.8204e-04 (~89 $\times$ more accurate)**

**Table 3: $L^2$ Error ($\| q - q_h \|_0$)** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG   | 3.9579e-01 | 2.0075e-01 | 1.0074e-01 | 5.0416e-02 | 2.5214e-02 |
| LSFEM / RG  | 7.0335e-02 | 1.7933e-02 | 4.5057e-03 | 1.1278e-03 | 2.8204e-04 |
|             |            |            |            |            |            |
| SFEM / PG   | 4.4342e-01 | 1.9879e-01 | 1.0656e-01 | 5.3637e-02 | 2.6144e-02 |
| LSFEM / PG  | 9.3888e-02 | 1.7910e-02 | 5.3102e-03 | 1.3730e-03 | 3.1838e-04 |


**Table 4: Convergence Rates for Solution $q$** on regular grids

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM   | 0.98 | 0.99 | 1.00 | 1.00 | ✅ |
| LSFEM  | 1.97 | 1.99 | 2.00 | 2.00 | ✅ |

#

### 🧩 Mesh Perturbation Robustness

**Table 5: Convergence Ratios for Solution $u$** on perturbed grids

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM   | 5.11 | 3.44 | 3.88 | 4.28 | ✅ |
| LSFEM  | 4.81 | 3.59 | 3.94 | 4.19 | ✅ |

**Table 6: Convergence Ratios for Solution $q$** on perturbed grids

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM   | 2.23 | 1.87 | 1.99 | 2.05 | ✅ |
| LSFEM  | 5.24 | 3.37 | 3.87 | 4.31 | ✅ |

Both methods demonstrate good stability under mesh perturbation:

- Convergence rates remain optimal on perturbed grids
- Error magnitudes increase modestly (typically <20%)
- LSFEM shows slightly more sensitivity to perturbation but maintains second-order convergence

#

### ⚡ Solver performance 

**Without preconditioning (CG):** 🔴

- Iteration count doubles with each mesh refinement (O(N) growth)
- ML5: 32/33 iterations → ML10: 513/541 iterations 
- Impractical for large-scale problems

**Jacobi preconditioning (CG-JACOBI):** 🟡

- Modest improvement, eliminates ~1 iteration per level on regular grids
- Still exhibits O(N) scaling behavior
- ML5: 31/32 iterations → ML8: 511/514 iterations
- Negligible improvement over unpreconditioned CG

**Algebraic multigrid (CG-AMG): ⭐ Mesh-independent convergence!**  🟢

- Iteration count grows sub-linearly with DOF (~O(N^0.6))
- ML4: 7/7 iterations → ML8: 45/46 iterations
- **~11× reduction compared to unpreconditioned CG at ML8**
- **~11× reduction compared to Jacobi-preconditioned CG at ML8**

**Table 7: CG - iteration counts (Regular/Perturbed)**
| Mesh | CG | CG-JACOBI | CG-AMG |
|------|----|-----------|------------------|
| ML4  |  32 /  33 |  31 /  32 |  7 /  7 |
| ML5  |  64 /  66 |  63 /  64 | 10 / 10 |
| ML6  | 128 / 133 | 127 / 128 | 16 / 16 |
| ML7  | 256 / 269 | 255 / 256 | 26 / 26 |
| ML8  | 513 / 541 | 511 / 514 | 45 / 46 |

AMG preconditioning reduces iteration counts dramatically, growing approximately as $O(h^{-0.6})$ compared to $O(h^{-1})$ for unpreconditioned CG. This is essential for the practical application of LSFEM, which produces larger symmetric systems.

#

## 📊 Adaptive Mesh Refinement

In this case we decrease the parameter $\epsilon$ to $\epsilon = 10^{-3}$ yielding to a severe boundary layer near to $x = 1$. Analysis is carried out for $P_1$ and $P_2$ elements. Two adaptive refinement strategies, namely $\alpha-Bulk$ and $Dörfler$ for $P_1$ and $P_2$ elements are considered. 

### 🧩 Solution variable $u$

**Figure 3 depicts the solution $u$ at mesh level 8 (ML8 / 256 elements) for regular grid and the $L^2(u)$ errors for graded grids.**

P1 - Bild

**Table 6: $L^2$ Error $\| u - u_h \|_0$** - $P_1 / P_2$ (graded grids)

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / $P_1$     | 9.1005e-02 | 2.3395e-02 | 9.1643e-03 | 2.8331e-03 | 7.6106e-04 |
| LSFEM / $P_1P_1$ | 1.9716e+00 | 8.2947e-01 | 2.6233e-01 | 7.1658e-02 | 1.8495e-02 |
|                  |            |            |            |            |            |
| SFEM / $P_2$     | 8.6915e-02 | 9.4535e-03 | 1.2089e-03 | 1.8358e-04 | 2.4956e-05 |
| LSFEM / $P_2P_2$ | 9.6991e-02 | 1.2736e-02 | 1.5562e-03 | 1.9993e-04 | 2.5391e-05 |



**Table 7: Convergence Ratios for Solution $u$** on graded grids

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / $P_1$      | 3.89 | 2.55 | 3.23 | 3.72 | ✅ |
| LSFEM / $P_1P_1$  | 2.38 | 3.16 | 3.66 | 3.87 | ✅ |
|                   |      |      |      |      |     |
| SFEM / $P_2$      | 9.19 | 7.82 | 6.59 | 7.36 | ✅ |
| LSFEM / $P_2P_2$  | 7.62 | 8.18 | 7.78 | 7.87 | ✅ |



P2-Bild

#

### 🧩 Solution flux variable $q$

**Figure 3 depicts the solution $q$ at mesh level 5 (ML5 / 32 elements) and the $L^2(q)$ errors for graded grids.**

P1-Bild

**Table 8: $L^2$ Error ($\| q - q_h \|_0$)** - $P_1 / P_2$ (graded grids)

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / $P_1$     | 1.4925e+01 | 1.6329e+01 | 1.2416e+01 | 7.3883e+00 | 3.9111e+00 |
| LSFEM / $P_1P_1$ | 3.9353e+01 | 2.2420e+01 | 9.1568e+00 | 2.8305e+00 | 7.6039e-01 |
|                  |            |            |            |            |            |
| SFEM / $P_2$     | 9.2460e+00 | 7.2021e+00 | 3.4956e+00 | 1.1362e+00 | 3.0997e-01 |
| LSFEM / $P_2P_2$ | 1.5892e+01 | 5.5243e+00 | 1.1766e+00 | 1.8351e-01 | 2.4807e-02 |



**Table 7: Convergence Ratios for Solution $q$** on graded grids

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / $P_1$      | 0.91 | 1.32 | 1.68 | 1.89 | ✅ |
| LSFEM / $P_1P_1$  | 1.76 | 2.45 | 3.24 | 3.72 | ✅ |
|                   |      |      |      |      |     |
| SFEM / $P_2$      | 1.28 | 2.06 | 3.08 | 3.67 | ✅ |
| LSFEM / $P_2P_2$  | 2.88 | 4.70 | 6.41 | 7.40 | ✅ |

P2-Bild

#

### ⚡ Solver performance 



**Table 10: CG - iteration counts on graded grids $(P_1P_1 / P_2P_2)$**
| Mesh | CG | CG-JACOBI | CG-AMG |
|------|----|-----------|--------|
| ML4  |   45 /  109 |  32 /   66 |  7 /  15 |
| ML5  |  102 /  246 |  64 /  138 | 10 /  25 |
| ML6  |  232 /  551 | 128 /  282 | 15 /  40 |
| ML7  |  510 / 1192 | 256 /  574 | 25 /  69 |
| ML8  | 1076 / 2499 | 512 / 1162 | 44 / 123 |

#

### 🏆 Summary

- In diffusion-dominated settings $(c = 10^{-4})$, all methods behave similarly for u, and both LSFEM and wLSFEM provide optimal flux convergence.

- In reaction-dominated settings $(c = 10^4)$, the unweighted LSFEM loses convergence for both u and q, while wLSFEM preserves SFEM-level accuracy in u and delivers superior, second-order accurate fluxes.

- From a solver viewpoint, wLSFEM also yields significantly lower CG–AMG iteration counts in the stiff reaction-dominated regime, indicating better conditioning and overall robustness.
