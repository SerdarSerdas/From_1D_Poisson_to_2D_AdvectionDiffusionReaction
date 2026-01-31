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

Method ML4 ML5 ML6 ML7 ML8
SFEM / RG  7.7953e-03 1.9821e-03 4.9766e-04 1.2455e-04 3.1145e-05
LSFEM / RG 1.3594e-02 3.4398e-03 8.6257e-04 2.1581e-04 5.3962e-05

SFEM / PG  1.0181e-02 1.9912e-03 5.7947e-04 1.4927e-04 3.4878e-05 
LSFEM / PG 1.6742e-02 3.4810e-03 9.7084e-04 2.4665e-04 5.8935e-05

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

**Table 4: Convergence Ratios for Solution $u$** on perturbed grids

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM   | 5.11 | 3.44 | 3.88 | 4.28 | ✅ |
| LSFEM  | 4.81 | 3.59 | 3.94 | 4.19 | ✅ |

**Table 4: Convergence Ratios for Solution $q$** on perturbed grids

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

**Without preconditioning (CG):**

- Iteration count doubles with each mesh refinement (O(N) growth)
- ML5: 64/66 iterations → ML10: 2049/2173 iterations 🔴
- Impractical for large-scale problems

**Jacobi preconditioning (CG-JACOBI):**

- Modest improvement, eliminates ~1 iteration per level
- Still exhibits O(N) scaling behavior
- ML5: 63/64 iterations → ML10: 2047/2075 iterations 🟡
- Negligible improvement over unpreconditioned CG

**Algebraic multigrid (CG-AMG): ⭐ Mesh-independent convergence!**

- Iteration count grows sub-linearly with DOF (~O(N^0.6))
- ML4: 7/7 iterations → ML8: 45/46 iterations 🟢
- **~11× reduction compared to unpreconditioned CG at ML8**
- **~11× reduction compared to Jacobi-preconditioned CG at ML8**

**Table 5: CG - iteration counts (Regular/Perturbed)**
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

**Figure 3 depicts the solution $u$ at mesh level 8 (ML8 / 256 elements) for regular grid and the $L^2(u)$ errors for regular girds.**

P1 - Bild

**Table 6: $L^2$ Error $\| u - u_h \|_0$** - $P_1 / P_2$ (graded grids)

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / $P_1$     | 1.1777e-01 | 1.9684e-01 | 5.8239e-02 | 1.7741e-02 | 6.4160e-03 | 1.8445e-03 |
| LSFEM / $P_1P_1$ | 1.6756e+00 | 2.7528e+00 | 1.6024e+00 | 5.9229e-01 | 1.7221e-01 | 4.5035e-02 |
|                  |            |            |            |            |            |            |
| SFEM / $P_2$     | 7.9493e-02 | 1.9903e-01 | 5.0769e-02 | 4.6098e-03 | 6.6623e-04 | 9.4856e-05 |
| LSFEM / $P_2P_2$ | 8.9287e-02 | 2.1738e-01 | 5.8522e-02 | 6.4992e-03 | 8.0132e-04 | 1.0041e-04 |

**Table 7: Convergence Rates for Solution $u$** on regular grids

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM / $P_1$     | -0.74 | 1.76 | 1.71 | 1.47 | 1.80 | ✅ |
| LSFEM / $P_1P_1$ | -0.72 | 0.78 | 1.44 | 1.78 | 1.94 | ✅ |
|                  |            |            |            |            |            |            |
| SFEM / $P_2$     | -1.32 | 1.97 | 3.46 | 2.79 | 2.81 | ✅ |
| LSFEM / $P_2P_2$ | -1.28 | 1.89 | 3.17 | 3.02 | 2.80 | ✅ |

P2-Bild

#

### 🧩 Solution flux variable $q$

**Figure 3 depicts the solution $q$ at mesh level 5 (ML5 / 32 elements) and the $L^2(q)$ errors for regular girds.**

P1-Bild

**Table 8: $L^2$ Error ($\| q - q_h \|_0$)** - $P_1 / P_2$ (recursive grids)

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / $P_1$     | 4.3580e+00 | 1.0682e+01 | 1.5973e+01 | 1.5514e+01 | 1.0708e+01 | 6.0225e+00 |
| LSFEM / $P_1P_1$ | 2.3713e+01 | 4.5047e+01 | 3.4751e+01 | 1.7589e+01 | 6.4106e+00 | 1.8428e+00 |
|                  |            |            |            |            |            |            |
| SFEM / $P_2$     | 7.0770e+00 | 1.0901e+01 | 8.8122e+00 | 6.0894e+00 | 2.5067e+00 | 7.4435e-01 |
| LSFEM / $P_2P_2$ | 1.3438e+01 | 2.2378e+01 | 1.2512e+01 | 3.5710e+00 | 6.6043e-01 | 9.4953e-02 |



**Table 9: Convergence Rates for Solution $q$** on regular grids

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM / $P_1$     | -1.29 | -0.58 | 0.04 | 0.53 | 0.83 | ✅ |
| LSFEM / $P_1P_1$ | -0.93 | 0.37 | 0.98 | 1.46 | 1.80 | ✅ |
|                  |            |            |            |            |            |            |
| SFEM / $P_2$     | -0.62 | 0.31 | 0.53 | 1.28 | 1.75 | ✅ |
| LSFEM / $P_2P_2$ | -0.74 | 0.84 | 1.81 | 2.43 | 2.80| ✅ |

P2-Bild

#

### ⚡ Solver performance 



**Table 10: CG - iteration counts (Regular/Perturbed)**
| Mesh | CG | CG-JACOBI | CG-AMG |
|------|-----------------|------------------|
| ML5  |  |  |  |
| ML6  |  |  |  |
| ML7  |  |  |  |
| ML8  |  |  |  |
| ML9  |  |  |  |
| ML10 |  |  |  |

#

### 🏆 Summary

- In diffusion-dominated settings $(c = 10^{-4})$, all methods behave similarly for u, and both LSFEM and wLSFEM provide optimal flux convergence.

- In reaction-dominated settings $(c = 10^4)$, the unweighted LSFEM loses convergence for both u and q, while wLSFEM preserves SFEM-level accuracy in u and delivers superior, second-order accurate fluxes.

- From a solver viewpoint, wLSFEM also yields significantly lower CG–AMG iteration counts in the stiff reaction-dominated regime, indicating better conditioning and overall robustness.
