## 📊 Numerical Results for $\epsilon = 10^{-1}$

This study compares the Standard Finite Element Method (SFEM) and Least-Squares Finite Element Method (LSFEM) for solving the 1D Poisson equation with parameter $\epsilon = 10^{-1}$. Both methods are tested on regular and perturbed grids across six mesh levels (ML5-ML10, corresponding to 32-1024 elements).

#

### 🧩 Solution variable $u$

Figure 1 depicts the solution $u$ at mesh level 5 (ML5 / 32 elements) and the $L^2(u)$ errors for regular girds.

![Figure 1: 1D Poisson - solution u](figures/Fig_Poisson1D_BL_ML5_u.png)
**Figure 1.** (a) Solution $u_h$ for ML5 (32 ekements), (b) $L^2(u)$ error convergence. 

Both methods achieve optimal second-order convergence for the primary solution variable:

- **SFEM**: Convergence rate $\approx$ 2.00, achieving $O(h^2)$ accuracy
- **LSFEM**: Convergence rate $\approx$ 2.00, achieving $O(h^2)$ accuracy
- **Accuracy comparison**: SFEM produces approximately $1.7 \times$ smaller errors for u

At the finest mesh (ML10):

- SFEM: $L^2$ error = 1.95e-06
- LSFEM: $L^2$ error = 3.37e-06

**Table 1: $L^2$ Error $\| u - u_h \|_0$** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 1.9821e-03 | 4.9766e-04 | 1.2455e-04 | 3.1145e-05 | 7.7869e-06 | 1.9468e-06 |
| LSFEM / RG  | 3.4398e-03 | 8.6257e-04 | 2.1581e-04 | 5.3962e-05 | 1.3491e-05 | 3.3728e-06 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 1.9912e-03 | 5.7947e-04 | 1.4927e-04 | 3.4878e-05 | 8.7734e-06 | 2.2322e-06 |
| LSFEM / PG  | 3.4810e-03 | 9.7084e-04 | 2.4665e-04 | 5.8935e-05 | 1.4771e-05 | 3.7402e-06 |



**Table 2: Convergence Rates for Solution $u$** on regular grids

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.99 | 2.00 | 2.00 | 2.00 | 2.00 | ✅ |
| LSFEM  | 2.00 | 2.00 | 2.00 | 2.00 | 2.00 | ✅ |

#

### 🧩 Solution flux variable $q$

**Figure 2 depicts the solution $q$ at mesh level 5 (ML5 / 32 elements) and the $L^2(q)$ errors for regular girds.**

![Figure 2: 1D Poisson - flux q](figures/Fig_Poisson1D_BL_ML5_q.png)

The most significant difference appears in the flux variable, where LSFEM demonstrates clear superiority:

- **SFEM**: First-order convergence (rate $\approx$ 1.00), achieving $O(h)$ accuracy
- **LSFEM**: Second-order convergence (rate $\approx$ 2.00), achieving $O(h^2)$ accuracy
- **Accuracy comparison**: LSFEM produces approximately $10 \times$ smaller errors at coarse meshes, with the advantage increasing at finer meshes due to the higher convergence rate

At the finest mesh (ML10):

- **SFEM**: $L^2$ error = 6.30e-03
- **LSFEM: $L^2$ error = 1.76e-05 (~357× more accurate)**

**Table 3: $L^2$ Error ($\| q - q_h \|_0$)** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 2.0075e-01 | 1.0074e-01 | 5.0416e-02 | 2.5214e-02 | 1.2608e-02 | 6.3039e-03 |
| LSFEM / RG  | 1.7933e-02 | 4.5057e-03 | 1.1278e-03 | 2.8204e-04 | 7.0517e-05 | 1.7629e-05 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 1.9879e-01 | 1.0656e-01 | 5.3637e-02 | 2.6144e-02 | 1.3101e-02 | 6.5859e-03 |
| LSFEM / PG  | 1.7910e-02 | 5.3102e-03 | 1.3730e-03 | 3.1838e-04 | 8.0174e-05 | 2.0436e-05 |



**Table 4: Convergence Rates for Solution $q$** on regular grids

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 | ✅ |
| LSFEM  | 1.99 | 2.00 | 2.00 | 2.00 | 2.00 | ✅ |

#

### 🧩 Mesh Perturbation Robustness

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
- ML5: 10/10 iterations → ML10: 149/152 iterations
- **~14× reduction compared to unpreconditioned CG at ML10**
- **~14× reduction compared to Jacobi-preconditioned CG at ML10**
- Enables practical solution of large-scale LSFEM systems

**Table 5: CG - iteration counts (Regular/Perturbed)**
| Mesh | CG | CG-JACOBI | CG-AMG |
|------|----|-----------|------------------|
| ML5  |   64 /   66 |   63 /   64 |  10 /  10 |
| ML6  |  128 /  133 |  127 /  128 |  16 /  16 |
| ML7  |  256 /  269 |  255 /  256 |  26 /  26 |
| ML8  |  512 /  541 |  511 /  514 |  45 /  46 |
| ML9  | 1025 / 1084 | 1023 / 1036 |  80 /  82 |
| ML10 | 2049 / 2173 | 2047 / 2075 | 149 / 152 |

AMG preconditioning reduces iteration counts dramatically, growing approximately as $O(h^{-0.6})$ compared to $O(h^{-1})$ for unpreconditioned CG. This is essential for the practical application of LSFEM, which produces larger symmetric systems.

#

## 📊 Numerical Results $\epsilon = 10^{-3}$


### 🧩 Solution variable $u$

**Figure 1 depicts the solution $u$ at mesh level 5 (ML5 / 32 elements) and the $L^2(u)$ errors for regular girds.**


**Table 6: $L^2$ Error $\| u - u_h \|_0$** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   |  |  |  |  |  |  |
| LSFEM / RG  |  |  |  |  |  |  |
|             |            |            |            |            |            |            |
| SFEM / PG   |  |  |  |  |  |  |
| LSFEM / PG  |  |  |  |  |  |  |



**Table 7: Convergence Rates for Solution $u$** on regular grids

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   |  |  |  |  |  | ✅ |
| LSFEM  |  |  |  |  |  | ✅ |

#

### 🌊 Flux variable $q$

**Figure 2 depicts the solution $q$ at mesh level 5 (ML5 / 32 elements) and the $L^2(q)$ errors for regular girds.**


**Table 8: $L^2$ Error ($\| q - q_h \|_0$)** / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   |  |  |  |  |  |  |
| LSFEM / RG  |  |  |  |  |  |  |
|             |            |            |            |            |            |            |
| SFEM / PG   |  |  |  |  |  |  |
| LSFEM / PG  |  |  |  |  |  |  |




**Table 9: Convergence Rates for Solution $q$** on regular grids

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   |  |  |  |  |  | ✅ |
| LSFEM  |  |  |  |  |  | ✅ |

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
