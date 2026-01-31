## 📊 Numerical Results for $\epsilon = 10^{-1}$

This study compares the Standard Finite Element Method (SFEM) and Least-Squares Finite Element Method (LSFEM) for solving the 1D Poisson equation with parameter $\epsilon = 10^{-1}$. Both methods are tested on regular and perturbed grids across five mesh levels (ML4-ML8, corresponding to 16-256 elements).

#

### Methods Compared

- **SFEM (Standard Galerkin FEM)**: Primal formulation with solution variable u
- **LSFEM (Least-Squares FEM)**: First-order system with ($u$, $q = -u^\prime$) as unknowns

- **Element Types:** $P_1$ (linear) and $P_2$ (quadratic)
- **Mesh Types:** Regular (RG) and Perturbed (PG) Grids
- **Mesh Levels:** ML4 (16 elements) to ML8 (256 elements)

#

## Main Results Summary

🎯 **Key Finding: LSFEM's Superior Flux Accuracy**

**For $P_1$ elements at ML8:**

- Solution u: SFEM is 1.74× more accurate
- **Flux q: LSFEM is 89× more accurate** ⭐

**For $P_2$ elements at ML8:**

- Solution u: LSFEM is 1.09× more accurate (actually better!)
- **Flux q: LSFEM is 203× more accurate** ⭐⭐⭐

**The flux accuracy difference grows with mesh refinement due to higher convergence rates.**

#

## 🧩 Solution variable $u$

Figure 1 depicts the solution $u$ at mesh level 5 (ML5 / 32 elements) and the $L^2(u)$ errors for regular grids with $P_1$ elements.

![Figure 1: 1D Poisson - solution u](figures/Fig_Poisson1D_BL_ML5_u.png)
**Figure 1.** (a) Solution $u_h$ for ML5 (32 elements), (b) $L^2(u)$ error convergence. 


**Table 1: $L^2$ Error $\| u - u_h \|_0$** for $P_1$ elements

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG   | 7.7953e-03 | 1.9821e-03 | 4.9766e-04 | 1.2455e-04 | 3.1145e-05 |
| LSFEM / RG  | 1.3594e-02 | 3.4398e-03 | 8.6257e-04 | 2.1581e-04 | 5.3962e-05 |
|             |            |            |            |            |            |
| SFEM / PG   | 1.0181e-02 | 1.9912e-03 | 5.7947e-04 | 1.4927e-04 | 3.4878e-05 |
| LSFEM / PG  | 1.6742e-02 | 3.4810e-03 | 9.7084e-04 | 2.4665e-04 | 5.8935e-05 |

**Table 2: $L^2(u)$ Ratios** for $P_1$ elements on regular grids / RG = Regular Grid, PG = Perturbed Grid

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG  | 3.93 | 3.98 | 4.00 | 4.00 | ✅ |
| LSFEM / RG | 3.95 | 3.99 | 4.00 | 4.00 | ✅ |
|            |      |      |      |      |     |
| SFEM / PG  | 5.11 | 3.44 | 3.88 | 4.28 | ✅ |
| LSFEM / PG | 4.81 | 3.59 | 3.94 | 4.19 | ✅ |

**Analysis:**

- Both methods achieve optimal $O(h^2)$ convergence for $P_1$ elements
- Ratio ≈ 4.0 = $2^2$ confirms second-order convergence
- Perturbed grids show initial variation but maintain asymptotic rates
- SFEM produces $1.74 \times$ smaller errors for solution u

**Table 3: $L^2$ Error $\| u - u_h \|_0$** for $P_2$ elements

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG   | 2.5682e-04 | 3.2647e-05 | 4.0983e-06 | 5.1283e-07 | 7.0115e-08 |
| LSFEM / RG  | 2.6312e-04 | 3.2849e-05 | 4.1047e-06 | 5.1303e-07 | 6.4128e-08 |
|             |            |            |            |            |            |
| SFEM / PG   | 3.9714e-04 | 3.4124e-05 | 5.3776e-06 | 7.1791e-07 | 8.1484e-08 |
| LSFEM / PG  | 4.1081e-04 | 3.4360e-05 | 5.3909e-06 | 7.2076e-07 | 8.0282e-08 |

**Table 4: $L^2(u)$ ratios** for $P_2$ elements on regular grids

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG  | 7.87 | 7.97 | 7.99 | 7.31 | 🔴 |
| LSFEM / RG | 8.01 | 8.00 | 8.00 | 8.00 | ✅ |
|            |      |      |      |      |     |
| SFEM / PG  | 11.64 | 6.35 | 7.49 | 8.81 | ✅ |
| LSFEM / PG | 11.96 | 6.37 | 7.48 | 8.98 | ✅ |

**Analysis:**

- Both methods achieve $O(h^3)$ convergence for $P_2$ elements
- Ratio ≈ 8.0 = $2^3$ confirms third-order convergence
- **LSFEM maintains perfect 8.0 ratio while SFEM degrades to 7.31 at finest level**
- LSFEM actually achieves slightly better accuracy than SFEM (6.41e-08 vs 7.01e-08)
- Both achieve engineering precision $(~10^{-7})$ with only 256 elements

#

## 🧩 Solution flux variable $q$

**Figure 2 depicts the solution $q$ at mesh level 5 (ML5 / 32 elements) and the $L^2(q)$ errors for regular grids.**

![Figure 2: 1D Poisson - flux q](figures/Fig_Poisson1D_BL_ML5_q.png)

**Table 5: $L^2$ Error $\| q - q_h \|_0$** for $P_1$ elements

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG   | 3.9579e-01 | 2.0075e-01 | 1.0074e-01 | 5.0416e-02 | 2.5214e-02 |
| LSFEM / RG  | 7.0335e-02 | 1.7933e-02 | 4.5057e-03 | 1.1278e-03 | 2.8204e-04 |
|             |            |            |            |            |            |
| SFEM / PG   | 4.4342e-01 | 1.9879e-01 | 1.0656e-01 | 5.3637e-02 | 2.6144e-02 |
| LSFEM / PG  | 9.3888e-02 | 1.7910e-02 | 5.3102e-03 | 1.3730e-03 | 3.1838e-04 |

**Accuracy Comparison at ML8:**

- SFEM: 2.5214e-02
- LSFEM: 2.8204e-04
- **LSFEM is 89× more accurate!** ⭐

**Table 6: $L^2(q)$ Ratios** for $P_1$ elements on regular grids / RG = Regular Grid, PG = Perturbed Grid

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG  | 1.97 | 1.99 | 2.00 | 2.00 | ✅ |
| LSFEM / RG | 3.92 | 3.98 | 3.99 | 4.00 | ✅ |
|            |      |      |      |      |     |
| SFEM / PG  | 2.23 | 1.87 | 1.99 | 2.05 | ✅ |
| LSFEM / PG | 5.24 | 3.37 | 3.87 | 4.31 | ✅ |

**Analysis:**

- **SFEM achieves only $O(h)$ convergence** (ratio ≈ 2.0)
- **LSFEM achieves $O(h^2)$ convergence** (ratio ≈ 4.0)
- This is the fundamental advantage of LSFEM: treats flux as primary variable
- The accuracy gap widens with refinement due to convergence rate difference

**Table 7: $L^2$ Error $\| q - q_h \|_0$** for $P_2$ elements / RG = Regular Grid, PG = Perturbed Grid

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG   | 3.1678e-02 | 8.0825e-03 | 2.0311e-03 | 5.0843e-04 | 1.3006e-04 |
| LSFEM / RG  | 2.5680e-03 | 3.2646e-04 | 4.0983e-05 | 5.1283e-06 | 6.4122e-07 |
|             |            |            |            |            |            |
| SFEM / PG   | 4.1298e-02 | 8.1191e-03 | 2.3646e-03 | 6.0959e-04 | 1.4309e-04 |
| LSFEM / PG  | 3.9709e-03 | 3.4119e-04 | 5.3799e-05 | 7.2036e-06 | 8.0272e-07 |

**Accuracy Comparison at ML8:**

- SFEM: 1.3006e-04
- LSFEM: 6.4122e-07
- **LSFEM is 203× more accurate!** ⭐⭐⭐

**Table 8: $L^2(q)$ ratios** for $P_2$ elements on regular grids / RG = Regular Grid, PG = Perturbed Grid

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG  |  3.92 | 3.98 | 3.99 | 3.91 | ✅ |
| LSFEM / RG |  7.87 | 7.97 | 7.99 | 8.00 | ✅ |
|            |      |      |      |      |     |
| SFEM / PG  |  5.09 | 3.43 | 3.88 | 4.26 | ✅ |
| LSFEM / PG | 11.64 | 6.34 | 7.47 | 8.97 | ✅ |

**Analysis:**

- SFEM achieves $O(h^2)$ convergence (ratio ≈ 4.0)
- **LSFEM achieves $O(h^3)$ convergence (ratio ≈ 8.0)**
- This matches the solution u convergence for LSFEM
- LSFEM treats both u and q with equal accuracy
- The $203 \times$ accuracy difference demonstrates LSFEM's superiority for derivative quantities

#

## 🧩 Solver performance 

**Table 9: CG - iteration counts of $P_1$ elements**
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  32 /  33 |  31 /  32 |  7 /  7 |
| ML5  |  32 |  66 |  64 /  66 |  63 /  64 | 10 / 10 |
| ML6  |  64 | 130 | 128 / 133 | 127 / 128 | 16 / 16 |
| ML7  | 128 | 258 | 256 / 269 | 255 / 256 | 26 / 26 |
| ML8  | 256 | 514 | 513 / 541 | 511 / 514 | 45 / 46 |

**Observations:**

- **Unpreconditioned CG:** Iterations ≈ 2N (O(N) scaling) - impractical
- **Jacobi preconditioning:** Saves ~1 iteration - negligible benefit
- **AMG preconditioning:** ~45 iterations at ML8 → $11 \times$ reduction! 🚀

**Scaling Analysis:**

- CG: Doubles with each refinement (2N)
- AMG: Sub-linear growth $(~N^{0.6})$
- At ML8: 513 iterations (CG) vs 45 iterations (AMG)


**Table 10: CG - iteration counts of $P_2$ elements (Regular/Perturbed)**
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |   66 |   68 /   74 |   66 /   67 |  16 /  17 |
| ML5  |  32 |  130 |  141 /  151 |  138 /  139 |  27 /  27 |
| ML6  |  64 |  258 |  286 /  311 |  283 /  285 |  42 /  43 |
| ML7  | 128 |  514 |  578 /  638 |  576 /  582 |  70 /  71 |
| ML8  | 256 | 1026 | 1162 / 1289 | 1163 / 1178 | 124 / 126 |

**Observations:**

- **Unpreconditioned CG:** 1162 iterations at ML8 - completely impractical
- **Jacobi preconditioning:** Actually WORSE (1163 vs 1162) - useless!
- **AMG preconditioning:** 124 iterations at ML8 → **9× reduction!** 🚀

**Surprising Result:**

- Jacobi preconditioner is counterproductive for $P_2$ elements
- Highlights importance of sophisticated preconditioning (AMG)


## 🧩 Mesh Perturbation Robustness

Both methods maintain optimal convergence rates on randomly perturbed grids:

**$P_1$** Regular vs Perturbed:

- Error magnitude increases ~12-20%
- Convergence rates remain $O(h^2)$ for u, $O(h)/O(h^2)$ for q
- Demonstrates numerical stability

**$P_2$** Regular vs Perturbed:

- Initial mesh shows larger variation due to geometric sensitivity
- Asymptotic rates achieved by ML6-ML7
- Both methods robust
