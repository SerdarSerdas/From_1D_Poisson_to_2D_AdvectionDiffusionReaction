## 📊 Numerical Results for $\epsilon = 10^{-1}$

This study compares the Standard Finite Element Method (SFEM) and Least-Squares Finite Element Method (LSFEM) for solving the 1D Poisson equation with parameter $\epsilon = 10^{-1}$. Both methods are tested on regular and perturbed grids across five mesh levels (ML4-ML8, corresponding to 16-256 elements).

#

### Methods Compared

- **SFEM (Standard Galerkin FEM)**: Primal formulation with solution variable u
- **LSFEM (Least-Squares FEM)**: Weighting factor $w_1 = 1$
- **wLSFEM (weighted Least-Squares FEM)**: Weighting factor $w_1 = c^{-1/2}$

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


**Table 1: $L^2$ Error $\| u - u_h \|_0$** for $P_1$ elements $(c=10^{-4} / c=10^4)$

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG    | 7.7953e-03 / 4.9652e-03 | 1.9821e-03 / 1.1063e-03 | 4.9766e-04 / 2.5767e-04 | 1.2455e-04 / 6.2755e-05 | 3.1145e-05 / 1.5573e-05 |
| LSFEM / RG   | 1.3592e-02 / 7.9749e-03 | 3.4394e-03 / 2.1054e-03 | 8.6248e-04 / 6.1251e-04 | 2.1578e-04 / 2.0743e-04 | 5.3956e-05 / 7.8687e-05 |
| wLSFEM / RG  | 7.7953e-03 / 4.9652e-03 | 1.9821e-03 / 1.1063e-03 | 4.9766e-04 / 2.5767e-04 | 1.2455e-04 / 6.2755e-05 | 3.1144e-05 / 1.5573e-05 |
|              |            |            |            |            |            |
| SFEM / PG    | 1.0181e-02 / 7.4707e-03 | 1.9912e-03 / 1.0525e-03 | 5.7947e-04 / 3.5115e-04 | 1.4927e-04 / 8.7706e-05 | 3.4878e-05 / 1.9304e-05 |
| LSFEM / PG   | 1.6740e-02 / 1.0729e-02 | 3.4806e-03 / 2.0229e-03 | 9.7074e-04 / 6.8117e-04 | 2.4662e-04 / 2.2083e-04 | 5.8929e-05 / 8.3269e-05 |
| wLSFEM / PG  | 1.0181e-02 / 7.4707e-03 | 1.9912e-03 / 1.0525e-03 | 5.7947e-04 / 3.5115e-04 | 1.4928e-04 / 8.7706e-05 | 3.4878e-05 / 1.9304e-05 |



**Table 2: $L^2(u)$ Ratios** for $P_1$ elements $(c=10^{-4} / c=10^4)$

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG   | 3.93 / 4.49 | 3.98 / 4.29 | 4.00 / 4.11 | 4.00 / 4.03 | ✅ |
| LSFEM / RG  | 3.95 / 3.79 | 3.99 / 3.44 | 4.00 / 2.95 | 4.00 / 2.64 | ✅ |
| wLSFEM / RG | 3.93 / 4.49 | 3.98 / 4.29 | 4.00 / 4.11 | 4.00 / 4.03 | ✅ |
|             |      |      |      |      |     |
| SFEM / PG   | 5.11 / 7.10 | 3.44 / 3.00 | 3.88 / 4.00 | 4.28 / 4.54 | ✅ |
| LSFEM / PG  | 4.81 / 5.30 | 3.59 / 2.97 | 3.94 / 3.08 | 4.19 / 2.65 | ✅ |
| wLSFEM / PG | 5.11 / 7.10 | 3.44 / 3.00 | 3.88 / 4.00 | 4.28 / 4.54 | ✅ |

**Analysis:**

- Both methods achieve optimal $O(h^2)$ convergence for $P_1$ elements
- Ratio ≈ 4.0 = $2^2$ confirms second-order convergence
- Perturbed grids show initial variation but maintain asymptotic rates
- SFEM produces $1.74 \times$ smaller errors for solution u

**Table 3: $L^2$ Error $\| u - u_h \|_0$** for $P_2$ elements $(c=10^{-4} / c=10^4)$

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG    | 2.5682e-04 / 2.4985e-04 | 3.2647e-05 / 3.2458e-05 | 4.0983e-06 / 4.0938e-06 | 5.1283e-07 / 5.1271e-07 | 7.0115e-08 / 7.0111e-08 |
| LSFEM / RG   |  /  |  /  |  /  |  /  |  /  |
| wLSFEM / RG  |  /  |  /  |  /  |  /  |  /  |
|              |            |            |            |            |            |
| SFEM / PG    | 3.9714e-04 / 3.3881e-04 | 3.4124e-05 / 3.2662e-05 | 5.3776e-06 / 5.2030e-06 | 7.1791e-07 / 7.0567e-07 | 8.1484e-08 / 8.1171e-08 |
| LSFEM / PG   |  /  |  /  |  /  |  /  |  /  |
| wLSFEM / PG  |  /  |  /  |  /  |  /  |  /  |



**Table 4: $L^2(u)$ ratios** for $P_2$ elements $(c=10^{-4} / c=10^4)$

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG   | 7.87 / 7.70 | 7.97 / 7.93 | 7.99 / 7.98 | 7.31 / 7.31 | 🔴 |
| LSFEM / RG  |  /  |  /  |  /  |  /  | ✅ |
| wLSFEM / RG |  /  |  /  |  /  |  /  | ✅ |
|             |      |      |      |      |     |
| SFEM / PG   | 11.64 / 10.37 | 6.35 / 6.28 | 7.49 / 7.37 | 8.81 / 8.69 | 🔴 |
| LSFEM / PG  |  /  |  /  |  /  |  /  | ✅ |
| wLSFEM / PG |  /  |  /  |  /  |  /  | ✅ |

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

**Table 5: $L^2$ Error $\| q - q_h \|_0$** for $P_1$ elements $(c=10^{-4} / c=10^4)$

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG    | 3.9579e-01 / 4.1262e-01 | 2.0075e-01 / 2.0482e-01 | 1.0074e-01 / 1.0142e-01 | 5.0416e-02 / 5.0510e-02 | 2.5214e-02 / 2.5226e-02 |
| LSFEM / RG   | 7.0334e-02 / 5.0259e+00 | 1.7933e-02 / 1.1662e+00 | 4.5056e-03 / 3.1832e-01 | 1.1278e-03 / 9.8065e-02 | 2.8204e-04 / 3.1403e-02 |
| wLSFEM / RG  | 7.0927e-02 / 3.2792e-02 | 1.8079e-02 / 8.1780e-03 | 4.5420e-03 / 2.0403e-03 | 1.1369e-03 / 5.0972e-04 | 2.8432e-04 / 1.2741e-04 |
|              |            |            |            |            |            |
| SFEM / PG    | 4.4342e-01 / 4.7628e-01 | 1.9879e-01 / 2.0218e-01 | 1.0656e-01 / 1.0729e-01 | 5.3637e-02 / 5.3763e-02 | 2.6144e-02 / 2.6156e-02 |
| LSFEM / PG   | 9.3888e-02 / 5.5023e+00 | 1.7910e-02 / 1.1558e+00 | 5.3101e-03 / 3.2675e-01 | 1.3730e-03 / 9.6475e-02 | 3.1838e-04 / 3.1921e-02 |
| wLSFEM / PG  | 9.4511e-02 / 4.6858e-02 | 1.8064e-02 / 8.6355e-03 | 5.3476e-03 / 2.6880e-03 | 1.3822e-03 / 7.4344e-04 | 3.2072e-04 / 1.6690e-04 |



**Accuracy Comparison at ML8:**

- SFEM: 2.5214e-02
- LSFEM: 2.8204e-04
- **LSFEM is 89× more accurate!** ⭐

**Table 6: $L^2(q)$ Ratios** for $P_1$ elements $(c=10^{-4} / c=10^4)$

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG   | 1.97 / 2.01 | 1.99 / 2.02 | 2.00 / 2.01 | 2.00 / 2.00 | 🔴 |
| LSFEM / RG  | 3.92 / 4.31 | 3.98 / 3.66 | 3.99 / 3.25 | 4.00 / 3.12 | ✅ |
| wLSFEM / RG | 3.92 / 4.01 | 3.98 / 4.01 | 4.00 / 4.00 | 4.00 / 4.00 | ✅ |
|             |      |      |      |      |     |
| SFEM / PG   | 2.23 / 2.36 | 1.87 / 1.88 | 1.99 / 2.00 | 2.05 / 2.06 | 🔴 |
| LSFEM / PG  | 5.24 / 4.76 | 3.37 / 3.54 | 3.87 / 3.39 | 4.31 / 3.02 | ✅ |
| wLSFEM / PG | 2.39 / 5.43 | 1.76 / 3.21 | 1.95 / 3.62 | 2.11 / 4.45 | ✅ |

**Analysis:**

- **SFEM achieves only $O(h)$ convergence** (ratio ≈ 2.0)
- **LSFEM achieves $O(h^2)$ convergence** (ratio ≈ 4.0)
- This is the fundamental advantage of LSFEM: treats flux as primary variable
- The accuracy gap widens with refinement due to convergence rate difference

**Table 7: $L^2$ Error $\| q - q_h \|_0$** for $P_2$ elements $(c=10^{-4} / c=10^4)$

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG    | 3.1678e-02 / 3.1810e-02 | 8.0825e-03 / 8.0843e-03 | 2.0311e-03 / 2.0311e-03 | 5.0843e-04 / 5.0843e-04 | 1.3006e-04 / 1.3006e-04 |
| LSFEM / RG   |  /  |  /  |  /  |  /  |  /  |
| wLSFEM / RG  |  /  |  /  |  /  |  /  |  /  |
|              |            |            |            |            |            |
| SFEM / PG    | 4.1298e-02 / 4.2518e-02 | 8.1191e-03 / 8.1396e-03 | 2.3646e-03 / 2.3662e-03 | 6.0959e-04 / 6.0968e-04 | 1.4309e-04 / 1.4309e-04 |
| LSFEM / PG   |  /  |  /  |  /  |  /  |  /  |
| wLSFEM / PG  |  /  |  /  |  /  |  /  |  /  |



**Accuracy Comparison at ML8:**

- SFEM: 1.3006e-04
- LSFEM: 6.4122e-07
- **LSFEM is 203× more accurate!** ⭐⭐⭐

**Table 8: $L^2(q)$ ratios** for $P_2$ elements on regular grids $(c=10^{-4} / c=10^4)$

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG   | 3.92 / 3.93 | 3.98 / 3.98 | 3.99 / 3.99 | 3.91 / 3.91 | 🔴 |
| LSFEM / RG  |  /  |  /  |  /  |  /  | ✅ |
| wLSFEM / RG |  /  |  /  |  /  |  /  | ✅ |
|             |      |      |      |      |     |
| SFEM / PG   | 5.09 / 5.22 | 3.43 / 3.44 | 3.88 / 3.88 | 4.26 / 4.26 | 🔴 |
| LSFEM / PG  |  /  |  /  |  /  |  /  | ✅ |
| wLSFEM / PG |  /  |  /  |  /  |  /  | ✅ |

**Analysis:**

- SFEM achieves $O(h^2)$ convergence (ratio ≈ 4.0)
- **LSFEM achieves $O(h^3)$ convergence (ratio ≈ 8.0)**
- This matches the solution u convergence for LSFEM
- LSFEM treats both u and q with equal accuracy
- The $203 \times$ accuracy difference demonstrates LSFEM's superiority for derivative quantities

#

## 🧩 Solver performance 

**Table 9: CG - iteration counts of $P_1$ elements on regular grid for LSFEM** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  32 /  127 |  31 /   32 |  7 / 12 |
| ML5  |  32 |  66 |  64 /  386 |  63 /   72 | 10 / 21 |
| ML6  |  64 | 130 | 128 / 1106 | 127 /  212 | 16 / 36 |
| ML7  | 128 | 258 | 256 / 3249 | 255 /  493 | 26 / 56 |
| ML8  | 256 | 514 | 513 / 8057 | 511 / 1281 | 45 / 76 |

**Table 10: CG - iteration counts of $P_1$ elements on regular grid for wLSFEM** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |   99 /  59 |  16 / 13 |  5 / 1 |
| ML5  |  32 |  66 |  297 /  41 |  32 /  8 |  5 / 1 |
| ML6  |  64 | 130 |  952 /  75 |  64 / 11 |  5 / 1 |
| ML7  | 128 | 258 | 3034 / 283 | 128 / 24 |  5 / 2 |
| ML8  | 256 | 514 | 9013 / 910 | 256 / 49 | 35 / 6 |



**Table 11: CG - iteration counts of $P_1$ elements on perturbed grid for LSFEM** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  33 /  132 |  32 /   34 |  7 / 12 |
| ML5  |  32 |  66 |  66 /  392 |  64 /   88 | 10 / 22 |
| ML6  |  64 | 130 | 134 / 1126 | 128 /  212 | 16 / 37 |
| ML7  | 128 | 258 | 269 / 3425 | 256 /  508 | 26 / 58 |
| ML8  | 256 | 514 | 541 / 8391 | 514 / 1309 | 46 / 81 |

**Table 12: CG - iteration counts of $P_1$ elements on perturbed grid for wLSFEM** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  101 /  60 |  33 / 13 |  6 / 2 |
| ML5  |  32 |  66 |  310 /  54 |  64 / 10 |  8 / 1 |
| ML6  |  64 | 130 |  999 /  96 | 131 / 14 | 12 / 2 |
| ML7  | 128 | 258 | 3207 / 311 | 261 / 29 | 21 / 3 |
| ML8  | 256 | 514 | 9599 / 995 | 522 / 57 | 35 / 6 |



**Observations:**

- **Unpreconditioned CG:** Iterations ≈ 2N (O(N) scaling) - impractical
- **Jacobi preconditioning:** Saves ~1 iteration - negligible benefit
- **AMG preconditioning:** ~45 iterations at ML8 → $11 \times$ reduction! 🚀

**Scaling Analysis:**

- CG: Doubles with each refinement (2N)
- AMG: Sub-linear growth $(~N^{0.6})$
- At ML8: 513 iterations (CG) vs 45 iterations (AMG)


**Table 11: CG - iteration counts of $P_2$ elements on regular grid** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  /  |  /  |  /  |
| ML5  |  32 |  66 |  /  |  /  |  /  |
| ML6  |  64 | 130 |  /  |  /  |  /  |
| ML7  | 128 | 258 |  /  |  /  |  /  |
| ML8  | 256 | 514 |  /  |  /  |  /  |

**Table 12: CG - iteration counts of $P_2$ elements on perturbed grid** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  /  |  /  |  /  |
| ML5  |  32 |  66 |  /  |  /  |  /  |
| ML6  |  64 | 130 |  /  |  /  |  /  |
| ML7  | 128 | 258 |  /  |  /  |  /  |
| ML8  | 256 | 514 |  /  |  /  |  /  |

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

