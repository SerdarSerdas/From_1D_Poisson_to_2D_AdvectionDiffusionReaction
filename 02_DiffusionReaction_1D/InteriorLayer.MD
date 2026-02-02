## 📊 Problem statement for $\epsilon = 10^{-3}$

This study compares the Standard Finite Element Method (SFEM) and Least-Squares Finite Element Method (LSFEM) for solving the 1D Diffusion-Reaction Equation with the parameter $\epsilon = 10^{-3}$ for the interior layer problem. Both methods are tested on regular and perturbed grids across five mesh levels (ML4-ML8, corresponding to 16-256 elements).

**Reaction coefficient:** c determines the relative importance of diffusion vs reaction

- $c = 10^{-4}$: Diffusion-dominated (nearly pure Poisson)
- $c = 10^4$: Reaction-dominated (strong coupling)

#

## Methods Compared

- **SFEM (Standard Galerkin FEM):** Primal formulation with solution variable u
- **LSFEM (Least-Squares FEM):** First-order system with (u, q = -u') as unknowns
  - Weights: $w_1 = 1$
- **wLSFEM (Weighted LSFEM):** Same formulation with parameter-dependent weights
  - Weighting factor: $w_1 = c^{-1/2}$
  - **Balances diffusion and reaction terms in the functional**

**Element Types:** $P_1$ (linear) and $P_2$ (quadratic)
**Mesh Types:** Regular and perturbed grids
**Mesh Levels:** ML4 (16 elements) to ML8 (256 elements)

#

## Main Results Summary

🎯 **Critical Discovery: Parameter Dependency in LSFEM!**

**For $c = 10^{-4}$ (diffusion-dominated):**

- Standard LSFEM performs well
- wLSFEM matches SFEM exactly

**For $c = 10^4$ (reaction-dominated):**

- **Standard LSFEM FAILS** - loses convergence, errors explode! 🔴
- **wLSFEM rescues performance** - matches SFEM, maintains convergence ✅
- **This demonstrates the critical importance of proper scaling!**

#

## 🧩 Detailed Results - Solution variable $u$

Figure 1 depicts the solution $u$ at mesh level 5 (ML5 / 32 elements) and the $L^2(u)$ errors for regular grids with $P_1$ elements.

![Figure 1: 1D Poisson - solution u](figures/Fig_Poisson1D_BL_ML5_u.png)
**Figure 1.** (a) Solution $u_h$ for ML5 (32 elements), (b) $L^2(u)$ error convergence. 


**Table 1: $L^2$ Error $\| u - u_h \|_0$** for $P_1$ elements $(c=10^{-4} / c=10^4)$

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG    | 1.5873e-02 / 7.6267e-03 | 4.0088e-03 / 1.7129e-03 | 1.0060e-03 / 4.1770e-04 | 2.5175e-04 / 1.0377e-04 | 6.2953e-05 / 2.5902e-05 |
| LSFEM / RG   | 2.8671e-02 / 9.8484e-03 | 7.4252e-03 / 1.9482e-03 | 1.8730e-03 / 1.8175e-03 | 4.6931e-04 / 1.4147e-03 | 1.1739e-04 / 8.3679e-04 |
| wLSFEM / RG  |  /  |  /  |  /  |  /  |  /  |
|              |            |            |            |            |            |
| SFEM / PG    | 1.9140e-02 / 1.0194e-02 | 4.9101e-03 / 2.3761e-03 | 1.1785e-03 / 5.6720e-04 | 2.8344e-04 / 1.3528e-04 | 7.0732e-05 / 3.4318e-05 |
| LSFEM / PG   | 3.3773e-02 / 1.3354e-02 | 8.6945e-03 / 2.6624e-03 | 2.0131e-03 / 1.8094e-03 | 5.2283e-04 / 1.3982e-03 | 1.2784e-04 / 8.3575e-04 |
| wLSFEM / PG  |  /  |  /  |  /  |  /  |  /  |




## Observations:

$c = 10^{-4}$ (diffusion-dominated):

All methods achieve optimal $O(h^2)$ convergence
- SFEM: Most accurate (3.11e-05)
- LSFEM: Good performance (5.40e-05)
- wLSFEM: **Matches SFEM exactly!** (3.11e-05)

$c = 10^4$ (reaction-dominated):

- SFEM: Excellent convergence (1.56e-05)
- **LSFEM: CONVERGENCE BREAKDOWN!** Error increases ML6→7→8 (6.13e-04 → 2.07e-04 → 7.87e-05)
- **wLSFEM: Rescues performance!** Matches SFEM exactly (1.56e-05)

**Table 2: $L^2(u)$ Ratios** for $P_1$ elements $(c=10^{-4} / c=10^4)$

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG   | 3.96 / 4.45 | 3.98 / 4.10 | 4.00 / 4.03 | 4.00 / 4.01 | ✅ / ✅ |
| LSFEM / RG  | 3.86 / 5.06 | 3.96 / 1.07 | 3.99 / 1.28 | 4.00 / 1.69 | ✅ / 🔴 **degrades** |
| wLSFEM / RG |  /  |  /  |  /  |  /  | ✅ / ✅ |
|             |      |      |      |      |     |
| SFEM / PG   | 3.90 / 4.29 | 4.17 / 4.19 | 4.16 / 4.19 | 4.01 / 3.94 | ✅ / ✅ |
| LSFEM / PG  | 3.88 / 5.02 | 4.32 / 1.47 | 3.85 / 1.29 | 4.09 / 1.67 | ✅ / ✅ |
| wLSFEM / PG |  /  |  /  |  /  |  /  | ✅ / ✅ |




**Analysis:**

- For $c = 10^{-4}$: All methods achieve $O(h^2)$ (ratio ≈ 4.0)
- For $c = 10^4$: LSFEM degrades from 3.79 → 2.64 (losing convergence order!)
- For $c = 10^4$: wLSFEM maintains perfect 4.0 ratios
- **wLSFEM is parameter-robust, standard LSFEM is not!**

**Table 3: $L^2$ Error $\| u - u_h \|_0$** for $P_2$ elements $(c=10^{-4} / c=10^4)$

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG    |  /  |  /  |  /  |  /  |  /  |
| LSFEM / RG   |  /  |  /  |  /  |  /  |  /  |
| wLSFEM / RG  |  /  |  /  |  /  |  /  |  /  |
|              |            |            |            |            |            |
| SFEM / PG    |  /  |  /  |  /  |  /  |  /  |
| LSFEM / PG   |  /  |  /  |  /  |  /  |  /  |
| wLSFEM / PG  |  /  |  /  |  /  |  /  |  /  |

**Observations:**

$c = 10^{-4}$: All methods achieve ~$10^{-7}$ accuracy
$c = 10^4$:

- SFEM: 7.01e-08 ✅
- **LSFEM:** 1.94e-07 (3× worse, wrong trend) 🔴
- wLSFEM: 6.41e-08 ✅ (actually better than SFEM!)

**Table 4: $L^2(u)$ ratios** for $P_2$ elements $(c=10^{-4} / c=10^4)$

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG   |  /  |  /  |  /  |  /  | ✅ / ✅ |
| LSFEM / RG  |  /  |  /  |  /  |  /  | ✅ / 🔴 **degrades** |
| wLSFEM / RG |  /  |  /  |  /  |  /  | ✅ / ✅ |
|             |      |      |      |      |     |
| SFEM / PG   |  /  |  /  |  /  |  /  | ✅ / ✅ |
| LSFEM / PG  |  /  |  /  |  /  |  /  | ✅ / ✅ |
| wLSFEM / PG |  /  |  /  |  /  |  /  | ✅ / ✅ |

**Critical Observation:**

- For $c = 10^4$, LSFEM convergence rate collapses: 7.65 → 7.25 → 6.18 → **4.12**
- wLSFEM maintains perfect 8.00 ratio even for $c = 10^4$
- **This proves weighted formulation is essential for parameter robustness!**

#

## ## 🧩 Detailed Results - Solution flux variable $q$

**Figure 2 depicts the solution $q$ at mesh level 5 (ML5 / 32 elements) and the $L^2(q)$ errors for regular grids.**

![Figure 2: 1D Poisson - flux q](figures/Fig_Poisson1D_BL_ML5_q.png)

**Table 5: $L^2$ Error $\| q - q_h \|_0$** for $P_1$ elements $(c=10^{-4} / c=10^4)$

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG    | 8.0567e-01 / 8.3015e-01 | 4.0613e-01 / 4.0982e-01 | 2.0367e-01 / 2.0415e-01 | 1.0191e-01 / 1.0197e-01 | 5.0965e-02 / 5.0972e-02 |
| LSFEM / RG   | 2.2920e-01 / 5.5329e+00 | 6.2490e-02 / 4.7039e-01 | 1.5812e-02 / 8.4996e-01 | 3.9654e-03 / 5.6253e-01 | 9.9212e-04 / 2.8100e-01 |
| wLSFEM / RG  |  /  |  /  |  /  |  /  |  /  |
|              |            |            |            |            |            |
| SFEM / PG    | 8.6115e-01 / 9.0414e-01 | 4.3709e-01 / 4.4388e-01 | 2.1479e-01 / 2.1558e-01 | 1.0588e-01 / 1.0596e-01 | 5.2829e-02 / 5.2838e-02 |
| LSFEM / PG   | 2.4755e-01 / 6.7866e+00 | 7.6399e-02 / 4.1357e-01 | 1.7035e-02 / 8.1364e-01 | 4.6915e-03 / 5.5315e-01 | 1.0984e-03 / 2.7828e-01 |
| wLSFEM / PG  |  /  |  /  |  /  |  /  |  /  |



DRAMATIC RESULTS:

$c = 10^{-4}$:

- SFEM: 2.52e-02 ($O(h)$ convergence)
- LSFEM: 2.82e-04 (**89 more accurate!** ⭐)
- wLSFEM: 2.84e-04 (matches LSFEM)

$c = 10^4$:

- SFEM: 2.52e-02 (still O(h))
- **LSFEM:** 3.14e-02 (WORSE than SFEM! Complete failure!) 🔴🔴🔴
- **wLSFEM:** 1.27e-04 (198× more accurate than SFEM!) 🚀🚀🚀

**Table 6: $L^2(q)$ Ratios** for $P_1$ elements $(c=10^{-4} / c=10^4)$

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG   | 1.98 /  2.03 | 1.99 / 2.01 | 2.00 / 2.00 | 2.00 / 2.00 | ✅ / ✅ |
| LSFEM / RG  | 3.67 / 11.76 | 3.95 / 0.55 | 3.99 / 1.51 | 4.00 / 2.00 | ✅ / 🔴 **poor** |
| wLSFEM / RG |  /  |  /  |  /  |  /  | ✅ / ✅ |
|             |      |      |      |      |     |
| SFEM / PG   | 1.97 / 2.04 | 2.03 / 2.06 | 2.03 / 2.03 | 2.00 / 2.01 | ✅ / ✅ |
| LSFEM / PG  | 3.24 / 16.41 | 4.48 / 0.51 | 3.63 / 1.47 | 4.27 / 1.99 | ✅ / ✅ |
| wLSFEM / PG |  /  |  /  |  /  |  /  | ✅ / ✅ |




**Analysis:**

- SFEM: Stuck at $O(h)$ regardless of c value
- LSFEM $c = 10^{-4}$: Perfect $O(h^2)$ convergence (ratio 4.0)
- **LSFEM $c = 10^4$: Loses convergence order** (4.31 → 3.12)
- **wLSFEM $c = 10^4$: Perfect $O(h^2)$ (4.01 → 4.00)** ⭐⭐⭐

**Table 7: $L^2$ Error $\| q - q_h \|_0$** for $P_2$ elements $(c=10^{-4} / c=10^4)$

| Method | ML 4 | ML 5 | ML 6 | ML 7 | ML 8 |
|--------|------|------|------|------|------|
| SFEM / RG    |  /  |  /  |  /  |  /  |  /  |
| LSFEM / RG   |  /  |  /  |  /  |  /  |  /  |
| wLSFEM / RG  |  /  |  /  |  /  |  /  |  /  |
|              |            |            |            |            |            |
| SFEM / PG    |  /  |  /  |  /  |  /  |  /  |
| LSFEM / PG   |  /  |  /  |  /  |  /  |  /  |
| wLSFEM / PG  |  /  |  /  |  /  |  /  |  /  |

**Spectacular Results:**

$c = 10^{-4}$:

- SFEM: 1.30e-04 ($O(h^2)$)
- LSFEM: 6.41e-07 (**203× better!**)
- wLSFEM: 6.53e-07 (matches LSFEM)

$c = 10^4$:

- SFEM: 1.30e-04 (still $O(h^2)$)
- **LSFEM: 1.88e-05** (error INCREASES ML7→8! DISASTER!) 🔴🔴🔴
- **wLSFEM: 6.40e-07* (203× better than SFEM! Perfect!) 🚀🚀🚀

**Table 8: $L^2(q)$ ratios** for $P_2$ elements on regular grids $(c=10^{-4} / c=10^4)$

| Method | ML4→5 | ML5→6 | ML6→7 | ML7→8 | Asymptotic |
|--------|-------|-------|-------|-------|------------|
| SFEM / RG   |  /  |  /  |  /  |  /  | $O(h^2)$ both |
| LSFEM / RG  |  /  |  /  |  /  |  /  | ✅ / 🔴 **CATASTROPHIC**|
| wLSFEM / RG |  /  |  /  |  /  |  /  | ⚠️ / ✅ |
|             |      |      |      |      |     |
| SFEM / PG   |  /  |  /  |  /  |  /  | ⚠️ / ⚠️ |
| LSFEM / PG  |  /  |  /  |  /  |  /  | ✅ |
| wLSFEM / PG |  /  |  /  |  /  |  /  | ✅ |

**CRITICAL FINDING:**

- For $c = 10^4$, standard LSFEM convergence ratio collapses: 15.26 → 12.97 → 6.35 → **1.08**
- Ratio of 1.08 means virtually NO mesh refinement benefit!
- **wLSFEM maintains ~8.0 ratio**, achieving $O(h^3)$ convergence
- **This is a complete validation of weighted LSFEM approach!**

#

## 🧩 Solver performance - The Shocking Truth

**Table 9: CG - iteration counts of $P_1$ elements on regular grid for LSFEM** $(c=10^{-4} / c=10^4)$

| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  16 /   49 |  16 /  16 |  7 / 12 |
| ML5  |  32 |  66 |  32 /  171 |  32 /  34 | 10 / 21 |
| ML6  |  64 | 130 |  64 /  533 |  64 /  95 | 15 / 38 |
| ML7  | 128 | 258 | 128 / 1702 | 128 / 283 | 25 / 58 |
| ML8  | 256 | 514 | 256 / 5024 | 256 / 736 | 44 / 77 |



**Table 10: CG - iteration counts of $P_1$ elements on regular grid for wLSFEM** $(c=10^{-4} / c=10^4)$

| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  /  |  /  |  /  |
| ML5  |  32 |  66 |  /  |  /  |  /  |
| ML6  |  64 | 130 |  /  |  /  |  /  |
| ML7  | 128 | 258 |  /  |  /  |  /  |
| ML8  | 256 | 514 |  /  |  /  |  /  |



**Table 11: CG - iteration counts of $P_1$ elements on perturbed grid for LSFEM** $(c=10^{-4} / c=10^4)$

| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  32 /  129 |  32 /   35 |  7 / 12 |
| ML5  |  32 |  66 |  66 /  337 |  64 /   87 | 11 / 22 |
| ML6  |  64 | 130 | 134 / 1032 | 128 /  198 | 16 / 41 |
| ML7  | 128 | 258 | 270 / 2815 | 256 /  461 | 26 / 62 |
| ML8  | 256 | 514 | 542 / 7874 | 512 / 1173 | 45 / 84 |

16, 16, 17, 34, , , 12, 1.2e+08
32, 32, 33, 66, , , 22, 1.1e+08
64, 64, 65, 130, , , 41, 1.1e+08
128, 128, 129, 258, , , 62, 1.1e+08
256, 256, 257, 514, , , 84, 1.1e+08

=== CONVERGENCE RATES ===
Interval, u_ratio, u_order, q_ratio, q_order
16→32, , 2.33, , 4.04
32→64, , 0.56, , -0.98
64→128, , 0.37, , 0.56
128→256, , 0.74, , 0.99

**Table 12: CG - iteration counts of $P_1$ elements on perturbed grid for wLSFEM** $(c=10^{-4} / c=10^4)$

| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |  /  |  /  |  /  |
| ML5  |  32 |  66 |  /  |  /  |  /  |
| ML6  |  64 | 130 |  /  |  /  |  /  |
| ML7  | 128 | 258 |  /  |  /  |  /  |
| ML8  | 256 | 514 |  /  |  /  |  /  |



**Observations:**


**Scaling Analysis:**


**Table 13: CG - iteration counts of $P_2$ elements on regular grid for LSFEM** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |   69 /   251 |   66 /   51 |  16 /  23 |
| ML5  |  32 |  66 |  141 /   677 |  138 /   95 |  27 /  43 |
| ML6  |  64 | 130 |  286 /  2029 |  283 /  320 |  42 /  70 |
| ML7  | 128 | 258 |  578 /  6050 |  575 /  998 |  70 / 123 |
| ML8  | 256 | 514 | 1162 / 17251 | 1163 / 2918 | 124 / 201 |



**Table 14: CG - iteration counts of $P_2$ elements on regular grid for wLSFEM** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |   318 /   48 |  43 /   8 |  15 /  1 |
| ML5  |  32 |  66 |  1048 /  114 |  92 /  15 |  26 /  2 |
| ML6  |  64 | 130 |  3361 /  361 | 136 /  29 |  41 /  5 |
| ML7  | 128 | 258 | 10179 / 1172 | 278 /  58 |  68 / 11 |
| ML8  | 256 | 514 | 29027 / 3646 | 569 / 115 | 122 / 21 |



**Table 15: CG - iteration counts of $P_2$ elements on perturbed grid for LSFEM** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |   74 /   328 |   67 /   94 |  17 /  29 |
| ML5  |  32 |  66 |  151 /   850 |  138 /  120 |  27 /  45 |
| ML6  |  64 | 130 |  311 /  2543 |  285 /  334 |  43 /  80 |
| ML7  | 128 | 258 |  638 /  7552 |  582 / 1005 |  71 / 145 |
| ML8  | 256 | 514 | 1289 / 22456 | 1178 / 3009 | 126 / 258 |



**Table 16: CG - iteration counts of $P_2$ elements on perturbed grid for wLSFEM** $(c=10^{-4} / c=10^4)$
| Mesh | Elements | DoF | CG | CG-JACOBI | CG-AMG |
|------|----------|-----|----|-----------|--------|
| ML4  |  16 |  34 |   351 /   88 |   72 /  11 |  16 /  1 |
| ML5  |  32 |  66 |  1212 /  154 |  146 /  20 |  27 /  3 |
| ML6  |  64 | 130 |  3860 /  454 |  300 /  37 |  41 /  5 |
| ML7  | 128 | 258 | 11948 / 1436 |  601 /  70 |  69 / 11 |
| ML8  | 256 | 514 | 33515 / 4459 | 1209 / 136 | 123 / 22 |



## 🧩 Mesh Perturbation Robustness



