## 1D Advection-Diffusion Equation

### Strong Form

$$
\begin{aligned}
-\nu u^{\prime\prime} + au^\prime &= f \quad &&\text{in } {\cal B} = (0,1), \\
u(0) = u(1) &= 0 \quad &&\text{on } \partial{\cal B}_D
\end{aligned}
$$

where $\nu \geq 0$ is the diffusion coefficient with the following parameters

**Parameters:**
- $\nu = 10^{-4}$ (**moderate diffusion**)
- $a = 1$ (**advection velocity**)

Analytical solution with $\epsilon = 10^{-4}$:

$$
u(x) = 4\left(\arctan\left(\frac{2(1/16-(x-1/2)^2)}{\pi\sqrt{\epsilon}}\right) + \frac{1}{2}\right)(1-x)x
$$

**Test Configuration:**
- Mesh levels: ML5 (32 elements) to ML9 (512 elements)
- Element types: $P_1$ (linear) and $P_2$ (quadratic)
- Mesh types: Regular and perturbed grids

---
## Methods Compared

### 1. Standard Galerkin FEM (SFEM)

The classical Galerkin finite element formulation seeks $u_h \in V_h \subset H^1_0({\cal B})$ such that

$$(\nu u^\prime_h, v^\prime_h)_0 + (a u^\prime_h, v_h)_0 = (f, v_h)_0 \quad \forall v_h \in V_h$$

### 2. Least Squares FEM $(LSFEM^{(d)})$:

**Flux:** $q = \nu u^\prime$  
**Balance:** $-q^\prime + a u^\prime = f$  
**Constitutive:** $q - \nu u^\prime = 0$

### Functional

$${\cal J}^{(d)}(u,q;f) = \frac{1}{2} \parallel (-q^\prime + a u^\prime - f) \parallel^2_0 + \frac{1}{2} \parallel w_2 (q - \nu u^\prime) \parallel^2_0$$

**Weighting analyzed:**
- **$(LSFEM^{(d)})$**: $w_2 = 1$ 
- **$(wLSFEM^{(d)})$**: $w_2 = \nu^{-1/2}$ 








### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}((u_h,q_h);(v_h,r_h)) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}((u_h,q_h);(v_h,r_h)) = w^2_1 (-q^\prime_h + c u_h , -r^\prime_h + c v_h)_0 + (q_h - u^\prime_h, r_h - v^\prime_h)_0$$

$${\cal F}(v_h,r_h) = w^2_1 (f, -r^\prime_h + c v_h)_0$$

---




## Least Squares FEM 

### First-Order System - Diffusive Flux 

### First-Order System - Total Flux $(LSFEM^{(t)})$ :

**Flux:** $q = \nu u^\prime - a u$  
**Balance:** $-q^\prime = f$  
**Constitutive:** $q - \nu u^\prime + a u = 0$

### Functionals



$${\cal J}^{(t)}(u,q;f) = \frac{1}{2} \parallel w_1 (-q^\prime - f) \parallel^2_0 + \frac{1}{2} \parallel w_2 (q - \nu u^\prime + a u) \parallel^2_0$$


### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}(u_h,q_h;v_h,r_h) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$ .

**Diffusive flux:**

$${\cal A}^{(d)}(u_h,q_h);(v_h,r_h)) = w^2_1 (-q^\prime_h + a u^\prime_h , -r^\prime_h + a v^\prime_h)_0 + w^2_2 (q_h - \nu u^\prime_h, r_h - \nu v^\prime_h)_0$$

$${\cal F}^{(d)}(v_h,r_h) = w_1 (f, -r^\prime_h + a v^\prime_h)_0$$

**Total flux:**

$${\cal A}^{(t)}(u_h,q_h);(v_h,r_h)) = w^2_1 (-q^\prime_h, -r^\prime_h)_0 + w^2_2 (q_h - \nu u^\prime_h + a u^\prime_h, r_h - \nu v^\prime_h + a v^\prime_h)_0$$

$${\cal F}^{(t)}(v_h,r_h) = w_1 (f, -r^\prime_h)_0$$


## Weighting factor
**Weighting analyzed:**
- **Unweighted**: $w_1 = 1$, $w_2 = 1$ 
- **Weighted**: $w_1 = 1$, $w_2 = \nu^{-1/2}$ 



## 🔍 Key Findings

All configurations are discretized using linear ($P_1$) finite elements with $\nu = 1\cdot10^{-3}$ and $a = 1$ (**moderate diffusion**). Both **regular** and **perturbed** grids are studied.


#

### 🧩 Solution variable $u$

1. SFEM suffers from advection-induced instabilities:

ML5 error (regular): 4.77×10−24.77 \times 10^{-2}
4.77×10−2
ML5 error (perturbed): 1.91×10−11.91 \times 10^{-1}
1.91×10−1 (
4× worse)
Convergence rate 4.28 indicates severe solution pollution on coarse meshes
Requires 4+ mesh refinements to reach asymptotic regime
Classic Galerkin instability for advection-dominated problems
2. LSFEM provides built-in stabilization:

ML5 error: 1.64×10−21.64 \times 10^{-2}
1.64×10−2 (
3× better than SFEM on regular grids)
ML5 error (perturbed): 1.86×10−21.86 \times 10^{-2}
1.86×10−2 (
10× better than SFEM)
Optimal second-order convergence from coarsest mesh ✅
No ad-hoc stabilization (SUPG, etc.) needed
3. Flux formulation choice has no impact on solution accuracy:

Diffusive and total flux formulations produce identical uu
u errors (to 4+ decimal places)

Both LSFEM variants maintain stability and optimal convergence
4. Weighting has minimal effect on solution accuracy:

wLSFEM ≈ LSFEM for both formulations
Differences appear only at 5th decimal place

**Table 1: $L^2$ Error $\| u - u_h \|_0$ / Diffusive Flux** for regular and perturbed grids

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.7704e-02 | 2.4613e-03 | 4.2058e-04 | 1.0339e-04 | 2.5744e-05 | 6.4296e-06 |
| LSFEM / RG  | 1.6371e-02 | 3.9346e-03 | 9.9363e-04 | 2.4912e-04 | 6.2323e-05 | 1.5585e-05 |
| wLSFEM / RG | 1.6372e-02 | 3.9349e-03 | 9.9372e-04 | 2.4914e-04 | 6.2333e-05 | 1.5588e-05 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 1.9123e-01 | 1.9898e-02 | 1.7512e-03 | 2.8413e-04 | 4.5773e-05 | 1.0628e-05 |
| LSFEM / PG  | 1.8624e-02 | 4.0468e-03 | 1.1457e-03 | 2.8020e-04 | 7.0392e-05 | 1.7625e-05 |
| wLSFEM / PG | 1.8638e-02 | 4.0470e-03 | 1.1458e-03 | 2.8022e-04 | 7.0400e-05 | 1.7628e-05 |



**Table 2: Convergence Rates for Solution $u$ / Diffusive Flux** (regular grids)

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 4.28 | 2.55 | 2.02 | 2.01 | 2.00 | ✅ |
| LSFEM  | 2.06 | 1.99 | 2.00 | 2.00 | 2.00 | ✅ |
| wLSFEM | 2.06 | 1.99 | 2.00 | 2.00 | 2.00 | ✅ |


**Table 3: $L^2$ Error ($\| q - q_h \|_0$) / Diffusive Flux**

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 3.1521e-04 | 8.5861e-05 | 4.0570e-05 | 2.0209e-05 | 1.0096e-05 | 5.0468e-06 |
| LSFEM / RG  | 8.8749e-05 | 2.2778e-05 | 5.4635e-06 | 1.3493e-06 | 3.3634e-07 | 8.4409e-08 |
| wLSFEM / RG | 8.8765e-05 | 2.2786e-05 | 5.4654e-06 | 1.3503e-06 | 3.3673e-07 | 8.4162e-08 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 9.3350e-04 | 1.9935e-04 | 8.6637e-05 | 3.1559e-05 | 1.2299e-05 | 5.8724e-06 |
| LSFEM / PG  | 1.0585e-04 | 3.0626e-05 | 1.3298e-05 | 5.0175e-06 | 2.4546e-06 | 1.1992e-06 |
| wLSFEM / PG | 1.0011e-04 | 3.0630e-05 | 1.3301e-05 | 5.0184e-06 | 2.4548e-06 | 1.1992e-06 |

**Table 4: Convergence Rates for Solution $q$ / Diffusive Flux**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.88 | 1.08 | 1.01 | 1.00 | 1.00 | ✅ |
| LSFEM  | 1.96 | 2.06 | 2.02 | 2.00 | 2.00 | ✅ |
| wLSFEM | 1.96 | 2.06 | 2.02 | 2.00 | 1.99 | ✅ |





**Table 7: CG-AMG - iteration counts / Diffusive Flux** (Regular / Perturbed)
| Mesh | LSFEM | wLSFEM |
|------|-----------------|------------------|
| ML5  |  43 /   45 |   6 /   6 |
| ML6  |  67 /   74 |   7 /   9 |
| ML7  | 119 /  128 |  16 /  19 |
| ML8  | 223 /  264 |  48 /  52 |
| ML9  | 435 /  533 |  96 / 105 |
| ML10 | 860 / 1084 | 183 / 201 |












**Table 3: $L^2$ Error $\| u - u_h \|_0$ / Total Flux**

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.7704e-02 | 2.4613e-03 | 4.2058e-04 | 1.0339e-04 | 2.5744e-05 | 6.4296e-06 |
| LSFEM / RG  | 1.6371e-02 | 3.9346e-03 | 9.9363e-04 | 2.4912e-04 | 6.2324e-05 | 1.5577e-05 |
| wLSFEM / RG | 1.6372e-02 | 3.9349e-03 | 9.9372e-04 | 2.4914e-04 | 6.2333e-05 | 1.5588e-05 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 1.9123e-01 | 1.9898e-02 | 1.7512e-03 | 2.8413e-04 | 4.5773e-05 | 1.0628e-05 |
| LSFEM / PG  | 1.8624e-02 | 4.0468e-03 | 1.1457e-03 | 2.8020e-04 | 7.0392e-05 | 1.7625e-05 |
| wLSFEM / PG | 1.8638e-02 | 4.0470e-03 | 1.1458e-03 | 2.8022e-04 | 7.0400e-05 | 1.7628e-05 |




**Table 4: Convergence Rates for Solution $u$ / Total Flux**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 4.28 | 2.55 | 2.02 | 2.01 | 2.00 | ✅ |
| LSFEM  | 2.06 | 1.99 | 2.00 | 2.00 | 2.00 | ✅ |
| wLSFEM | 2.06 | 1.99 | 2.00 | 2.00 | 2.00 | ✅ |



**Table 7: $L^2$ Error ($\| q - q_h \|_0$) / Total Flux**

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.7705e-02 | 2.4627e-03 | 4.2253e-04 | 1.0534e-04 | 2.7652e-05 | 8.1737e-06 |
| LSFEM / RG  | 1.6371e-02 | 3.9345e-03 | 9.9362e-04 | 2.4911e-04 | 6.2323e-05 | 1.5577e-05 |
| wLSFEM / RG | 1.6372e-02 | 3.9349e-03 | 9.9371e-04 | 2.4914e-04 | 6.2332e-05 | 1.5588e-05 |
|             |            |            |            |            |            |            |
| SFEM / PG   | 1.9123e-01 | 1.9899e-02 | 1.7534e-03 | 2.8588e-04 | 4.7397e-05 | 1.2142e-05 |
| LSFEM / PG  | 1.8607e-02 | 4.0475e-03 | 1.1452e-03 | 2.8014e-04 | 7.0341e-05 | 1.7587e-05 |
| wLSFEM / PG | 1.8620e-02 | 4.0477e-03 | 1.1453e-03 | 2.8016e-04 | 7.0349e-05 | 1.7589e-05 |





**Table 8: Convergence Rates for Solution $q$ / Total Flux**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 4.28 | 2.54 | 2.00 | 1.93 | 1.76 | ❌ |
| LSFEM  | 2.06 | 1.99 | 2.00 | 2.00 | 2.00 | ✅ |
| wLSFEM | 2.06 | 1.99 | 2.00 | 2.00 | 2.00 | ✅ |






**Table 8: CG-AMG - iteration counts / Total Flux** (Regular / Perturbed)
| Mesh | LSFEM | wLSFEM |
|------|-----------------|------------------|
| ML5  | 5 / 5 | 15 / 15 |
| ML6  | 5 / 5 | 11 / 12 |
| ML7  | 5 / 5 |  8 /  8 |
| ML8  | 5 / 6 |  6 /  7 |
| ML9  | 5 / 6 |  5 /  6 |
| ML10 | 5 / 6 |  5 /  6 |


#

### 🌊 Flux variable $q$

**For the diffusion-dominated regime $c = 10^{-4}$, the behavior matches the pure Poisson case: SFEM delivers only first-order convergence for the flux $q = u^\prime$, while both LSFEM and wLSFEM reach optimal second-order convergence with very similar error levels. The weighting has only a marginal influence here because diffusion dominates the operator and the residual scaling is already well balanced (Table 4 - Table 6).**

**The reaction-dominated regime $c = 10^{4}$, reveals the main advantage of wLSFEM. SFEM remains locked at first-order flux convergence, and unweighted LSFEM not only fails to improve the flux but actually shows erratic and suboptimal convergence behavior for q. In contrast, wLSFEM achieves clean second-order convergence with flux errors several orders of magnitude smaller than both SFEM and LSFEM on the finest meshes (Table 5, Table 6). This demonstrates that the choice $w_1 = c^{-1/2}$ successfully rebalances the reaction and flux residuals in the least-squares functional.**







#

### ⚡ Solver performance (CG–AMG)

**The CG–AMG iteration counts show that unweighted LSFEM becomes increasingly difficult to solve as c grows. For $c=10^{-4}$, both LSFEM and wLSFEM exhibit mesh-dependent but still moderate iteration counts, with wLSFEM sometimes requiring fewer iterations on regular grids but more on strongly perturbed meshes (Table 7).**

**For $c=10^4$, the contrast becomes much stronger. LSFEM requires between 136 and 222 CG–AMG iterations at ML10, whereas wLSFEM converges in fewer than 25 iterations even on perturbed meshes (Table 8). The weighted formulation yields a much better conditioned system in the reaction-dominated regime, making it clearly preferable from a solver efficiency perspective.**





#

### 🏆 Summary

- In diffusion-dominated settings $(c = 10^{-4})$, all methods behave similarly for u, and both LSFEM and wLSFEM provide optimal flux convergence.

- In reaction-dominated settings $(c = 10^4)$, the unweighted LSFEM loses convergence for both u and q, while wLSFEM preserves SFEM-level accuracy in u and delivers superior, second-order accurate fluxes.

- From a solver viewpoint, wLSFEM also yields significantly lower CG–AMG iteration counts in the stiff reaction-dominated regime, indicating better conditioning and overall robustness.



