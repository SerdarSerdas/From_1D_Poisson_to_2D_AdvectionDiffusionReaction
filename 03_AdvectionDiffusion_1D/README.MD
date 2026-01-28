## 1D Advection-Diffusion Equation

### Strong Form

$$
\begin{aligned}
-\nu u^{\prime\prime} + au^\prime &= f \quad &&\text{in } {\cal B} = (0,1), \\
u(0) = u(1) &= 0 \quad &&\text{on } \partial{\cal B}_D
\end{aligned}
$$

**Parameters:**
- $\nu = 10^{-4}$ (**moderate diffusion**)
- $a = 1$ (**advection velocity**)

## Standard Galerkin FEM (SFEM)

The classical Galerkin finite element formulation seeks $u_h \in V_h \subset H^1_0({\cal B})$ such that

$$(\nu u^\prime_h, v^\prime_h)_0 + (a u^\prime_h, v_h)_0 = (f, v_h)_0 \quad \forall v_h \in V_h$$


## Least Squares FEM 

### First-Order System - Diffusive Flux $(LSFEM^{(d)})$ :

**Flux:** $q = \nu u^\prime$  
**Balance:** $-q^\prime + a u^\prime = f$  
**Constitutive:** $q - \nu u^\prime = 0$

### First-Order System - Total Flux $(LSFEM^{(t)})$ :

**Flux:** $q = \nu u^\prime - a u$  
**Balance:** $-q^\prime = f$  
**Constitutive:** $q - \nu u^\prime + a u = 0$

### Functionals

$${\cal J}^{(d)}(u,q;f) = \frac{1}{2} \parallel w_1 (-q^\prime + a u^\prime - f) \parallel^2_0 + \frac{1}{2} \parallel w_2 (q - \nu u^\prime) \parallel^2_0$$

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


**Table 1: $L^2$ Error $\| u - u_h \|_0$ / Diffusive Flux** for regular and perturbed grids

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.7704e-02 | 2.4613e-03 | 4.2058e-04 | 1.0339e-04 | 2.5744e-05 | 6.4296e-06 |
| LSFEM / RG  | 1.6371e-02 | 3.9346e-03 | 9.9363e-04 | 2.4912e-04 | 6.2323e-05 | 1.5583e-05 |
| wLSFEM / RG | 1.6372e-02 | 3.9349e-03 | 9.9372e-04 | 2.4914e-04 | 6.2333e-05 | 1.5588e-05 |
|--------|------|------|------|------|------|-------|
| SFEM / PG   | 1.9123e-01 | 1.9898e-02 | 1.7512e-03 | 2.8413e-04 | 4.5773e-05 | 1.0628e-05 |
| LSFEM / PG  |  |  |  |  |  |  |
| wLSFEM / PG |  |  |  |  |  |  |



**Table 2: Convergence Rates for Solution $u$ / Diffusive Flux**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 4.28 | 2.55 | 2.02 | 2.01 | 2.00 | ✅ |
| LSFEM  | 2.06 | 1.99 | 2.00 | 2.00 | 2.00 | ✅ |
| wLSFEM | 2.06 | 1.99 | 2.00 | 2.00 | 2.00 | ✅ |








**Table 3: $L^2$ Error $\| u - u_h \|_0$ / Total Flux**

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   | 4.7704e-02 | 2.4613e-03 | 4.2058e-04 | 1.0339e-04 | 2.5744e-05 | 6.4296e-06 |
| LSFEM / RG  |  |  |  |  |  |  |
| wLSFEM / RG |  |  |  |  |  |  |
|--------|------|------|------|------|------|-------|
| SFEM / PG   |  |  |  |  |  |  |
| LSFEM / PG  |  |  |  |  |  |  |
| wLSFEM / PG |  |  |  |  |  |  |



**Table 4: Convergence Rates for Solution $u$ / Total Flux**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 4.28 | 2.55 | 2.02 | 2.01 | 2.00 | ✅ |
| LSFEM  |  |  |  |  |  | ❌ |
| wLSFEM |  |  |  |  |  | ✅ |


#

### 🌊 Flux variable $q$

**For the diffusion-dominated regime $c = 10^{-4}$, the behavior matches the pure Poisson case: SFEM delivers only first-order convergence for the flux $q = u^\prime$, while both LSFEM and wLSFEM reach optimal second-order convergence with very similar error levels. The weighting has only a marginal influence here because diffusion dominates the operator and the residual scaling is already well balanced (Table 4 - Table 6).**

**The reaction-dominated regime $c = 10^{4}$, reveals the main advantage of wLSFEM. SFEM remains locked at first-order flux convergence, and unweighted LSFEM not only fails to improve the flux but actually shows erratic and suboptimal convergence behavior for q. In contrast, wLSFEM achieves clean second-order convergence with flux errors several orders of magnitude smaller than both SFEM and LSFEM on the finest meshes (Table 5, Table 6). This demonstrates that the choice $w_1 = c^{-1/2}$ successfully rebalances the reaction and flux residuals in the least-squares functional.**

**Table 5: $L^2$ Error ($\| q - q_h \|_0$) / Diffusive Flux**

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM / RG   |  |  |  |  |  |  |
| LSFEM / RG  |  |  |  |  |  |  |
| wLSFEM / RG |  |  |  |  |  |  |
|                               |
| SFEM / PG   |  |  |  |  |  |  |
| LSFEM / PG  |  |  |  |  |  |  |
| wLSFEM / PG |  |  |  |  |  |  |



**Table 6: Convergence Rates for Solution $q$ / Diffusive Flux**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 1.88 | 1.08 | 1.01 | 1.00 | 1.00 | ✅ |
| LSFEM  | 1.96 | 2.06 | 2.02 | 2.00 | 2.00 | ✅ |
| wLSFEM | 1.96 | 2.06 | 2.02 | 2.00 | 1.99 | ✅ |


**Table 7: $L^2$ Error ($\| q - q_h \|_0$) / Total Flux**

| Method | ML 5 | ML 6 | ML 7 | ML 8 | ML 9 | ML 10 |
|--------|------|------|------|------|------|-------|
| SFEM   | 4.7705e-02 | 2.4627e-03 | 4.2253e-04 | 1.0534e-04 | 2.7652e-05 | 8.1737e-06 |
| LSFEM  |  |  |  |  |  |  |
| wLSFEM |  |  |  |  |  |  |





**Table 8: Convergence Rates for Solution $q$ / Total Flux**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| SFEM   | 4.28 | 2.54 | 2.00 | 1.93 | 1.76 | ❌ |
| LSFEM  |  |  |  |  |  | ❌ |
| wLSFEM |  |  |  |  |  | ✅ |



#

### ⚡ Solver performance (CG–AMG)

**The CG–AMG iteration counts show that unweighted LSFEM becomes increasingly difficult to solve as c grows. For $c=10^{-4}$, both LSFEM and wLSFEM exhibit mesh-dependent but still moderate iteration counts, with wLSFEM sometimes requiring fewer iterations on regular grids but more on strongly perturbed meshes (Table 7).**

**For $c=10^4$, the contrast becomes much stronger. LSFEM requires between 136 and 222 CG–AMG iterations at ML10, whereas wLSFEM converges in fewer than 25 iterations even on perturbed meshes (Table 8). The weighted formulation yields a much better conditioned system in the reaction-dominated regime, making it clearly preferable from a solver efficiency perspective.**

**Table 7: CG-AMG - iteration counts / Diffusive Flux** (Regular / Perturbed)
| Mesh | LSFEM | wLSFEM |
|------|-----------------|------------------|
| ML5  |  43 /   45 |   6 /   6 |
| ML6  |  67 /   77 |   7 /   9 |
| ML7  | 124 /  128 |  16 /  19 |
| ML8  | 226 /  257 |  48 /  52 |
| ML9  | 432 /  538 |  96 / 105 |
| ML10 | 859 / 1085 | 183 / 200 |




**Table 8: CG-AMG - iteration counts / Total Flux** (Regular / Perturbed)
| Mesh | LSFEM | wLSFEM |
|------|-----------------|------------------|
| ML5  |  /  |  /  |
| ML6  |  /  |  /  |
| ML7  |  /  |  /  |
| ML8  |  /  |  /  |
| ML9  |  /  |  /  |
| ML10 |  /  |  /  |

#

### 🏆 Summary

- In diffusion-dominated settings $(c = 10^{-4})$, all methods behave similarly for u, and both LSFEM and wLSFEM provide optimal flux convergence.

- In reaction-dominated settings $(c = 10^4)$, the unweighted LSFEM loses convergence for both u and q, while wLSFEM preserves SFEM-level accuracy in u and delivers superior, second-order accurate fluxes.

- From a solver viewpoint, wLSFEM also yields significantly lower CG–AMG iteration counts in the stiff reaction-dominated regime, indicating better conditioning and overall robustness.



