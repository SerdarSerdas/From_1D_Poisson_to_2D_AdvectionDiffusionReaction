## 📊 Numerical Results

![Figure 1: 1D Poisson method comparison](Fig_Poisson1D_BL_SolutionConverg.png)
**Figure 1.** (a) Solution $u_h$, (b) $u$ error convergence; (c) flux $q_h$, (d) $q$ error convergence. ML5 solutions (32 elements); log-log convergence (regular meshes). wLSFEM achieves SFEM-level $u$ accuracy with LSFEM-level flux recovery.

### 🧩 Solution variable $u$

Both **SFEM** and **wLSFEM** demonstrate excellent accuracy and optimal second-order convergence in the $L^2$ norm of $u$. In contrast, **LSFEM** produces consistently larger absolute errors across all refinement levels.

**SFEM** yields the smallest $u$-errors overall, from $1.5\times10^{-2}$ at ML5 to $1.55\times10^{-5}$ at ML10, reflecting textbook second-order convergence for $P_1$ elements.

**LSFEM** ($w_2 = 1$) exhibits larger solution errors — roughly **five times higher** than SFEM at every mesh level — suggesting an imperfect balance between the PDE and constitutive residuals in the unweighted least-squares functional.

**wLSFEM** ($w_2 = \epsilon^{-1/2}$), however, **recovers SFEM-level performance in $u$**. Its $L^2$-errors nearly coincide with SFEM values.

### 🌊 Flux variable $q = u'$

A more pronounced contrast appears in the flux errors **$\|q - q_h\|_0$**:

**SFEM** achieves only first-order convergence, consistent with derivative-recovery limitations. Flux errors remain large — **1.53** at ML5 and **$5.0\times10^{-2}$** at ML10.

**LSFEM** dramatically improves flux recovery, achieving second-order convergence and reducing errors by ~**five orders of magnitude** relative to SFEM at ML10.

**wLSFEM** further refines this: second-order convergence with **consistently smaller flux errors** than LSFEM, especially on perturbed meshes.



### 📊 Comparative summary

| Method | $\parallel u - u_h\parallel_0$ (solution) | $\parallel q - q_h\parallel_0$ (flux) | Convergence | Key characteristic |
|--------|----------------------------|-------------------------|-------------|-------------------|
| **SFEM**   | Excellent (smallest errors) | Poor (1st order) | $u:O(h^2), q:O(h)$ | Highly accurate $u$, weak flux recovery |
| **LSFEM**  | Moderate ($5\times$ larger than SFEM) | Excellent (2nd order) | $u,q:O(h^2)$ | Robust fluxes, suboptimal $u$ |
| **wLSFEM** | Excellent (≈ SFEM) | Very good (better than LSFEM) | $u,q:O(h^2)$ | Balanced performance and robustness |

### 🌊 Key observation
The log-log slopes confirm: **wLSFEM = best of both worlds** 💎

