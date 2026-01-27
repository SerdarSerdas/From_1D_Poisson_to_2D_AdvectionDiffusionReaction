## 📊 Numerical Results

### 📊 Comparative summary

| Method | $\parallel u − u_h \parallel_0$ (solution) | $\parallel q − q_h \parallel_0$ (flux) | Convergence                                  | Key characteristic                           |
| ------ | -------------------------------------- | ---------------------------------- | -------------------------------------------- | -------------------------------------------- |
| SFEM   | Excellent (smallest errors) | Poor (1st order) | $u:O(h^2), q:O(h)$ | Highly accurate u, weak flux recovery      |
| LSFEM  | Moderate ($5 \times$ larger than SFEM) | Excellent (2nd order) | $u,q:O(h^2)$ | Robust fluxes, suboptimal u |
| wLSFEM | Excellent (≈ SFEM) | Very good (better than LSFEM) | $u,q:O(h^2)$ | Balanced performance and improved robustness |

### 🧩 Solution variable u

Both SFEM and wLSFEM demonstrate excellent accuracy and optimal second‑order convergence in the $L^2$ norm of u. In contrast, LSFEM produces consistently larger absolute errors across all refinement levels. SFEM yields the smallest u-errors overall, from $1.5 \times 10^{−2}$ at ML5 to $1.55 \times 10^{-5}$ at ML10, reflecting textbook second‑order convergence for $P_1$ elements.

LSFEM $(w_2 = 1)$ exhibits larger solution errors — roughly five times higher than SFEM at every mesh level — suggesting an imperfect balance between the PDE and constitutive residuals in the unweighted least‑squares functional.

wLSFEM $(w_2 = \epsilon^{-1/2})$, however, recovers SFEM‑level performance in u. Its $L^2$ norm-errors nearly coincide with SFEM values, indicating that the applied weighting successfully restores the balance lost in the standard LSFEM formulation.




### 🌊 Flux variable $q = u^\prime$
 
A more pronounced contrast appears in the flux errors $\parallel q - q _h \parallel_0$:

SFEM achieves only first‑order convergence, consistent with the expected derivative‑recovery limitation of the classical Galerkin approach. Flux errors remain large — 1.53 at ML5 and $5.0 \times 10^{-2}$ at ML10 — showing weak gradient accuracy even on refined meshes.

LSFEM dramatically improves flux recovery, achieving second‑order convergence and reducing errors by approximately five orders of magnitude relative to SFEM at the finest mesh. This highlights the benefit of treating the flux as a primary unknown in the least‑squares framework.

wLSFEM further refines this behavior: it maintains second‑order convergence while yielding consistently smaller flux errors than unweighted LSFEM on both regular and perturbed meshes. The weighting enhances stability and provides the overall best compromise between solution and flux accuracy.
