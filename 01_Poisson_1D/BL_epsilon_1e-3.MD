## 📊 Numerical Results

### 🧩 Solution variable u

Both SFEM and wLSFEM demonstrate excellent accuracy and optimal second‑order convergence in the $L^2$ norm of u. In contrast, LSFEM produces consistently larger absolute errors across all refinement levels. SFEM yields the smallest u-errors overall, from $1.5 \times 10^{−2} at ML5 to $1.55 \times 10^{-5}$ at ML10, reflecting textbook second‑order convergence for $P_1$ elements.

LSFEM $(w_2 = 1)$ exhibits larger solution errors — roughly five times higher than SFEM at every mesh level — suggesting an imperfect balance between the PDE and constitutive residuals in the unweighted least‑squares functional.

wLSFEM, however, recovers SFEM‑level performance in u. Its $L^2$ norm-errors nearly coincide with SFEM values, indicating that the applied weighting successfully restores the balance lost in the standard LSFEM formulation.
