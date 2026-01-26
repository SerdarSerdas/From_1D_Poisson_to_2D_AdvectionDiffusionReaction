## 1D Diffusion-Reaction Equation

### Strong Form

$$
\begin{aligned}
-u^{\prime\prime} + cu &= f \quad &&\text{in } {\cal B} = (0,1), \\
u(0) = u(1) &= 0 \quad &&\text{on } \partial{\cal B}_D
\end{aligned}
$$

where $c \geq 0$ is the reaction coefficient.

## Standard Galerkin FEM (SFEM)

The classical Galerkin finite element formulation seeks $u_h \in V_h \subset H^1_0({\cal B})$ such that

$$(u^\prime_h, v^\prime_h)_0 + c(u_h, v_h)_0 = (f, v_h)_0 \quad \forall v_h \in V_h$$

## Least Squares FEM (LSFEM)

### First-Order System

Define the diffusive flux $q = u^\prime$:

$$-q^\prime + cu = f, \qquad q - u^\prime = 0$$

### Functional

$${\cal J}(u,q;f) = \frac{1}{2} \parallel w_1 (-q^\prime + cu - f) \parallel^2_0 + \frac{1}{2} \parallel q - u^\prime \parallel^2_0$$

**Weighting analyzed:**
- **LSFEM**: $w_1 = 1$ (unweighted)
- **wLSFEM**: $w_1 = c^{-1/2}$ (balanced scaling)

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}((u_h,q_h);(v_h,r_h)) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}((u,q);(v,r)) = w^2_1 (-q^\prime + c u , -r^\prime + c v)_0 + (q - u^\prime, r - v^\prime)_0$$

$${\cal F}(v,r) = w^2_1 (f, -r^\prime + c v)_0$$


## 📊 Numerical Results

### 🔧 Test Problem Configuration

We consider the **boundary layer problem** and the **interior layer problem** with $\epsilon = 10^{-3}$ for two reaction regimes:

| Regime | Reaction coefficient |
|--------|---------------------|
| **Diffusion-dominated** | `c = 10^{-4}` (weak reaction) |
| **Reaction-dominated** | `c = 10^{4}` (strong reaction) |

**Both problems**: $P_1$ (linear) elements on **regular** and **perturbed** meshes.


---

### 🎯 Solution Accuracy: `$L^2$(u)` Errors


### Convergence Rate Analysis

The convergence rates are computed as the logarithmic ratio of successive $L^2$ errors between consecutive mesh refinement levels.


**Solution Variable $u$ — Boundary Layer**

Diffusion-dominated regime $(c=10^{-4})$: 
All three methods achieve optimal second-order convergence $({\cal O}(h^2))$, with rates asymptotically approaching 2.00 as expected for $P_1$ elements.

Reaction-dominated regime $(c = 10^{4})$:
A critical distinction emerges. SFEM and wLSFEM maintain optimal convergence with rates consistently near 2.00. However, standard LSFEM experiences severe degradation, with rates collapsing to 0.80–1.18. This demonstrates that the unweighted least-squares functional fails when reaction terms dominate. The weighted formulation (wLSFEM) successfully restores optimal convergence across all regimes.

**Table 3: Convergence Rates for Solution $u$**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| **Diffusion-dominated ($c=10^{-4}$)** |
| SFEM   | 1.94  | 1.98  | 2.00  | 2.00  | 2.00   | ✅ |
| LSFEM  | 1.98  | 1.99  | 2.00  | 2.00  | 2.00   | ✅ |
| wLSFEM | 1.94  | 1.99  | 2.00  | 2.00  | 2.00   | ✅ |
| **Reaction-dominated ($c=10^{4}$)** |
| SFEM   | 2.06  | 2.03  | 2.01  | 2.00  | 2.00   | ✅ |
| LSFEM  | 0.97  | 0.80  | 0.84  | 0.98  | 1.18   | ❌ (Degraded) |
| wLSFEM | 2.06  | 2.03  | 2.01  | 2.00  | 2.00   | ✅ |

> **🔍 Observation**: LSFEM errors **stop decreasing** at high reaction (`c=10⁴`), while **wLSFEM maintains optimal convergence**.

**Flux Variable $q$ — Boundary Layer** 

Diffusion-dominated regime $(c = 10^{-4})$:

SFEM exhibits standard first-order convergence $(\cal{O}(h), rate \approx 1.00)$ for flux recovery. Both LSFEM and wLSFEM demonstrate superconvergence, achieving second-order rates $(\approx 2.00)$ — a full order of magnitude improvement over SFEM.

Reaction-dominated regime $(c = 10^{4})$:
SFEM maintains first-order convergence. Standard LSFEM shows irregular and deteriorating behavior (rates 1.12–1.51), losing its superconvergent property. wLSFEM preserves second-order superconvergence $(rate \approx 2.00)$, demonstrating robust flux accuracy independent of the reaction coefficient magnitude. 

*Key Finding:* Weighted LSFEM is essential for maintaining optimal convergence and flux superconvergence in reaction-dominated problems, while standard LSFEM is only suitable for diffusion-dominated regimes.

**Table 4: Convergence Rates for Flux $q$**

| Method | ML5→6 | ML6→7 | ML7→8 | ML8→9 | ML9→10 | Asymptotic |
|--------|-------|-------|-------|-------|--------|------------|
| **Diffusion-dominated ($c=10^{-4}$)** |
| SFEM   | 0.95  | 0.99  | 1.00  | 1.00  | 1.00   | ✅ |
| LSFEM  | 1.94  | 1.98  | 2.00  | 2.00  | 2.00   | ✅ |
| wLSFEM | 1.94  | 1.98  | 2.00  | 2.00  | 2.00   | ✅ |
| **Reaction-dominated ($c=10^{4}$)** |
| SFEM   | 0.98  | 1.00  | 1.00  | 1.00  | 1.00   | ✅ |
| LSFEM  | 1.19  | 1.12  | 1.20  | 1.33  | 1.51   | ❌ (Degraded) |
| wLSFEM | 1.97  | 2.00  | 2.00  | 2.00  | 2.00   | ✅ |



---

### 📈 Convergence Rates: Solution `u`

**Table 2: Convergence Rates for `u` (Regular Grids)**


Solver Performance
The symmetric positive definite structure of LSFEM systems enables efficient solution via the conjugate gradient (CG) method with algebraic multigrid (AMG) preconditioning. Table 5 presents iteration counts for LSFEM and wLSFEM across different preconditioners.
Diffusion-dominated regime (c=10−4c = 10^{-4}
c=10−4):

Without preconditioning, both methods exhibit iteration counts that scale linearly with problem size (CG: 32→1024 for LSFEM, 153→32402 for wLSFEM from ML5→ML10). AMG preconditioning dramatically improves efficiency: LSFEM requires 10→147 iterations, while wLSFEM achieves near mesh-independent convergence with only 5→25 iterations—representing up to a 1000× speedup over unpreconditioned CG at the finest mesh level.
Reaction-dominated regime (c=104c = 10^{4}
c=104):

The performance characteristics reverse. Standard LSFEM becomes poorly conditioned, requiring 173→24864 unpreconditioned iterations. In contrast, wLSFEM demonstrates superior conditioning with only 39→10801 iterations—already 2× better without preconditioning. With AMG, wLSFEM achieves remarkable efficiency: merely 1→18 iterations across all mesh levels. LSFEM-AMG requires 21→150 iterations, showing wLSFEM is 8× more efficient at the finest level.
Conclusion: The weighted formulation not only restores convergence accuracy but also significantly improves system conditioning, making wLSFEM the preferred choice for general-purpose applications, particularly in reaction-dominated regimes.
