# 1D Poisson Equation: SFEM vs LSFEM

### Strong Form

$$
\begin{aligned}
-u^{\prime\prime} &= f \quad &&\text{in } {\cal B} = (0,1) \\
\end{aligned}
$$

with

$$
\begin{aligned}
u &= 0 \quad &&\text{on } \partial{\cal B}_D .
\end{aligned}
$$

and $\epsilon = 10^{-3}$ for the analytical solution

$$
u(x) = 4\left(\arctan\left(\frac{2(1/16-(x-1/2)^2)}{\pi\sqrt{\epsilon}}\right) + \frac{1}{2}\right)(1-x)x
$$

- Mesh levels: ML5 (32 elements) to ML9 (512 elements)
- Element types: $P_1$ (linear) and $P_2$ (quadratic)


## Standard Galerkin FEM (SFEM)

The classical Galerkin finite element formulation seeks $u \in V \subset H^1_0({\cal B})$ such that

$$(u^\prime_h, v^\prime_h)_0 = (f, v_h)_0 \quad \forall v_h \in V_h$$

## Least Squares FEM (LSFEM)

### First-Order System

Define the diffusive flux $q = u^\prime$:

$$-q^\prime = f, \qquad q - u^\prime = 0$$

### Functional

$${\cal J}(u,q;f) = \frac{1}{2} \parallel w_1 (-q^\prime - f) \parallel^2_0 + \frac{1}{2} \parallel w_2 (q - u^\prime) \parallel^2_0$$

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}(u_h,q_h;v_h,r_h) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}(u_h,q_h;v_h,r_h) = (-q^\prime_h , -r^\prime_h)_0 + (q_h - \nu u^\prime_h, r_h - v^\prime_h)_0$$

$${\cal F}(v_h,r_h) = (f, -r^\prime_h )_0$$


#

## 🎯 Accuracy and Convergence

### $P_1$ Elements

Both methods demonstrate excellent convergence behavior:

- **SFEM:** Achieves optimal convergence rates with ratio ≈ 4.0 for $\parallel u-u_h \parallel_0$ and ratio ≈ 2.0 for $\parallel q-q_h \parallel_0$
- **LSFEM:** Similarly achieves ratio ≈ 4.0 for $\parallel u-u_h \parallel_0$ and ratio ≈ 4.0 for $\parallel q-q_h \parallel_0$

LSFEM shows notably better accuracy for the flux variable q, with errors approximately 6-7 times smaller than SFEM at comparable mesh refinements.

![Figure 1: 1D Poisson - solution u](Fig_Poisson1D_Convergence.png)

### $P_2$ Elements

Both methods show higher-order convergence:

- **SFEM:** Generally achieves ratio ≈ 8.0 for $\parallel u-u_h \parallel_0$ and ratio ≈ 4.0 for $\parallel q-q_h \parallel_0$
- **LSFEM:** Achieves ratio ≈ 8.0 for both $\parallel u-u_h \parallel_0$ and $\parallel q-q_h \parallel_0$

**Important Note:** SFEM shows a convergence rate drop-off at ML9 (ratio drops to 3.77 for u and 2.78 for q), likely due to numerical precision effects. Interestingly, this occurs despite SFEM having a better condition number than LSFEM at ML9:

- SFEM: 5.7×10⁵
- LSFEM: 1.4×10⁶ (≈2.5× worse)

This suggests LSFEM's formulation is more robust to numerical precision issues at fine mesh levels, even with worse conditioning.

#

## 🎯 Mesh Perturbation

Both methods demonstrate robustness to mesh perturbations, maintaining similar convergence rates and error levels on perturbed grids compared to regular grids.

#

## 🎯 Computational Efficiency

### Iteration Counts

LSFEM requires significantly more iterations than would be needed for SFEM (not shown in tables), but responds very well to preconditioning:

**P1 Elements (ML9, Regular mesh):**

- CG: 512 iterations
- CG-JACOBI: 512 iterations
- CG-AMG: 79 iterations (85% reduction)

**P2 Elements (ML9, Regular mesh):**

- CG: 1163 iterations
- CG-JACOBI: 1166 iterations
- CG-AMG: 231 iterations (80% reduction)

Algebraic multigrid (AMG) preconditioning proves highly effective, reducing iteration counts by approximately 80-85% across all refinement levels. This is particularly important for LSFEM given its higher condition numbers.

![Figure 2: 1D Poisson - Solver](Fig_Poisson1D_Solver.png)

## 💡 Key Takeaways

- **Flux Accuracy:** LSFEM provides superior flux approximation compared to SFEM
- **Numerical Robustness:** LSFEM maintains optimal convergence rates even at fine meshes where SFEM begins to degrade, despite having worse conditioning
- **Mesh Perturbation:** Both methods handle mesh perturbations well with minimal impact on accuracy
- **Conditioning vs. Convergence:** Higher condition numbers don't necessarily imply worse convergence behavior—the problem formulation matters
- **Preconditioning:** AMG is essential for LSFEM efficiency, dramatically reducing iteration counts and mitigating conditioning issues

