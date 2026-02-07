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

and $\epsilon = 10^{-4}$ for the analytical solution

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

$${\cal J}(u,q;f) = \frac{1}{2} \parallel (-q^\prime - f) \parallel^2_0 + \frac{1}{2} \parallel (q - u^\prime) \parallel^2_0$$

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}(u_h,q_h;v_h,r_h) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}(u_h,q_h;v_h,r_h) = (-q^\prime_h , -r^\prime_h)_0 + (q_h - \nu u^\prime_h, r_h - v^\prime_h)_0$$

$${\cal F}(v_h,r_h) = (f, -r^\prime_h )_0$$


#

## 🎯 Accuracy and Solver Performance for $$\epsilon = 10^{-4}$$

**P1 Elements (Linear)**

**Figure 1** compares Standard FEM (SFEM) and Least-Squares FEM (LSFEM) using P1 elements:

- Solution accuracy ($L^2$ norm of u): Both methods achieve second-order convergence $(O(h^2))$, with SFEM being 1.9× more accurate at the finest mesh level (ML9).
- Flux accuracy ($L^2$ norm of q): LSFEM significantly outperforms SFEM, achieving second-order convergence $(O(h^2))$ compared to SFEM's first-order convergence $(O(h))$. LSFEM is 30× more accurate for flux computation.
- AMG preconditioning: LSFEM demonstrates excellent preconditioning efficiency with 92.5% iteration reduction using AMG, achieving a 13.4× speedup on perturbed grids (PG) compared to regular grids (RG).

![Figure 1: 1D Poisson - solution](figures/Fig_Poisson1D_P1.png)

**P2 Elements (Quadratic)**

**Figure 2** compares the same methods using P2 elements:

- Solution accuracy (($L^2$ norm of u): Both methods achieve nearly identical third-order convergence $(O(h^3))$, with comparable accuracy at ML9.
- Flux accuracy (($L^2$ norm of q): LSFEM again excels, maintaining third-order convergence $(O(h^3))$ while SFEM drops to second-order $(O(h^2))$. LSFEM is 45.9× more accurate for flux computation.
- AMG preconditioning: LSFEM shows 91.0% iteration reduction with AMG and an 11.1× speedup on perturbed grids.

![Figure 2: 1D Poisson - solution](figures/Fig_Poisson1D_P2.png)


## 💡 Key Takeaways

- LSFEM provides **superior** flux accuracy for both element types, maintaining one order higher convergence than SFEM
- Solution accuracy is comparable between methods
- LSFEM exhibits excellent AMG preconditioning properties, especially on distorted meshes
