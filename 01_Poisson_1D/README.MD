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

## 🎯 Accuracy and Convergence

LSFEM shows notably better accuracy for the flux variable q, with errors approximately 6-7 times smaller than SFEM at comparable mesh refinements.

![Figure 1: 1D Poisson - solution u](Fig_Poisson1D_Convergence.png)

LSFEM provides superior flux accuracy: **47× more accurate** than SFEM for $P_2$ elements at ML9.

#

## 🎯 Solver Performance

![Figure 2: 1D Poisson - Solver](Fig_Poisson1D_Solver.png)

AMG preconditioning reduces LSFEM iterations by **80-93%** compared to unpreconditioned CG.

## 💡 Key Takeaways

LSFEM achieves higher-order flux convergence $(O(h^3) vs O(h^2))$ while maintaining excellent solver efficiency with AMG.

