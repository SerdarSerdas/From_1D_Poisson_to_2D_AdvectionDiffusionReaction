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

$${\cal J}(u,q;f) = \frac{1}{2} \parallel w_1 (-q^\prime + cu - f) \parallel^2_0 + \frac{1}{2} \parallel (q^\prime - u^\prime) \parallel^2_0$$

where $w_1$ is the weighting factor. We perform the computations with $w_1 = 1$ (LSFEM) and $w_1 = c^{-1/2}$ (wLSFEM).

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}((u_h,q_h);(v_h,r_h)) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}((u,q);(v,r)) = (q', r')_0 + (q - u', r - v')_0 + c^2(u, v)_0$$

$${\cal F}(v,r) = -(f, r')_0 + c(f, v)_0$$


## 4. Numerical Results

### Test Problem Configuration

We consider two configurations with $\epsilon = 10^{-3}$:
- **Boundary Layer**: Sharp gradients near domain boundaries
- **Interior Layer**: Sharp gradients in the domain interior

Both problems are discretized using linear ($P_1$) elements on regular and perturbed meshes.
