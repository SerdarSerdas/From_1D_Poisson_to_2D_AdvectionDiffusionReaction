# 1D Poisson Equation: SFEM vs LSFEM

## 1. Problem Statement

**Strong form ($\nu = 1$):**

$$
\begin{aligned}
-u^{\prime\prime} &= f \quad &&\text{in } {\cal B}, \\
u &= 0 \quad &&\text{on } \partial{\cal B}_D
\end{aligned}
$$

## 2. Standard Galerkin FEM (SFEM)

The classical Galerkin finite element formulation seeks $u_h \in V_h \subset H^1_0({\cal B})$ such that

$$(u^\prime_h, v^\prime_h)_0 = (f, v_h)_0 \quad \forall v_h \in V_h$$

### Bilinear and Linear Forms

$$a(u_h, v_h) = (u^\prime_h, v^\prime_h)_0 = \int_{\cal B} u^\prime_h v^\prime_h dx$$

$$\ell(v_h) = (f, v_h)_0 = \int_{\cal B} f v_h dx$$

## 3. Least Squares FEM (LSFEM)

### First-Order System

Define the diffusive flux $q = u'$. The strong form is reformulated as:

$$-q^\prime  = f, \qquad q - u^\prime = 0$$

### Least-Squares Functional

$${\cal J}(u,q;f) = \|-q^\prime - f\|_0^2 + \|q - u'\|_0^2$$

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}((u_h,q_h);(v_h,r_h)) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}((u,q);(v,r)) = (-q^\prime, -r^\prime)_0 + (q - u', r - v')_0$$

$${\cal F}(v,r) = (f, r')_0$$


## 4. Numerical Results
