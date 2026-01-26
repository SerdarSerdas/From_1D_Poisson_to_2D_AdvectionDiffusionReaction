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



> **🔍 Observation**: LSFEM errors **stop decreasing** at high reaction (`c=10⁴`), while **wLSFEM maintains optimal convergence**.

---

### 📈 Convergence Rates: Solution `u`

**Table 2: Convergence Rates for `u` (Regular Grids)**


