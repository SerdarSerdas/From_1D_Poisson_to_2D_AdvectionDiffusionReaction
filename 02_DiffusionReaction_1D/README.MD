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

**Table 1: `||u - uₕ||₀` Convergence (Boundary Layer)**

| Method | Mesh | `c=10⁻⁴` Reg | `c=10⁻⁴` Pert | `c=10⁴` Reg | `c=10⁴` Pert |
|--------|------|--------------|---------------|-------------|--------------|
| **SFEM** | ML5 | 1.50e-02 | 1.35e-02 | 1.08e-02 | 8.70e-03 |
| | ML6 | 3.91e-03 | 4.46e-03 | 2.60e-03 | 3.07e-03 |
| | ML7 | 9.87e-04 | 1.19e-03 | 6.35e-04 | 8.09e-04 |
| | ML8 | 2.47e-04 | 2.85e-04 | 1.58e-04 | 1.90e-04 |
| | ML9 | 6.19e-05 | 7.27e-05 | 3.94e-05 | 4.90e-05 |
| | ML10| 1.55e-05 | 1.82e-05 | 9.83e-06 | 1.25e-05 |
| **LSFEM**| ML5 | 7.47e-02 | 7.15e-02 | 3.09e-02 | 2.92e-02 |
| | ML6 | 1.90e-02 | 2.10e-02 | 1.58e-02 | 1.59e-02 |
| | ML7 | 4.76e-03 | 5.44e-03 | 9.04e-03 | 9.01e-03 |
| | ML8 | 1.19e-03 | 1.32e-03 | 5.04e-03 | 5.01e-03 |
| | ML9 | 2.98e-04 | 3.30e-04 | 2.56e-03 | 2.51e-03 |
| | ML10| 7.45e-05 | 8.22e-05 | 1.13e-03 | 1.11e-03 |
| **wLSFEM** | ML5 | 1.50e-02 | 1.35e-02 | 1.08e-02 | 8.70e-03 |
| | ML6 | 3.91e-03 | 4.46e-03 | 2.60e-03 | 3.07e-03 |
| | ML7 | 9.87e-04 | 1.19e-03 | 6.35e-04 | 8.09e-04 |
| | ML8 | 2.47e-04 | 2.85e-04 | 1.58e-04 | 1.90e-04 |
| | ML9 | 6.19e-05 | 7.27e-05 | 3.94e-05 | 4.90e-05 |
| | ML10| 1.55e-05 | 1.82e-05 | 9.83e-06 | 1.25e-05 |

> **🔍 Observation**: LSFEM errors **stop decreasing** at high reaction (`c=10⁴`), while **wLSFEM maintains optimal convergence**.

---

### 📈 Convergence Rates: Solution `u`

**Table 2: Convergence Rates for `u` (Regular Grids)**


