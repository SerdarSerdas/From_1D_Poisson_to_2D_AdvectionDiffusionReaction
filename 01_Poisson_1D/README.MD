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

$${\cal A}((u_h,q_h);(v_h,r_h)) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}((u_h,q_h);(v_h,r_h)) = w^2_1(-q^\prime_h , -r^\prime_h)_0 + w^2_2 (q_h - \nu u^\prime_h, r_h - v^\prime_h)_0$$

$${\cal F}(v_h,r_h) = w^2_1 (f, -r^\prime_h )_0$$


## 🔍 Summary

This repository presents a comprehensive comparison of **Standard Galerkin FEM (SFEM)** and **Least-Squares FEM (LSFEM)** for solving the 1D Poisson equation. Two challenging test cases with sharp solution features are analyzed to validate convergence rates and demonstrate LSFEM's superior flux accuracy.

### 📊 Test Cases

**1. Boundary Layer Problem**

$$
u(x) = x - (\exp((x - 1) / \epsilon) - \exp(-1 / \epsilon)) / (1 - \exp(-1 / \epsilon))
$$

with $\epsilon = 10^{-1}$.

- Mesh levels: ML4 (16 elements) to ML8 (256 elements)
- Element types: $P_1$ (linear) and $P_2$ (quadratic)


**2. Interior Layer Problem**

$$
u(x) = 4\left(\arctan\left(\frac{2(1/16-(x-1/2)^2)}{\pi\sqrt{\epsilon}}\right) + \frac{1}{2}\right)(1-x)x
$$

#

### 🎯 Key Results

**Convergence Rates (Theory vs Practice)**

### Summary of Convergence Rates

| Method    | Variable | Theoretical | Observed Ratio    | Status              |
|-----------|----------|-------------|--------------|---------------------|
| SFEM $P_1$   | Solution u | $O(h^2)$     | ✅ 4.00     | Perfect            |
| SFEM $P_1$   | Flux q     | $O(h)$      | ✅ 2.00     | Perfect            |
| LSFEM $P_1$  | Solution u | $O(h^2)     | ✅ 4.00     | Perfect            |
| **LSFEM $P_1$** | **Flux q** | **$O(h^2)$** | **✅ 4.00** | **Perfect ⭐** |
| SFEM $P_2$   | Solution u | $O(h^3)$     | ⚠️ 7.31-7.69| Degrades          |
| SFEM $P_2$   | Flux q     | $O(h^2)$     | ⚠️ 3.91-3.95| Degrades          |
| **LSFEM $P_2$** | **Solution u** | **$O(h^3)$** | **✅ 8.00-8.01** | **Perfect ⭐** |
| **LSFEM $P_2$** | **Flux q** | **$O(h^3)$** | **✅ 7.99-8.00** | **Perfect ⭐⭐** |

#

**Key Insight**: LSFEM achieves one order higher convergence for flux (derivative) quantities!

### Flux Accuracy Comparison (at Finest Mesh Level)

| Problem Type | Element | SFEM Error | LSFEM Error | LSFEM Advantage |
|--------------|---------|------------|-------------|-----------------|
| Boundary Layer | $P_1$ | 2.52e-02 | 2.82e-04 | **89× more accurate** 🚀    |
| Boundary Layer | $P_2$ | 1.30e-04 | 6.41e-07 | **203× more accurate** 🚀🚀 |
| Interior Layer | $P_1$ | 5.10e-02 | 9.92e-04 | **51× more accurate** 🚀    |
| Interior Layer | $P_2$ | 4.08e-04 | 5.72e-06 | **71× more accurate** 🚀🚀 |



**Takeaway:** LSFEM provides **50-200× better flux accuracy** depending on element order and problem type.

#

### Solution Accuracy Comparison

| Problem Type | Element | SFEM Error | LSFEM Error | LSFEM Advantage |
|--------------|---------|------------|-------------|-----------------|
| Boundary Layer | $P_1$ | 3.11e-05 | 5.40e-05 | SFEM 1.74× better    |
| Boundary Layer | $P_2$ | 7.01e-08 | 6.41e-08 | **LSFEM 1.09× better** ✅ |
| Interior Layer | $P_1$ | 6.30e-05 | 1.17e-04 | SFEM 1.86× better    |
| Interior Layer | $P_2$ | 2.12e-07 | 2.04e-07 | **LSFEM 1.04× better** ✅ |

**Takeaway:** For $P_1$, SFEM slightly more accurate for solution. For $P_2$, **LSFEM equals or exceeds SFEM** for solution while maintaining massive flux advantage!

#

## Robustness: SFEM vs LSFEM ($P_2$ Elements)

**Convergence Rate Stability at Finest Levels:**

| Problem Type | Variable | SFEM Ratio (ML7→8) | LSFEM Ratio (ML7→8) | Winner |
|--------------|---------|------------|-------------|-----------------|
| Boundary Layer | $u$ | 7.31 🔴 | 8.00 ✅ | **LSFEM** |
| Boundary Layer | $q$ | 3.91 🔴 | 8.00 ✅ | **LSFEM** |
| Interior Layer | $q$ | 7.69 🔴 | 8.01 ✅ | **LSFEM** |
| Interior Layer | $q$ | 3.95 🔴 | 7.99 ✅ | **LSFEM** |

**Observation:** SFEM shows degradation at finest mesh levels, while **LSFEM maintains perfect convergence rates** (target: 8.00 for $O(h^3)$).


### ⚡ Computational Efficiency

**Solver Performance (Iteration Counts at ML8)**

**$P_1$ Elements:**

| Problem Type | Preconditioner | Iterations | Speedup vs. CG|
|--------------|---------|------------|-------------|
| Boundary Layer | None (CG) | 513 | 1.0×   |
| Boundary Layer | JACOBI    | 511 | 1.0× (negligible) |
| Boundary Layer | AMG       |  45 | **11.4×** 🚀 |
| Interior Layer | None (CG) | 256 | 1.0×   |
| Interior Layer | JACOBI    | 256 | 1.0× (zero improvement) |
| Interior Layer | AMG       |  44 | **5.8×** 🚀 |

**$P_2$ Elements:**

| Problem Type | Preconditioner | Iterations | Speedup vs. CG |
|--------------|---------|------------|-------------|
| Boundary Layer | None (CG) | 1162 | 1.0×          |
| Boundary Layer | JACOBI    | 1163 | 0.99× (worse) |
| Boundary Layer | AMG       |  45 | **9.4×** 🚀    |
| Interior Layer | None (CG) | 578 | 1.0×           |
| Interior Layer | JACOBI    | 577 | 1.0× (useless) |
| Interior Layer | AMG       | 123 | **4.7×** 🚀   |

**Critical Finding:**

- ❌ **Jacobi preconditioning is ineffective** (or even counterproductive)
- ✅ **AMG preconditioning is essential** for practical LSFEM computations
- 🚀 **5-11× iteration reduction** with AMG

#
