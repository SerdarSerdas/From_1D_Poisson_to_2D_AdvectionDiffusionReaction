# 1D Poisson Equation: Mathematical Formulation

> **Complete derivation of SFEM and LSFEM formulations** for the 1D Poisson equation with Dirichlet boundary conditions. This document provides the mathematical foundation for all numerical results presented in [README.md](README.md).

---

## Strong Form

The one-dimensional Poisson equation with homogeneous Dirichlet boundary conditions is:

**Find** $u: \mathcal{B} = (0,1) \to \mathbb{R}$ **such that:**

$$-u'' = f \quad \text{in } \mathcal{B} = (0,1)$$

$$u = 0 \quad \text{on } \partial\mathcal{B}_D = \{0, 1\}$$

where:
- $u(x)$ is the unknown scalar field (temperature, displacement, potential, etc.)
- $f(x)$ is the prescribed source term
- Prime notation: $u' = \frac{du}{dx}$, $u'' = \frac{d^2u}{dx^2}$

### Physical Interpretation

The Poisson equation models steady-state diffusion processes:
- **Heat conduction:** $u$ = temperature, $f$ = heat source density
- **Electrostatics:** $u$ = electric potential, $f$ = charge density
- **Elasticity:** $u$ = displacement, $f$ = body force

---

## Test Problems

### 1. Boundary Layer Configuration

**Governing equation:**

$$-u'' = f, \quad u(0) = 0, \quad u(1) = 0$$

**Analytical solution:**

$$u(x) = 1 - \frac{\sinh(x/\sqrt{\epsilon}) + \sinh((1-x)/\sqrt{\epsilon})}{\sinh(1/\sqrt{\epsilon})}$$

**Layer structure:**
- Boundary layer at $x=0$ and $x=1$
- $f = -u''$ computed from analytical solution

**Perturbation parameter:** ε = 10⁻³, 5×10⁻⁶ (moderate to extreme)

### 2. Interior Layer Configuration

**Governing equation:**

$$-u'' = f(x)$$

**Analytical solution:**

$$u(x) = 4\left(\arctan\left(\frac{2(1/16-(x-1/2)^2)}{\pi\sqrt{\epsilon}}\right) + \frac{1}{2}\right)(1-x)x$$

**Layer structure:**
- Interior shock at $x = 0.25$ and $x=0.75$ 
- Source term: $f = -u''$ computed from analytical solution

**Perturbation parameter:** ε = 10⁻³, 10⁻⁷ (moderate to extreme)

---

## Standard Galerkin FEM (SFEM)

### Weak Formulation

Multiply the strong form by a test function $v \in H^1_0(\mathcal{B})$ and integrate by parts:

$$\int_\mathcal{B} u' v' dv = \int_\mathcal{B} f v dv$$

**Variational problem:** Find $u \in H^1_0(\mathcal{B})$ such that

$$(u', v') = (f, v) \quad \forall v \in H^1_0(\mathcal{B})$$

where the inner product is defined as:
$$(w, v) = \int_\mathcal{B} w(x) v(x) dv$$

### Discrete Formulation

Introduce a finite-dimensional subspace $V_h \subset H^1_0(\mathcal{B})$ spanned by piecewise polynomial basis functions $\{\phi_i(x)\}_{i=1}^N$:

$$u_h(x) = \sum_{i=1}^N u_i \phi_i(x)$$

**Galerkin approximation:** Find $u_h \in V_h$ such that

$$(u'_h, v'_h) = (f, v_h) \quad \forall v_h \in V_h$$

Substituting the discrete representation and testing with each basis function yields the **linear system:**

```
Ku = f
```

where:
- **Stiffness matrix:** $K_{ij} = (\phi'_j, \phi'_i) = \int_0^1 \phi'_j(x) \phi'_i(x) \, dx$
- **Load vector:** $f_i = (f, \phi_i) = \int_0^1 f(x) \phi_i(x) \, dx$
- **Solution vector:** $\mathbf{u} = [u_1, u_2, \ldots, u_N]^T$

### Properties

**Matrix structure:**
- **Symmetric:** $K_{ij} = K_{ji}$ (follows from inner product symmetry)
- **Positive definite:** $\mathbf{u}^T K \mathbf{u} = \|u'_h\|^2_0 > 0$ for $\mathbf{u} \neq \mathbf{0}$
- **Sparse:** Tridiagonal for 1D uniform mesh with P₁ elements
- **Band structure:** Bandwidth = $2p+1$ for P_p elements

**Flux recovery:**
The flux $q = -u'$ is obtained through **post-processing:**

$$q_h(x) = -u'_h(x) = -\sum_{i=1}^N u_i \phi'_i(x)$$

This derivative is:
- **Discontinuous** across element boundaries (for C⁰ elements)
- **Lower accuracy:** Achieves only O(h^p) convergence for P_p elements
- **Not variationally consistent:** No residual orthogonality property

---

## Least-Squares FEM (LSFEM)

### First-Order System Reformulation

Introduce the **flux as an independent variable:**

$$q = -u'$$

The second-order Poisson equation becomes a **first-order system:**

$$q + u' = 0$$

$$-q' = f$$

Equivalently:
$$-q' - f = 0, \quad q - u' = 0$$

### Least-Squares Functional

Define the **L² least-squares functional:**

$$\mathcal{J}(u, q; f) = \frac{1}{2}\|-q' - f\|_0^2 + \frac{1}{2}\|q - u'\|_0^2$$

Expanding the norms:

$$\mathcal{J}(u, q; f) = \frac{1}{2}\int_0^1 (-q' - f)^2 \, dx + \frac{1}{2}\int_0^1 (q - u')^2 \, dx$$

**Physical interpretation:**
- First term: Residual of equilibrium equation ($-q' = f$)
- Second term: Residual of constitutive equation ($q = -u'$)
- Minimization: Seeks best simultaneous satisfaction of both equations in L² sense

### Variational Formulation

The **continuous minimization problem** is:

**Find** $(u, q) \in H^1_0(\mathcal{B}) \times L^2(\mathcal{B})$ **such that:**

$$\mathcal{J}(u, q; f) = \min_{(v, r) \in H^1_0 \times L^2} \mathcal{J}(v, r; f)$$

Taking the first variation (Euler-Lagrange equations):

$$\mathcal{A}((u, q); (v, r)) = \mathcal{F}(v, r) \quad \forall (v, r) \in H^1_0 \times L^2$$

where the **bilinear form** is:

$$\mathcal{A}((u, q); (v, r)) = (-q', -r') + (q - u', r - v')$$

and the **linear functional** is:

$$\mathcal{F}(v, r) = (f, -r')$$

Expanding:

$$\mathcal{A}((u, q); (v, r)) = \int_0^1 q' r' \, dx + \int_0^1 (q - u')(r - v') \, dx$$

$$\mathcal{F}(v, r) = -\int_0^1 f r' \, dx$$

### Discrete Formulation

Introduce finite element spaces:
- $V_h \subset H^1_0$ for the scalar field $u$
- $Q_h \subset L^2$ for the flux field $q$

For **equal-order approximation** (P₁-P₁ or P₂-P₂):

$$u_h(x) = \sum_{i=1}^N u_i \phi_i(x), \quad q_h(x) = \sum_{i=1}^N q_i \phi_i(x)$$

**Discrete minimization:** Find $(u_h, q_h) \in V_h \times Q_h$ such that

$$\mathcal{A}((u_h, q_h); (v_h, r_h)) = \mathcal{F}(v_h, r_h) \quad \forall (v_h, r_h) \in V_h \times Q_h$$

Testing with each basis function yields the **block linear system:**

```
┌              ┐ ┌   ┐   ┌      ┐
│  M    -K^T   │ │ q │   │  Mf  │
│              │ │   │ = │      │
│ -K     M     │ │ u │   │  0   │
└              ┘ └   ┘   └      ┘
```

where:
- **Mass matrix:** $M_{ij} = (\phi_i, \phi_j) = \int_0^1 \phi_i(x) \phi_j(x) \, dx$
- **Stiffness matrix:** $K_{ij} = (\phi'_i, \phi'_j) = \int_0^1 \phi'_i(x) \phi'_j(x) \, dx$
- **Load vector:** $(Mf)_i = \int_0^1 f(x) \phi_i(x) \, dx$
- **Solution vectors:** $\mathbf{q} = [q_1, \ldots, q_N]^T$, $\mathbf{u} = [u_1, \ldots, u_N]^T$

### Properties

**System matrix structure:**
- **Symmetric:** Block structure with $M^T = M$, $K^T = K$
- **Positive definite:** $\mathcal{A}((u,q); (u,q)) = \|q'\|^2_0 + \|q-u'\|^2_0 > 0$
- **Doubled dimension:** 2N unknowns vs N for SFEM
- **No saddle-point:** Unlike mixed methods, no zero diagonal blocks

**Advantages over SFEM:**
1. **Flux is primary variable:** $q_h$ solved variationally, not post-processed
2. **Superconvergence:** $\|q - q_h\|_0 = O(h^{p+1})$ vs SFEM's $O(h^p)$
3. **Always SPD:** No loss of positive definiteness for any problem parameters
4. **Natural error estimator:** $\eta^2_K = \|q'_h + f\|^2_{0,K} + \|q_h - u'_h\|^2_{0,K}$

**Trade-offs:**
- 2× DOFs (but CG-AMG keeps iterations bounded ~150)
- More complex implementation (block system, flux boundary conditions)

---

## Finite Element Spaces

### P₁ Elements (Linear)

**Basis functions on reference element** $\hat{K} = [0, h]$:

$$\hat{\phi}_1(\xi) = 1 - \frac{\xi}{h}, \quad \hat{\phi}_2(\xi) = \frac{\xi}{h}$$

**Properties:**
- C⁰ continuity (continuous, but derivatives discontinuous)
- 2 DOFs per element
- Piecewise linear approximation

**Convergence:**
- Solution: $\|u - u_h\|_0 = O(h^2)$, $\|u - u_h\|_1 = O(h)$
- SFEM flux: $\|u' - u'_h\|_0 = O(h)$
- LSFEM flux: $\|q - q_h\|_0 = O(h^2)$ ✨ (superconvergence!)

### P₂ Elements (Quadratic)

**Basis functions on reference element** $\hat{K} = [0, h]$ with midpoint node at $h/2$:

$$\hat{\phi}_1(\xi) = \frac{(h-\xi)(h-2\xi)}{h^2}, \quad \hat{\phi}_2(\xi) = \frac{4\xi(h-\xi)}{h^2}, \quad \hat{\phi}_3(\xi) = \frac{\xi(2\xi-h)}{h^2}$$

**Properties:**
- C⁰ continuity
- 3 DOFs per element (2 vertices + 1 midpoint)
- Piecewise quadratic approximation

**Convergence:**
- Solution: $\|u - u_h\|_0 = O(h^3)$, $\|u - u_h\|_1 = O(h^2)$
- SFEM flux: $\|u' - u'_h\|_0 = O(h^2)$
- LSFEM flux: $\|q - q_h\|_0 = O(h^3)$ ✨ (superconvergence!)

**Advantage for thin layers:**
- Extended asymptotic regime: P₂ achieves cubic convergence where P₁ shows pre-asymptotic chaos
- Better resolution of sharp gradients with fewer DOFs

---

## Boundary Conditions

### Homogeneous Dirichlet (Essential)

**Strong form:** $u(0) = 0$, $u(1) = 0$

**Implementation:**
- **SFEM:** Eliminate boundary DOFs from system (nodes at $x=0$ and $x=1$)
- **LSFEM:** Eliminate boundary DOFs for $u$ only; flux $q$ has natural boundary conditions

**Function space:** $H^1_0(\mathcal{B}) = \{v \in H^1 : v(0) = v(1) = 0\}$

### Inhomogeneous Dirichlet

**Strong form:** $u(0) = g_0$, $u(1) = g_1$

**Implementation:**
- **Lift function:** $u = u_0 + \tilde{u}$ where $u_0$ satisfies BCs, $\tilde{u} \in H^1_0$
- **Example lift:** $u_0(x) = g_0(1-x) + g_1 x$ (linear interpolation)
- Solve for $\tilde{u}$ with homogeneous BCs and modified RHS

### Neumann (Natural)

**Strong form:** $-u'(0) = g_N^0$, $-u'(1) = g_N^1$

**SFEM weak form:**
$$(u', v') = (f, v) + g_N^1 v(1) - g_N^0 v(0)$$

**LSFEM treatment:**
- Flux $q$ has essential BC: $q(0) = g_N^0$, $q(1) = g_N^1$
- More natural in first-order formulation

---

## Error Norms

**L² norm (solution accuracy):**
$$\|u - u_h\|_0 = \left(\int_0^1 (u - u_h)^2 \, dx\right)^{1/2}$$

**H¹ seminorm (derivative accuracy):**
$$|u - u_h|_1 = \left(\int_0^1 (u' - u'_h)^2 \, dx\right)^{1/2}$$

**Flux L² norm (LSFEM primary metric):**
$$\|q - q_h\|_0 = \left(\int_0^1 (q - q_h)^2 \, dx\right)^{1/2}$$

**Numerical integration:**
All integrals computed using Gaussian quadrature with $p+1$ points for P_p elements (exact integration of polynomial integrands).

---

## Convergence Rate Calculation

**Empirical convergence rate** between consecutive mesh levels:

$$r = \log_2\left(\frac{E_{\text{ML}(k)}}{E_{\text{ML}(k+1)}}\right)$$

where $E$ is the error norm ($\|u-u_h\|_0$ or $\|q-q_h\|_0$).

**Interpretation:**
- $r \approx 2$: O(h²) convergence (optimal for P₁)
- $r \approx 3$: O(h³) convergence (optimal for P₂)
- $r < 0$: Error growth (pre-asymptotic regime!)
- Erratic $r$: Pre-asymptotic chaos, mesh too coarse

**Asymptotic regime criteria:**
- Convergence rates stabilize within ±0.05 over 3+ consecutive refinements
- Typically requires $h/\lambda < 0.1$ for singular perturbations

---

## Implementation Notes

**Assembly:**
- Element-by-element loop over mesh
- Map to reference element, apply quadrature
- Add contributions to global matrices

**SFEM system size:** N × N (N = number of interior nodes)

**LSFEM system size:** 2N × 2N (block structure)

**Memory:** SFEM < LSFEM, but both sparse (nnz ~ O(N))

**Solver choice:**
- SFEM: CG with incomplete Cholesky often sufficient
- LSFEM: CG-AMG mandatory for production (7-14× speedup)

---

## References

- **Brenner & Scott (2008).** *The Mathematical Theory of Finite Element Methods.* [DOI:10.1007/978-0-387-75934-0](https://doi.org/10.1007/978-0-387-75934-0)
- **Bochev & Gunzburger (2009).** *Least-Squares Finite Element Methods.* [DOI:10.1007/b13382](https://doi.org/10.1007/b13382)

---

**See also:**
- [Main Results](README.md) — Comprehensive numerical analysis
- [Theory Documentation](../01_Theory/README.md) — General ADR framework

**Last updated:** January 2025
