# 1D Poisson Equation

---

## Strong Form

The one-dimensional Poisson equation with homogeneous Dirichlet boundary conditions is:

**Find** $u: \mathcal{B} = (0,1) \to \mathbb{R}$ **such that:**

$$-u'' = f \quad \text{in } \mathcal{B} = (0,1)$$

$$u = 0 \quad \text{on } \partial\mathcal{B}_D = \{0, 1\}$$

where:
- $u(x)$ is the unknown scalar field (temperature, displacement, potential, etc.)
- $f(x)$ is the prescribed source term
- Prime notation: $u' = \frac{du}{dx}$ , $u'' = \frac{d^2u}{dx^2}$
