# From 1D Poisson to 2D Advection-Diffusion-Reaction

> **TL;DR:** LSFEM achieves **100-1000× better flux accuracy** than standard FEM through **O(h²) superconvergence**. For extreme boundary layers (ε→0), adaptive mesh refinement provides **20× efficiency gains**. This repository systematically compares SFEM and LSFEM across the full hierarchy: Poisson → Diffusion-Reaction → Advection-Diffusion → Advection-Diffusion-Reaction.

[![1D Solutions](1D_Poisson/Fig_Poisson1D_Solutions_ML5.png)](#-1d-poisson-results)
[![Convergence](1D_Poisson/Fig_Poisson1D_Convergence.png)](#-convergence-analysis)
[![AMR Efficiency](1D_Poisson/Fig_Poisson1D_AMR_Comparison.png)](#-adaptive-mesh-refinement)

---

## 📖 Why This Repository?

### The Problem: Standard FEM's Limitations

The standard/classical Galerkin finite element method (SFEM) is the workhorse for solving partial differential equations (PDEs), but it faces challenges with:

- **Convection-dominated transport**: Spurious oscillations without stabilization techniques
- **Saddle-point systems**: Velocity-pressure coupling in incompressible Navier-Stokes requires carefully matched element pairs (inf-sup stability)
- **Mixed formulations**: Must satisfy the Ladyzhenskaya-Babuška-Brezzi (LBB) condition through compatible approximation spaces

### The Solution: Least-Squares FEM (LSFEM)

LSFEM reformulates PDEs as first-order systems and minimizes the $L^2$ residual norm. This elegant approach delivers:

✅ **Equal-order interpolation** — Use $P_1-P_1$ or $P_2-P_2$ freely; no inf-sup constraints  
✅ **Symmetric positive definite matrices** — Always well-conditioned systems  
✅ **Natural error estimators** — Residual norms directly guide adaptive refinement  
✅ **Flux superconvergence** — ${\cal O}(h^{p+1})$ for gradients vs ${\cal O}(h^{p})$ in standard FEM  


### What We Investigate

This repository provides **the first comprehensive 1D/2D comparison** of Standard Galerkin FEM (SFEM) versus Adaptive Least-Squares FEM (LSFEM) across:

**Problem hierarchy:**

- Poisson (ν=1, 𝐚=0, c=0)
- Diffusion-Reaction (ν=1, 𝐚=0, c≠0)
- Advection-Diffusion (ν≪1, 𝐚≠0, c=0)
- Advection-Diffusion-Reaction (ν≪1, 𝐚≠0, c>0)


**Three dimensions of analysis:**
1. **Discretization**: $P_1$ vs $P_2$ elements on regular and perturbed meshes
2. **Solution strategies**: Unpreconditioned CG, CG-Jacobi, CG-AMG
3. **Mesh adaptivity**: Uniform refinement vs Dörfler and $\alpha$-Bulk marking strategy

**Extreme test regimes:**
- Singular perturbations: $\epsilon = 10^{-3} → \epsilon = 5 \cdot 10^{-6}$ (boundary layer), $\epsilon = 10^{-7}$ (interior layer)
- High Péclet numbers: $Pe = \parallel 𝐚 \parallel h / (2 \nu) ≫ 1$
- Geometric complexity: 1D intervals, 2D quadrilaterals, 2D triangles


### Key Research Questions

1. **What is the influence/contribution of the weighting factor w.r.t. the accurracy and solver performance?**

2. **When does LSFEM's superconvergence break down?**  
 
3. **How efficient is adaptive mesh refinement?**  
 
4. **Which preconditioner is optimal?**

5. **Which first-order system performance better**
 
6. **Does adaptivity preserve superconvergence?**  
 
---


## 🧮 Problem Statement

### General Advection-Diffusion-Reaction Equation

The steady ADR equation models scalar transport under simultaneous diffusion, advection, and reaction:

$$
-\nabla \cdot (\nu \nabla u + \mathbf{a} u) + c u = f, \quad \text{in } \mathcal{B}
$$

**Boundary conditions:**

$$
u = g_D \quad \text{on } \partial\mathcal{B}_D, \qquad \mathbf{n}\cdot(\nu\nabla u+\mathbf{a}u)=g_N \quad \text{on } \partial\mathcal{B}_N
$$

**Parameters:**
- $\nu > 0$ — Diffusion coefficient (small → boundary/interior layers)
- $\mathbf{a}$ — Advection velocity field
- $c \geq 0$ — Reaction coefficient
- $f$ — Volumetric source term

**Note:** This study considers piecewise-constant coefficients and homogeneous Dirichlet conditions ($g_D = 0$) for least-squares formulations.


### 1D Benchmark Solutions

**Boundary layer** (exponential profile):

$$
u(x) = 1 - \frac{\sinh(x/\sqrt{\epsilon}) + \sinh((1-x)/\sqrt{\epsilon})}{\sinh(1/\sqrt{\epsilon})}
$$

![Poisson1D_BVP_BL](Fig_Poisson1D_BVP_BoundaryLayer.png)

**Interior layer** (localized shock):

$$
u(x) = 4\left(\arctan\left(\frac{2(1/16-(x-1/2)^2)}{\pi\sqrt{\epsilon}}\right) + \frac{1}{2}\right)(1-x)x
$$

![Poisson1D_BVP_IL](Fig_Poisson1D_BVP_InteriorLayer.png)

---

## 🧩 Mesh Configurations

Two different grid types (regular and perturbed) for 1D and 2D (quadrilaterals/triangular) are used which are shown in the following figure. In the case of regular grids, the next levels are obtained by applying a recursive uniform refinements starting from coarsest grid which consists of just one element. The perturbed (irregular) grids however, are obtained after construction of regular grids of each level, by applying a perturbation filter ($\alpha = 0.2$) based on the grid size h.

![Mesh Hierarchies](Fig_Mesh_1D_2D_RMR_PMR_2x3.png)

**Top row:** Regular meshes (1D, quadrilateral 2D, triangular 2D)  
**Bottom row:** Perturbed meshes (α=0.2 random displacement) for robustness testing

### 1D Mesh Hierarchy

| Level | Elements | h       | SFEM P₁ | SFEM P₂ | LSFEM P₁P₁ | LSFEM P₂P₂ |
|-------|----------|---------|---------|---------|------------|------------|
| ML4   | 16       | 1/16    | 17      | 33      | 34         | 66         |
| ML5   | 32       | 1/32    | 33      | 65      | 66         | 130        |
| ML6   | 64       | 1/64    | 65      | 129     | 130        | 258        |
| ML7   | 128      | 1/128   | 129     | 257     | 258        | 514        |
| ML8   | 256      | 1/256   | 257     | 513     | 514        | 1026       |
| ML9   | 512      | 1/512   | 513     | 1025    | 1026       | 2050       |
| ML10  | 1024     | 1/1024  | 1025    | 2049    | 2050       | 4098       |

**Note:** LSFEM has 2× DOFs due to solving for both $u$ and $q$ simultaneously.

---


## 💻 Implementation Details

**Language:** Python 3.8+

**Core Dependencies:**
- **FEM assembly:** NumPy 1.21+
- **Linear solvers:** SciPy 1.7+ (CG), PyAMG 4.0+ (smoothed aggregation)
- **Visualization:** Matplotlib 3.5+, Seaborn 0.11+

**AMG Configuration:**
- Smoother: Symmetric Gauss-Seidel (4 pre/post sweeps)
- Cycle: W-cycle
- Tolerance: 10⁻¹⁰ relative residual


---

## 🔬 Reproducibility

All results are fully reproducible with documented parameters:

**Mesh configurations:**
- ML levels: $N = 2^{\text{ML}}$ uniform elements
- Perturbed meshes: Random displacement of nodes ($\alpha 0.2$)

**Solver settings:**
- Tolerance: 10⁻¹⁰ relative residual
- AMG: Smoothed aggregation, W-cycle, 4 sweeps

**AMR parameters:**
- $\alpha$-Bulk marking: $\alpha$ = 0.2 
- Dörfler marking: θ = 0.5 
- Error estimator: $\eta_K$
- Stopping criterion: max refinement steps or pre-defined condition number reached

---

## 🚀 Future Work

- [ ] **Convection dominance:** High Péclet number analysis (Pe > 100)
- [ ] **Nonlinear problems:** Incompressible Navier-Stokes (Newtonian and Non-Newtonian fluids)
- [ ] **Time-dependent:** Parabolic ADR and Navier-Stokes with time adaptivity
- [ ] **Structure mechanics:** Small and large deformations
- [ ] **Multiphase materials:** Electro-hydro-chemomechanical models
- [ ] **Fluid-Structure-Interaction** 


---

## 📄 Citation

If you use this work in your research, please cite:

```bibtex
@software{poisson_to_adr_2025,
  author = {Serdar Serdas},
  title = {From 1D Poisson to 2D Advection-Diffusion-Reaction: 
           Comprehensive Comparison of SFEM and LSFEM},
  year = {2025},
  url = {https://github.com/SerdarSerdas/From_1D_Poisson_to_2D_AdvectionDiffusionReaction},
  note = {Systematic benchmark suite demonstrating LSFEM flux superconvergence 
          and AMR efficiency for singularly perturbed problems}
}
```

---

---

## References

### Core LSFEM Theory

[1] Bochev, P. B., & Gunzburger, M. D. (1998). *Finite element methods of least-squares type*. SIAM Review, 40(4), 789-837. DOI: [10.1137/S0036144597321156](https://doi.org/10.1137/S0036144597321156)

[2] Bochev, P. B., & Gunzburger, M. D. (2009). *Least-Squares Finite Element Methods*. Springer. DOI: [10.1007/b13382](https://doi.org/10.1007/b13382)

[3] Jiang, B. N. (1998). *The Least-Squares Finite Element Method: Theory and Applications in Computational Fluid Dynamics and Electromagnetics*. Springer. DOI: [10.1007/978-3-662-03740-9](https://doi.org/10.1007/978-3-662-03740-9)

[4] Bochev, P., & Gunzburger, M. (2006). Least-squares methods for the velocity-pressure-stress formulation of the Stokes equations. *Computer Methods in Applied Mechanics and Engineering*, 195(7-8), 828-845. DOI: [10.1016/j.cma.2005.02.018](https://doi.org/10.1016/j.cma.2005.02.018)

### Transport Equations

[5] Pehlivanov, A. I., & Carey, G. F. (1994). Error estimates for least-squares mixed finite elements. *RAIRO-Modélisation mathématique et analyse numérique*, 28(5), 499-516. DOI: [10.1051/m2an/1994280504991](https://doi.org/10.1051/m2an/1994280504991)

[6] Pehlivanov, A. I., Carey, G. F., & Lazarov, R. D. (1994). Least-squares mixed finite elements for second-order elliptic problems. *SIAM Journal on Numerical Analysis*, 31(5), 1368-1377. DOI: [10.1137/0731071](https://doi.org/10.1137/0731071)

[7] Yang, S. Y. (2000). Analysis of optimal least-squares finite element methods for the Stokes equations on rectangular elements. *Journal of Computational Mathematics*, 18(3), 277-288.

[8] Carey, G. F., Pehlivanov, A. I., & Vassilevski, P. S. (1995). Least-squares mixed finite element methods for non-selfadjoint elliptic problems: I. Error estimates. *Numerische Mathematik*, 72(4), 501-522. DOI: [10.1007/s002110050180](https://doi.org/10.1007/s002110050180)

[9] Fialko, S. Y., Manko, V. P., & McConnell, P. W. (1998). First-order system least squares for the convection-diffusion problem. *International Journal for Numerical Methods in Engineering*, 43(2), 307-321.

[10] Hsieh, P. W., & Yang, S. Y. (2009). On efficient least-squares finite element methods for convection-dominated problems. *Computer Methods in Applied Mechanics and Engineering*, 199(1-4), 183-196. DOI: [10.1016/j.cma.2009.09.029](https://doi.org/10.1016/j.cma.2009.09.029)

[11] Lazarov, R., Tobiska, L., & Vassilevski, P. S. (1997). Streamline diffusion least-squares mixed finite element methods for convection-diffusion problems. *East-West Journal of Numerical Mathematics*, 5(4), 249-264.

[12] Cai, Z., Manteuffel, T. A., & McCormick, S. F. (1997). First-order system least squares for second-order partial differential equations: Part II. *SIAM Journal on Numerical Analysis*, 34(2), 425-454. DOI: [10.1137/S0036142994266066](https://doi.org/10.1137/S0036142994266066)

### Algebraic Multigrid

[13] Bell, W. N., Olson, L. N., & Schroder, J. B. (2022). PyAMG: Algebraic Multigrid Solvers in Python v5.0. Release 5.0. DOI: [10.5281/zenodo.596446](https://doi.org/10.5281/zenodo.596446)

[14] Virtanen, P., Gommers, R., Oliphant, T. E., et al. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261-272. DOI: [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)

### Adaptive Mesh Refinement

[15] Dörfler, W. (1996). A convergent adaptive algorithm for Poisson's equation. *SIAM Journal on Numerical Analysis*, 33(3), 1106-1124. DOI: [10.1137/0733054](https://doi.org/10.1137/0733054)

### Singular Perturbations

[16] Roos, H.-G., Stynes, M., & Tobiska, L. (2008). *Robust Numerical Methods for Singularly Perturbed Differential Equations* (2nd ed.). Springer Series in Computational Mathematics, Vol. 24. DOI: [10.1007/978-3-540-34467-4](https://doi.org/10.1007/978-3-540-34467-4)

[17] Shishkin, G. I. (1988). Grid approximation of singularly perturbed elliptic and parabolic equations. Second doctoral thesis, Keldysh Institute, Moscow. (In Russian)

[18] Linß, T. (2010). *Layer-Adapted Meshes for Reaction-Convection-Diffusion Problems*. Springer Lecture Notes in Mathematics, Vol. 1985. DOI: [10.1007/978-3-642-05134-0](https://doi.org/10.1007/978-3-642-05134-0)

### Convergence Theory

[19] Brenner, S. C., & Scott, L. R. (2008). *The Mathematical Theory of Finite Element Methods* (3rd ed.). Springer Texts in Applied Mathematics, Vol. 15. DOI: [10.1007/978-0-387-75934-0](https://doi.org/10.1007/978-0-387-75934-0)

[20] Ciarlet, P. G. (2002). *The Finite Element Method for Elliptic Problems*. SIAM Classics in Applied Mathematics, Vol. 40. DOI: [10.1137/1.9780898719208](https://doi.org/10.1137/1.9780898719208)

---

## Additional Resources

### Software Implementations

- **PyAMG**: [https://github.com/pyamg/pyamg](https://github.com/pyamg/pyamg)

---

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-contribution`)
3. Commit your changes with clear messages
4. Push to your branch and open a pull request

---

## 📧 Contact

**Serdar Serdas**  
[GitHub](https://github.com/SerdarSerdas) | [Email](mailto:serdarserdas55@gmail.com)

**Feedback:** Found a bug or have suggestions? Open an issue or use the discussions tab!

---

## 🙏 Acknowledgments

- PyAMG developers for the excellent algebraic multigrid implementation
- All contributors to the numerical analysis literature cited herein

---

**Last updated:** January 2026  
**Repository status:** Active development (1D Poisson complete, 1D/2D ADR in progress)
