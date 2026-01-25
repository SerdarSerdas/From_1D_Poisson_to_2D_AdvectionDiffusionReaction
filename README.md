# 🚀 From 1D Poisson to 2D Advection-Diffusion-Reaction

**SFEM vs LSFEM** comparative study across progressively complex PDEs (1D→2D).

## 📋 Contents

| Section | Problem | Dimensions |
|---------|---------|------------|
| [00] | [ADR Equation](00_ProblemDefinition/) | 1D/2D |
| [01] | Poisson | 1D ⏳ |
| [02] | Diffusion-Reaction | 1D ⏳ |
| [03] | Advection-Diffusion | 1D ⏳ |
| [04] | Advection-Diffusion-Reaction | 1D ⏳ |
| [05] | Poisson | 2D ⏳ |
| [06] | Diffusion-Reaction | 2D ⏳ |
| [07] | Advection-Diffusion | 2D ⏳ |
| [08] | Advection-Diffusion-Reaction | 2D ⏳ |

## 🔬 Methods Compared

| Method | Formulation | Key Features |
|--------|-------------|--------------|
| **SFEM** | Standard Galerkin FEM | Uniform refinement only |
| **LSFEM** | First-order least-squares FEM | **Dörfler + α-bulk adaptive refinement**<br/>+ robust on distorted meshes |  

**Discretizations:** P₁/P₂ (1D), Q₁/Q₂ & P₁/P₂ (2D)  
**Meshes:** Regular + perturbed (α=0.2) + LSFEM adaptive
