## 1D Advection-Diffusion Equation

### Strong Form

$$
\begin{aligned}
-\nu u^{\prime\prime} + au^\prime &= f \quad &&\text{in } {\cal B} = (0,1), \\
u(0) = u(1) &= 0 \quad &&\text{on } \partial{\cal B}_D
\end{aligned}
$$

where $\nu \geq 0$ is the diffusion coefficient and $a$ is the advection velocity with the following parameters

**Parameters:**
- $\nu = \epsilon = 10^{-3}$ / $\nu = \epsilon = 10^{-5}$
- $a = 1$ 

Analytical solution:

$$
u(x) = 4\left(\arctan\left(\frac{2(1/16-(x-1/2)^2)}{\pi\sqrt{\epsilon}}\right) + \frac{1}{2}\right)(1-x)x
$$

**Test Configuration:**
- Mesh levels: ML5 (32 elements) to ML9 (512 elements)
- Element types: $P_1$ (linear) and $P_2$ (quadratic)
- Mesh types: Regular and perturbed grids

---
## Methods Compared

### 1. Standard Galerkin FEM (SFEM)

The classical Galerkin finite element formulation seeks $u_h \in V_h \subset H^1_0({\cal B})$ such that

$$(\nu u^\prime_h, v^\prime_h)_0 + (a u^\prime_h, v_h)_0 = (f, v_h)_0 \quad \forall v_h \in V_h$$

#

### 2. Least Squares FEM - Diffusive Flux $(LSFEM^{(d)})$:

**Flux:** $q = \nu u^\prime$  
**Balance:** $-q^\prime + a u^\prime = f$  
**Constitutive:** $q - \nu u^\prime = 0$

### Functional

$${\cal J}^{(d)}(u,q;f) = \frac{1}{2} \parallel (-q^\prime + a u^\prime - f) \parallel^2_0 + \frac{1}{2} \parallel w_2 (q - \nu u^\prime) \parallel^2_0$$

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}(u_h,q_h;v_h,r_h) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}^{(d)}(u_h,q_h;v_h,r_h) = (-q^\prime_h + a u^\prime_h , -r^\prime_h + a v^\prime_h)_0 + w^2_2 (q_h - \nu u^\prime_h, r_h - \nu v^\prime_h)_0$$

$${\cal F}^{(d)}(v_h,r_h) = (f, -r^\prime_h + a v^\prime_h)_0$$

**Weighting analyzed:**
- **$LSFEM^{(d)}$**: $w_2 = 1$ 
- **$wLSFEM^{(d)}$**: $w_2 = \nu^{-1/2}$ 

#

### 3. Least Squares FEM - Total Flux $(LSFEM^{(t)})$:

**Flux:** $q = \nu u^\prime - a u$  
**Balance:** $-q^\prime = f$  
**Constitutive:** $q - \nu u^\prime + a u = 0$

### Functionals

$${\cal J}^{(t)}(u,q;f) = \frac{1}{2} \parallel (-q^\prime - f) \parallel^2_0 + \frac{1}{2} \parallel w_2 (q - \nu u^\prime + a u) \parallel^2_0$$

### Variational Formulation

The discrete variational form seeks $(u_h, q_h) \in V_h \times Q_h$ satisfying

$${\cal A}(u_h,q_h;v_h,r_h) = {\cal F}(v_h,r_h), \quad \forall (v_h,r_h) \in V_h \times Q_h$$

with

$${\cal A}^{(t)}(u_h,q_h;v_h,r_h) = (-q^\prime_h, -r^\prime_h)_0 + w^2_2 (q_h - \nu u^\prime_h + a u^\prime_h, r_h - \nu v^\prime_h + a v^\prime_h)_0$$

$${\cal F}^{(t)}(v_h,r_h) = (f, -r^\prime_h)_0$$


**Weighting analyzed:**
- **$LSFEM^{(t)}$**: $w_2 = 1$ 
- **$wLSFEM^{(t)}$**: $w_2 = \nu^{-1/2}$ 

---






