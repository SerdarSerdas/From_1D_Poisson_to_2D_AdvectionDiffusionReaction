# 1D Poisson Equation

## Strong Form ($\nu = 1$)

$$
\begin{aligned}
-u'' &= f &&\text{in } \mathcal{B} = (0,1), \\
u &= 0 &&\text{on } \partial\mathcal{B}_D = \{0,1\}.
\end{aligned}
$$

## SFEM (Standard Galerkin)

Find 
$u_h \in V_h \subset H^1_0(\mathcal{B})$ such that:

$$
(u_h', v_h') = (f, v_h) \quad \forall v_h \in V_h
$$

## LSFEM (Least-Squares)

**First-order system:** $\mathbf{q} = u'$

$$
\begin{aligned}
-q' &= f, \\
q - u' &= 0,
\end{aligned}
\qquad \text{in } \mathcal{B}
$$

**Least-squares functional:**
$$
\mathcal{J}(u,q;f) = \|q' + f\|_0^2 + \|q - u'\|_0^2
$$

**Discrete variational form:** Find $(u_h,\mathbf{q}_h) \in V_h \times \mathbf{Q}_h$

$$
\mathcal{A}((u_h,\mathbf{q}_h);(v_h,\mathbf{r}_h)) = \mathcal{F}(v_h,\mathbf{r}_h)
$$

where
$$
\begin{aligned}
\mathcal{A}((u,q);(v,r)) &= (q',r')_0 + (q - u',r - v')_0, \\
\mathcal{F}(v,r) &= (f,r')_0.
\end{aligned}
$$
