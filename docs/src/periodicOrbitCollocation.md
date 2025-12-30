# Periodic orbits based on orthogonal collocation

```@contents
Pages = ["periodicOrbitCollocation.md"]
Depth = 3
```

We compute `Ntst` time slices of a periodic orbit using orthogonal collocation. This is implemented in the structure `BifurcationKit.PeriodicOrbitOCollProblem`.

!!! warning "Large scale"
    The current implementation is not yet optimized for large scale problems. This will be improved in the future.

!!! warning "Floquet coefficients"
    The current implementation does not yet allow for computing stability of periodic orbits. This will be improved in the future.  

The general method is explained in [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbitCollocation/).

We recall the basics for completeness.

## Introduction

The general method is very well exposed in [^Dankowicz],[^Doedel] and we adopt the notations of [^Dankowicz]. However our implementation is based on [^Doedel] because it is more economical (less equations) when it enforces the continuity of the solution.

We look for periodic orbits as solutions $(x(0), T)$ of

$$\dot x(t) = T\cdot F(x(t), x(t-\tau/T)),\ x(0)=x(1)\in\mathbb R^n.$$

We focus on the differential equality and consider a partition of the time domain

$$0=\tau_{1}<\cdots<\tau_{j}<\cdots<\tau_{N_{tst}+1}=1$$

where the points are referred to as **mesh points**. On each mesh interval $[\tau_j,\tau_{j+1}]$ for $j=1,\cdots,N_{tst}$, we define the affine transformation

$$\tau=\tau^{(j)}(\sigma):=\tau_{j}+\frac{(1+\sigma)}{2}\left(\tau_{j+1}-\tau_{j}\right), \sigma \in[-1,1].$$

The functions $x^{(j)}$ defined on $[-1,1]$ by $x^{(j)}(\sigma) \equiv x(\tau_j(\sigma))$ satisfies the following equation on $[-1,1]$:

$$\dot x^{(j)} = T\frac{\tau_{j+1}-\tau_j}{2}\cdot F(x^{(j)})\tag{$E_j$}$$

with the continuity equation $x^{(j+1)}(-1) = x^{(j)}(1)$.

We now aim at  solving $(E_j)$ by using an approximation with a polynomial of degree $m$. Following [^Dankowicz], we define a (uniform) partition:

$$-1=\sigma_{1}<\cdots<\sigma_{i}<\cdots<\sigma_{m+1}=1.$$

The points $\tau_{i,j} = \tau^{(i)}(\sigma_j)$ are called the **base points**: they serve as collocation points.

The associated $m+1$ Lagrange polynomials of degree $m$ are:

$$\mathcal{L}_{i}(\sigma):=\prod_{k=1, k \neq i}^{m+1} \frac{\sigma-\sigma_{k}}{\sigma_{i}-\sigma_{k}}, i=1, \ldots, m+1.$$

We then introduce the approximation $p_j$ of $x^{(j)}$:

$$\mathcal p_j(\sigma)\equiv \sum\limits_{k=1}^{m+1}\mathcal L_k(\sigma)x_{j,k}$$

and the problem to be solved at the **nodes** $z_l$, $l=1,\cdots,m$:

$$\forall 1\leq l\leq m,\quad 1\leq j\leq N_{tst},\quad \dot p_j(z_l) = T\frac{\tau_{j+1}-\tau_j}{2}\cdot F(p_j(z_l), p_{j_0}(``t_l-\tau''))\tag{$E_j^2$}.$$

The **nodes** $(z_l)$ are associated with a Gauss–Legendre quadrature.


## Mesh adaptation

Mesh adaptation can be turned on like in the case of [ODEs](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/periodicOrbitCollocation/#Mesh-adaptation).
