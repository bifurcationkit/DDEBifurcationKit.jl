# Hutchinson with Diffusion

```@contents
Pages = ["HutchinsonDiff.md"]
Depth = 3
```
Consider the Hutchinson equation with diffusion

$$\begin{aligned}
& \frac{\partial u(t, x)}{\partial t}=d \frac{\partial^2 u(t, x)}{\partial x^2}-a u(t-1, x)[1+u(t, x)], \quad t>0, x \in(0, \pi), \\
& \frac{\partial u(t, x)}{\partial x}=0, x=0, \pi
\end{aligned}$$

where $a>0, d>0$. 

## Problem discretization

We start by discretizing the above PDE based on finite differences.

```@example TUTHut
using Revise, DDEBifurcationKit, Plots, SparseArrays
using BifurcationKit
const BK = BifurcationKit

function Hutchinson(u, ud, p)
   (; a, d, Δ) = p
   d .* (Δ * u) .- a .* ud[1] .* (1 .+ u)
end

delaysF(par) = [1.]
nothing #hide
```

## Bifurcation analysis
We can now instantiate the model

```@example TUTHut
# discretisation
Nx = 200; Lx = pi/2;
X = -Lx .+ 2Lx/Nx*(0:Nx-1) |> collect
h = 2Lx/Nx
Δ = spdiagm(0 => -2ones(Nx), 1 => ones(Nx-1), -1 => ones(Nx-1) ) / h^2; Δ[1,1]=Δ[end,end]=-1/h^2
nothing #hide
```

We are now ready to compute the bifurcation of the trivial (constant in space) solution:

```@example TUTHut
# bifurcation problem
pars = (a = 0.5, d = 1, τ = 1.0, Δ = Δ, N = Nx)
x0 = zeros(Nx)

prob = ConstantDDEBifProblem(Hutchinson, delaysF, x0, pars, (@optic _.a))

optn = NewtonPar(eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 10., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 0, plot = false, normC = norminf)
br
```

We note that the first Hopf point is close to the theoretical value $a=\frac\pi 2$. This can be improved by increasing `opts.n_inversion`.

We can now plot the branch

```@example TUTHut
scene = plot(br)
```

## Performance improvements
The previous implementation being simple, it leaves a lot performance on the table. For example, the jacobian is dense because it is computed with automatic differentiation without sparsity detection. 

We show how to specify the jacobian and speed up the code a lot.

```@example TUTHut
# analytical jacobian
function JacHutchinson(u, p)
   (;a, d, Δ) = p
   # we compute the jacobian at the steady state
   J0 = d * Δ .- a .* Diagonal(u)
   J1 = -a .* Diagonal(1 .+ u)
   return J0, [J1]
end

prob2 = ConstantDDEBifProblem(Hutchinson, delaysF, x0, pars, (@optic _.a); J = JacHutchinson)

optn = NewtonPar(eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 10., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 4)
br = continuation(prob2, PALC(), opts; verbosity = 1, plot = true, normC = norminf)
br
```


## References
