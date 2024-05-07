# Bifurcation Problem


## Constant delays (DDE)

Consider the system of delay differential equations with constant delays (DDEs)

$$\frac{\mathrm{d}}{\mathrm{d} t} x(t)=\mathbf F\left(x(t), x\left(t-\tau_1\right), \ldots, x\left(t-\tau_m\right); p\right)$$

where the delays $\tau_i>0$ are constant and $p$ is a set of parameters. In order to specify this, we need to provide the vector field and the delays. The delays are provided using a delay function which must return an `AbstractVecor`

```julia
function mydelays(pars)
	[1, pars.tau1]
end
```

where `pars` are some user defined variables. The vector field is then specified as follows

```julia
function myF(x, xd, pars)
	[
		x[1] + xd[2][1]^2,
		x[2] + xd[3][2]^2,
	]
end
```

where `xd` is a vector holding `[x(t-d[1]), x(t-d[2])]` where `d = mydelays(pars)`. Some simple examples can be found in the tutorials.

The structure [`ConstantDDEBifProblem`](@ref) encapsulates the bifurcation problem.


## State-dependent delays (SD-DDE)

Consider the system of delay differential equations with state-dependent delays.

$$\frac{\mathrm{d}}{\mathrm{d} t} x(t)=\mathbf F\left(x(t), x\left(t-\tau_1(x(t))\right), \ldots, x\left(t-\tau_m(x(t))\right); p\right)$$


where the delays $\tau_i>0$ are functions of $x(t)$ and $p$ is a set of parameters. The only difference with the previous case is the specification of the delay function which now depends on `x`

```julia
function mydelays(x, pars)
	[
		1 + x[1]^2,
		2 + x[2]^2
	]
end
```

The structure [`SDDDEBifProblem`](@ref) encapsulates the bifurcation problem.

> A more elaborate problem would be to allow $\tau_i$ to depend on the history of $x(\theta+t)$ for $\theta\in[-\tau_{max},0]$ and not just on the current value of $x(t)$. This is not implemented yet.
