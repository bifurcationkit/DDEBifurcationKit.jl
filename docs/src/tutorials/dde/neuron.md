# Neuron model (codim 2, periodic orbits)

```@contents
Pages = ["neuron.md"]
Depth = 3
```
Consider the neuron model

$$\left\{\begin{array}{l}
\dot{x}_1(t)=-\kappa x_1(t)+\beta \tanh \left(x_1\left(t-\tau_s\right)\right)+a_{12} \tanh \left(x_2\left(t-\tau_2\right)\right) \\
\dot{x_2}(t)=-\kappa x_2(t)+\beta \tanh \left(x_2\left(t-\tau_s\right)\right)+a_{21} \tanh \left(x_1\left(t-\tau_1\right)\right)
\end{array}\right.$$


## Continuation and codim 1 bifurcations

We first instantiate the model

```@example TUTneuron
using Revise, DDEBifurcationKit, Parameters, Setfield, LinearAlgebra, Plots
using BifurcationKit
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

function neuronVF(x, xd, p)
   @unpack κ, β, a12, a21, τs, τ1, τ2 = p
   [
      -κ * x[1] + β * tanh(xd[3][1]) + a12 * tanh(xd[2][2]),
      -κ * x[2] + β * tanh(xd[3][2]) + a21 * tanh(xd[1][1])
   ]
end

delaysF(par) = [par.τ1, par.τ2, par.τs]

pars = (κ = 0.5, β = -1, a12 = 1, a21 = 0.5, τ1 = 0.2, τ2 = 0.2, τs = 1.5)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@lens _.τs))

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 13., pMin = 0., newtonOptions = optn, ds = -0.01, detectBifurcation = 3, nev = 5, dsmax = 0.2, nInversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = true, normC = norminf)
```

We then plot the branch

```@example TUTneuron
scene = plot(br)
```

## Normal forms computation

As in [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl), it is straightforward to compute the normal forms.

```@example TUTneuron
hopfpt = BK.getNormalForm(br, 2)
```

## Continuation of Hopf points
We follow the Hopf points in the parameter plane $(a_{21},\tau_s)$. We tell the solver to consider br.specialpoint[3] and continue it.

```@example TUTneuron
# continuation of the first Hopf point
brhopf = continuation(br, 3, (@lens _.a21),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.04, maxSteps = 230, pMax = 15., pMin = -1.,ds = -0.02);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 2,
         # bothside = true,
         startWithEigen = true)

# continuation of the second Hopf point
brhopf2 = continuation(br, 2, (@lens _.a21),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.1, maxSteps = 56, pMax = 15., pMin = -1.,ds = -0.01, nInversion = 4);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 2,
         startWithEigen = true,
         bothside=true)

scene = plot(brhopf, brhopf2, legend = :top)
```

## Branch of periodic orbits

We change the continuation parameter and study the bifurcations as function of $a_{21}$.

```@example TUTneuron
prob2 = ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@lens _.a21))
br2 = BK.continuation(prob2, PALC(), setproperties(opts, ds = 0.1, pMax = 3., nInversion=8); verbosity = 0, plot = false, normC = norminf)
```

We then compute the branch of periodic orbits from the Hopf bifurcation points using orthogonal collocation.

```@example TUTneuron
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, pMax = 10., pMin=-0., maxSteps = 120, detectBifurcation = 0, saveSolEveryStep=1)
@set! opts_po_cont.newtonOptions.tol = 1e-8
@set! opts_po_cont.newtonOptions.verbose = false

# arguments for periodic orbits
args_po = (	recordFromSolution = (x, p) -> begin
			xtt = BK.getPeriodicOrbit(p.prob, x, nothing)
			return (max = maximum(xtt[1,:]),
					min = minimum(xtt[1,:]),
					period = getPeriod(p.prob, x, nothing))
		end,
		plotSolution = (x, p; k...) -> begin
			xtt = BK.getPeriodicOrbit(p.prob, x, nothing)
			plot!(xtt.t, xtt[1,:]; label = "V1", k...)
			plot!(xtt.t, xtt[2,:]; label = "V2", k...)
			plot!(br2; subplot = 1, putspecialptlegend = false)
			end,
		normC = norminf)

probpo = PeriodicOrbitOCollProblem(60, 4; N = 2)
	br_pocoll = @time continuation(
		br2, 1, opts_po_cont,
		probpo;
		verbosity = 0,	plot = false,
		args_po...,
		ampfactor = 1/0.3593 * 0.0610*2.2,
		δp = 0.001,
		normC = norminf,
		)
scene = plot(br2, br_pocoll)		
```

We can plot the periodic orbit as they approach the homoclinic point.

```@example TUTneuron
scene = plot(layout = 2)
for ii = 1:10:110
	solpo = BK.getPeriodicOrbit(br_pocoll.γ.prob.prob, br_pocoll.sol[ii].x, nothing)
	plot!(scene, solpo.t ./ solpo.t[end], solpo.u[1,:], label = "", subplot = 1)
end
xlabel!(scene, "t / period", subplot = 1)
plot!(scene, br_pocoll, vars = (:param, :period), subplot = 2, xlims=(2.2,2.4))
scene
```

## References
