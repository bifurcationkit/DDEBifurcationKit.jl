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
using Revise, DDEBifurcationKit, Plots
using BifurcationKit
const BK = BifurcationKit

function neuronVF(x, xd, p)
   (; κ, β, a12, a21, τs, τ1, τ2) = p
   [
      -κ * x[1] + β * tanh(xd[3][1]) + a12 * tanh(xd[2][2]),
      -κ * x[2] + β * tanh(xd[3][2]) + a21 * tanh(xd[1][1])
   ]
end

delaysF(par) = [par.τ1, par.τ2, par.τs]

pars = (κ = 0.5, β = -1, a12 = 1, a21 = 0.5, τ1 = 0.2, τ2 = 0.2, τs = 1.5)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@optic _.τs))

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 13., p_min = 0., newton_options = optn, ds = -0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = true, normC = norminf)
```

We then plot the branch

```@example TUTneuron
scene = plot(br)
```

## Normal forms computation

As in [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl), it is straightforward to compute the normal forms.

```@example TUTneuron
hopfpt = BK.get_normal_form(br, 2)
```

## Continuation of Hopf points
We follow the Hopf points in the parameter plane $(a_{21},\tau_s)$. We tell the solver to consider br.specialpoint[3] and continue it.

```@example TUTneuron
# continuation of the first Hopf point
brhopf = continuation(br, 3, (@optic _.a21),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.04, max_steps = 230, p_max = 15., p_min = -1.,ds = -0.02);
         detect_codim2_bifurcation = 2,
         # bothside = true,
         start_with_eigen = true)

# continuation of the second Hopf point
brhopf2 = continuation(br, 2, (@optic _.a21),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.1, max_steps = 56, p_max = 15., p_min = -1.,ds = -0.01, n_inversion = 4);
         detect_codim2_bifurcation = 2,
         start_with_eigen = true,
         bothside=true)

scene = plot(brhopf, brhopf2, legend = :top)
```

## Branch of periodic orbits

We change the continuation parameter and study the bifurcations as function of $a_{21}$.

```@example TUTneuron
prob2 = ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@optic _.a21))
br2 = BK.continuation(prob2, PALC(), ContinuationPar(opts, ds = 0.1, p_max = 3., n_inversion = 8); verbosity = 0, plot = false, normC = norminf)
```

We then compute the branch of periodic orbits from the Hopf bifurcation points using orthogonal collocation.

```@example TUTneuron
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, p_max = 10., p_min=-0., max_steps = 120, detect_bifurcation = 0, save_sol_every_step=1)
@reset opts_po_cont.newton_options.tol = 1e-8
@reset opts_po_cont.newton_options.verbose = false

# arguments for periodic orbits
args_po = (	record_from_solution = (x, p; k...) -> begin
			xtt = BK.get_periodic_orbit(p.prob, x, nothing)
			return (max = maximum(xtt[1,:]),
					min = minimum(xtt[1,:]),
					period = getperiod(p.prob, x, nothing))
		end,
		plot_solution = (x, p; k...) -> begin
			xtt = BK.get_periodic_orbit(p.prob, x, nothing)
			plot!(xtt.t, xtt[1,:]; label = "V1", k...)
			plot!(xtt.t, xtt[2,:]; label = "V2", k...)
			plot!(br2; subplot = 1, putspecialptlegend = false)
			end,
		normC = norminf)

probpo = PeriodicOrbitOCollProblem(60, 4; N = 2, jacobian = BK.AutoDiffDense())
br_pocoll = @time continuation(
	br2, 1, opts_po_cont,
	probpo;
	verbosity = 0,	plot = false,
	args_po...,
	δp = 0.003,
	normC = norminf,
	)
scene = plot(br2, br_pocoll)
```

We can plot the periodic orbit as they approach the homoclinic point.

```@example TUTneuron
scene = plot(layout = 2)
for ii = 1:10:110
	solpo = BK.get_periodic_orbit(br_pocoll.γ.prob.prob, br_pocoll.sol[ii].x, nothing)
	plot!(scene, solpo.t ./ solpo.t[end], solpo.u[1,:], label = "", subplot = 1)
end
xlabel!(scene, "t / period", subplot = 1)
plot!(scene, br_pocoll, vars = (:param, :period), subplot = 2, xlims=(2.2,2.4))
scene
```

## References
