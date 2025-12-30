# (Another) Neuron model (codim 2, periodic orbits)

```@contents
Pages = ["neuronV2.md"]
Depth = 3
```

Consider the neuron model

$$\begin{aligned}
& \dot{x}_1(t)=-x_1(t)-a g\left(b x_1\left(t-\tau_1\right)\right)+c g\left(d x_2\left(t-\tau_2\right)\right) \\
& \dot{x}_2(t)=-x_2(t)-a g\left(b x_2\left(t-\tau_1\right)\right)+c g\left(d x_1\left(t-\tau_2\right)\right)
\end{aligned}$$

where $g(z)=[\tanh (z-1)+\tanh (1)] \cosh (1)^2$. 

## Continuation and codim 1 bifurcations

We first instantiate the model

```@example TUTneuron2
using Revise, DDEBifurcationKit, Plots
using BifurcationKit
const BK = BifurcationKit

g(z) = (tanh(z − 1) + tanh(1)) * cosh(1)^2
function neuron2VF(x, xd, p)
   (; a,b,c,d) = p
   [
      -x[1] - a * g(b * xd.u[1][1]) + c * g(d * xd.u[2][2]),
      -x[2] - a * g(b * xd.u[1][2]) + c * g(d * xd.u[2][1])
   ]
end

delaysF(par) = [par.τ1, par.τ2]

pars = (a = 0.25, b = 2., c = 15/29, d = 1.2, τ1 = 12.7, τ2 = 20.2)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuron2VF, delaysF, x0, pars, (@optic _.a))

optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit=100))
opts = ContinuationPar(p_max = 1., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 9, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(tangent=Bordered()), opts)
```

We then plot the branch

```@example TUTneuron2
scene = plot(br)
```

## Normal forms computation

As in [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl), it is straightforward to compute the normal forms.

```@example TUTneuron2
hopfpt = BK.get_normal_form(br, 2)
```

## Continuation of Hopf points
We follow the Hopf points in the parameter plane $(a,c)$. We tell the solver to consider br.specialpoint[1] and continue it.

```@example TUTneuron2
# continuation of the first Hopf point
brhopf = continuation(br, 1, (@optic _.c),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 1.1, p_min = -0.1,ds = 0.01, n_inversion = 2);
         verbosity = 0,
         detect_codim2_bifurcation = 2,
         bothside = true,
         start_with_eigen = true)

brhopf2 = continuation(br, 2, (@optic _.c),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 1.1, p_min = -0.1,ds = -0.01);
         verbosity = 0,
         detect_codim2_bifurcation = 2,
         bothside = true,
         start_with_eigen = true)

scene = plot(brhopf, vars = (:a, :c), xlims = (0,0.7), ylims = (0, 1))
plot!(scene, brhopf2, vars = (:a, :c), xlims = (-0,0.7), ylims = (-0.1, 1))
scene
```

## Continuation of Fold points
We follow the Fold points in the parameter plane $(a, c)$. We tell the solver to consider br2.specialpoint[3] and continue it.

```@example TUTneuron2
prob2 = ConstantDDEBifProblem(neuron2VF, delaysF, x0, (pars..., a = 0.12), (@optic _.c))
br2 = continuation(prob2, PALC(), ContinuationPar(opts, p_max = 1.22);)


# change tolerance for avoiding error computation of the EV
opts_fold = br.contparams
@reset opts_fold.newton_options.eigsolver.σ = 1e-7

brfold = continuation(br2, 3, (@optic _.a),
         ContinuationPar(opts_fold; detect_bifurcation = 1, dsmax = 0.01, max_steps = 70, p_max = 0.6, p_min = -0.6,ds = -0.01, n_inversion = 2, tol_stability = 1e-6);
         verbosity = 1, plot = true,
         detect_codim2_bifurcation = 2,
         update_minaug_every_step = 1,
         bothside = false,
         start_with_eigen = true)

scene = plot(brfold, vars = (:a, :c), branchlabel = "Fold")
plot!(scene, brhopf, vars = (:a, :c), branchlabel = "Hopf")
plot!(scene, brhopf2, vars = (:a, :c), branchlabel = "Hopf")
scene
```
