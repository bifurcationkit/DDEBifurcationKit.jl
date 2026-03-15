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
opts = ContinuationPar(p_max = 1., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 9, dsmax = 0.1, n_inversion = 4)
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

scene = plot(brhopf,  vars = (:a, :c))
plot!(scene, brhopf2, vars = (:a, :c), xlims = (0, 0.7), ylims = (-0.1, 1))
scene
```

## Continuation of branch points
We follow the branch points in the parameter plane $(a, c)$.

```@example TUTneuron2
prob2 = ConstantDDEBifProblem(neuron2VF, delaysF, x0, (pars..., a = 0.12), (@optic _.c))
br2 = continuation(prob2, PALC(), ContinuationPar(opts, p_max = 1.22);)

# change tolerance for avoiding error computation of the eigenvalues
opts_bp = br.contparams
@reset opts_bp.newton_options.eigsolver.σ = 1e-7

# index of the first branch point
index_bp = findfirst(x -> x.type == :bp, br2.specialpoint)
brbp = continuation(br2, index_bp, (@optic _.a),
         ContinuationPar(opts_bp; detect_bifurcation = 1, dsmax = 0.01, max_steps = 70, p_max = 0.6, p_min = -0.6, ds = -0.01, n_inversion = 2, tol_stability = 1e-6);
         verbosity = 0,
         plot = true,
         detect_codim2_bifurcation = 2,
         bothside = false,
         start_with_eigen = true)

scene = plot(brbp, vars = (:a, :c), branchlabel = "Branch points")
plot!(scene, brhopf, vars = (:a, :c), branchlabel = "Hopf")
plot!(scene, brhopf2, vars = (:a, :c), branchlabel = "Hopf")
scene
```

## Branch of periodic orbits

We compute the branch of periodic orbits from a Hopf bifurcation.

```@example TUTneuron2
pars = (a = 0.069, b = 2., c = 0.6, d = 1.2, τ1 = 11.6, τ2 = 20.3)
delaysF(par) = [par.τ1, par.τ2]
prob3 = ConstantDDEBifProblem(neuron2VF, delaysF, zeros(2), pars, (@optic _.c))
br3 = continuation(prob3, PALC(), ContinuationPar(opts, p_max = 1.0);)
```

```@example TUTneuron2
opts_po_cont = ContinuationPar(dsmax = 0.02, ds = -0.0001, max_steps = 100, nev = 10, tol_stability = 1e-4)

probpo = PeriodicOrbitOCollProblem(20, 5; N = 2, 
            jacobian = BK.AutoDiffDense(),)

br_po = @time continuation(
            br3, 1, ContinuationPar(opts_po_cont; max_steps = 185, detect_bifurcation = 3),
            probpo;
            plot = true,
            normC = norminf,
            eigsolver = BK.FloquetGEV(DDE_DefaultEig(maxit=100, tol = 1e-12, σ = 1e-3)),
            )
plot(br_po, vars = (:param, :amplitude))
```

```@example TUTneuron2
plot(br3, br_po)
```