using Revise, DDEBifurcationKit
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

g(x,p) = (exp(x)-1)/(p.c1 * exp(x) + p.c2)

function vdp(x, xd, p)
   (;ϵ, τ) = p
   x1, x2 = x
   [
        p.τ * x2,
        p.τ * (ϵ * g(xd.u[1][1], p) - ϵ * (x1^2 - 1) * x2 - x1)
   ]
end

delaysF(par) = [1.0]

pars = (c1 = 1/4, c2 = 1/2, ϵ = 0.6, τ = 0.751)

prob = ConstantDDEBifProblem(vdp, delaysF, [1.01, 0.001], pars, (@optic _.ϵ); record_from_solution = (x,p; k...) -> (x1 = x[1], x2 = x[2]))

optn = NewtonPar(eigsolver = DDE_DefaultEig(tol = 1e-8))
opts = ContinuationPar(p_max = 2., p_min = 0.2, newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 2, plot = true, bothside = true, normC = norminf)
plot(br)

br2 = continuation(br, 2, bothside = true, normC = norminf)
plot(br,br2)
################################################################################
error("ca devrait trouver une BT!!")
brhopf = continuation(br2, 2, (@optic _.τ),
         ContinuationPar(br2.contparams; detect_bifurcation = 0, dsmax = 0.01, max_steps = 70, p_max = 0.9, p_min = 0.6,ds = 0.01, n_inversion = 2, tol_stability = 1e-6, nev = 5);
         verbosity = 1, plot = true,
         detect_codim2_bifurcation = 2,
         update_minaug_every_step = 1,
         bothside = true,
         start_with_eigen = false)

plot(brhopf, vars=(:ϵ, :τ))

brfold = continuation(br2, 4, (@optic _.τ),
         ContinuationPar(br2.contparams; detect_bifurcation = 0, dsmax = 0.01, max_steps = 70, p_max = 0.9, p_min = 0.1,ds = 0.01, n_inversion = 2, tol_stability = 1e-6, nev = 5);
         verbosity = 1, plot = true,
         detect_codim2_bifurcation = 2,
         update_minaug_every_step = 1,
         bothside = true,
         start_with_eigen = false)

plot(brfold, vars=(:ϵ, :τ)); plot!(brhopf, vars=(:ϵ, :τ), plotspecialpoints= false)

plot(brfold.param, brfold.BT)

plot(brfold)