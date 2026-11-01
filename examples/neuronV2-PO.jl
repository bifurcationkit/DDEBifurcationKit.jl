cd(@__DIR__)
using Pkg, LinearAlgebra, Test
pkg"activate ."
using Revise, BifurcationKit, DDEBifurcationKit, Plots
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit
# https://www.dropbox.com/scl/fi/7ndu8tideuskogzdj1l5z/DDE_seminar_tutorial_IV.pdf?rlkey=7h789sjohazen8kswaj9uhg8b&e=1&dl=0
g(z) = (tanh(z − 1) + tanh(1)) * cosh(1)^2

function neuron2VF(x, xd, p)
   (; a,b,c,d) = p
   [
      -x[1] - a * g(b * xd.u[1][1]) + c * g(d * xd.u[2][2]),
      -x[2] - a * g(b * xd.u[1][2]) + c * g(d * xd.u[2][1])
   ]
end

delaysF(par) = [par.τ1, par.τ2]

pars = (a = 0.069, b = 2., c = 0.6, d = 1.2, τ1 = 11.6, τ2 = 20.3)
prob = ConstantDDEBifProblem(neuron2VF, delaysF, zeros(2), pars, (@optic _.c))

optn = NewtonPar(eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 1.0, p_min = 0.25, newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 9, dsmax = 0.05, n_inversion = 8)
br = continuation(prob, PALC(), opts; bothside = Val(true))
plot(br)

hpnf = BK.get_normal_form(br, 2)
################################################################################
# computation periodic orbit
@views function plot_solution(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, nothing)
    plot!(xtt.t, xtt[1,:]; label = "V1", marker = :d, markersize = 2, k...)
    plot!(xtt.t, xtt[2,:]; label = "V2", k...)
    # plot!(br; subplot = 1, putspecialptlegend = false)
end

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.02, ds = -0.0001, p_max = 1.2, p_min = 0.100661, max_steps = 200, nev = 10, tol_stability = 1e-4, plot_every_step = 10, n_inversion = 6)

# arguments for periodic orbits
args_po = (  ;
    plot_solution,
    normC = norminf)

probpo = Collocation(30, 5, 
            # jacobian = BK.AutoDiffDense(),
            # jacobian = BK.DenseAnalytical(),
            jacobian = BK.DenseAnalyticalInplace(),
            # meshadapt = true,
            )
br_pocoll = @time continuation(
            br, 2, ContinuationPar(opts_po_cont; max_steps = 200, detect_bifurcation = 3),
            probpo;
            # verbosity = 2,
            # plot = true,
            args_po...,
            )

plot(br, br_pocoll)

# ┌────┬────────────┬─────────────┬────────────────────────────┬───────────┬─────────────┬──────────┬──────┐
# │ Nb │    type    │      c      │          interval          │ precision │   status    │    δ     │ step │
# ├────┼────────────┼─────────────┼────────────────────────────┼───────────┼─────────────┼──────────┼──────┤
# │  1 │         ns │ +0.70997815 │ (+0.70997815, +0.70997822) │     7e-08 │   converged │ ( 2,  2) │   18 │
# │  2 │         pd │ +0.65042660 │ (+0.65042660, +0.65049435) │     7e-05 │   converged │ (-1, -1) │   22 │
# │  3 │         bp │ +0.46191911 │ (+0.46191911, +0.46191911) │     2e-09 │   converged │ (-1,  0) │   43 │
# │  4 │         pd │ +0.46533676 │ (+0.46533668, +0.46533676) │     8e-08 │   converged │ (-1, -1) │   49 │
# │  5 │         pd │ +0.59672929 │ (+0.59672926, +0.59672929) │     3e-08 │   converged │ ( 1,  1) │   91 │
# │  6 │         bp │ +0.61515065 │ (+0.61515058, +0.61515065) │     7e-08 │   converged │ ( 1,  0) │  112 │
# │  7 │         pd │ +0.52173920 │ (+0.52173920, +0.52173924) │     4e-08 │   converged │ ( 1,  1) │  159 │
# │  8 │         pd │ +0.52174009 │ (+0.52173776, +0.52174009) │     2e-06 │   converged │ (-1, -1) │  163 │
# │  9 │   endpoint │ +0.61314530 │                            │           │             │          │  201 │
# └────┴────────────┴─────────────┴────────────────────────────┴───────────┴─────────────┴──────────┴──────┘

BK.get_normal_form(deepcopy(br_pocoll), 2; detailed = Val(false))

br_pd = @time continuation(
            deepcopy(br_pocoll), 2, ContinuationPar(BK.getcontparams(br_pocoll); max_steps = 12, detect_bifurcation = 0, ds = 0.001,);
            verbosity = 2,
            plot = true,
            args_po...,
            normal_form_kwargs = (detailed = Val(false),),
            # δp = 0.001, ampfactor = 0.25, # for #4
            # δp = 0.001, ampfactor = 0.25, # for #7
            δp = 0.005, ampfactor = 0.125, bothside = Val(true), # for #2
            )

plot(br_pocoll);plot!(br_pd, color = :red)
plot(br);plot!(br_pocoll, br_pd)
