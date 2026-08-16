cd(@__DIR__)
using Pkg, LinearAlgebra, Test
pkg"activate ."

# Humphries et al ( A. R. Humphries, O. A. DeMasi, F. M. G. Magpantay, F. Upham (2012), Dynamics of a delay differential equation with multiple state-dependent delays, Discrete and Continuous Dynamical Systems 32(8) pp. 2701-2727 http://dx.doi.org/10.3934/dcds.2012.32.2701)
using Revise, DDEBifurcationKit, Plots
using BifurcationKit
const DDEBK = DDEBifurcationKit

function humpriesVF(x, xd, p)
   (; κ1, κ2, γ) = p
   [
      -γ * x[1] - κ1 * xd.u[1][1] - κ2 * xd.u[2][1]
   ]
end

function delaysF(x, par)
   [
      par.a1 + par.c * x[1],
      par.a2 + par.c * x[1],
   ]
end

function plotSolution(x, p; k...)
    xtt = DDEBK.get_periodic_orbit(p.prob, x, nothing)
    plot!(xtt.t, xtt[1,:]; label = "x", marker = :d, markersize = 1, k...)
    plot!(br; subplot = 1, putspecialptlegend = false)    
end


pars = (κ1=0., κ2=2.3, a1=1.3, a2=6, γ=4.75, c=1.)
prob = SDDDEBifProblem(humpriesVF, delaysF, zeros(1), pars, (@optic _.κ1))
optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 13., p_min = 0., newton_options = optn, nev = 25, )
br = continuation(prob, PALC(), opts)
plot(br)

get_normal_form(br, 2)
################################################################################
brh = continuation(BifurcationKit.re_make(prob; params = @set pars.κ2=3.), PALC(), opts; plot = true)
brhopf = [continuation(brh, ii, (@optic _.κ2),
         ContinuationPar(br.contparams, detect_bifurcation = 3, dsmax = 0.04, max_steps = 230, p_max = 5., p_min = -1.,ds = -0.02);
        #  verbosity = 2, 
         plot = true,
         detect_codim2_bifurcation = 2,
         bothside = true,
         start_with_eigen = true) for ii in eachindex(brh.specialpoint) if (brh.specialpoint[ii].type == :hopf) ]

plot(brhopf..., vars = (:κ1, :κ2))
################################################################################
# computation periodic orbit
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, p_max = 12., p_min=-5., max_steps = 3000,
nev = 20, tol_stability = 1e-3, plot_every_step = 4, newton_options = NewtonPar(tol = 1e-10, verbose = true))

# arguments for periodic orbits
args_po = (    
    plot_solution = plotSolution,
    normC = norminf)

probpo = Collocation(200, 2; N = 1, jacobian = DDEBK.BifurcationKit.AutoDiffDense())
br_pocoll = @time continuation(
    br, 1, ContinuationPar(opts_po_cont; detect_bifurcation = 2),
    probpo;
    alg = PALC(tangent = Bordered()),
    # regular continuation options
    verbosity = 1, plot = true,
    args_po...,
    δp = 0.01,
    )

plot(br);plot!(br_pocoll, plotfold=false, ylabel = "amplitude")
################################################################################
import DelayDiffEq as DDE

function humpriesVF_DE2(x,h,p,t)
    (;κ1,κ2,γ,a1,a2,c) = p
   -γ * x - κ1 * h(p, t-(a1 + c * x)) - κ2 * h(p, t-(a2 + c * x))
end

function h0(p, t)
     t ≤ 0 || error("history function is only implemented for t ≤ 0")
     0 .+ 0.03sin(t)
end
prob_de = DDE.DDEProblem(humpriesVF_DE2,h0,(0.,10200.), (DDEBK.@set pars.κ1 = br.specialpoint[2].param + 0.01); dependent_lags=((x,par,t)->par.a1 + par.c * x, (x,par,t)->par.a2 + par.c * x))
alg = DDE.MethodOfSteps(DDE.Rosenbrock23())
sol = DDE.solve(prob_de,alg)
plot(plot(sol, xlims = (sol.t[end]-30,sol.t[end])), plot(sol))

