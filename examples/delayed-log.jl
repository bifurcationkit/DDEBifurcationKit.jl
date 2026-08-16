cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."
using Revise, DDEBifurcationKit, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

function delayedlogVF(x, xd, p)
   (;λ) = p
   y = xd.u[1][1]
   [
      (λ - y) * x[1]
   ]
end

function delaysF(par)
   [
      1.0,
   ]
end

pars = (λ=1.1,b=0.)
x0 = [1.]

prob = ConstantDDEBifProblem(delayedlogVF, delaysF, x0, pars, (@optic _.λ), record_from_solution=(x,p;k...)-> (x=x[1], _x=1))

optn = NewtonPar(eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 2., p_min = 0., newton_options = optn, ds = 0.01, nev = 6, n_inversion = 6 )
br = BK.continuation(prob, PALC(), opts)
plot(br)
################################################################################
# computation periodic orbit
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= 0.001, dsmin = 1e-4, p_max = 10., p_min=-5., max_steps = 130,
    nev = 5, tol_stability = 1e-8, detect_bifurcation = 2, plot_every_step = 2)
@reset opts_po_cont.newton_options.tol = 1e-8
@reset opts_po_cont.newton_options.verbose = true

# arguments for periodic orbits
args_po = (    record_from_solution = (x, p;k...) -> begin
		xtt = BK.get_periodic_orbit(p.prob, x, nothing)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getperiod(p.prob, x, nothing))
	end,
	plot_solution = (x, p; k...) -> begin
		xtt = BK.get_periodic_orbit(p.prob, x, nothing)
		plot!(xtt.t, xtt[1,:]; label = "V1", k...)
		plot!(br; subplot = 1, putspecialptlegend = false)
		end,
	normC = norminf)

probpo = Collocation(20, 4; N = 1, jacobian = BifurcationKit.DenseAnalyticalInplace())
br_pocoll = @time continuation(
		br, 1, opts_po_cont,
		probpo;
		verbosity = 2,	plot = true,
		eigsolver = BK.FloquetGEV(DDE_DefaultEig(maxit=200, tol = 1e-10, σ = 1e-4), length(probpo), 1),
		args_po...,
		ampfactor = 2,
		δp = 0.01,
		normC = norminf,
		)

#####
################################################################################
using DifferentialEquations

function delayedlog_DE(du,x,h,p,t)
    @unpack λ = p
    du[1] = (λ - h(p,t-1)[1]) * x[1]
end

h(p, t) = -1*zeros(1) .+ 1. *cos(t/4)
prob_de = DDEProblem(delayedlog_DE,h(pars, 0),h,(0.,1320.),setproperties(pars, λ = pi/2 + .121); constant_lags=delaysF(pars))
alg = MethodOfSteps(BS3())
sol = solve(prob_de,alg; progress = true)
plot(plot(sol, xlims = (sol.t[end]-100,sol.t[end])), plot(sol),title = "p=$(sol.prob.p.λ)")

