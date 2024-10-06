cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."

using Revise, DDEBifurcationKit
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

function ikedaVF(x, xd, p)
   (; Λ) = p
   y = xd[1][1]
   [
      -pi/2 + Λ/2 * y^2;
   ]
end

delaysF(par) = [1.0]

pars = (Λ=0.1,b=0.)
x0 = [-sqrt(pi)]

prob = ConstantDDEBifProblem(ikedaVF, delaysF, x0, pars, (@optic _.Λ), record_from_solution=(x,p;k...)-> (x=x[1], _x=1))

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 2., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 4, n_inversion = 12 )
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = false)
plot(br)

BK.get_normal_form(br, 1) # l1= -0.0591623057, b = 0.09293196762669392
################################################################################
# computation periodic orbit

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.2, ds= -0.001, dsmin = 1e-4, p_max = 10., p_min=-5., max_steps = 40,
    nev = 3, tol_stability = 1e-8, detect_bifurcation = 0, plot_every_step = 1, save_sol_every_step=1)
@reset opts_po_cont.newton_options.tol = 1e-7
@reset opts_po_cont.newton_options.verbose = true

# arguments for periodic orbits
args_po = (    record_from_solution = (x, p; k...) -> begin
		xtt = BK.get_periodic_orbit(p.prob, x, nothing)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getperiod(p.prob, x, nothing))
	end,
	plot_solution = (x, p; k...) -> begin
		xtt = BK.get_periodic_orbit(p.prob, x, nothing)
		plot!(xtt.t, xtt[1,:]; label = "x", k...)
		plot!(br; subplot = 1, putspecialptlegend = false)
		end,
	normC = norminf)

probpo = PeriodicOrbitTrapProblem(M = 100, jacobian = :DenseAD, N = 1)
probpo = PeriodicOrbitOCollProblem(90, 4; N = 1, jacobian = BK.AutoDiffDense())
br_pocoll = @time continuation(
    br, 1, opts_po_cont,
    # we want to use the Collocation method to locate PO, with polynomial degree 5
    probpo;
    # regular continuation options
    verbosity = 2,    plot = true,
    args_po...,
    # ampfactor = 1/0.467829783456199 ,#* 0.014,
    δp = 0.01,
    callback_newton = (state; k...) -> begin
        xtt = BK.get_periodic_orbit(probpo,state.x,nothing)
        # plot(xtt.t, xtt[1,:], title = "it = $(state.it)") |> display
        printstyled(color=:red, "amp = ", BK.amplitude(xtt[:,:],1),"\n")
        @show state.x[end]
        state.step < 18
    end
    )

plot(br, br_pocoll)

plot(br_pocoll, vars = (:param, :period))

################################################################################
using DifferentialEquations

function ikedaVF_DE(du,u,h,p,t)
    (; Λ) = p
   du[1] = -pi/2 + Λ/2 * h(p,t-1)[1]^2;
end

u0 = -2ones(1)
h(p, t) = -2ones(1) .+ 0.001cos(t/4)
h(p,t) = br_pocoll.orbit(t)
prob_de = DDEProblem(ikedaVF_DE,u0,h,(0.,1240.),setproperties(pars, Λ = br.specialpoint[1].param + 0.01); constant_lags=delaysF(pars))
alg = MethodOfSteps(Rosenbrock23())
sol = solve(prob_de,alg)
plot(plot(sol, xlims = (1000,1020)), plot(sol))

