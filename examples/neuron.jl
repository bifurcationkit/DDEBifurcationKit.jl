using Revise, DDEBifurcationKit, Parameters, Setfield, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

using Plots

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

plot(br)
################################################################################
prob2 = ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@lens _.a21))
br2 = BK.continuation(prob2, PALC(), setproperties(opts, ds = 0.1, pMax = 3., nInversion=8); verbosity = 3, plot = true, bothside = false, normC = norminf)

# @set! br2.contparams.newtonOptions.eigsolver.σ = 1e-5
BK.getNormalForm(br2, 2)
#Hopf l1 ≈ −0.0601.
################################################################################
brhopf = continuation(br, 3, (@lens _.a21),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.04, maxSteps = 230, pMax = 15., pMin = -1.,ds = -0.02);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 2,
         # bothside = true,
         startWithEigen = true)

plot(brhopf, vars = (:a21, :τs))
plot(brhopf, vars = (:τs, :ω))

brhopf2 = continuation(br, 2, (@lens _.a21),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.1, maxSteps = 56, pMax = 15., pMin = -1.,ds = -0.01, nInversion = 4);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 2,
         startWithEigen = true,
         bothside=true)

plot(brhopf, brhopf2, legend = :top)


################################################################################
# computation periodic orbit

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, pMax = 10., pMin=-5., maxSteps = 130,
	nev = 3, tolStability = 1e-8, detectBifurcation = 0, plotEveryStep = 2, saveSolEveryStep=1)
	@set! opts_po_cont.newtonOptions.tol = 1e-8
	@set! opts_po_cont.newtonOptions.verbose = true

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
	# probpo = PeriodicOrbitTrapProblem(M = 2000, jacobian = :DenseAD, N = 2)
	br_pocoll = @time continuation(
		br2, 1, opts_po_cont,
		probpo;
		verbosity = 2,	plot = true,
		args_po...,
		# ampfactor = 1/0.3593 * 0.0610*2.2,
		δp = 0.001,
		normC = norminf,
		)

plot(br2, br_pocoll)
plot(br_pocoll, vars = (:param, :period))

# plot the periodic orbit
plot(layout = 2)
	for ii = 1:10:110
		solpo = BK.getPeriodicOrbit(br_pocoll.γ.prob.prob, br_pocoll.sol[ii].x, nothing)
		plot!(solpo.t ./ solpo.t[end], solpo.u[1,:], label = "", subplot = 1)
	end
	xlabel!("t / period", subplot = 1)
	plot!(br_pocoll, vars = (:param, :period), subplot = 2, xlims=(2.2,2.4))

################################################################################
using  DifferentialEquations

function neuron_DE(du,u,h,p,t)
	@unpack κ, β, a12, a21, τs, τ1, τ2 = p
	du[1] = -κ * u[1] + β * tanh(h(p, t-τs)[1]) + a12 * tanh(h(p, t-τ2)[2])
	du[2] = -κ * u[2] + β * tanh(h(p, t-τs)[2]) + a21 * tanh(h(p, t-τ1)[1])
end

h(p, t) = -0*ones(2) .+ 0.25sin(t/4)
	prob_de = DDEProblem(neuron_DE,h(pars, 0),h,(0.,20000.),setproperties(pars, a21 = br2.specialpoint[3].param + 0.001); constant_lags=delaysF(pars))
	alg = MethodOfSteps(BS3())
	sol = solve(prob_de,alg)
	plot(plot(sol, xlims = (sol.t[end]-100,sol.t[end])), plot(sol),title = "a21=$(sol.prob.p.a21)")

