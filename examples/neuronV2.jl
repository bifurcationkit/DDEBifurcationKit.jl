cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."
using Revise, DDEBifurcationKit, Parameters, Setfield, LinearAlgebra, Plots
using BifurcationKit
const BK = BifurcationKit

g(z) = (tanh(z − 1) + tanh(1))*cosh(1)^2
norminf(x) = norm(x, Inf)
function neuron2VF(x, xd, p)
   @unpack a,b,c,d = p
   [
      -x[1] - a * g(b*xd[1][1]) + c * g(d*xd[2][2]),
      -x[2] - a * g(b*xd[1][2]) + c * g(d*xd[2][1])
   ]
end

function delaysF(par)
   [par.τ1, par.τ2]
end

pars = (a = 0.25, b = 2., c = 15/29, d = 1.2, τ1 = 12.7, τ2 = 20.2)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuron2VF, delaysF, x0, pars, (@lens _.a))

optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit=100))
opts = ContinuationPar(pMax = 0.4, pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 9, dsmax = 0.2, nInversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 0, plot = true, bothside = false)

plot(br)

hpnf = BK.getNormalForm(br, 1)
################################################################################
brhopf = continuation(br, 1, (@lens _.c),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.01, maxSteps = 100, pMax = 1.1, pMin = -0.1,ds = 0.01, nInversion = 2);
         verbosity = 0, plot = true,
         detectCodim2Bifurcation = 2,
         bothside = true,
         startWithEigen = true)

brhopf2 = continuation(br, 2, (@lens _.c),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.01, maxSteps = 100, pMax = 1.1, pMin = -0.1,ds = -0.01);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 2,
         bothside = true,
         startWithEigen = true)

plot(brhopf, vars = (:a, :c), xlims = (0,0.7), ylims = (0,1))
   plot!(brhopf2, vars = (:a, :c), xlims = (-0,0.7), ylims = (-0.1,1))

################################################################################
prob2 = ConstantDDEBifProblem(neuron2VF, delaysF, x0, (@set pars.a = 0.12), (@lens _.c))
br2 = continuation(prob2, PALC(), setproperties(opts, pMax = 1.22); verbosity = 1, plot = true, bothside = false)

plot(br2)

# change tolerance for avoiding error computation of the EV
opts_fold = br.contparams
@set! opts_fold.newtonOptions.eigsolver.σ = 1e-7

brfold = continuation(br2, 3, (@lens _.a),
         setproperties(opts_fold; detectBifurcation = 1, dsmax = 0.01, maxSteps = 100, pMax = 0.6, pMin = -0.5,ds = -0.01, nInversion = 2, tolStability = 1e-6);
         verbosity = 1, plot = true,
         detectCodim2Bifurcation = 2,
         bothside = false,
         startWithEigen = true)

plot(brfold)

plot(brfold, vars = (:a, :c), branchlabel = "Fold")
   plot!(brhopf, vars = (:a, :c), branchlabel = "Hopf")
   plot!(brhopf2, vars = (:a, :c), branchlabel = "Hopf")
################################################################################
# computation periodic orbit

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, pMax = 10., pMin=-5., maxSteps = 30, nev = 3, tolStability = 1e-8, detectBifurcation = 0, plotEveryStep = 2, saveSolEveryStep = 1)
	@set! opts_po_cont.newtonOptions.tol = 1e-9
	@set! opts_po_cont.newtonOptions.verbose = true
	@set! opts_po_cont.newtonOptions.maxIter = 8

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
		plot!(br; subplot = 1, putspecialptlegend = false)
		end,
	normC = norminf)

probpo = PeriodicOrbitOCollProblem(100, 3; N = 2)
# probpo = PeriodicOrbitTrapProblem(M = 2000, jacobian = :DenseAD, N = 2)
br_pocoll = @time continuation(
	br, 1, opts_po_cont,
	# PeriodicOrbitOCollProblem(100, 4);
	probpo;
	verbosity = 2,	plot = true,
	args_po...,
	# ampfactor = 1/0.24391300209895822 * 0.1,
	ampfactor = 1.42,
	δp = 0.001,
	normC = norminf,
	callbackN = (state; k...) -> begin
		xtt = BK.getPeriodicOrbit(probpo,state.x,nothing)
		# plot(xtt.t, xtt[1,:], title = "it = $(state.it)") |> display
		printstyled(color=:red, "amp = ", BK.amplitude(xtt[:,:],1),"\n")
		# @show state.x[end]
		# @show state.f[end]
		state.step < 16
	end
	)
################################################################################
using  DifferentialEquations

function neuronV2_DE(du,x,h,p,t)
	@unpack a,b,c,d,τ1,τ2 = p
   du[1] = -x[1] - a * g(b*h(p, t-τ1)[1]) + c * g(d*h(p, t-τ2)[2])
   du[2] = -x[2] - a * g(b*h(p, t-τ1)[2]) + c * g(d*h(p, t-τ2)[1])
end

u0 = -2ones(2)
	h(p, t) = -0*ones(2) .+ 0.01cos(t/4)
	# h(p,t) = br_pocoll.orbit(t)
	prob_de = DDEProblem(neuronV2_DE,h(pars,0),h,(0.,54240.),setproperties(pars, a = br.specialpoint[1].param + 0.001); constant_lags=delaysF(pars))
	alg = MethodOfSteps(Rosenbrock23())
	sol = solve(prob_de,alg)
	plot(plot(sol, xlims = (sol.t[end]-100,sol.t[end])), plot(sol))

