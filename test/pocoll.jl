using Test, DDEBifurcationKit
using Parameters, Setfield, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)
g(z) = (tanh(z − 1) + tanh(1))*cosh(1)^2
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
br = continuation(prob, PALC(), opts, bothside = false)

#######################################################
# hopf aBS
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, pMax = 10., pMin=-5., maxSteps = 5, nev = 3, tolStability = 1e-8, detectBifurcation = 0, plotEveryStep = 2, saveSolEveryStep = 1)
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
		# verbosity = 2,	plot = true,
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