# using Revise, Plots
using DDEBifurcationKit, Parameters, Setfield, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

function humpriesVF(x, xd, p)
   @unpack κ1,κ2,γ,a1,a2,c = p
   [
      -γ * x[1] - κ1 * xd[1][1] - κ2 * xd[2][1]
   ]
end

function delaysF(x, par)
   [
      par.a1 + par.c * x[1],
      par.a2 + par.c * x[1],
   ]
end


pars = (κ1=0.,κ2=2.3,a1=1.3,a2=6,γ=4.75,c=1.)
x0 = zeros(1)

prob = SDDDEBifProblem(humpriesVF, delaysF, x0, pars, (@lens _.κ1))

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 13., pMin = 0., newtonOptions = optn, ds = -0.01, detectBifurcation = 3, nev = 3, )

alg = PALC()
br = continuation(prob, alg, opts; verbosity = 0, plot = false, bothside = true)

# plot(br)
################################################################################
# computation of periodic orbit
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, pMax = 12., pMin=-5., maxSteps = 3,
	nev = 3, tolStability = 1e-8, detectBifurcation = 0, plotEveryStep = 20, saveSolEveryStep=1)
	@set! opts_po_cont.newtonOptions.tol = 1e-9
	@set! opts_po_cont.newtonOptions.verbose = true

	# arguments for periodic orbits
	args_po = (	recordFromSolution = (x, p) -> begin
			xtt = BK.getPeriodicOrbit(p.prob, x, nothing)
			_max = maximum(xtt[1,:])
			_min = minimum(xtt[1,:])
			return (amp = _max - _min,
					max = _max,
					min = _min,
					period = getPeriod(p.prob, x, nothing))
		end,
		plotSolution = (x, p; k...) -> begin
			xtt = BK.getPeriodicOrbit(p.prob, x, nothing)
			plot!(xtt.t, xtt[1,:]; label = "x", k...)
			plot!(br; subplot = 1, putspecialptlegend = false)
			end,
		normC = norminf)

probpo = PeriodicOrbitOCollProblem(200, 2; N = 1)
	br_pocoll = @time continuation(
	br, 2, opts_po_cont,
	probpo;
	alg = PALC(tangent = Bordered()),
	# regular continuation options
	# verbosity = 2,	plot = true,
	args_po...,
	ampfactor = 1/0.467829783456199 * 0.1,
	δp = 0.01,
	callbackN = (state; k...) -> begin
		xtt = BK.getPeriodicOrbit(probpo,state.x,nothing)
		# plot(xtt.t, xtt[1,:], title = "it = $(state.it)") |> display
		printstyled(color=:red, "amp = ", BK.amplitude(xtt[:,:],1),"\n")
		printstyled(color=:green, "T = ", (state.x[end]),"\n")
		@show state.x[end]
		state.step < 15 && BK.cbMaxNorm(10.0)(state; k...)
	end
	)
