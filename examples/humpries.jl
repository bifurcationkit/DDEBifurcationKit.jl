cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."

# https://ddebiftool.sourceforge.net/demos/neuron/html/demo1_stst.html
using Revise, DDEBifurcationKit, Parameters, Setfield, LinearAlgebra, Plots
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

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 13., pMin = 0., newtonOptions = optn, ds = -0.01, detectBifurcation = 3, nev = 3, )
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = true)

plot(br)
################################################################################
brhopf = continuation(br, 2, (@lens _.κ2),
         setproperties(br.contparams, detectBifurcation = 2, dsmax = 0.04, maxSteps = 230, pMax = 5., pMin = -1.,ds = -0.02);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 0,
         bothside = true,
         startWithEigen = true)

plot(brhopf, vars = (:κ1, :κ2))
################################################################################
# computation periodic orbit

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, pMax = 12., pMin=-5., maxSteps = 3000,
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
	verbosity = 2,	plot = true,
	args_po...,
	ampfactor = 1/0.467829783456199 * 0.1,
	δp = 0.01,
	callbackN = (state; k...) -> begin
		xtt = BK.getPeriodicOrbit(probpo,state.x,nothing)
		# plot(xtt.t, xtt[1,:], title = "it = $(state.it)") |> display
		printstyled(color=:red, "amp = ", BK.amplitude(xtt[:,:],1),"\n")
		printstyled(color=:green, "T = ", (state.x[end]),"\n")
		@show state.x[end]
		state.it < 15 && BK.cbMaxNorm(10.0)(state; k...)
	end
	)

plot(br);plot!(br_pocoll, plotfold=false, ylabel = "amplitude")
################################################################################
using  DifferentialEquations

function humpriesVF_DE2(x,h,p,t)
	@unpack κ1,κ2,γ,a1,a2,c = p
   -γ * x - κ1 * h(p, t-(a1 + c * x)) - κ2 * h(p, t-(a2 + c * x))
end

function h0(p, t)
	 t ≤ 0 || error("history function is only implemented for t ≤ 0")
	 0 .+ 0.03sin(t)
 end
	prob_de = DDEProblem(humpriesVF_DE2,h0,(0.,10200.),setproperties(pars, κ1 = br.specialpoint[2].param + 0.01); dependent_lags=((x,par,t)->par.a1 + par.c * x, (x,par,t)->par.a2 + par.c * x))
	alg = MethodOfSteps(Rosenbrock23())
	sol = solve(prob_de,alg)
	plot(plot(sol, xlims = (sol.t[end]-30,sol.t[end])), plot(sol))

