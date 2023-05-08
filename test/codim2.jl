# using Revise
using Test, DDEBifurcationKit
using Parameters, Setfield
using BifurcationKit
const BK = BifurcationKit

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

#####
show(prob)
BK.isInplace(prob)
BK.isSymmetric(prob)
BK.getVectorType(prob)
#####

optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit=100))
opts = ContinuationPar(pMax = 0.4, pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 9, dsmax = 0.2, nInversion = 4)
br = continuation(prob, PALC(), opts, bothside = false)

hpnf = BK.getNormalForm(br, 1)
################################################################################
brhopf = continuation(br, 1, (@lens _.c),
			setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.01, maxSteps = 100, pMax = 1.1, pMin = -0.1,ds = 0.01, nInversion = 2);
			verbosity = 0, plot = false,
			detectCodim2Bifurcation = 2,
			bothside = true,
			startWithEigen = true)

brhopf2 = continuation(br, 2, (@lens _.c),
			setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.01, maxSteps = 100, pMax = 1.1, pMin = -0.1,ds = -0.01);
			verbosity = 0, plot = false,
			detectCodim2Bifurcation = 2,
			bothside = true,
			startWithEigen = true)

@test length(brhopf.specialpoint) == 10
@test brhopf.specialpoint[2].type == :hh
@test brhopf.specialpoint[3].type == :gh

@test length(brhopf2.specialpoint) == 10
@test brhopf2.specialpoint[2].type == :hh
@test brhopf2.specialpoint[3].type == :hh
################################################################################
prob2 = ConstantDDEBifProblem(neuron2VF, delaysF, x0, (@set pars.a = 0.12), (@lens _.c))
br2 = continuation(prob2, PALC(), setproperties(opts, pMax = 1.22); verbosity = 0, plot = false, bothside = false)

# change tolerance for avoiding error computation of the EV
opts_fold = br.contparams
@set! opts_fold.newtonOptions.eigsolver.σ = 1e-7

brfold = continuation(br2, 3, (@lens _.a),
			setproperties(opts_fold; detectBifurcation = 1, dsmax = 0.01, maxSteps = 100, pMax = 0.6, pMin = -0.6,ds = -0.01, nInversion = 2, tolStability = 1e-6);
			verbosity = 0, plot = false,
			detectCodim2Bifurcation = 2,
			bothside = false,
			startWithEigen = true)

@test brfold.specialpoint[2].type == :zh
@test brfold.specialpoint[3].type == :zh
################################################################################
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
#####
show(prob)
BK.isInplace(prob)
BK.isSymmetric(prob)
BK.getVectorType(prob)
#####