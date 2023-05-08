# using Revise
using DDEBifurcationKit, Parameters, Setfield, BifurcationKit
using Test
const BK = BifurcationKit

function wrightVF(x, xd, p)
	@unpack a = p
	y = xd[1][1]
	[
		-a * y * (1 + x[1])
	]
end

delaysF(par) = [1.0]

pars = (a=0.1,b=0.)
x0 = [0.]

prob = ConstantDDEBifProblem(wrightVF, delaysF, x0, pars, (@lens _.a))

optn = NewtonPar(eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 9., pMin = -1., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 6, nInversion = 12 )
br = continuation(prob, PALC(), opts)
# hopf bifurcations
function NFcoeff(N)
	# Faria 2006
	ωN = pi/2 + N*pi
	aN = ωN*(-1)^N
	K1 = aN / (1+aN^2)
	K2 = ωN / (5(1+aN^2))*((-1)^N-3ωN)
	return (a=K1, b=K2)
end

h1 = BifurcationKit.getNormalForm(br, 1)
@test isapprox(real(h1.nf.a), NFcoeff(0).a; rtol = 1e-5)
@test isapprox(real(h1.nf.b), NFcoeff(0).b; rtol = 1e-5)

h2 = BifurcationKit.getNormalForm(br, 2)
@test isapprox(real(h2.nf.a), NFcoeff(2).a; rtol = 1e-5)
@test isapprox(real(h2.nf.b), NFcoeff(2).b; rtol = 1e-5)

# static bifurcations
pars = (a=-0.1,b=0.)
prob = ConstantDDEBifProblem(wrightVF, delaysF, x0, pars, (@lens _.a))
optn = NewtonPar(eigsolver = DDE_DefaultEig(γ = .01))
opts = ContinuationPar(pMax = 9., pMin = -1., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 7, nInversion = 4 )
br = continuation(prob, PALC(), opts; verbosity = 0)
BifurcationKit.getNormalForm(br, 1)

BK.isInplace(prob)
BK.isSymmetric(prob)
BK.getVectorType(prob)
show(prob)
##########################################################################################
# SDDE, test dummy Hopf normal form
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
br = continuation(prob, PALC(), opts; verbosity = 0, bothside = true)
BifurcationKit.getNormalForm(br, 2)

BK.isInplace(prob)
BK.isSymmetric(prob)
BK.getVectorType(prob)
BK.getLens(prob)
BK.hasAdjoint(prob)
BK.getDelta(prob)
show(prob)