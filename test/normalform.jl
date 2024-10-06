# using Revise

using DDEBifurcationKit, BifurcationKit
using Test
const BK = BifurcationKit

function wrightVF(x, xd, p)
    (;a) = p
    y = xd[1][1]
    [
        -a * y * (1 + x[1])
    ]
end

delaysF(par) = [1.0]

pars = (a=0.1, b=0.)
x0 = [0.]

prob = ConstantDDEBifProblem(wrightVF, delaysF, x0, pars, (@optic _.a))

optn = NewtonPar(eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 9., p_min = -1., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 6, n_inversion = 12 )
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

h1 = BifurcationKit.get_normal_form(br, 1)
@test isapprox(real(h1.nf.a), NFcoeff(0).a; rtol = 1e-5)
@test isapprox(real(h1.nf.b), NFcoeff(0).b; rtol = 1e-5)

h2 = BifurcationKit.get_normal_form(br, 2)
@test isapprox(real(h2.nf.a), NFcoeff(2).a; rtol = 1e-5)
@test isapprox(real(h2.nf.b), NFcoeff(2).b; rtol = 1e-5)

# static bifurcations
pars = (a = -0.1, b = 0.)
prob = ConstantDDEBifProblem(wrightVF, delaysF, x0, pars, (@optic _.a))
optn = NewtonPar(eigsolver = DDE_DefaultEig(γ = .01))
opts = ContinuationPar(p_max = 9., p_min = -1., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 7, n_inversion = 4 )
br = continuation(prob, PALC(), opts; verbosity = 0)
BifurcationKit.get_normal_form(br, 1)

BK.isinplace(prob)
BK.is_symmetric(prob)
BK._getvectortype(prob)
show(prob)

DDEBifurcationKit.newton_hopf(br, 2)
##########################################################################################
# SDDE, test dummy Hopf normal form
function humpriesVF(x, xd, p)
    (;κ1,κ2,γ,a1,a2,c) = p
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
prob = SDDDEBifProblem(humpriesVF, delaysF, x0, pars, (@optic _.κ1))
optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 13., p_min = 0., newton_options = optn, ds = -0.01, detect_bifurcation = 3, nev = 3, )
br = continuation(prob, PALC(), opts; verbosity = 0, bothside = true)
BifurcationKit.get_normal_form(br, 2)

BK.isinplace(prob)
BK.is_symmetric(prob)
BK._getvectortype(prob)
BK.getlens(prob)
BK.has_adjoint(prob)
BK.getdelta(prob)
show(prob)