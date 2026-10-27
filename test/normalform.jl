# using Revise

using DDEBifurcationKit, BifurcationKit
using Test
const BK = BifurcationKit

function wrightVF(x, xd, p)
    (;a) = p
    y = xd.u[1][1]
    [
        -a * y * (1 + x[1])
    ]
end

delaysF(par) = [1.0]

let
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
end
##########################################################################################
# SDDE, test dummy Hopf normal form
function humpriesVF(x, xd, p)
    (;κ1,κ2,γ,a1,a2,c) = p
    [
       -γ * x[1] - κ1 * xd.u[1][1] - κ2 * xd.u[2][1]
    ]
end

function delaysF(x, par)
[
    par.a1 + par.c * x[1],
    par.a2 + par.c * x[1],
]
end
let
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
end
##########################################################################################
# Hopf normal form, coefficient a = dμ/dp when the parameter moves the equilibrium (R01 ≠ 0)
# The equilibrium x* = (ρ/(1-ρ), ρ) depends on ρ but J(x, p) does not depend on ρ at fixed x.
function movingEqVF(x, xd, p)
    (;a, b, ρ) = p
    [
        x[2] - x[1] / (1 + x[1]),
        a * (ρ - x[2]) - b * (xd.u[1][1] - ρ / (1 - ρ))
    ]
end
movingEqDelays(par) = [par.τ0]
movingEqDelays(x, par) = [par.τ0 + par.χ * (1 + x[1])]

# characteristic function f(μ) = (μ + c)(μ + a) + b⋅exp(-μτ) with c = (1-ρ)² and τ = τ0 + χ/(1-ρ)
# at the equilibrium, hence dμ/dρ = -∂ρf / ∂μf
function movingEq_dμdρ(par, μ)
    (;a, b, ρ, τ0, χ) = par
    c = (1 - ρ)^2
    τ = τ0 + χ / (1 - ρ)
    ∂μf = 2μ + a + c - b * τ * exp(-μ * τ)
    ∂ρf = -2(1 - ρ) * (μ + a) - b * μ * exp(-μ * τ) * χ / (1 - ρ)^2
    return -∂ρf / ∂μf
end

let
optn = NewtonPar(eigsolver = DDE_DefaultEig(), tol = 1e-11)
opts = ContinuationPar(p_min = 0.3, p_max = 0.75, ds = 0.002, dsmax = 0.01, newton_options = optn, detect_bifurcation = 3, nev = 6, n_inversion = 8)
for (pars, Prob) in (((a = 0.1, b = 0.1, ρ = 0.40, τ0 = 5.0, χ = 0.0), ConstantDDEBifProblem),
                     ((a = 0.1, b = 0.1, ρ = 0.55, τ0 = 0.0, χ = 1.0), SDDDEBifProblem))
    x0 = [pars.ρ / (1 - pars.ρ), pars.ρ]
    prob = Prob(movingEqVF, movingEqDelays, x0, pars, (@optic _.ρ))
    br = continuation(prob, PALC(), opts)
    hp = BifurcationKit.get_normal_form(br, findfirst(sp -> sp.type == :hopf, br.specialpoint))
    @test isapprox(hp.nf.a, movingEq_dμdρ(merge(pars, (ρ = hp.p,)), Complex(0, hp.ω)); rtol = 1e-4)
end
end
