# using Revise
using Test, DDEBifurcationKit
using BifurcationKit
const BK = BifurcationKit

g(z) = (tanh(z − 1) + tanh(1))*cosh(1)^2
function neuron2VF(x, xd, p)
    (;a, b, c, d) = p
    [
        -x[1] - a * g(b*xd.u[1][1]) + c * g(d*xd.u[2][2]),
        -x[2] - a * g(b*xd.u[1][2]) + c * g(d*xd.u[2][1])
    ]
end

function delaysF(par)
    [par.τ1, par.τ2]
end

pars = (a = 0.25, b = 2., c = 15/29, d = 1.2, τ1 = 12.7, τ2 = 20.2)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuron2VF, delaysF, x0, pars, (@optic _.a))

#####
show(prob)
BK.isinplace(prob)
BK.is_symmetric(prob)
BK._getvectortype(prob)
#####

optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit=100))
opts = ContinuationPar(p_max = 0.4, p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 9, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts, bothside = false)

hpnf = BK.get_normal_form(br, 2)
################################################################################
brhopf = continuation(br, 1, (@optic _.c),
            setproperties(br.contparams, detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 1.1, p_min = -0.1,ds = 0.01, n_inversion = 2);
            verbosity = 0, plot = false,
            # update_minaug_every_step = 1,
            detect_codim2_bifurcation = 2,
            bothside = true,
            start_with_eigen = true)

BK.isinplace(brhopf.prob)
BK.is_symmetric(brhopf.prob)
BK.has_adjoint(brhopf.prob.prob)
BK.has_adjoint_MF(brhopf.prob.prob)
BK.getlens(brhopf.prob.prob)

brhopf2 = continuation(br, 2, (@optic _.c),
            setproperties(br.contparams, detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 1.1, p_min = -0.1,ds = -0.01);
            verbosity = 0, plot = false,
            detect_codim2_bifurcation = 2,
            bothside = true,
            start_with_eigen = true)

@test length(brhopf.specialpoint) == 11
@test brhopf.specialpoint[2].type == :hh
@test brhopf.specialpoint[3].type == :gh

@test length(brhopf2.specialpoint) == 9
@test brhopf2.specialpoint[2].type == :hh
@test brhopf2.specialpoint[3].type == :hh
################################################################################
prob2 = ConstantDDEBifProblem(neuron2VF, delaysF, x0, (BK.@set pars.a = 0.12), (@optic _.c))
br2 = continuation(prob2, PALC(), setproperties(opts, p_max = 1.22); verbosity = 0, plot = false, bothside = false)

# change tolerance for avoiding error computation of the EV
opts_fold = br.contparams
@reset opts_fold.newton_options.eigsolver.σ = 1e-7

brfold = continuation(br2, 3, (@optic _.a),
            setproperties(opts_fold; detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 0.6, p_min = -0.6,ds = -0.01, n_inversion = 2, tol_stability = 1e-6);
            verbosity = 0, plot = false,
            detect_codim2_bifurcation = 1,
            jacobian_ma = BK.AutoDiff(),
            bothside = false,
            start_with_eigen = true)

@test brfold.specialpoint[2].type == :zh
@test brfold.specialpoint[3].type == :zh
################################################################################
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
#####
show(prob)
BK.isinplace(prob)
BK.is_symmetric(prob)
BK._getvectortype(prob)
#####
################################################################################
# first Lyapunov coefficient along a Hopf curve, compared with the Hopf normal form
# computed at the same parameters (the delay Jacobians are not symmetric here)
function movingEqVF(x, xd, p)
    (;a, b, ρ) = p
    [
        x[2] - x[1] / (1 + x[1]),
        a * (ρ - x[2]) - b * (xd.u[1][1] - ρ / (1 - ρ))
    ]
end
movingEqDelays(par) = [par.τ0]
movingEqDelays(x, par) = [par.τ0 + par.χ * (1 + x[1])]

let
pars = (a = 0.5, b = 1., ρ = 0.5, τ0 = 0.3, χ = 0.1)
x0 = [pars.ρ / (1 - pars.ρ), pars.ρ]
optn = NewtonPar(eigsolver = DDE_DefaultEig(), tol = 1e-11)
opts = ContinuationPar(p_min = 0., p_max = 2., ds = 0.01, dsmax = 0.05, newton_options = optn, detect_bifurcation = 3, nev = 6, n_inversion = 8)

# Hopf point on the branch of equilibria continued in τ0
function hopf_point(Prob, par)
    prob = Prob(movingEqVF, movingEqDelays, x0, par, (@optic _.τ0))
    br = continuation(prob, PALC(), opts)
    return br, findfirst(sp -> sp.type == :hopf, br.specialpoint)
end

for (Prob, lens2, p_min, p_max) in ((ConstantDDEBifProblem, (@optic _.b), 0.5, 1.),
                                    (SDDDEBifProblem,       (@optic _.χ), 0.,  0.35))
    br, ind = hopf_point(Prob, pars)
    brhopf = continuation(br, ind, lens2,
                setproperties(opts, p_min = p_min, p_max = p_max, ds = -0.01, dsmax = 0.02, max_steps = 8, detect_bifurcation = 1);
                detect_codim2_bifurcation = 2,
                update_minaug_every_step = 1,
                start_with_eigen = true)
    for pt in (brhopf.branch[1], brhopf.branch[end])
        hp = BK.get_normal_form(hopf_point(Prob, BK.set(pars, lens2, lens2(pt)))...)
        @test pt.τ0 ≈ hp.p rtol = 1e-4
        @test real(pt.l1) ≈ real(hp.nf.b) rtol = 1e-4
    end
end
end
