# using Revise
using Test, DDEBifurcationKit
using BifurcationKit
const BK = BifurcationKit

g(z) = (tanh(z − 1) + tanh(1))*cosh(1)^2
function neuron2VF(x, xd, p)
    (;a, b, c, d) = p
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

@test length(brhopf.specialpoint) == 10
@test brhopf.specialpoint[2].type == :hh
@test brhopf.specialpoint[3].type == :gh

@test length(brhopf2.specialpoint) == 10
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