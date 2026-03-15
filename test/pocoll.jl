using Test, DDEBifurcationKit
using LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

function neuron2VF(x, xd, p)
    (;a,b,c,d) = p
    [
        -x[1] - a * g(b*xd.u[1][1]) + c * g(d*xd.u[2][2]),
        -x[2] - a * g(b*xd.u[1][2]) + c * g(d*xd.u[2][1])
    ]
end
g(z) = (tanh(z − 1) + tanh(1))*cosh(1)^2
delaysF(par) = [par.τ1, par.τ2]

pars = (a = 0.25, b = 2., c = 15/29, d = 1.2, τ1 = 12.7, τ2 = 20.2)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuron2VF, delaysF, x0, pars, (@optic _.a))
optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit = 100))
opts = ContinuationPar(p_max = 0.4, p_min = 0., newton_options = optn, ds = 0.01, nev = 10, dsmax = 0.2, n_inversion = 8)
br = continuation(prob, PALC(), opts, bothside = false)
#######################################################
# hopf aBS
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, p_max = 10., p_min=-5., max_steps = 5, nev = 10, tol_stability = 1e-8, detect_bifurcation = 0)

for m in (3,4,5), Ntst in (99, 100)
    probpo = PeriodicOrbitOCollProblem(Ntst, m; N = 2, jacobian = BK.AutoDiffDense())
    br_pocoll = @time continuation(
            br, 1, opts_po_cont,
            probpo;
            δp = 0.001,
            normC = norminf,
            )
    # test anaytical jacobian
    ind_po = 5
    _po = br_pocoll.sol[ind_po].x
    _pars = BK.setparam(br,br_pocoll.sol[ind_po].p)
    _J = @time BK.jacobian(br_pocoll.prob, _po, _pars)
    @test ((_J-_J2)[1:end-1,1:end-1] |> norminf) ≈ 0 atol = 1e-14
    _J2 = @time DDEBifurcationKit.analytical_jacobian_dde_cst(br_pocoll.prob.prob, _po, _pars)
    _J2 = @time DDEBifurcationKit.analytical_jacobian_dde_cst_floquetcoll(br_pocoll.prob.prob, _po, _pars)
    @test (_J2.J0 + _J2.Jd -_J)[1:end-1,1:end-1] |> norminf ≈ 0 atol = 1e-14
    _J2 = @time DDEBifurcationKit.analytical_jacobian_dde_cst_floquetgev(br_pocoll.prob.prob, _po, _pars)
    @test (_J2.J0 + sum(_J2.Jd) -_J)[1:end-1,1:end-1] |> norminf ≈ 0 atol = 1e-14
end