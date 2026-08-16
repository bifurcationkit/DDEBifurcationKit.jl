# using Revise, Plots
using Test
using DDEBifurcationKit, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

function humpriesVF(x, xd, p)
    (;κ1,κ2,γ,a1,a2) = p
    [
        -γ * x[1] - κ1 * xd.u[1][1] - κ2 * xd.u[2][1]
    ]
end

function delaysF(x, par)
    [
        par.a1 + par.c1 * x[1],
        par.a2 + par.c2 * x[1],
    ]
end

let
pars = (κ1=0.,κ2=2.3,a1=1.3,a2=6,γ=4.75,c1=1.,c2=1.234)
prob = SDDDEBifProblem(humpriesVF, delaysF, zeros(1), pars, (@optic _.κ1))
optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 13., p_min = 0., newton_options = optn, ds = 0.01, nev = 20, )

alg = PALC()
br = continuation(prob, alg, opts)
@test br.specialpoint[1].type == :hopf
BK.get_normal_form(br, 1)
################################################################################
# computation of periodic orbit
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= 0.1, dsmin = 1e-4, p_max = 12., p_min=-5., max_steps = 10,
     nev = 10, tol_stability = 1e-8, detect_bifurcation = 0)
@reset opts_po_cont.newton_options.verbose = true

for ma in (false)
    probpo = Collocation(100, 4; N = 1, jacobian = BK.AutoDiffDense(), meshadapt = ma, K = 100)
    br_pocoll = continuation(
        br, 1, opts_po_cont,
        probpo;
        alg = PALC(tangent = Bordered()),
        # regular continuation options
        verbosity = 3,
        ampfactor = 1.0,
        δp = 0.2,
        normC = norminf,
        )

    ind_po = 9
    _po_saved = br_pocoll.sol[ind_po].x
    _pars = @set pars.κ1 = br_pocoll.sol[ind_po].p
    _wrap = BK.getprob(br_pocoll)
    BK.restore_problem!(_wrap, _po_saved, _pars)
    _po = BK.saved_solution(_po_saved)
   #  @test norminf(BK.residual(_wrap, _po, _pars)) < 1e-8

    # jacobian of the PO functional
    _J = BK.jacobian(_wrap, _po, _pars);
    # heatmap(iszero.(_J) , yflip=true, color = :viridis)

    _J2 = DDEBK.analytical_jacobian_dde_cst(BK.get_discretization(_wrap), _po, _pars)
    # @test norm(_J - _J2, Inf) < 1e-10

    # SD-DDE Floquet collocation jacobian: must match the ForwardDiff jacobian of the residual,
    # i.e. it correctly captures the state-dependent delays and the delay-derivative coupling.
    # (the Floquet jacobian has no phase-condition row, so we compare the collocation rows only)
    _J3 = DDEBK.analytical_jacobian_dde_floquet(BK.get_discretization(_wrap), _po, _pars)
    @test norm((_J3.J0 + _J3.Jd - _J)[1:end-1, 1:end-1], Inf) < 1e-9

    # DenseAnalytical uses the analytical jacobian (same as the Floquet one, full matrix)
    J_da = BK._jacobian_po(br_pocoll.prob, BK.DenseAnalytical(), _po, _pars)
    @test norm(J_da - _J, Inf) < 1e-8

    # Floquet exponents (extended monodromy)
    vals = DDEBK.__floquet_coll(BK.FloquetColl(), BK.get_discretization(_wrap), _po, _pars, 15)[1]
    # the phase mode (trivial Floquet exponent) must be ≈ 0: this only holds once the PO has
    # a real amplitude (amp ≈ 0.4 here); for a near-equilibrium orbit the monodromy is
    # ill-conditioned and the trivial exponent drifts away from 0.
    @test minimum(abs.(real.(vals))) < 1e-7
end # for
end