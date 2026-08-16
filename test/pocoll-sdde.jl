# using Revise, Plots
using Test
using DDEBifurcationKit, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

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

alg = PALC()
br = continuation(prob, alg, opts; verbosity = 0, plot = false, bothside = true)
BK.get_normal_form(br, 2)

# plot(br)
################################################################################
# computation of periodic orbit
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, p_max = 12., p_min=-5., max_steps = 3,
    nev = 3, tol_stability = 1e-8, detect_bifurcation = 0, plot_every_step = 20,)
# @reset opts_po_cont.newton_options.tol = 1e-9
@reset opts_po_cont.newton_options.verbose = true

# arguments for periodic orbits
args_po = ( 
        plot_solution = (x, p; k...) -> begin
            xtt = BK.get_periodic_orbit(p.prob, x, nothing)
            plot!(xtt.t, xtt[1,:]; label = "x", k...)
            plot!(br; subplot = 1, putspecialptlegend = false)
            end,
        normC = norminf)

probpo = Collocation(50, 4; N = 1, jacobian = BK.AutoDiffDense())
br_pocoll = @time continuation(
    br, 2, opts_po_cont,
    probpo;
    alg = PALC(tangent = Bordered()),
    # regular continuation options
    # verbosity = 2,    plot = true,
    args_po...,
    ampfactor = 0.2,
    δp = 0.01,
    )

ind_po = 3
br_pocoll.sol[ind_po].p
period = br_pocoll.sol[ind_po].x[end]
_pars = BK.setparam(br_pocoll,br_pocoll.sol[ind_po].p)
_po = br_pocoll.sol[ind_po].x

# jacobian of the PO functional
_J = BK.jacobian(br_pocoll.prob, _po, _pars);
# heatmap(iszero.(_J) , yflip=true, color = :viridis)

_J2 = DDEBK.analytical_jacobian_dde_cst(BK.get_discretization(BK.getprob(br_pocoll)), _po, _pars)
# @test norm(_J - _J2, Inf) < 1e-10

# SD-DDE Floquet collocation jacobian: must match the ForwardDiff jacobian of the residual,
# i.e. it correctly captures the state-dependent delays and the delay-derivative coupling.
_J3 = DDEBK.analytical_jacobian_dde_floquet(BK.get_discretization(BK.getprob(br_pocoll)), _po, _pars)
@test norm(_J3.J0 + _J3.Jd - _J, Inf) < 1e-9

# DenseAnalytical uses the analytical jacobian (same as the Floquet one, full matrix)
J_da = BK._jacobian_po(br_pocoll.prob, BK.DenseAnalytical(), _po, _pars)
@test norm(J_da - _J, Inf) < 1e-8

# Floquet exponents (Verheyden-Lust monodromy)
vals = DDEBK.__floquet_coll(BK.FloquetColl(), BK.get_discretization(BK.getprob(br_pocoll)), _po, _pars, 6)[1]
@test length(vals) == 6
end