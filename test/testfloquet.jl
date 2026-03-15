using Test, DDEBifurcationKit
using LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

# using Plots
# plotH(x) = heatmap(x, yflip=true, color = :viridis)

function sinusvf(X, xd, p)
    (;A, ω, r) = p
    xτ = xd.u[1][1]
    x = X[1]
    [
        -xτ + (x^2 + xτ^2 - A^2)^2 * x
    ]
end
delaysF(par) = [pi/2]

pars = (A = 0.5, ω = 1.0, r = -0.1) 
x0 = [0.01,]

prob = ConstantDDEBifProblem(sinusvf, delaysF, x0, pars, (@optic _.A))
args_po = (
    plot_solution = (x, p; k...) -> begin
        xtt = BK.get_periodic_orbit(p.prob, x, nothing)
        plot!(xtt.t, xtt[1,:]; label = "V1", k...)
        end,
    normC = norminf)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.01, dsmin = 1e-4, p_max = 0.9, max_steps = 120, nev = 15, tol_stability = 1e-5, plot_every_step = 1, newton_options = NewtonPar(tol = 1e-12, verbose = false))

# build the po functional
probpo = PeriodicOrbitOCollProblem(20, 5; N = 1, 
    jacobian = BK.AutoDiffDense(),
    # jacobian = BK.DenseAnalytical(), 
    prob_vf = prob, xπ = zeros(1), ϕ = zeros(2))
@reset probpo.xπ = zeros(length(probpo))
@reset probpo.ϕ = zeros(length(probpo))
ci = BK.generate_solution(probpo, t->[(pars.A)*cos(t)], 2pi);
BK.updatesection!(probpo, ci, nothing)
BK.residual(probpo, ci, pars) |> BK.norminf

br_pocoll = @time continuation(
            probpo, ci, BK.PALC(), ContinuationPar(opts_po_cont; detect_bifurcation = 0);
            # verbosity = 2,
            # plot = true,
            args_po...,
            # eigsolver = BK.FloquetGEV(DDE_DefaultEig(maxit=200, tol = 1e-10, σ = 1e-3)),
            normC = norminf,
            )

# plot(br_pocoll)

# variational equation
# ∂p = (a(t)- λ)⋅p(t) + exp(-λ⋅τ)⋅b(t)⋅p(t-τ)
# we find a = 0 and b = -ω

# the floquet exponents are analytical
# using LambertW#, SparseArrays, RecursiveArrayTools
# λs = [2*pars.ω/pi * lambertw(complex(-pi/2,0), k) for k in -0:7]
# μs = exp.(λs*2pi)

# log.(μs) gives
# 8-element Vector{ComplexF64}:
#  2.999519565323715e-32 - 2.4492935982947064e-16im
#     -6.417163653792045 - 0.8271574313995605im
#     -8.793370519927757 - 0.623834539656118im
#    -10.266815687546782 - 0.503218647834753im
#    -11.339508661193713 - 0.42473619559681874im
#    -12.183953392697111 - 0.3693428703906798im
#     -12.88060245239543 - 0.32795127295702187im
#    -13.473627512100485 - 0.2957194071943831im

ind_po = 3
br_pocoll.sol[ind_po].p
period = br_pocoll.sol[ind_po].x[end]
_pars = BK.setparam(br_pocoll,br_pocoll.sol[ind_po].p)
_po = br_pocoll.sol[ind_po].x
_sol = BK.get_periodic_orbit(br_pocoll, ind_po)

# jacobian of the PO functional
_J = BK.jacobian(br_pocoll.prob, BK.AutoDiffDense(), _po, _pars);
# plotH(iszero.(_J))

_J2 = DDEBK.analytical_jacobian_dde_cst(br_pocoll.prob.prob, _po, _pars)
# @test norm(_J - _J2, Inf) < 1e-10


_J2 = DDEBK.analytical_jacobian_dde_cst_floquetgev(br_pocoll.prob.prob, _po, _pars)
# plotH(iszero.(_J2.J0))
# plotH(iszero.(_J2.Jd[1]))
# (_J2.J0 + _J2.Jd[1] -_J)[1:end-1,1:end-1] |> norminf

_J2 = DDEBK.analytical_jacobian_dde_cst_floquetcoll(br_pocoll.prob.prob, _po, _pars)
# plotH(iszero.(_J2.J0))
# plotH(iszero.(_J2.Jd))
# (_J2.J0 + _J2.Jd -_J)[1:end-1,1:end-1] |> norminf

# computation of Floquet exponents based in GEV: it works!
res = @time DDEBK.__floquet_coll_gev(BK.FloquetGEV(DDE_DefaultEig(maxit=300, tol = 1e-12, σ = 1e-2)), br_pocoll.prob, _po, _pars, 15)[1]

# computation of Floquet exponents based Verheyden, Lust 2005
# it does not work !!
vals = @time DDEBK.__floquet_coll(BK.FloquetColl(), br_pocoll.prob.prob, _po, _pars, 15)[1]
