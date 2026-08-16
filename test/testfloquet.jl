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
delaysF(par) = [par.T]

function plot_solution_po(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, nothing)
    plot!(xtt.t, xtt[1,:]; label = "V1", k...)
end

# case delays < period
let
pars = (A = 0.5, ω = 1.0, r = -0.1, T = pi/2) 
x0 = [0.01,]

prob = ConstantDDEBifProblem(sinusvf, delaysF, x0, pars, (@optic _.A))

args_po = (
    plot_solution = plot_solution_po,
    normC = norminf)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.01, dsmin = 1e-4, p_max = 0.9, max_steps = 120, nev = 15, tol_stability = 1e-5, plot_every_step = 1, newton_options = NewtonPar(tol = 1e-12, verbose = false))

# build the po functional
probpo = Collocation(20, 5; N = 1, 
    jacobian = BK.AutoDiffDense(), 
    prob_vf = prob, xπ = zeros(1), ϕ = zeros(2))
@reset probpo.xπ = zeros(length(probpo))
@reset probpo.ϕ = zeros(length(probpo))
ci = BK.generate_solution(probpo, t->[(pars.A)*cos(t)], 2pi);
BK.updatesection!(probpo, ci, nothing)
BK.po_residual(probpo, ci, pars) |> BK.norminf

br_pocoll = @time continuation(
            probpo, ci, BK.PALC(), ContinuationPar(opts_po_cont; detect_bifurcation = 0);
            args_po...,
            normC = norminf,
            )

# variational equation
# ∂p = (a(t)- λ)⋅p(t) + exp(-λ⋅τ)⋅b(t)⋅p(t-τ)
# we find a = 0 and b = -ω

# the floquet exponents are analytical
# using LambertW#, SparseArrays, RecursiveArrayTools
# λs = [pars.ω/pars.T * lambertw(complex(-pars.T,0), k) for k in -0:7]
# μs = exp.(λs*2pi)

# log.(μs) # gives
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
_J = BK._jacobian_po(br_pocoll.prob, BK.AutoDiffDense(), _po, _pars);
# plotH(iszero.(_J))

_J2 = DDEBK.analytical_jacobian_dde_cst(BK.get_discretization(br_pocoll.prob), _po, _pars)
# @test norm(_J - _J2, Inf) < 1e-10


_J2 = DDEBK.analytical_jacobian_dde_cst_floquetgev(BK.get_discretization(br_pocoll.prob), _po, _pars)
# plotH(iszero.(_J2.J0))
# plotH(iszero.(_J2.Jd[1]))
# (_J2.J0 + _J2.Jd[1] -_J)[1:end-1,1:end-1] |> norminf

_J2 = DDEBK.analytical_jacobian_dde_floquet(BK.get_discretization(br_pocoll.prob), _po, _pars)
# plotH(iszero.(_J2.J0))
# plotH(iszero.(_J2.Jd))
# (_J2.J0 + _J2.Jd -_J)[1:end-1,1:end-1] |> norminf

# computation of Floquet exponents based in GEV: it works!
@time DDEBK.__floquet_coll_gev(BK.FloquetGEV(DDE_DefaultEig(maxit=300, tol = 1e-12, σ = 1e-2)), br_pocoll.prob, _po, _pars, 15)[1]

# computation of Floquet exponents based Verheyden, Lust 2005
# it should be close to log.(μs), the analytical exponents of ẏ = -y(t-π/2)
vals = @time DDEBK.__floquet_coll(BK.FloquetColl(), BK.get_discretization(br_pocoll.prob), _po, _pars, 15)[1]
Is = sortperm(real.(vals), rev = true)
# analytical exponents: roots of λ e^{λπ/2} = -1, principal log of e^{λ 2π}
@test vals[Is[3]] ≈ (-6.417163653792045 - 0.8271574313995605im) atol = 1e-4
@test vals[Is[5]] ≈ (-8.793370519927757 - 0.623834539656118im) atol = 1e-4
@test vals[Is[7]] ≈ (-10.266815687546782 - 0.503218647834753im) atol = 1e-3
@test vals[Is[9]] ≈ (-11.339508661193713 - 0.42473619559681874im) atol = 1e-2
@test vals[Is[11]] ≈ (-12.183953392697111 - 0.3693428703906798im) atol = 5e-2
@test vals[Is[13]] ≈ (-12.88060245239543 - 0.32795127295702187im) atol = 1e-1
end

# case delays > period
let
pars = (A = 0.5, ω = 1.0, r = -0.1, T = 5pi/2)
x0 = [0.01,]

prob = ConstantDDEBifProblem(sinusvf, delaysF, x0, pars, (@optic _.A))

args_po = (
    plot_solution = plot_solution_po,
    normC = norminf)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.01, dsmin = 1e-4, p_max = 0.9, max_steps = 120, nev = 15, tol_stability = 1e-5, plot_every_step = 1, newton_options = NewtonPar(tol = 1e-12, verbose = false))

# build the po functional
probpo = Collocation(20, 5; N = 1, 
    jacobian = BK.AutoDiffDense(), 
    prob_vf = prob, xπ = zeros(1), ϕ = zeros(2))
@reset probpo.xπ = zeros(length(probpo))
@reset probpo.ϕ = zeros(length(probpo))
ci = BK.generate_solution(probpo, t->[(pars.A)*cos(t)], 2pi);
BK.updatesection!(probpo, ci, nothing)
BK.po_residual(probpo, ci, pars) |> BK.norminf

br_pocoll = @time continuation(
            probpo, ci, BK.PALC(), ContinuationPar(opts_po_cont; detect_bifurcation = 0);
            args_po...,
            normC = norminf,
            )

# the floquet exponents are analytical
# using LambertW#, SparseArrays, RecursiveArrayTools
# λs = [pars.ω/pars.T * lambertw(complex(-pars.T,0), k) for k in -0:7]
# μs = exp.(λs*2pi)
# log.(μs) # gives
# 8-element Vector{ComplexF64}:
#     0.9483095575323278 + 1.6698278968351687im
#  2.999519565323715e-32 - 2.4492935982947064e-16im
#   -0.46856588619946243 - 1.289859855250175im
#    -0.7634499910877753 - 2.5507193661000267im
#    -0.9784867467771906 + 2.4765942352387134im
#    -1.1477693415341552 + 1.2218182055154072im
#    -1.2873872201583652 - 0.03279891887379571im
#     -1.406204295114303 - 1.2875174628123405im

ind_po = 3
br_pocoll.sol[ind_po].p
period = br_pocoll.sol[ind_po].x[end]
_pars = BK.setparam(br_pocoll,br_pocoll.sol[ind_po].p)
_po = br_pocoll.sol[ind_po].x
_sol = BK.get_periodic_orbit(br_pocoll, ind_po)

vals = @time DDEBK.__floquet_coll(BK.FloquetColl(), BK.get_discretization(br_pocoll.prob), _po, _pars, 15)[1]
Is = sortperm(real.(vals), rev = true)
close_z(v, z, tol) = abs(v - z) < tol || abs(v - conj(z)) < tol
# analytical exponents: roots of λ e^{5πλ/2} = -1 (Lambert W), principal log of e^{λ 2π}
@test close_z(vals[Is[1]],  0.9483095575323278 + 1.6698278968351687im, 1e-4)
@test abs(vals[Is[3]]) < 1e-4 && abs(vals[Is[4]]) < 1e-4   # trivial μ = 1 (phase + e^{±it})
@test close_z(vals[Is[5]],  -0.46856588619946243 - 1.289859855250175im, 1e-4)
@test close_z(vals[Is[7]],  -0.7634499910877753 - 2.5507193661000267im, 1e-3)
@test close_z(vals[Is[9]],  -0.9784867467771906 + 2.4765942352387134im, 1e-2)
@test close_z(vals[Is[11]], -1.1477693415341552 + 1.2218182055154072im, 5e-2)
@test close_z(vals[Is[13]], -1.2873872201583652 - 0.03279891887379571im, 1e-1)
@test close_z(vals[Is[15]], -1.406204295114303 - 1.2875174628123405im, 1e-1)

end