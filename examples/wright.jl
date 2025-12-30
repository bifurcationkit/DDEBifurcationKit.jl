cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."

using Revise, DDEBifurcationKit
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

wrightVF(x, xd, p) = [-p.a * xd.u[1][1] * (1 + x[1])]
delaysF(par) = [1.0]
pars = (a=0.1, b=0.)
x0 = [0.]

prob = ConstantDDEBifProblem(wrightVF, delaysF, x0, pars, (@optic _.a), record_from_solution=(x,p;k...)-> (x=x[1], _x=1))

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 9., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 4, n_inversion = 12 )
br = continuation(prob, PALC(), opts; verbosity = 0, plot = false, bothside = false)
plot(br)

get_normal_form(br, 1)

function NFcoeff(N)
   # Faria 2006
   ωN = pi/2 + N*pi
   aN = ωN*(-1)^N
   K1 = aN / (1+aN^2)
   K2 = ωN / (5(1+aN^2))*((-1)^N-3ωN)
   return (a = K1, b = K2)
end

h1 = get_normal_form(br, 1)
@test isapprox(real(h1.nf.a), NFcoeff(0).a; rtol = 1e-5)
@test isapprox(real(h1.nf.b), NFcoeff(0).b; rtol = 1e-5)

h2 = get_normal_form(br, 2)
@test isapprox(real(h2.nf.a), NFcoeff(2).a; rtol = 1e-5)
@test isapprox(real(h2.nf.b), NFcoeff(2).b; rtol = 1e-5)

# arguments for periodic orbits
args_po = (    record_from_solution = (x, p; k...) -> begin
        xtt = get_periodic_orbit(p.prob, x, nothing)
        mi, ma = extrema(xtt[1,:])
        return (max = ma,
                min = mi,
                amp = ma-mi,
                period = getperiod(p.prob, x, nothing))
    end,
    plot_solution = (x, p; k...) -> begin
        xtt = get_periodic_orbit(p.prob, x, nothing)
        plot!(xtt.t, xtt[1,:]; label = "V1", marker = :d, k...)
        plot!(br; subplot = 1, putspecialptlegend = false)
        end,
    normC = norminf)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.01, dsmin = 1e-4, p_max = 10., p_min=-5., max_steps = 150, nev = 10, tol_stability = 1e-5, detect_bifurcation = 0, plot_every_step = 10, newton_options = NewtonPar(tol = 1e-10, verbose = true))

probpo = PeriodicOrbitOCollProblem(30, 5; N = 1, jacobian = DDEBK.BifurcationKit.AutoDiffDense(), meshadapt=false, K = 100)
br_pocoll = @time continuation(
            br, 1, ContinuationPar(opts_po_cont;detect_bifurcation = 0, nev = 10),
            probpo;
            verbosity = 2,
            plot = true,
            args_po...,
            normC = norminf,
            )

plot(br, br_pocoll)
plot(br_pocoll.param, br_pocoll.period)
plot(br_pocoll.param, delaysF(0) ./ br_pocoll.period)
