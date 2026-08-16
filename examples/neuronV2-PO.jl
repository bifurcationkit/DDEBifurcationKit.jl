cd(@__DIR__)
using Pkg, LinearAlgebra, Test
pkg"activate ."
using Revise, BifurcationKit, DDEBifurcationKit, Plots
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit
# https://www.dropbox.com/scl/fi/7ndu8tideuskogzdj1l5z/DDE_seminar_tutorial_IV.pdf?rlkey=7h789sjohazen8kswaj9uhg8b&e=1&dl=0
g(z) = (tanh(z − 1) + tanh(1)) * cosh(1)^2

function neuron2VF(x, xd, p)
   (; a,b,c,d) = p
   [
      -x[1] - a * g(b * xd.u[1][1]) + c * g(d * xd.u[2][2]),
      -x[2] - a * g(b * xd.u[1][2]) + c * g(d * xd.u[2][1])
   ]
end

delaysF(par) = [par.τ1, par.τ2]

pars = (a = 0.069, b = 2., c = 0.6, d = 1.2, τ1 = 11.6, τ2 = 20.3)
prob = ConstantDDEBifProblem(neuron2VF, delaysF, zeros(2), pars, (@optic _.c))

optn = NewtonPar(eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 1.0, p_min = 0.25, newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 9, dsmax = 0.05, n_inversion = 8)
br = continuation(prob, PALC(), opts; bothside = true)
plot(br)

hpnf = BK.get_normal_form(br, 2)
################################################################################
# computation periodic orbit
@views function plot_solution(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, nothing)
    plot!(xtt.t, xtt[1,:]; label = "V1", marker = :d, markersize = 2, k...)
    plot!(xtt.t, xtt[2,:]; label = "V2", k...)
    plot!(br; subplot = 1, putspecialptlegend = false)
end

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.02, ds = -0.0001, p_max = 1.2, p_min = 0.100661, max_steps = 200, nev = 10, tol_stability = 1e-4, plot_every_step = 10)

# arguments for periodic orbits
args_po = (  ;
    plot_solution,
    normC = norminf)

probpo = Collocation(20, 4, 
            # jacobian = BK.AutoDiffDense(),
            # jacobian = BK.DenseAnalytical(),
            jacobian = BK.DenseAnalyticalInplace(),
            # meshadapt = true,
            )
br_pocoll = @time continuation(
            br, 2, ContinuationPar(opts_po_cont; max_steps = 200, detect_bifurcation = 3),
            probpo;
            # verbosity = 2,
            plot = true,
            args_po...,
            )

plot(br, br_pocoll)

# - #  1,       ns at c ≈ +0.70955033 ∈ (+0.70955033, +0.71008605), |δp|=5e-04, δ = ( 2,  2), step =  18
# - #  2,       pd at c ≈ +0.64986821 ∈ (+0.64986821, +0.65094888), |δp|=1e-03, δ = (-1, -1), step =  22
# - #  3,       bp at c ≈ +0.46191900 ∈ (+0.46191900, +0.46191900), |δp|=2e-09, δ = (-1,  0), step =  44
# - #  4,       pd at c ≈ +0.46534971 ∈ (+0.46533014, +0.46534971), |δp|=2e-05, δ = (-1, -1), step =  50
# - #  5,       pd at c ≈ +0.59697038 ∈ (+0.59640882, +0.59697038), |δp|=6e-04, δ = ( 1,  1), step =  92
# - #  6,       bp at c ≈ +0.61515927 ∈ (+0.61515927, +0.61515953), |δp|=3e-07, δ = ( 1,  0), step = 114
# - #  7,       pd at c ≈ +0.52169354 ∈ (+0.52169354, +0.52176700), |δp|=7e-05, δ = ( 1,  1), step = 161
# - #  8,       pd at c ≈ +0.52174066 ∈ (+0.52173603, +0.52174066), |δp|=5e-06, δ = (-1, -1), step = 165
# - #  9, endpoint at c ≈ +0.61178495, 