cd(@__DIR__)
# cd("..")
using Pkg, LinearAlgebra, Test
pkg"activate ."

# https://ddebiftool.sourceforge.net/demos/neuron/html/demo1_stst.html
using Revise, DDEBifurcationKit, Plots
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

function humpriesVF(x, xd, p)
   (; κ1, κ2, γ, a1, a2, c) = p
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


pars = (κ1=0., κ2=2.3, a1=1.3, a2=6, γ=4.75, c=1.)
x0 = zeros(1)

prob = DDEBK.SDDDEBifProblem(humpriesVF, delaysF, x0, pars, (@optic _.κ1))

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 13., p_min = 0., newton_options = optn, ds = -0.01, detect_bifurcation = 3, nev = 3, )

br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = true)

plot(br)
################################################################################
brhopf = continuation(br, 2, (@optic _.κ2),
         ContinuationPar(br.contparams, detect_bifurcation = 2, dsmax = 0.04, max_steps = 230, p_max = 5., p_min = -1.,ds = -0.02);
         verbosity = 2, plot = true,
         detect_codim2_bifurcation = 0,
         bothside = true,
         start_with_eigen = true)

plot(brhopf, vars = (:κ1, :κ2))
################################################################################
# computation periodic orbit

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, p_max = 12., p_min=-5., max_steps = 3000,
nev = 3, tol_stability = 1e-8, detect_bifurcation = 0, plot_every_step = 20)
@reset opts_po_cont.newton_options.tol = 1e-9
@reset opts_po_cont.newton_options.verbose = true

# arguments for periodic orbits
args_po = (    record_from_solution = (x, p; k...) -> begin
        xtt = BK.get_periodic_orbit(p.prob, x, nothing)
        _max = maximum(xtt[1,:])
        _min = minimum(xtt[1,:])
        return (amp = _max - _min,
                max = _max,
                min = _min,
                period = getperiod(p.prob, x, nothing))
    end,
    plot_solution = (x, p; k...) -> begin
        xtt = BK.get_periodic_orbit(p.prob, x, nothing)
        plot!(xtt.t, xtt[1,:]; label = "x", k...)
        plot!(br; subplot = 1, putspecialptlegend = false)
        end,
    normC = norminf)

probpo = PeriodicOrbitOCollProblem(200, 2; N = 1, jacobian = BK.AutoDiffDense())
br_pocoll = @time continuation(
    br, 2, opts_po_cont,
    probpo;
    alg = PALC(tangent = Bordered()),
    # regular continuation options
    verbosity = 2, plot = true,
    args_po...,
    ampfactor = 0.2,
    use_normal_form = false,
    δp = 0.01,
    callback_newton = (state; k...) -> begin
        xtt = BK.get_periodic_orbit(probpo,state.x,nothing)
        m1,m2 = extrema(xtt[:,:])
        printstyled(color=:red, "amp = ", m2-m1,"\n")
        printstyled(color=:green, "T = ", (state.x[end]),"\n")
        @show state.x[end]
        state.step < 15 && BK.cbMaxNorm(10.0)(state; k...)
    end
    )

plot(br);plot!(br_pocoll, plotfold=false, ylabel = "amplitude")
################################################################################
using DifferentialEquations

function humpriesVF_DE2(x,h,p,t)
    (;κ1,κ2,γ,a1,a2,c) = p
   -γ * x - κ1 * h(p, t-(a1 + c * x)) - κ2 * h(p, t-(a2 + c * x))
end

function h0(p, t)
     t ≤ 0 || error("history function is only implemented for t ≤ 0")
     0 .+ 0.03sin(t)
 end
prob_de = DDEProblem(humpriesVF_DE2,h0,(0.,10200.),ContinuationPar(pars, κ1 = br.specialpoint[2].param + 0.01); dependent_lags=((x,par,t)->par.a1 + par.c * x, (x,par,t)->par.a2 + par.c * x))
alg = MethodOfSteps(Rosenbrock23())
sol = solve(prob_de,alg)
plot(plot(sol, xlims = (sol.t[end]-30,sol.t[end])), plot(sol))

