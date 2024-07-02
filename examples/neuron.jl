using Revise, DDEBifurcationKit, Parameters, Setfield, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

function neuronVF(x, xd, p)
   (;κ, β, a12, a21, τs, τ1, τ2) = p
   [
      -κ * x[1] + β * tanh(xd[3][1]) + a12 * tanh(xd[2][2]),
      -κ * x[2] + β * tanh(xd[3][2]) + a21 * tanh(xd[1][1])
   ]
end

delaysF(par) = [par.τ1, par.τ2, par.τs]

pars = (κ = 0.5, β = -1, a12 = 1, a21 = 0.5, τ1 = 0.2, τ2 = 0.2, τs = 1.5)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@lens _.τs); record_from_solution = (x,p)->(x1 = x[1], x2 = x[2]))

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 5., p_min = 0., newton_options = optn, ds = -0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 0, plot = true, bothside = true, normC = norminf)

plot(br)
################################################################################
prob2 = ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@lens _.a21); record_from_solution = prob.recordFromSolution)
br2 = BK.continuation(prob2, PALC(), ContinuationPar(opts, ds = 0.1, p_max = 4., n_inversion=8); verbosity = 0, plot = true, normC = norminf)

# @set! br2.contparams.newton_options.eigsolver.σ = 1e-5
BK.get_normal_form(br2, 1)
#Hopf l1 ≈ −0.0601.
BK.get_normal_form(br2, 2)

opts_pitch = ContinuationPar(br2.contparams, ds = 0.005, dsmax = 0.03, n_inversion = 4)
@set! opts_pitch.newton_options.eigsolver.γ = 0.001
br_pitchfork = continuation(br2, 2, opts_pitch)
plot(br2, br_pitchfork)
################################################################################
brhopf = continuation(br, 2, (@lens _.a21),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.04, max_steps = 230, p_max = 15., p_min = -1.,ds = -0.02);
         verbosity = 2, plot = true,
         detect_codim2_bifurcation = 2,
         # bothside = true,
         start_with_eigen = true)

plot(brhopf, vars = (:a21, :τs))
plot(brhopf, vars = (:τs, :ω))

brhopf2 = continuation(br, 2, (@lens _.a21),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.1, max_steps = 56, p_max = 1.5, p_min = -1.,ds = -0.01, n_inversion = 4);
         verbosity = 2, plot = true,
         detect_codim2_bifurcation = 2,
         start_with_eigen = true,
         bothside=true)

plot(brhopf, brhopf2, legend = :top)


################################################################################
# computation periodic orbit
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= 0.0001, dsmin = 1e-4, p_max = 10., p_min=-5., max_steps = 20,
    nev = 3, tol_stability = 1e-8, detect_bifurcation = 0, plot_every_step = 2, save_sol_every_step=1, detect_fold = true)
@set! opts_po_cont.newton_options.tol = 1e-8
@set! opts_po_cont.newton_options.verbose = true

# arguments for periodic orbits
args_po = (    record_from_solution = (x, p) -> begin
        _par = BK.getparams(p.prob)
        xtt = BK.get_periodic_orbit(p.prob, x, @set _par.a21=p.p)
        return (max = maximum(xtt[1,:]),
                min = minimum(xtt[1,:]),
                period = getperiod(p.prob, x, @set _par.a21=p.p))
    end,
    plot_solution = (x, p; k...) -> begin
        _par = BK.getparams(p.prob)
        xtt = BK.get_periodic_orbit(p.prob, x, @set _par.a21=p.p)
        plot!(xtt.t, xtt[1,:]; label = "V1", k...)
        plot!(xtt.t, xtt[2,:]; label = "V2", k...)
        plot!(br2; subplot = 1, putspecialptlegend = false)
        end,
    normC = norminf)

probpo = PeriodicOrbitOCollProblem(50, 3; N = 2, jacobian = BK.AutoDiffDense())
# probpo = PeriodicOrbitTrapProblem(M = 150, N = 2; jacobian = BK.AutoDiffDense())
br_pocoll = @time continuation(
    br2, 1, opts_po_cont,
    probpo;
    alg = PALC(tangent = Bordered()),
    verbosity = 2,    plot = true,
    args_po...,
    # eigsolver = BK.FloquetCollGEV(DefaultEig(), length(probpo), probpo.N),
    δp = 0.001,
    normC = norminf,
    callback_newton = (state; k...) -> begin
        xtt = BK.get_periodic_orbit(probpo,state.x,nothing)
        # plot(xtt.t, xtt[1,:], title = "it = $(state.it)") |> display
        printstyled(color=:red, "amp = ", BK.amplitude(xtt[:,:],1),"\n")
        # @show state.x[end]
        # @show state.f[end]
        state.step < 6
    end
    )

plot(br2, br_pocoll)
plot(br_pocoll, vars = (:param, :period))

plot(br2, br_pocoll, br_pitchfork);plot!(br_pocoll, vars = (:param,:min))

# plot the periodic orbit
plot(layout = 2)
    for ii = 1:10:110
        solpo = BK.get_periodic_orbit(br_pocoll.γ.prob.prob, br_pocoll.sol[ii].x, 1)
        plot!(solpo.t ./ solpo.t[end], solpo.u[1,:], label = "", subplot = 1)
    end
    xlabel!("t / period", subplot = 1)
    plot!(br_pocoll, vars = (:param, :period), subplot = 2, xlims=(2.2,2.4))

plot(br2, br_pocoll, br_pitchfork);plot!(br_pocoll, vars = (:param,:min))
############
################################################################################
using  DifferentialEquations

function neuron_DE(du,u,h,p,t)
    @unpack κ, β, a12, a21, τs, τ1, τ2 = p
    du[1] = -κ * u[1] + β * tanh(h(p, t-τs)[1]) + a12 * tanh(h(p, t-τ2)[2])
    du[2] = -κ * u[2] + β * tanh(h(p, t-τs)[2]) + a21 * tanh(h(p, t-τ1)[1])
end

h(p, t) = -0*ones(2) .+ 0.25sin(t/4)
    prob_de = DDEProblem(neuron_DE,h(pars, 0),h,(0.,20000.),ContinuationPar(pars, a21 = br2.specialpoint[1].param + 0.01); constant_lags=delaysF(pars), abstol = 1e-10, reltol = 1e-9)
    alg = MethodOfSteps(Tsit5())
sol = solve(prob_de,alg)
    plot(plot(sol, xlims = (sol.t[end]-100,sol.t[end]), vars=(:t,:1),ylims=(-0.1,0.1)), plot(sol.t,sol[1,:]),title = "a21=$(sol.prob.p.a21)")

