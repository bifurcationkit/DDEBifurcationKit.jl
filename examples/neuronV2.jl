cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."
using Revise, DDEBifurcationKit, LinearAlgebra, Plots, Accessors
using BifurcationKit
const BK = BifurcationKit

g(z) = (tanh(z − 1) + tanh(1))*cosh(1)^2
function neuron2VF(x, xd, p)
   (; a,b,c,d) = p
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

optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit=100))
opts = ContinuationPar(p_max = 0.4, p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 9, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 0, plot = true, bothside = false)

plot(br)

hpnf = BK.get_normal_form(br, 1)
################################################################################
brhopf = continuation(br, 1, (@optic _.c),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 1.1, p_min = -0.1,ds = 0.01, n_inversion = 2);
         verbosity = 0, plot = true,
         detect_codim2_bifurcation = 2,
         bothside = true,
         start_with_eigen = true)

brhopf2 = continuation(br, 2, (@optic _.c),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 1.1, p_min = -0.1,ds = -0.01);
         verbosity = 2, plot = true,
         detect_codim2_bifurcation = 2,
         bothside = true,
         start_with_eigen = true)

plot(brhopf, vars = (:a, :c), xlims = (0,0.7), ylims = (0,1))
plot!(brhopf2, vars = (:a, :c), xlims = (-0,0.7), ylims = (-0.1,1))

################################################################################
prob2 = ConstantDDEBifProblem(neuron2VF, delaysF, x0, (@set pars.a = 0.12), (@optic _.c))
br2 = continuation(prob2, PALC(), ContinuationPar(opts, p_max = 1.22); verbosity = 1, plot = true, bothside = false)

plot(br2)

# change tolerance for avoiding error computation of the EV
opts_fold = br.contparams
@reset opts_fold.newton_options.eigsolver.σ = 1e-7

brfold = continuation(br2, 3, (@optic _.a),
         ContinuationPar(opts_fold; detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 0.6, p_min = -0.5,ds = -0.01, n_inversion = 2, tol_stability = 1e-6);
         verbosity = 1, plot = true,
         detect_codim2_bifurcation = 2,
         bothside = false,
         start_with_eigen = true)

plot(brfold)

plot(brfold, vars = (:a, :c), branchlabel = "Fold")
plot!(brhopf, vars = (:a, :c), branchlabel = "Hopf")
plot!(brhopf2, vars = (:a, :c), branchlabel = "Hopf")
################################################################################
# computation periodic orbit

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, p_max = 10., p_min=-5., max_steps = 30, nev = 3, tol_stability = 1e-8, detect_bifurcation = 0, plot_every_step = 2, save_sol_every_step = 1)
@reset opts_po_cont.newton_options.tol = 1e-9
@reset opts_po_cont.newton_options.verbose = true
@reset opts_po_cont.newton_options.max_iterations = 8

# arguments for periodic orbits
args_po = (    record_from_solution = (x, p; k...) -> begin
        xtt = BK.get_periodic_orbit(p.prob, x, nothing)
        return (max = maximum(xtt[1,:]),
                min = minimum(xtt[1,:]),
                period = getperiod(p.prob, x, nothing))
    end,
    plot_solution = (x, p; k...) -> begin
        xtt = BK.get_periodic_orbit(p.prob, x, nothing)
        plot!(xtt.t, xtt[1,:]; label = "V1", k...)
        plot!(xtt.t, xtt[2,:]; label = "V2", k...)
        plot!(br; subplot = 1, putspecialptlegend = false)
        end,
    normC = norminf)

probpo = PeriodicOrbitOCollProblem(100, 3; N = 2, jacobian = BK.AutoDiffDense())
br_pocoll = @time continuation(
    br, 1, ContinuationPar(opts_po_cont; detect_bifurcation = 0, tol_stability = 1e-5),
    probpo;
    verbosity = 2,
    plot = true,
    args_po...,
    # ampfactor = 1/0.24391300209895822 * 0.1,
    ampfactor = 1.42,
    δp = 0.001,
    normC = norminf,
    # eigsolver = BK.FloquetCollGEV(DefaultEig(), 602, 2),
    callback_newton = (state; k...) -> begin
        xtt = BK.get_periodic_orbit(probpo,state.x,nothing)
        m1,m2 = extrema(xtt[:,:])
        # plot(xtt.t, xtt[1,:], title = "it = 0") |> display
        printstyled(color=:red, "amp = ", m2-m1,"\n")
        # @show state.x[end]
        # @show state.f[end]
        state.step < 16
    end
    )
################################################################################
using DifferentialEquations

function neuronV2_DE(du,x,h,p,t)
    (; a,b,c,d,τ1,τ2) = p
   du[1] = -x[1] - a * g(b*h(p, t-τ1)[1]) + c * g(d*h(p, t-τ2)[2])
   du[2] = -x[2] - a * g(b*h(p, t-τ1)[2]) + c * g(d*h(p, t-τ2)[1])
end

u0 = -2ones(2)
h0(p, t) = -0*ones(2) .+ 0.01cos(t/4)# h(p,t) = br_pocoll.orbit(t)prob_de = DDEProblem(neuronV2_DE,h0(pars,0),h0,(0.,54240.),(pars..., a = br.specialpoint[1].param + 0.001); constant_lags=delaysF(pars))alg = MethodOfSteps(Rosenbrock23())sol = solve(prob_de,alg)plot(plot(sol, xlims = (sol.t[end]-100,sol.t[end])), plot(sol))