using Test, DDEBifurcationKit
using LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

g(z) = (tanh(z − 1) + tanh(1))*cosh(1)^2
function neuron2VF(x, xd, p)
    (;a,b,c,d) = p
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
optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit = 100))
opts = ContinuationPar(p_max = 0.4, p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 9, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts, bothside = false)

#######################################################
# hopf aBS
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, p_max = 10., p_min=-5., max_steps = 5, nev = 3, tol_stability = 1e-8, detect_bifurcation = 0, plot_every_step = 2, save_sol_every_step = 1)
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
        br, 1, opts_po_cont,
        # PeriodicOrbitOCollProblem(100, 4);
        probpo;
        # verbosity = 2,    plot = true,
        args_po...,
        ampfactor = 1.42,
        δp = 0.001,
        normC = norminf,
        callback_newton = (state; k...) -> begin
            state.step < 16
        end
        )