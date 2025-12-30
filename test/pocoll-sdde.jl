# using Revise, Plots
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


pars = (κ1=0.,κ2=2.3,a1=1.3,a2=6,γ=4.75,c=1.)
x0 = zeros(1)

prob = SDDDEBifProblem(humpriesVF, delaysF, x0, pars, (@optic _.κ1))

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 13., p_min = 0., newton_options = optn, ds = -0.01, detect_bifurcation = 3, nev = 3, )

alg = PALC()
br = continuation(prob, alg, opts; verbosity = 0, plot = false, bothside = true)

# plot(br)
################################################################################
# computation of periodic orbit
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, p_max = 12., p_min=-5., max_steps = 3,
    nev = 3, tol_stability = 1e-8, detect_bifurcation = 0, plot_every_step = 20, save_sol_every_step=1)
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
    # verbosity = 2,    plot = true,
    args_po...,
    ampfactor = 1/0.467829783456199 * 0.1,
    δp = 0.01,
    callback_newton = (state; k...) -> begin
        xtt = BK.get_periodic_orbit(probpo,state.x,nothing)
        # plot(xtt.t, xtt[1,:], title = "it = $(state.it)") |> display
        # printstyled(color=:red, "amp = ", BK.amplitude(xtt[:,:],1),"\n")
        printstyled(color=:green, "T = ", (state.x[end]),"\n")
        @show state.x[end]
        state.step < 15 && BK.cbMaxNorm(10.0)(state; k...)
    end
    )
