using DDEBifurcationKit, LinearAlgebra, SparseArrays
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

function Hutchinson(u, ud, p)
   (;a,d,Δ) = p
   d .* (Δ*u) .- a .* ud.u[1] .* (1 .+ u)
end

delaysF(par) = [1.]

# discretisation
Nx = 100; Lx = pi/2;
X = -Lx .+ 2Lx/Nx*(0:Nx-1) |> collect
h = 2Lx/Nx
Δ = spdiagm(0 => -2ones(Nx), 1 => ones(Nx-1), -1 => ones(Nx-1) ) / h^2; Δ[1,1]=Δ[end,end]=-1/h^2

pars = (a = 0.5, d = 1, τ = 1.0, Δ = Δ, N = Nx)
x0 = zeros(Nx)

################################################################################
# case where we specify the jacobian
function JacHutchinson(u, p)
	(;a,d,Δ) = p
   # we compute the jacobian at the steady state
   J0 = d * Δ .- a .* Diagonal(u)
   J1 = -a .* Diagonal(1 .+ u)
   return J0, [J1]
end

prob = ConstantDDEBifProblem(Hutchinson, delaysF, x0, pars, (@optic _.a); J = JacHutchinson)
optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 10., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 1, plot = false, normC = norminf)
################################################################################
# periodic orbits
using DifferentiationInterface, SparseConnectivityTracer, SparseMatrixColorings

br_pocoll = @time continuation(
            br, 3, ContinuationPar(br.contparams; detect_bifurcation = 0, max_steps = 3, newton_options = NewtonPar(eigsolver = DDE_DefaultEig(), verbose = true), plot_every_step = 1),
            PeriodicOrbitOCollProblem(20, 4; jacobian = DDEBK.AutoSparseDI());
            verbosity = 2,
            # plot = true,
            normC = norminf,
            )
