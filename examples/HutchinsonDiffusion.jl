# Hutchinson with diffusion
cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."
using Revise, DDEBifurcationKit, Parameters, Setfield, LinearAlgebra, Plots, SparseArrays
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

using DiffEqOperators

function Hutchinson(u, ud, p)
   @unpack a,d,Δ = p
   d .* (Δ*u) .- a .* ud[1] .* (1 .+ u)
end

delaysF(par) = [1.]

# discretisation
Nx = 200; Lx = pi/2;
X = -Lx .+ 2Lx/Nx*(0:Nx-1) |> collect

# boundary condition
Q = Neumann0BC(X[2]-X[1])
Dxx = sparse(CenteredDifference(2, 2, X[2]-X[1], Nx) * Q)[1]

pars = (a = 0.5, d = 1, τ = 1.0, Δ = Dxx, N = Nx)
x0 = zeros(Nx)

prob = ConstantDDEBifProblem(Hutchinson, delaysF, x0, pars, (@lens _.a))

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 10., pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 5, dsmax = 0.2, nInversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, normC = norminf)

plot(br)
hopfpt = BK.getNormalForm(br, 1)
################################################################################
# case where we specify the jacobian
function JacHutchinson(u, p)
   @unpack a,d,Δ = p
   # we compute the jacobian at the steady state
   J0 = d * Δ .- a .* Diagonal(u)
   J1 = -a .* Diagonal(1 .+ u)
   return J0, [J1]
end

prob2 = ConstantDDEBifProblem(Hutchinson, delaysF, x0, pars, (@lens _.a); J = JacHutchinson)

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 10., pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 5, dsmax = 0.2, nInversion = 4)
br = continuation(prob2, PALC(), opts; verbosity = 1, plot = true, normC = norminf)
