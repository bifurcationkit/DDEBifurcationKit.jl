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

function TW(u, ud, p)
   @unpack a,Δ = p
   (Δ*u) .+ u .- a .* ud[1] .* (1 .+ u)
end

delaysF(par) = [1.]

# discretisation
Nx = 200; Lx = pi/2;
X = -Lx .+ 2Lx/Nx*(0:Nx-1) |> collect

# boundary condition
Q = Dirichlet0BC(X[2]-X[1]|>typeof)
Dxx = sparse(CenteredDifference(2, 2, X[2]-X[1], Nx) * Q)[1]

pars = (a = 0.5, τ = 1.0, Δ = Dxx, N = Nx)
x0 = zeros(Nx)

prob = ConstantDDEBifProblem(TW, delaysF, x0, pars, (@lens _.a))

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 10., pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 5, dsmax = 0.2, nInversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, normC = norminf)

plot(br)
################################################################################
hopfpt = BK.getNormalForm(br, 1)
plot(hopfpt.ζ |> real)
################################################################################
# case where we specify the jacobian
function JacTW(u, p)
   @unpack a,Δ = p
   # we compute the jacobian at the steady state
   J0 = Δ + I .- a .* Diagonal(u)
   J1 = -a .* Diagonal(1 .+ u)
   return J0, [J1]
end

prob2 = ConstantDDEBifProblem(TW, delaysF, x0, pars, (@lens _.a); J = JacTW)

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 10., pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 5, dsmax = 0.2, nInversion = 4)
br = continuation(prob2, PALC(), opts; verbosity = 1, plot = true, normC = norminf)
################################################################################
using  DifferentialEquations

function TW_DE(du,u,h,p,t)
	@unpack a,Δ = p
   du .= (Δ*u) .+ u .- a .* h(p,t-1) .* (1 .+ u)
end

h(p, t) = real(hopfpt.ζ) .*(1+ 0.01cos(t/4))
	prob_de = DDEProblem(TW_DE,h(pars,0),h,(0.,120.),setproperties(pars, a = br.specialpoint[1].param + 0.1); constant_lags=delaysF(pars))
	alg = MethodOfSteps(Rosenbrock23())
	sol = solve(prob_de,alg, progress=true)
	plot(sol.t, sol[1,:])

heatmap(sol.t,1:Nx,reduce(hcat,sol.u))
