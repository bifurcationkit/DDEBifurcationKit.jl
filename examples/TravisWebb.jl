# Hutchinson with diffusion
cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."
using Revise, DDEBifurcationKit, LinearAlgebra, Plots, SparseArrays
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

function TW(u, ud, p)
   (;a,Δ) = p
   (Δ*u) .+ u .- a .* ud[1] .* (1 .+ u)
end

delaysF(par) = [1.]

# discretisation
Nx = 500; Lx = pi/2;
X = -Lx .+ 2Lx/Nx*(0:Nx-1) |> collect

# boundary condition
h = 2Lx/Nx
Δ = spdiagm(0 => -2ones(Nx), 1 => ones(Nx-1), -1 => ones(Nx-1) ) / h^2; Δ[1,1]=Δ

pars = (a = 0.5, τ = 1.0, Δ = Δ, N = Nx)
x0 = zeros(Nx)

prob = ConstantDDEBifProblem(TW, delaysF, x0, pars, (@optic _.a))

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 10., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, normC = norminf)

plot(br)
################################################################################
hopfpt = BK.get_normal_form(br, 1)
plot(hopfpt.ζ |> real)
################################################################################
# case where we specify the jacobian
function JacTW(u, p)
   (;a,Δ) = p
   # we compute the jacobian at the steady state
   J0 = Δ + I .- a .* Diagonal(u)
   J1 = -a .* Diagonal(1 .+ u)
   return J0, [J1]
end

prob2 = ConstantDDEBifProblem(TW, delaysF, x0, pars, (@optic _.a); J = JacTW)

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 2., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 6)
br = continuation(prob2, PALC(), opts; verbosity = 1, plot = true, normC = norminf)
################################################################################
using DifferentialEquations

function TW_DE(du,u,h,p,t)
	@unpack a,Δ = p
   du .= (Δ*u) .+ u .- a .* h(p,t-1) .* (1 .+ u)
end

h(p, t) = real(hopfpt.ζ) .*(1+ 0.01cos(t/4))


heatmap(sol.t,1:Nx,reduce(hcat,sol.u))
