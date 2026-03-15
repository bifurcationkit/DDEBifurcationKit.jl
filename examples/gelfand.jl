using Plots, SparseArrays
using Revise, Plots, SparseArrays
import LinearAlgebra: I, norm
using DDEBifurcationKit, BifurcationKit, LinearAlgebra, Plots
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

# define the sup norm and a L2 norm
normbratu(x) = norm(x .* w) / sqrt(length(x)) # the weight w is defined below

# some plotting functions to simplify our life
plotsol!(x, nx = Nx, ny = Ny; kwargs...) = heatmap!(reshape(x, nx, ny); color = :viridis, kwargs...)
plotsol(x, nx = Nx, ny = Ny; kwargs...) = (plot();plotsol!(x, nx, ny; kwargs...))

function Laplacian2D(Nx, Ny, lx, ly)
    hx = 2lx/Nx
    hy = 2ly/Ny
    D2x = spdiagm(0 => -2ones(Nx), 1 => ones(Nx-1), -1 => ones(Nx-1) ) / hx^2
    D2y = spdiagm(0 => -2ones(Ny), 1 => ones(Ny-1), -1 => ones(Ny-1) ) / hy^2

    D2x[1,1] = -1/hx^2
    D2x[end,end] = -1/hx^2

    D2y[1,1] = -1/hy^2
    D2y[end,end] = -1/hy^2

    D2xsp = sparse(D2x)
    D2ysp = sparse(D2y)
    A = kron(sparse(I, Ny, Ny), D2xsp) + kron(D2ysp, sparse(I, Nx, Nx))
    return A, D2x
end

ϕ(u, λ)  = -10(u-λ*exp(u))
dϕ(u, λ) = -10(1-λ*exp(u))

function Fmit(u::AbstractArray{T}, ud, p) where {T}
    f = similar(u, promote_type(T, typeof(p.λ)))
    mul!(f, p.Δ, u)
    f .= f .+ ϕ.(ud.u[1], p.λ)
    return f
end

function JFmit(x, p)
    J = p.Δ
    dg = dϕ.(x, p.λ)
    return J + spdiagm(0 => 0 .* dg), [spdiagm(0 => dg)]
end

Nx = 30; Ny = 30
lx = 0.5; ly = 0.5
# weight for the weighted norm
const w = (lx .+ LinRange(-lx,lx,Nx)) * (LinRange(-ly,ly,Ny))' |> vec

Δ, = Laplacian2D(Nx, Ny, lx, ly)
par_mit = (λ = .05, Δ = Δ, τ = 1.)

# initial guess f for newton
sol0 = zeros(Nx, Ny) |> vec

# Bifurcation Problem
delayF(par) = [par.τ]

prob = ConstantDDEBifProblem(Fmit, delayF, sol0, par_mit, (@optic _.λ),; J = JFmit,
  record_from_solution = (x, p; k...) -> (x = normbratu(x), n2 = norm(x), n∞ = norminf(x)),
  plot_solution = (x, p; k...) -> plotsol!(x ; k...))

# eigensolver
eigls = DDE_DefaultEig(σ = 1e-1, maxit=400, tol = 1e-8)

# options for Newton solver, we pass the eigen solver
opt_newton = BK.NewtonPar(tol = 1e-8, eigsolver = eigls, max_iterations = 20)

# options for continuation
opts_br = ContinuationPar(p_max = 3.5, p_min = 0.025,
    # for a good looking curve
    dsmin = 0.001, dsmax = 0.05, ds = 0.01,
    # number of eigenvalues to compute
    nev = 10,
    newton_options = (@set opt_newton.verbose = false),
    tol_stability = 1e-6,
    # detect codim 1 bifurcations
    detect_bifurcation = 3,
    plot_every_step = 3,
    # Optional: bisection options for locating bifurcations
    n_inversion = 2)

# optional arguments for continuation
kwargsC = (verbosity = 2, plot = true, normC = norminf)

br = continuation(prob, PALC(), opts_br; kwargsC...)

get_normal_form(br, 1)
