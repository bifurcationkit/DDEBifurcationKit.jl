abstract type AbstractDDEEigenSolver <: BifurcationKit.AbstractEigenSolver end

@with_kw mutable struct DDE_DefaultEig{T, Tw, Tv} <: AbstractDDEEigenSolver
    maxit::Int = 100
    which::Tw = real
    σ::T = 0.
    γ::T = 1.
    tol::T = 1e-10
    logger::Int = 0
    check_error_every::Int = 1
    v::Tv = nothing
end

function (eig::DDE_DefaultEig)(J::JacobianConstantDDE, nev; kwargs...)
    dep = NonlinearEigenproblems.DEP([J.J0, J.Jd...], [0, J.delays...])
    λ,V = NonlinearEigenproblems.iar_chebyshev(dep;
                    maxit = eig.maxit,
                    neigs = nev + 2,
                    tol = eig.tol,
                    v = isnothing(eig.v) ? rand(size(dep,1),1) : eig.v,
                    σ = eig.σ)
    @assert length(λ) >= nev
    I = sortperm(λ, by = eig.which, rev = true)
    return λ[I], V[:,I], true, 1
end
