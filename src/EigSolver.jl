abstract type AbstractDDEEigenSolver <: BifurcationKit.AbstractEigenSolver end

@with_kw struct DDE_NLEVEigSolver{T} <: AbstractDDEEigenSolver
    maxit::Int = 100
    which::T = real
end

function (eig::DDE_NLEVEigSolver)(J::JacobianConstantDDE, nev; kwargs...)
    dep = NonlinearEigenproblems.DEP([J.J0, J.Jd...] , [0, J.delays...])
    λ,V = NonlinearEigenproblems.iar_chebyshev(dep,maxit=eig.maxit, neigs = nev+2)
    @assert length(λ) >= nev
    I = sortperm(λ, by = eig.which, rev = true)
    return λ[I],V[:,I],true,1
end
