abstract type AbstractDDEEigenSolver <: BifurcationKit.AbstractEigenSolver end

_build_dep(J0, Jd, delays) = NLE.DEP([J0, Jd...], [0, delays...])

"""
$(TYPEDEF)

Default eigen solver for DDEBifurcationKit based on the julia package `NonlinearEigenproblems.jl`. More precisely, we rely on `NonlinearEigenproblems.iar_chebyshev` for the computation of eigenvalues.

## Fields

$(TYPEDFIELDS)

## Constructors

- `DDE_DefaultEig(; kwargs...)` and `kwargs` are the fields above.

!!! tip
    If it fails, it is very likely that the jacobian of the delayed part is almost zero.
"""
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

function (eig::DDE_DefaultEig)(J::JacobianDDE, nev; kwargs...)
    return _eig_defaulteig(eig, J.J0, J.Jd, J.delays, nev; kwargs...)
end

function _eig_defaulteig(eig::DDE_DefaultEig, J0, Jd, delays, nev; kwargs...)
    # we have to be careful: Jd matrices must be non-zero
    dep = _build_dep(J0, Jd, delays)
    λ, V = NLE.iar_chebyshev(dep;
                    maxit = eig.maxit,
                    logger = eig.logger,
                    neigs = nev + 2,
                    tol = eig.tol,
                    v = isnothing(eig.v) ? rand(size(dep, 1), 1) : eig.v,
                    σ = eig.σ
                    )
    @assert length(λ) >= nev
    I = sortperm(λ, by = eig.which, rev = true)
    return λ[I], V[:, I], true, 1
end
