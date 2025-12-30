import LinearAlgebra as LA
import SparseArrays as SA

# compute the Floquet exponents based on GEV. Seem online documentation.
function BK.compute_eigenvalues(eig::FloquetGEV{ <: AbstractDDEEigenSolver}, 
                                iter::BK.ContIterable{Tkind}, 
                                state, 
                                u0, 
                                par, 
                                nev = iter.contparams.nev; k...) where {Tkind <: BK.AbstractContinuationKind}
    return __floquet_coll_gev(eig, BK.get_wrap_po(iter), u0, par, nev)
end

function __floquet_coll_gev(eig::FloquetGEV{ <: AbstractDDEEigenSolver},
                            wrapcoll,
                            u0,
                            par,
                            nev = 2
                        )
    coll = wrapcoll.prob
    n, m, Ntst = size(coll)
    period = BK.getperiod(coll, u0, nothing)
    J = analytical_jacobian_dde_cst_floquetgev(wrapcoll, u0, par)
    @assert coll.prob_vf isa ConstantDDEBifProblem
    _delays = delays(coll.prob_vf, nothing, par)

    # λ⋅B * p + D * p - J0 * p - exp(-λ⋅τ) * Jd1 * p = 0
    # λ⋅B * p + J.J0 + exp(-λ⋅τ) * J.Jd[1] * p = 0
    fs = Function[λ -> λ, λ -> one(λ)]
    for τ in _delays
        push!(fs, λ -> exp(-λ*τ))
    end

    # B is the identity matrix for the collocation problem
    B = analytical_jacobian_dde_cst(wrapcoll, u0, par; ρD = 0, ρF = 0, ρI = -1)[1:end-1, 1:end-1]
    for i=1:n
        B[end-n+i,end-n+i] = 0
        B[end-n+i,i] = 0
    end
    B = SA.sparse(B)
        mats = [B, SA.sparse(J.J0[1:end-1, 1:end-1])]
    for i in eachindex(J.Jd)
        push!(mats, J.Jd[i][1:end-1, 1:end-1] |> SA.sparse)
    end

        nep = SPMF_NEP(mats, fs)
    v0 = isnothing(eig.eigsolver.v) ? rand(size(nep, 1)) : eig.eigsolver.v
    v0[1:n] .= v0[end-n+1:end]
    λ, V = NonlinearEigenproblems.iar_chebyshev(nep;
                    maxit = eig.eigsolver.maxit,
                    neigs = nev + 2,
                    tol = eig.eigsolver.tol,
                    v = v0,
                    σ = eig.eigsolver.σ
                    )
    λ = @. log(complex(exp(λ * period)))
    I = sortperm(λ, by = real, rev = true)
    return λ[I], nothing, true, 1
end
# A Newton-Picard Collocation Method for Periodic Solutions of Delay Differential Equations,
# author = Verheyden, Koen and Lust, Kurt,
	
# select the last K component of vector of length L
Itilde(K, L) = [zeros(K, L-K)  LA.I(K)]
isnonzero(x) = !iszero(x)

# compute the Floquet multipliers based on monodromy. See online documentation.
function BK.compute_eigenvalues(eig::FloquetColl, 
                                iter::BK.ContIterable{BK.PeriodicOrbitCont, <: BK.WrapPOColl{ <: BK.PeriodicOrbitOCollProblem{Tprob}}}, 
                                state, 
                                u0, 
                                par, 
                                nev = iter.contparams.nev; k...) where {Tkind <: BK.AbstractContinuationKind, Tprob <: AbstractDDEBifurcationProblem}
    wrapcoll = BK.get_wrap_po(iter)
    return __floquet_coll(eig, BK.get_wrap_po(iter), u0, par, nev)
end

function __floquet_coll(eig::FloquetColl,
                            wrapcoll,
                            u0,
                            par,
                            nev = 2
                        )
    coll = wrapcoll.prob
    n, = size(coll)
    # n = 0 # if we put this, we obtain the zero eigenvalue
    J = analytical_jacobian_dde_cst_floquetcoll(wrapcoll, u0, par)
    # let's find the effective of Jd, ie the number of mesh points in [-tau_max, 0]
    Jdnz = vec(sum(isnonzero, J.Jd; dims = 1))
    a_left = findfirst(isnonzero, vec(Jdnz))

    # remove the phase condition, periodicity condition and derivative wrt period
    B = @views J.J0[1:end-1-n, 1:end-1]
    A = @views J.Jd[1:end-1-n, a_left:end-1]
    @assert size(A, 1) == size(B, 1)

    Mₜ = -B\A
    Nₜ, N = size(Mₜ)

    if N <= Nₜ
        It = Itilde(N, Nₜ)
        M = It * Mₜ
    else
        It = Itilde(N - Nₜ, N)
        M = vcat(It, Mₜ)
    end

    vals = LA.eigvals(M)
    logvals = log.(complex.(vals))
    I = sortperm(logvals, by = real, rev = true)[1:nev]
    # floquet exponents
    σ = logvals[I]
    # remove the trivial multiplier
    # deleteat!(σ, findmin(abs, σ)[2])
    return σ, nothing, true, 1
end