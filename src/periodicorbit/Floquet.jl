import LinearAlgebra as LA
import SparseArrays as SA

# dispatch on the following constructor to pass options specific to DDE
BK.FloquetGEV(eig::AbstractDDEEigenSolver; ntot = 0, n = 0, k...) = BK.FloquetGEV(eig, ntot, n; k...)
########################################################################################
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
                            nev = 3
                        )
    coll = wrapcoll.prob
    n, m, Ntst = size(coll)
    period = BK.getperiod(coll, u0, nothing)
    J = analytical_jacobian_dde_cst_floquetgev(coll, u0, par)
    @assert coll.prob_vf isa ConstantDDEBifProblem
    _delays = delays(coll.prob_vf, nothing, par)

    # λ⋅B * p + D * p - J0 * p - exp(-λ⋅τ) * Jd1 * p = 0
    # λ⋅B * p + J.J0 + exp(-λ⋅τ) * J.Jd[1] * p = 0

    # Icoll is the identity matrix for the collocation problem
    Icoll = analytical_jacobian_dde_cst(coll, u0, par; ρD = 0, ρF = 0, ρI = -1)[1:end-1, 1:end-1] #remove phase condition
    # remove periodic boundary condition
    for i = 1:n
        Icoll[end-n+i, end-n+i] = 0
        Icoll[end-n+i, i] = 0
    end
    B = SA.sparse(Icoll)

    USENEP = true

    if USENEP
        mats = [B, SA.sparse(J.J0[1:end-1, 1:end-1])]
    else
        mats = [SA.sparse(J.J0[1:end-1, 1:end-1])]
    end

    for i in eachindex(J.Jd)
        push!(mats, J.Jd[i][1:end-1, 1:end-1] |> SA.sparse)
    end

    if USENEP == false
        dep = NLE.DEP(mats, [0, _delays...]) # M(λ) = -λI + Σ_i A_i exp(-τ_i λ)
        pep = NLE.PEP([SA.spzeros(size(B)), LA.I + B])
        nep = NLE.SumNEP(pep, dep)
    else
        fs = Function[λ -> λ, λ -> one(λ)]
        for τ in _delays
            push!(fs, λ -> exp(-λ * τ))
        end
        nep = NLE.SPMF_NEP(mats, fs)
    end

    v0 = isnothing(eig.eigsolver.v) ? rand(size(nep, 1)) : eig.eigsolver.v
    v0[1:n] .= v0[end-n+1:end]
    
    args_nep = (maxit = eig.eigsolver.maxit,
                    neigs = nev + 2,
                    tol = eig.eigsolver.tol,
                    logger = eig.eigsolver.logger,
                    v = v0,)
    λ, V = NLE.iar_chebyshev(nep; args_nep...,
                    σ = eig.eigsolver.σ,
                    )
    λ2, V = NLE.iar_chebyshev(nep; args_nep...,
                    σ = eig.eigsolver.σ + pi*im/period,
                    )
    λ3, V = NLE.iar_chebyshev(nep; args_nep...,
                    σ = eig.eigsolver.σ - pi*im/period,
                    )
    append!(λ, λ2)
    append!(λ, λ3)
    

    # λ = @. log(complex(exp(λ * period)))
    I = sortperm(λ, by = real, rev = true)
    λ = λ[I] .* period
    # we filter the eigenvalues with large imaginary part
    # this must be done only if the translated version is in the spectrum...
    # indeed, if λ is in the spectrum, so is λ + 2πZ
    mytol = 1e-5
    λ = filter(x -> -pi + mytol < imag(x) < pi + mytol, λ)
    λ =  unique(round.(λ; digits = abs(Int(log10(mytol)))) .+ (0+0im) ) # this trick is for -0 ≈ 0

    return λ, nothing, true, 1
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
    return __floquet_coll(eig, BK.get_wrap_po(iter).prob, u0, par, nev)
end

function __floquet_coll(eig::FloquetColl,
                            coll,
                            u0::AbstractVector{𝒯},
                            par,
                            nev = 3
                        ) where {𝒯}
    n, m, Ntst = size(coll)
    period = BK.getperiod(coll, u0, par)
    J = analytical_jacobian_dde_cst_floquetcoll(coll, u0, par)

    # let's find the effective of Jd, ie the number of mesh points in [-tau_max, 0]
    # we then find the closest coarse mesh point
    Jdnz = vec(sum(isnonzero, J.Jd; dims = 1))
    a_left = findfirst(isnonzero, Jdnz)

    # remove the phase/periodicity condition and derivative wrt period
    B = @views J.J0[1:end-1-n, 1:end-1]
    A = @views J.Jd[1:end-1-n, a_left:end-1]

    # the part corresponding to t = 0 must be in A
    ii = n # ii = m * n
    A = hcat(A[:, 1:end-ii], B[:, 1:ii]) # replace the last columns of A with the first ones of B
    B = B[:, ii+1:end]

    @assert size(A, 1) == size(B, 1)

    # C   = [A B; I 0]
    # Cbc = [0 0; 0 I]
    if true
        # using ii does not change the following
        C = hcat(A, B)
        dn = size(C, 2) - size(C, 1)+1
        dn = size(A, 2)
        # matrix to ensure periodicity boundary condition
        Cbc = vcat(zero(C), zeros(𝒯, dn, size(C, 2)))
        C = vcat(C, zeros(𝒯, dn, size(C, 2)))
        s2 = size(B, 2)
        for i in 1:dn
            Cbc[s2+i, s2+i] = 1
            C[s2+i, i] = 1
        end
        # return Cbc, C
        vals  = LA.eigvals(Cbc, C) # (Cbc - λC)x = 0
        logvals = log.(complex.(vals))
        return sort(logvals, by = real, rev = true)[1:nev], nothing, true, 1
    end
    Mₜ = -B\A
    Nₜ, N = size(Mₜ)

    if N <= Nₜ
        It = Itilde(N, Nₜ)
        M = It * Mₜ # same as Mₜ[end-N+1:end, :]
    else
        It = Itilde(N - Nₜ, N)
        M = vcat(It, Mₜ)
    end

    vals = LA.eigvals(M)
    logvals = log.(complex.(vals))
    I = sortperm(logvals, by = real, rev = true)[1:min(nev, length(logvals))]
    # floquet exponents
    σ = logvals[I]
    # remove the trivial multiplier
    # deleteat!(σ, findmin(abs, σ)[2])
    return σ, nothing, true, 1
end