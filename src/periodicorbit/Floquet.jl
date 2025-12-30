import LinearAlgebra as LA
# A Newton-Picard Collocation Method for Periodic Solutions of Delay Differential Equations,
# author = Verheyden, Koen and Lust, Kurt,
	
# select the last K component of vector of length L
Itilde(K, L) = [zeros(K, L-K)  LA.I(K)]
isnonzero(x) = !iszero(x)

function BK.compute_eigenvalues(eig::FloquetColl, 
                                iter::BK.ContIterable{BK.PeriodicOrbitCont, <: BK.WrapPOColl{ <: BK.PeriodicOrbitOCollProblem{Tprob}}}, 
                                state, 
                                u0, 
                                par, 
                                nev = iter.contparams.nev; k...) where {Tkind <: BK.AbstractContinuationKind, Tprob <: AbstractDDEBifurcationProblem}
    wrapcoll = BK.get_wrap_po(iter)
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