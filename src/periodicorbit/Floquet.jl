import LinearAlgebra as LA
import SparseArrays as SA

# dispatch on the following constructor to pass options specific to DDE
BK.FloquetGEV(eig::AbstractDDEEigenSolver; ntot = 0, n = 0, k...) = BK.FloquetGEV(eig, ntot, n; k...)
########################################################################################
# Compute the Floquet exponents based on Generalized eigenvalue Problem. See online documentation.
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
    coll = BK.get_discretization(wrapcoll)
    n, _, _ = size(coll)
    period = BK.getperiod(coll, u0, nothing)
    J = analytical_jacobian_dde_cst_floquetgev(coll, u0, par)
    if !(coll.prob_vf isa ConstantDDEBifProblem)
        error("FloquetGEV (Generalized Eigenvalue / NEP) requires constant delays.\nUse `FloquetColl` (Verheyden-Lust monodromy) for state-dependent delays.")
    end
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
    λ,  = NLE.iar_chebyshev(nep; args_nep..., σ = eig.eigsolver.σ)
    λ2, = NLE.iar_chebyshev(nep; args_nep..., σ = eig.eigsolver.σ + pi * im / period)
    λ3, = NLE.iar_chebyshev(nep; args_nep..., σ = eig.eigsolver.σ - pi * im / period)
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
    λ =  unique(round.(λ; digits = abs(Int(log10(mytol)))) .+ (0 + 0im) ) # this trick is for -0 ≈ 0
    return λ, nothing, true, 1
end
########################################################################################
# Compute the Floquet exponents based on
# A Newton-Picard Collocation Method for Periodic Solutions of Delay Differential Equations,
# author = Verheyden, Koen and Lust, Kurt,
	
# compute the Floquet multipliers based on monodromy. See online documentation.
function BK.compute_eigenvalues(eig::FloquetColl, 
                                iter::BK.ContIterable{BK.PeriodicOrbitCont, <: BK.PeriodicOrbitFunctionalColl{ <: BK.Collocation{Tprob}}}, 
                                state, 
                                u0, 
                                par, 
                                nev = iter.contparams.nev; k...) where {Tprob <: AbstractDDEBifurcationProblem}
    wrapcoll = BK.get_wrap_po(iter)
    return __floquet_coll(eig, BK.get_discretization(wrapcoll), u0, par, nev)
end

function __floquet_coll(::FloquetColl,
                        coll,
                        u::AbstractVector{𝒯},
                        par,
                        nev = 3
                    ) where {𝒯}
    # Floquet multipliers via the collocation Jacobian of the DDE linearized around a T-periodic orbit.
    # The delayed contributions are placed at their true delayed times (which may span
    # several periods) by _dde_coll_jac_extended, so delays larger than the period are
    # supported; the generalized pencil (Cbc, C) then yields the multipliers μ.
    C, n_history, period = _dde_coll_jac_extended(coll, u, par)
    n, m, Ntst = size(coll)
    N = m * Ntst
    # segment length (history + u₀), in components: seg·n
    seg = n_history + 1
    dn = seg * n
    # current period length (u₁,...,u_N)
    s2 = N * n
    Cbc = vcat(zero(C), zeros(eltype(C), dn, size(C, 2)))
    C = vcat(C, zeros(eltype(C), dn, size(C, 2)))
    for i in 1:dn
        Cbc[s2+i, s2+i] = 1
        C[s2+i, i] = 1
    end
    # TODO: this is super-slow. Use EV instead of GEV
    vals = LA.eigvals(Matrix(Cbc), Matrix(C))

    # computation of eigenvalues
    logvals = log.(complex.(vals))
    I = sortperm(logvals, by = real, rev = true)
    return logvals[I[1:min(nev, length(I))]], nothing, true, 1
end