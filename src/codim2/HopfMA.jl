import LinearAlgebra as LA

# The Hopf condition is Δ(iω, params)⋅q = 0
# where Δ(λ, params) ≡ λI − A0 - Aⱼ exp(−λτⱼ)

# Transpose / adjoint of JacobianDDE (needed by the MA formulation for JAd)
function LA.transpose(J::JacobianDDE)
    J0t = LA.transpose(J.J0)
    Jdt = [LA.transpose(A) for A in J.Jd]
    return JacobianDDE(J.prob, J0t, Jdt, J.delays)
end

function LA.adjoint(J::JacobianDDE)
    J0a = LA.adjoint(J.J0)
    Jda = [LA.adjoint(A) for A in J.Jd]
    return JacobianDDE(J.prob, J0a, Jda, J.delays)
end

# Dispatch hopf_ma_test for JacobianDDE.
# The standard MA formulation from BifurcationKit solves (J - iωI)v + aσ = 0
# For DDEs, we build Δ(iω) = iωI - J₀ - Σ exp(-iωτⱼ)·Jⱼ and solve Δ(iω)v + aσ = 0
function BK.hopf_ma_test(𝐇, J::JacobianDDE, a, b, J22, _zero, n, ω::𝒯) where 𝒯
    # -(iωI - J₀ - Σ exp(-iωτⱼ)Jⱼ) = J₀ + Σ exp(-iωτⱼ)Jⱼ - iωI = A(iω) - iωI,
    # i.e. the DDE analogue of the ODE operator (J - iωI).
    Δω = -Δ(J, Complex{𝒯}(0, ω))
    return 𝐇.linbdsolver(Δω, a, b, J22, _zero, n)
end

# Apply a JacobianDDE (the equilibrium Jacobian Jall = J₀ + Σ Jⱼ) to a vector.
# Needed for apply(J, x) calls in the MA framework.
function BK.apply(J::JacobianDDE, x)
    return J.Jall * x
end

###################################################################################################
# Compute bordered vectors using Δ(iω) for DDEs.
# Overrides _compute_bordered_vectors from BifurcationKit when J_at_xp is a JacobianDDE.
# The standard ODE version passes shift = ±iω to the linear solver, which doesn't work
# for JacobianDDE. Instead, we build Δ(iω) explicitly and solve without shift.
# This method handles both HopfMinimallyAugmentedFormulation (from BK) and
# HopfDDEFormulation (from DDEBifurcationKit) which share the same field structure.
# function BK.__compute_bordered_vectors(𝐇::BK.HopfMinimallyAugmentedFormulation, J_at_xp::JacobianDDE, JAd_at_xp, ω::𝒯) where 𝒯
function BK.__compute_bordered_vectors(linbdsolver, linbdsolver_adjoint, J_at_xp::JacobianDDE, JAd_at_xp, ω::𝒯, a, b, _zero) where {𝒯}
    # consistent with BK.hopf_ma_test: solve (A(iω) - iωI)v + aσ = 0
    Δω = -Δ(J_at_xp, Complex{𝒯}(0, ω))
    Δω_adj = JAd_at_xp isa JacobianDDE ? -Δ(JAd_at_xp, Complex{𝒯}(0, ω)) : JAd_at_xp
    v, _, cv1, itv = linbdsolver(Δω, a, b, zero(𝒯), _zero, one(𝒯))
    ~cv1 && @debug "Bordered linear solver for Δ(iω) did not converge."
    w, _, cv2, itw = linbdsolver_adjoint(Δω_adj, b, a, zero(𝒯), _zero, one(𝒯))
    ~cv2 && @debug "Bordered linear solver for Δ(iω)' did not converge."
    return (; v, w, itv, itw)
end
###################################################################################################
# Matrix-Based Jacobian for the DDE Hopf MA problem with MinAugMatrixBased.
# Builds the full Jacobian of the MA formulation:
#       [ ∂F/∂x    ∂F/∂p      0    ]
# J =   [ Re(σx)   Re(σp)  Re(σω)  ]
#       [ Im(σx)   Im(σp)  Im(σω)  ]
# where the state is X = (x; p; ω) and σ is the bordered variable of the functional.
# The σ rows are computed by finite differences of the functional σ (hopf_ma_test),
# which is guaranteed to be consistent with the functional used by the continuation.
function BK.jacobian(pdpb::BK.HopfMAProblem{<:BK.HopfMinimallyAugmentedFormulation{<:AbstractDDEBifurcationProblem}, BK.MinAugMatrixBased},
                     X::AbstractVector{𝒯}, par) where 𝒯
    𝐇 = BK.get_formulation(pdpb)
    lens = BK.getlens(𝐇)
    n = length(X) - 2
    x = @view X[begin:n]
    p = X[n+1]
    ω = X[n+2]

    δ = BK.getdelta(𝐇.prob_vf)
    ϵ = 𝒯(δ)

    # top-left block: ∂F/∂x = Jall, ∂F/∂p via finite differences
    J_at_xp = BK.jacobian(𝐇.prob_vf, x, BK.set(par, lens, p))
    Jmat = J_at_xp isa JacobianDDE ? J_at_xp.Jall : J_at_xp
    dₚF = (BK.residual(𝐇.prob_vf, x, BK.set(par, lens, p + ϵ)) -
           BK.residual(𝐇.prob_vf, x, BK.set(par, lens, p - ϵ))) / (2ϵ)

    # σ(x, p, ω): bordered variable of the functional, with the current vectors a, b
    σf(xl, pl, ωl) = begin
        Jl = BK.jacobian(𝐇.prob_vf, xl, BK.set(par, lens, pl))
        _, σl, _, _ = BK.hopf_ma_test(𝐇, Jl, 𝐇.a, 𝐇.b, zero(𝒯), 𝐇.zero, one(𝒯), ωl)
        σl
    end

    # ∂σ/∂x, ∂σ/∂p, ∂σ/∂ω by finite differences of the functional
    σx = zeros(Complex{𝒯}, n)
    for j in 1:n
        ej = zeros(𝒯, n); ej[j] = ϵ
        σx[j] = (σf(x .+ ej, p, ω) - σf(x .- ej, p, ω)) / (2ϵ)
    end
    σp = (σf(x, p + ϵ, ω) - σf(x, p - ϵ, ω)) / (2ϵ)
    σω = (σf(x, p, ω + ϵ) - σf(x, p, ω - ϵ)) / (2ϵ)

    # assemble the (n+2)×(n+2) Jacobian
    z = BK.VI.zerovector(dₚF)
    Jhopf = hcat(Jmat, dₚF, z)
    Jhopf = vcat(Jhopf, vcat(real(σx), real(σp), real(σω))')
    Jhopf = vcat(Jhopf, vcat(imag(σx), imag(σp), imag(σω))')
    return Jhopf
end

function BK.compute_eigenvalues(eig::BK.HopfEig{P, S}, 
                             iter::ContIterable,
                             state,
                             u0,
                             par,
                             nev = getcontparams(iter).nev;
                             kwargs...) where {P, S <: AbstractDDEEigenSolver}
    𝐏𝐛 = BK.getprob(iter)
    𝐇 = BK.get_formulation(𝐏𝐛)
    𝒯 = eltype(𝐇)
    zu = BK.getx(state)
    x = BK.getvec(zu, 𝐇) # hopf point
    newpar = BK.getparams(iter, state)
    J_at_xp = BK.jacobian(𝐇.prob_vf, x, newpar)
    return eig.eigsolver(J_at_xp, nev; iter, state, kwargs...)
end