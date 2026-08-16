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
# The standard MA formulation from BifurcationKit solves (J - iωI)v + aσ = 0.
# For DDEs, we build Δ(iω) = iωI - J₀ - Σ exp(-iωτⱼ)·Jⱼ and solve Δ(iω)v + aσ = 0.
function BK.hopf_ma_test(𝐇, J::JacobianDDE, a, b, J22, _zero, n, ω::𝒯) where 𝒯
    Δω = -Δ(J, Complex{𝒯}(0, ω))
    return 𝐇.linbdsolver(Δω, a, b, J22, _zero, n)
end

# Apply a JacobianDDE (the equilibrium Jacobian Jall = J₀ + Σ Jⱼ) to a vector.
# Needed for apply(J, x) calls in the MA framework.
function BK.apply(J::JacobianDDE, x)
    return J.Jall * x
end

###################################################################################################
# MA formulation call operator for DDEs (Δ(iω)-based)
# This dispatches (𝐇)(x, p, ω, params) when the prob_vf is a DDE problem.
# function (𝐇::BK.HopfMinimallyAugmentedFormulation{<:AbstractDDEBifurcationProblem})(x, p::𝒯, ω::𝒯, params) where 𝒯
#     a = 𝐇.a
#     b = 𝐇.b
#     par = BK.set(params, BK.getlens(𝐇), p)
#     J = BK.jacobian(𝐇.prob_vf, x, par)
#     _, σ1, cv = BK.hopf_ma_test(𝐇, J, a, b, zero(𝒯), 𝐇.zero, one(𝒯), ω)
#     ~cv && @debug "[Hopf DDE MA] Linear solver for Δ(iω) did not converge."
#     return BK.residual(𝐇.prob_vf, x, par), real(σ1), imag(σ1)
# end
###################################################################################################
# Compute bordered vectors using Δ(iω) for DDEs.
# Overrides _compute_bordered_vectors from BifurcationKit when J_at_xp is a JacobianDDE.
# The standard ODE version passes shift = ±iω to the linear solver, which doesn't work
# for JacobianDDE. Instead, we build Δ(iω) explicitly and solve without shift.
# This method handles both HopfMinimallyAugmentedFormulation (from BK) and
# HopfDDEFormulation (from DDEBifurcationKit) which share the same field structure.
# function BK.__compute_bordered_vectors(𝐇::BK.HopfMinimallyAugmentedFormulation, J_at_xp::JacobianDDE, JAd_at_xp, ω::𝒯) where 𝒯
function BK.__compute_bordered_vectors(linbdsolver, linbdsolver_adjoint, J_at_xp::JacobianDDE, JAd_at_xp, ω::𝒯, a, b, _zero) where {𝒯}
    Δω = Δ(J_at_xp, Complex{𝒯}(0, ω))
    Δω_adj = JAd_at_xp isa JacobianDDE ? Δ(JAd_at_xp, Complex{𝒯}(0, ω)) : JAd_at_xp
    v, _, cv1, itv = linbdsolver(Δω, a, b, zero(𝒯), _zero, one(𝒯))
    ~cv1 && @debug "Bordered linear solver for Δ(iω) did not converge."
    w, _, cv2, itw = linbdsolver_adjoint(Δω_adj, b, a, zero(𝒯), _zero, one(𝒯))
    ~cv2 && @debug "Bordered linear solver for Δ(iω)' did not converge."
    return (; v, w, itv, itw)
end

function BK._compute_bordered_vectors(𝐇::HopfDDEFormulation, J_at_xp::JacobianDDE, JAd_at_xp, ω::𝒯) where 𝒯
    @assert false
    Δω = Δ(J_at_xp, Complex{𝒯}(0, ω))
    Δω_adj = JAd_at_xp isa JacobianDDE ? Δ(JAd_at_xp, Complex{𝒯}(0, ω)) : JAd_at_xp
    v, _, cv1, itv = 𝐇.linbdsolver(Δω, 𝐇.a, 𝐇.b, zero(𝒯), 𝐇.zero, one(𝒯))
    ~cv1 && @debug "Bordered linear solver for Δ(iω) did not converge."
    w, _, cv2, itw = 𝐇.linbdsolverAdjoint(Δω_adj, 𝐇.b, 𝐇.a, zero(𝒯), 𝐇.zero, one(𝒯))
    ~cv2 && @debug "Bordered linear solver for Δ(iω)' did not converge."
    return (; v, w, itv, itw)
end
###################################################################################################
# Matrix-Based Jacobian for the DDE Hopf MA problem with MinAugMatrixBased.
# Builds the full Jacobian of the MA formulation:
#       [ ∂F/∂x    ∂F/∂p      0    ]
# J =   [ Re(σx)   Re(σp)  Re(σω)  ]
#       [ Im(σx)   Im(σp)  Im(σω)  ]
# where the state is X = (x; p; ω)
function BK.jacobian(pdpb::BK.HopfMAProblem{<:AbstractDDEBifurcationProblem, BK.MinAugMatrixBased},
                     X::AbstractVector{𝒯}, par) where 𝒯
    @assert false
    𝐇 = BK.get_formulation(pdpb)
    @assert false
    lens = BK.getlens(𝐇)
    n = length(X) - 2
    x = @view X[begin:n]
    p = X[n+1]
    ω = X[n+2]

    par0 = BK.set(par, lens, p)
    δ = BK.getdelta(𝐇.prob_vf)
    ϵ = 𝒯(δ)

    # Jacobians at the equilibrium (both may be JacobianDDE)
    J_at_xp  = BK.jacobian(𝐇.prob_vf, x, par0)
    JAd_at_xp = BK.has_adjoint(𝐇) ? BK.jacobian_adjoint(𝐇.prob_vf, x, par0) : LA.adjoint(J_at_xp)

    # full matrix form of the equilibrium Jacobian (for the top-left block of Jhopf)
    Jmat = J_at_xp isa JacobianDDE ? J_at_xp.Jall : J_at_xp

    # Bordered vectors: solve Δ(iω)·v + a·σ = 0, b'·v = 1
    Δω = Δ(J_at_xp, Complex{𝒯}(0, ω))
    v, σ1, cv1, itv = 𝐇.linbdsolver(Δω, 𝐇.a, 𝐇.b, zero(𝒯), 𝐇.zero, one(𝒯))
    ~cv1 && @debug "[DDE Hopf MA Jac] Bordered solver for Δ(iω) did not converge."

    # Adjoint: solve Δ(iω)'·w + b·σ̃ = 0, a'·w = 1
    Δω_adj = JAd_at_xp isa JacobianDDE ? Δ(JAd_at_xp, Complex{𝒯}(0, ω)) : JAd_at_xp
    w, σ2, cv2, itw = 𝐇.linbdsolverAdjoint(Δω_adj, 𝐇.b, 𝐇.a, zero(𝒯), 𝐇.zero, one(𝒯))
    ~cv2 && @debug "[DDE Hopf MA Jac] Bordered solver for Δ(iω)' did not converge."

    cw = LA.conj(w)
    vr = real(v); vi = imag(v)

    # ∂F/∂p via finite differences
    Fp_plus  = BK.residual(𝐇.prob_vf, x, BK.set(par, lens, p + ϵ))
    Fp_minus = BK.residual(𝐇.prob_vf, x, BK.set(par, lens, p - ϵ))
    dₚF = (Fp_plus - Fp_minus) / (2ϵ)

    # σₚ = -⟨w, ∂(Jv)/∂p⟩ via finite differences
    Jp_plus  = BK.jacobian(𝐇.prob_vf, x, BK.set(par, lens, p + ϵ))
    Jp_minus = BK.jacobian(𝐇.prob_vf, x, BK.set(par, lens, p - ϵ))
    dₚJv_num = (BK.apply(Jp_plus, v) - BK.apply(Jp_minus, v)) / (2ϵ)
    σₚ = -BK.VI.inner(w, dₚJv_num)

    # σω = i·⟨w, v⟩
    σω = Complex{𝒯}(0, 1) * BK.VI.inner(w, v)

    # σₓ via finite differences on the adjoint Jacobian
    u1r = BK.apply_jacobian(𝐇.prob_vf, x + ϵ * vr, par0, cw, true)
    u1i = BK.apply_jacobian(𝐇.prob_vf, x + ϵ * vi, par0, cw, true)
    u2r = BK.apply(JAd_at_xp, cw)
    σxv2r = @. -(u1r - u2r) / ϵ
    σxv2i = @. -(u1i - u2r) / ϵ
    σₓ = @. σxv2r + Complex{𝒯}(0, 1) * σxv2i

    # Assemble the full 3×3-block Jacobian
    z = BK.VI.zerovector(dₚF)
    Jhopf = hcat(Jmat, dₚF, z)
    Jhopf = vcat(Jhopf, vcat(real(σₓ), real(σₚ), real(σω))')
    Jhopf = vcat(Jhopf, vcat(imag(σₓ), imag(σₚ), imag(σω))')
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