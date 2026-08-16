function BK.get_normal_form1d(prob::ConstantDDEBifProblem, 
                                br::BK.AbstractBranchResult,
                                ind_bif::Int,
                                Teigvec::Type{𝒯eigvec} = BK._getvectortype(br);
                                nev::Int = length(BK.eigenvalsfrombif(br, ind_bif)),
                                verbose::Bool = false,
                                lens = BK.getlens(br),
                                tol_fold = 1e-3,
                                scaleζ = LA.norm,

                                ζ = nothing,
                                ζ_ad = nothing,

                                autodiff::Bool = true,
                                detailed::Val{detailed_type} = Val(true),

                                bls = BK.MatrixBLS(),
                                ) where {𝒯eigvec, detailed_type}
    # parameters for normal form
    kwargs_nf = (;nev, verbose, lens, scaleζ)
    @warn "Computation of normal form based on a little hack ;)"
    Fode = (x,p) -> prob.VF.F(x, VectorOfArray([x for _ in eachindex(prob.delays0)]), p)
    prob_ode = BK.BifurcationProblem(Fode, prob.u0, prob.params, prob.lens; record_from_solution = prob.recordFromSolution)
    br_ode = @set br.contparams.newton_options.eigsolver = BK.DefaultEig()
    BK.get_normal_form1d(prob_ode, br_ode, ind_bif; kwargs_nf...)
end

function BK.__hopf_normal_form(prob::AbstractDDEBifurcationProblem, 
                            pt::BK.Hopf, 
                            ls::BK.AbstractLinearSolver; # for dispatch from BK 
                            verbose::Bool = false,
                            L = nothing)
    x0 = pt.x0
    p = pt.p
    lens = pt.lens
    parbif = set(pt.params, lens, p)
    ω = pt.ω
    λ0 = Complex(0, ω)
    δ = BK.getdelta(prob)

    # jacobian at the bifurcation point
    # do not recompute it if passed
    if isnothing(L)
        L = BK.jacobian(prob, x0, parbif)
    end
    Δ0  = Δ(L, 0λ0)
    Δ2ω = Δ(L, 2λ0)

    ζ = pt.ζ
    cζ = conj.(pt.ζ)
    ζ★ = copy(pt.ζ★)
    ζ★ ./= conj(LA.dot(ζ★, Δ(Val(:der), L, ζ, λ0)))
    # test the normalisation
    if ~isapprox(LA.dot(ζ★, Δ(Val(:der), L, ζ, λ0)), 1; rtol = 1e-3)
        @warn "We found instead $(LA.dot(ζ★, Δ(Val(:der), L, ζ, λ0)))"
    end

    x0c = VectorOfArray([copy(x0) for _ in 1:length(prob.delays0)+1])

    ζθ = expθ(L, ζ, λ0)
    ζθc = conj.(ζθ)

    # we use BilinearMap to be able to call on complex valued arrays
    R2 = BK.BilinearMap( (dx1, dx2)      -> BK.d2F(prob, x0c, parbif, dx1, dx2) ./2)
    R3 = BK.TrilinearMap((dx1, dx2, dx3) -> BK.d3F(prob, x0c, parbif, dx1, dx2, dx3) ./6 )

    # −LΨ001 = R01
    r01 = R01(prob, x0, parbif)
    Ψ001, cv, it = ls(Δ0, -r01)
    ~cv && @debug "[Hopf Ψ001] Linear solver for J did not converge. it = $it"
    Ψ001θ = Complex.(expθ(L, Ψ001, 0))

    # (2iω−L)Ψ200 = R20(ζ, ζ)
    r20 = R2(ζθ, ζθ)
    Ψ200, cv, it = ls(Δ2ω, r20)
    ~cv && @debug "[Hopf Ψ200] Linear solver for J did not converge. it = $it"
    Ψ200θ = expθ(L, Ψ200, 2λ0)
    # @assert Ψ200 ≈ (Complex(0, 2ω)*I - L) \ R20

    # −LΨ110 = 2R20(ζ,cζ)
    r20 = 2 .* R2(ζθ, ζθc)
    Ψ110, cv, it = ls(Δ0, r20)
    ~cv && @debug "[Hopf Ψ110] Linear solver for J did not converge. it = $it"
    Ψ110θ = Complex.(expθ(L, Ψ110, 0))

    # a = ⟨R11(ζ) + 2R20(ζ,Ψ001), ζ∗⟩
    _Jp = BK.jacobian(prob, x0, set(parbif, lens, p + δ))
    _Jm = BK.jacobian(prob, x0, set(parbif, lens, p - δ))
    r11 = (A(_Jp, ζ, λ0) .- A(_Jm, ζ, λ0)) ./ (2δ)
    av = r11 .+ 2 .* R2(ζθ, Ψ001θ)
    a = LA.dot(ζ★, av)

    # b = ⟨2R20(ζ,Ψ110) + 2R20(cζ,Ψ200) + 3R30(ζ,ζ,cζ), ζ∗⟩)
    bv = 2 .* R2(ζθ, Ψ110θ) .+ 2 .* R2(ζθc, Ψ200θ) .+ 3 .* R3(ζθ, ζθ, ζθc)
    b = LA.dot(ζ★, bv)

    # @error "info" b real(b)/ω/2 parbif δ Ψ110 Ψ200 2λ0

    verbose && println((a = a, b = b))

    # we set this type of normal form coefficients because the second order
    # hopf predictor does not work otherwise.
    @reset pt.nf = (;a, b, 
                    Ψ110_dde = Ψ110,
                    Ψ001_dde = Ψ001,
                    Ψ200_dde = Ψ200,
                    Ψ110 = zero(x0),
                    Ψ001 = zero(x0),
                    Ψ200 = zero(x0))
    if real(b) < 0
        pt.type = :SuperCritical
    elseif real(b) > 0
        pt.type = :SubCritical
    else
        pt.type = :Singular
    end
    verbose && printstyled(color = :red, "──▶ Hopf bifurcation point is: ", pt.type, "\n")
    return pt
end

# bilinear / trilinear contractions without conjugation, used for the delay-derivative
# corrections (the Taylor expansion of the delays is analytic, hence bilinear).
_lin(α, v) = sum(α[a] * v[a] for a in eachindex(α))
_bi(q, u, v) = sum(u[a] * q[a, b] * v[b] for a in axes(q, 1), b in axes(q, 2))

# Second derivative of the evaluation map E(φ) = (φ(0), φ(-τ₁(φ(0))), ..., φ(-τ_m(φ(0))))
# at the equilibrium: E''ⱼ(φ, ψ) = -(cⱼ·ψ(0)) φ'(-τⱼ*) - (cⱼ·φ(0)) ψ'(-τⱼ*).
# The perturbation functions are exponentials θ ↦ Ψ e^{λθ}, hence φ'(-τⱼ*) = λφ·φ.u[j+1].
function _sdde_E2(m, dτ, φ, λφ, ψ, λψ)
    out = VectorOfArray([zero(φ.u[1]) for _ in 1:m+1])
    for j in 1:m
        cj = dτ[j]
        out.u[j+1] .= -_lin(cj, ψ.u[1]) .* (λφ .* φ.u[j+1]) .-
                        _lin(cj, φ.u[1]) .* (λψ .* ψ.u[j+1])
    end
    return out
end

# delay-derivative correction to the second derivative G''(x*)·(φ, ψ):
# -Σⱼ (cⱼ·ψ(0)) Aⱼ φ'(-τⱼ*) - Σⱼ (cⱼ·φ(0)) Aⱼ ψ'(-τⱼ*)
function _sdde_corr2(Jd, dτ, φ, λφ, ψ, λψ)
    res = zero(φ.u[1])
    for j in eachindex(Jd)
        cj = dτ[j]
        res .-= _lin(cj, ψ.u[1]) .* (Jd[j] * (λφ .* φ.u[j+1]))
        res .-= _lin(cj, φ.u[1]) .* (Jd[j] * (λψ .* ψ.u[j+1]))
    end
    return res
end

# third derivative of the evaluation map at the equilibrium, slot j:
# E'''ⱼ(φ,ψ,χ) = φ''(cψ)(cχ) + ψ''(cχ)(cφ) + χ''(cφ)(cψ)
#               - φ'(ψqχ) - ψ'(φqχ) - χ'(φqψ)
function _sdde_E3j(dτ, dτ2, j, φ, λφ, ψ, λψ, χ, λχ)
    cj = dτ[j]; qj = dτ2[j]
    φ0 = φ.u[1]; ψ0 = ψ.u[1]; χ0 = χ.u[1]
    φj = φ.u[j+1]; ψj = ψ.u[j+1]; χj = χ.u[j+1]
    res = (λφ^2 .* φj) .* (_lin(cj, ψ0) * _lin(cj, χ0))
    res .+= (λψ^2 .* ψj) .* (_lin(cj, χ0) * _lin(cj, φ0))
    res .+= (λχ^2 .* χj) .* (_lin(cj, φ0) * _lin(cj, ψ0))
    res .-= (λφ .* φj) .* _bi(qj, ψ0, χ0)
    res .-= (λψ .* ψj) .* _bi(qj, φ0, χ0)
    res .-= (λχ .* χj) .* _bi(qj, φ0, ψ0)
    return res
end

# delay-derivative correction to the third derivative G'''(x*)·(φ, ψ, χ):
# F''(E''(φ,ψ), E'χ) + F''(E''(φ,χ), E'ψ) + F''(E''(ψ,χ), E'φ) + Σⱼ Aⱼ E'''ⱼ
function _sdde_corr3(m, R2_int, Jd, dτ, dτ2, φ, λφ, ψ, λψ, χ, λχ)
    res = zero(φ.u[1])
    res .+= 2 .* R2_int(_sdde_E2(m, dτ, φ, λφ, ψ, λψ), χ)
    res .+= 2 .* R2_int(_sdde_E2(m, dτ, φ, λφ, χ, λχ), ψ)
    res .+= 2 .* R2_int(_sdde_E2(m, dτ, ψ, λψ, χ, λχ), φ)
    for j in 1:m
        res .+= Jd[j] * _sdde_E3j(dτ, dτ2, j, φ, λφ, ψ, λψ, χ, λχ)
    end
    return res
end

function BK.__hopf_normal_form(prob::SDDDEBifProblem, 
                         pt::BK.Hopf, 
                        ls::BK.AbstractLinearSolver; # for dispatch from BK 
                        verbose::Bool = false,
                        L = nothing)
    x0 = pt.x0
    p = pt.p
    lens = pt.lens
    parbif = set(pt.params, lens, p)
    ω = pt.ω
    λ0 = Complex(0, ω)
    δ = BK.getdelta(prob)
    m = length(prob.delays0)

    # jacobian at the bifurcation point, do not recompute it if passed
    if isnothing(L)
        L = BK.jacobian(prob, x0, parbif)
    end
    Δ0  = Δ(L, 0λ0)
    Δ2ω = Δ(L, 2λ0)

    ζ = pt.ζ
    cζ = conj.(pt.ζ)
    ζ★ = copy(pt.ζ★)
    ζ★ ./= conj(LA.dot(ζ★, Δ(Val(:der), L, ζ, λ0)))
    # test the normalisation
    if ~isapprox(LA.dot(ζ★, Δ(Val(:der), L, ζ, λ0)), 1; rtol = 1e-3)
        @warn "We found instead $(LA.dot(ζ★, Δ(Val(:der), L, ζ, λ0)))"
    end

    # derivatives of the delays w.r.t. the state at the equilibrium
    dτ  = [ForwardDiff.gradient(x -> delays(prob, x, parbif)[j], x0) for j in eachindex(prob.delays0)]
    dτ2 = [ForwardDiff.hessian(x -> delays(prob, x, parbif)[j], x0) for j in eachindex(prob.delays0)]

    x0c = VectorOfArray([copy(x0) for _ in 1:length(prob.delays0)+1])

    ζθ = expθ(L, ζ, λ0)
    ζθc = conj.(ζθ)

    # intrinsic (delay-frozen) bilinear / trilinear forms
    R2_int = BK.BilinearMap( (dx1, dx2)      -> BK.d2F(prob, x0c, parbif, dx1, dx2) ./2)
    R3_int = BK.TrilinearMap((dx1, dx2, dx3) -> BK.d3F(prob, x0c, parbif, dx1, dx2, dx3) ./6 )

    # SD-DDE aware maps including the delay-derivative corrections
    _R2(φ, λφ, ψ, λψ) = R2_int(φ, ψ) .+ _sdde_corr2(L.Jd, dτ, φ, λφ, ψ, λψ) ./ 2
    _R3(φ, λφ, ψ, λψ, χ, λχ) = R3_int(φ, ψ, χ) .+ _sdde_corr3(m, R2_int, L.Jd, dτ, dτ2, φ, λφ, ψ, λψ, χ, λχ) ./ 6

    # −LΨ001 = R01
    r01 = R01(prob, x0, parbif)
    Ψ001, cv, it = ls(Δ0, -r01)
    ~cv && @debug "[Hopf Ψ001] Linear solver for J did not converge. it = $it"
    Ψ001θ = Complex.(expθ(L, Ψ001, 0))

    # a = ⟨R11(ζ) + 2R20(ζ, Ψ001), ζ∗⟩
    _Jp = BK.jacobian(prob, x0, set(parbif, lens, p + δ))
    _Jm = BK.jacobian(prob, x0, set(parbif, lens, p - δ))
    r11 = (A(_Jp, ζ, λ0) .- A(_Jm, ζ, λ0)) ./ (2δ)
    av = r11 .+ 2 .* _R2(ζθ, λ0, Ψ001θ, 0)
    a = LA.dot(ζ★, av)

    # (2iω−L)Ψ200 = R20(ζ, ζ)
    r20 = _R2(ζθ, λ0, ζθ, λ0)
    Ψ200, cv, it = ls(Δ2ω, r20)
    ~cv && @debug "[Hopf Ψ200] Linear solver for J did not converge. it = $it"
    Ψ200θ = expθ(L, Ψ200, 2λ0)

    # −LΨ110 = 2R20(ζ, cζ)
    r20 = 2 .* _R2(ζθ, λ0, ζθc, conj(λ0))
    Ψ110, cv, it = ls(Δ0, r20)
    ~cv && @debug "[Hopf Ψ110] Linear solver for J did not converge. it = $it"
    Ψ110θ = Complex.(expθ(L, Ψ110, 0))

    # b = ⟨2R20(ζ, Ψ110) + 2R20(cζ, Ψ200) + 3R30(ζ, ζ, cζ), ζ∗⟩)
    bv = 2 .* _R2(ζθ, λ0, Ψ110θ, 0) .+ 2 .* _R2(ζθc, conj(λ0), Ψ200θ, 2λ0) .+ 3 .* _R3(ζθ, λ0, ζθ, λ0, ζθc, conj(λ0))
    b = LA.dot(ζ★, bv)

    verbose && println((a = a, b = b))

    @reset pt.nf = (;a, b, 
                    Ψ110_dde = Ψ110,
                    Ψ001_dde = Ψ001,
                    Ψ200_dde = Ψ200,
                    Ψ110 = zero(x0),
                    Ψ001 = zero(x0),
                    Ψ200 = zero(x0))
    if real(b) < 0
        pt.type = :SuperCritical
    elseif real(b) > 0
        pt.type = :SubCritical
    else
        pt.type = :Singular
    end
    verbose && printstyled(color = :red, "──▶ Hopf bifurcation point is: ", pt.type, "\n")
    return pt
end

function BK.hopf_normal_form(prob::AbstractDDEBifurcationProblem,
                             br::BK.AbstractBranchResult, 
                             ind_hopf::Int,
                             Teigvec::Type{𝒯eigvec} = BK._getvectortype(br);
                             nev::Int = length(BK.eigenvalsfrombif(br, ind_hopf)),
                             verbose::Bool = false,
                             lens = BK.getlens(br),
                             detailed::Val{detailed_type} = Val(true),
                             start_with_eigen::Val{start_with_eigen_type} = Val(true),
                             scaleζ = LA.norm,
                             bls = BK.MatrixBLS(),
                             bls_adjoint = bls) where {detailed_type, 𝒯eigvec, start_with_eigen_type}
    # the kwargs detailed is only here to allow to extend BK.hopf_normal_form
    @assert br.specialpoint[ind_hopf].type == :hopf "The provided index does not refer to a Hopf Point"
    verbose && println("#"^53*"\n──▶ Hopf Normal form computation")

    options = br.contparams.newton_options

    # bifurcation point
    bifpt = br.specialpoint[ind_hopf]
    eigRes = br.eig

    # eigenvalue
    λ = eigRes[bifpt.idx].eigenvals[bifpt.ind_ev]
    ω = imag(λ)
    λ0 = Complex(0, ω)

    # parameter for vector field
    p = bifpt.param
    parbif = set(getparams(br), lens, p)
    L = BK.jacobian(prob, convert(Teigvec, bifpt.x), parbif)

    # right eigenvector
    if BK.haseigenvector(br) == false
        # we recompute the eigen-elements if there were not saved during the computation of the branch
        _λ, _ev, _ = options.eigsolver(L, bifpt.ind_ev + 2)
        @assert _λ[bifpt.ind_ev] ≈ λ "We did not find the correct eigenvalue $λ. We found $(_λ)"
        ζ = geteigenvector(options.eigsolver, _ev, bifpt.ind_ev)
    else
        ζ = copy(geteigenvector(options.eigsolver ,br.eig[bifpt.idx].eigenvecs, bifpt.ind_ev))
    end
    ζ ./= scaleζ(ζ)

    # left eigen-elements
    _Jt = BK.has_adjoint(prob) ? BK.jacobian_adjoint(prob, convert(Teigvec, bifpt.x), parbif) : adjoint(L)
    ζ★, λ★ = BK.get_adjoint_basis(_Jt, conj(λ), options.eigsolver; nev = nev, verbose = verbose)

    # check that λ★ ≈ conj(λ)
    abs(λ + λ★) > 1e-2 && @warn "We did not find the left eigenvalue for the Hopf point to be very close to the imaginary part:\nλ ≈ $λ,\nλ★ ≈ $λ★?\n You can perhaps increase the number of computed eigenvalues, the number is nev = $nev"

    # ζ, ζ★ = get_null_vectors(Δ(L, λ0))

    # normalise left eigenvector
    ζ★ ./= LA.dot(ζ, ζ★)
    @assert LA.dot(ζ, ζ★) ≈ 1

    hopfpt = BK.Hopf(bifpt.x, bifpt.τ, bifpt.param,
        ω,
        parbif, lens,
        ζ, ζ★,
        (a = missing, b = missing ),
        Symbol("?")
    )
    if ~detailed_type
        return hopfpt
    end
    return BK.__hopf_normal_form(prob, hopfpt, options.linsolver ; verbose)
end
