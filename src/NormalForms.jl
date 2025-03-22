function BK.get_normal_form1d(prob::ConstantDDEBifProblem, 
                                br::ContResult, 
                                ind_bif::Int; 
                                kwargs_nf...)
    @warn "Computation of normal form based on a little hack ;)"
    Fode = (x,p) -> prob.VF.F(x, VectorOfArray([x for _ in eachindex(prob.delays0)]),p)
    prob_ode = BK.BifurcationProblem(Fode, prob.u0, prob.params, prob.lens; record_from_solution = prob.recordFromSolution)
    br_ode = @set br.contparams.newton_options.eigsolver = BK.DefaultEig()
    BK.get_normal_form1d(prob_ode, br_ode, ind_bif; kwargs_nf...)
end

function hopf_normal_form(prob::ConstantDDEBifProblem, 
                            pt::BK.Hopf, 
                            ls; 
                            autodiff = false,
                            verbose::Bool = false)
    x0 = pt.x0
    p = pt.p
    lens = pt.lens
    parbif = set(pt.params, lens, p)
    ω = pt.ω
    λ0 = Complex(0, ω)
    δ = BK.getdelta(prob)

    # jacobian at the bifurcation point
    L = BK.jacobian(prob, x0, parbif)
    Δ0 = Δ(L, 0λ0)
    Δ2ω = Δ(L, 2λ0)

    ζ = pt.ζ
    cζ = conj.(pt.ζ)
    ζstar = copy(pt.ζ★)
    ζstar ./= conj(dot(ζstar, Δ(Val(:der), L, ζ, λ0)))
    # test the normalisation
    if ~isapprox(dot(ζstar, Δ(Val(:der), L, ζ, λ0)), 1; rtol = 1e-2)
        @warn "We found instead $(dot(ζstar, Δ(Val(:der), L, ζ, λ0)))"
    end

    x0c = VectorOfArray([copy(x0) for _ in 1:length(prob.delays0)+1])

    ζθ = expθ(L, ζ, λ0)
    ζθc = conj.(ζθ)

    # we use BilinearMap to be able to call on complex valued arrays
    R2 = BK.BilinearMap( (dx1, dx2)      -> BK.d2F(prob, x0c, parbif, dx1, dx2) ./2)
    R3 = BK.TrilinearMap((dx1, dx2, dx3) -> BK.d3F(prob, x0c, parbif, dx1, dx2, dx3) ./6 )

    # −LΨ001 = R01
    R01 = (BK.residual(prob, x0, set(parbif, lens, p + δ)) .- 
           BK.residual(prob, x0, set(parbif, lens, p - δ))) ./ (2δ)
    Ψ001 = Complex.(expθ(L, ls(Δ0, -R01)[1], 0))

    # (2iω−L)Ψ200 = R20(ζ,ζ)
    R20 = R2(ζθ, ζθ)
    Ψ200 = expθ(L, ls(Δ2ω, R20)[1], 2λ0)
    # @assert Ψ200 ≈ (Complex(0, 2ω)*I - L) \ R20

    # −LΨ110 = 2R20(ζ,cζ).
    R20 = 2 .* R2(ζθ, ζθc)
    Ψ110 = Complex.(expθ(L, ls(Δ0, -R20)[1], 0))

    # a = ⟨R11(ζ) + 2R20(ζ,Ψ001), ζ∗⟩
    _Jp = BK.jacobian(prob, x0, set(parbif, lens, p + δ))
    _Jm = BK.jacobian(prob, x0, set(parbif, lens, p - δ))
    av = (A(_Jp, ζ, λ0) .- A(_Jm,  ζ, λ0)) ./ (2δ)
    av .+= 2 .* R2(ζθ, Ψ001)
    a = dot(ζstar, av)

    # b = ⟨2R20(ζ,Ψ110) + 2R20(cζ,Ψ200) + 3R30(ζ,ζ,cζ), ζ∗⟩)
    bv = 2 .* R2(ζθ, Ψ110) .+ 2 .* R2(ζθc, Ψ200) .+ 3 .* R3(ζθ, ζθ, ζθc)
    b = dot(ζstar, bv)

    # @info "info" b real(b)/ω/2 parbif δ Ψ110 Ψ200 2λ0

    # return coefficients of the normal form
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
    if real(a) * real(b) < 0
        pt.type = :SuperCritical
    elseif real(a) * real(b) > 0
        pt.type = :SubCritical
    else
        pt.type = :Singular
    end
    verbose && printstyled(color = :red, "--> Hopf bifurcation point is: ", pt.type, "\n")
    return pt
end

function hopf_normal_form(prob::SDDDEBifProblem, 
                        pt::BK.Hopf, 
                        ls; 
                        verbose::Bool = false)
    @error "Hopf normal form computation for SD-DDE is not implemented"
    a = Complex{eltype(pt.x0)}(1, 0)
    b = Complex{eltype(pt.x0)}(1, 0)
    x0 = pt.x0
    @reset pt.nf = (a = a, b = b,
                    Ψ110 = zero(x0),
                    Ψ001 = zero(x0),
                    Ψ200 = zero(x0))
    if real(a) * real(b) < 0
        pt.type = :SuperCritical
    elseif real(a) * real(b) > 0
        pt.type = :SubCritical
    else
        pt.type = :Singular
    end
    verbose && printstyled(color = :red,"--> Hopf bifurcation point is: ", pt.type, "\n")
    return pt
end

function BK.hopf_normal_form(prob::AbstractDDEBifurcationProblem,
                             br::BK.AbstractBranchResult, 
                             ind_hopf::Int;
                             nev = length(BK.eigenvalsfrombif(br, id_bif)),
                             verbose::Bool = false,
                             lens = BK.getlens(br),
                             Teigvec = BK._getvectortype(br),
                             scaleζ = norm,
                             autodiff = false,
                             detailed = false)
    # the kwargs detailed is only here to allow to extend BK.hopf_normal_form
    @assert br.specialpoint[ind_hopf].type == :hopf "The provided index does not refer to a Hopf Point"
    verbose && println("#"^53*"\n--> Hopf Normal form computation")

    options = br.contparams.newton_options

    # bifurcation point
    bifpt = br.specialpoint[ind_hopf]
    eigRes = br.eig

    # eigenvalue
    λ = eigRes[bifpt.idx].eigenvals[bifpt.ind_ev]
    ω = imag(λ)

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
    _Jt = BK.has_adjoint(prob) ? BK.jad(prob, convert(Teigvec, bifpt.x), parbif) : adjoint(L)
    ζstar, λstar = BK.get_adjoint_basis(_Jt, conj(λ), options.eigsolver; nev = nev, verbose = verbose)

    # check that λstar ≈ conj(λ)
    abs(λ + λstar) > 1e-2 && @warn "We did not find the left eigenvalue for the Hopf point to be very close to the imaginary part:\nλ ≈ $λ,\nλstar ≈ $λstar?\n You can perhaps increase the number of computed eigenvalues, the number is nev = $nev"

    # normalise left eigenvector
    ζstar ./= dot(ζ, ζstar)
    @assert dot(ζ, ζstar) ≈ 1

    hopfpt = BK.Hopf(bifpt.x, bifpt.τ, bifpt.param,
        ω,
        parbif, lens,
        ζ, ζstar,
        (a = zero(Complex{eltype(bifpt.x)}), b = zero(Complex{eltype(bifpt.x)}) ),
        :SuperCritical
    )
    return hopf_normal_form(prob, hopfpt, options.linsolver ; verbose, autodiff)
end
