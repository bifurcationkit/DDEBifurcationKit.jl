BK.jacobian(hopfpb::BK.FoldMAProblem{Tprob, BK.AutoDiff, Tu0, Tp, Tl, Tplot, Trecord}, x, p) where {Tprob <: BK.FoldProblemMinimallyAugmented{ <: ConstantDDEBifProblem}, Tu0, Tp, Tl <: Union{BK.AllOpticTypes, Nothing}, Tplot, Trecord} = JacobianCodim2DDE(hopfpb, ForwardDiff.jacobian(z -> hopfpb.prob(z, p), x), x, p)


function (eig::BK.FoldEig)(Jdde::JacobianCodim2DDE, nev; kwargs...)
    xh = BK.getvec(Jdde.x)          # hopf point
    p1 = BK.getp(Jdde.x)            # first parameter
    newpar = set(Jdde.p, BK.getlens(Jdde.prob.prob), p1)
    J = BK.jacobian(Jdde.prob.prob.prob_vf, xh, newpar)
    eigenelts = eig.eigsolver(J, nev; kwargs...)
    return eigenelts
end

# Bogdanov-Takens / Cusp test function for the Fold functional
function BK.test_bt_cusp(iter::BK.ContIterable{BK.FoldCont, <: BK.FoldMAProblem{ <: BK.FoldProblemMinimallyAugmented{Tprob}} }, state) where {Tprob <: AbstractDDEBifurcationProblem}
    probma = BK.getprob(iter)
    lens1, lens2 = BK.get_lenses(probma)

    z = BK.getx(state)
    x = BK.getvec(z)    # fold point
    p1 = BK.getp(z)     # first parameter
    p2 = BK.getp(state) # second parameter
    par = BK.getparams(probma)
    newpar = set(par, lens1, p1)
    newpar = set(newpar, lens2, p2)

    𝐅 = probma.prob
    𝒯 = eltype(𝐅)

    # expression of the jacobian
    J_at_xp = BK.jacobian(𝐅.prob_vf, x, newpar)
    JAd_at_xp = BK.has_adjoint(𝐅) ? BK.jacobian_adjoint(𝐅, x, newpar) : transpose(J_at_xp)

    bd_vec = BK._compute_bordered_vectors(𝐅, J_at_xp, JAd_at_xp)

    # compute new b
    ζ = bd_vec.v
    BK.VI.scale!(ζ, 1 / 𝐅.norm(ζ))

    # compute new a
    ζstar = bd_vec.w
    BK.VI.scale!(ζstar, 1 / 𝐅.norm(ζstar))

    # test function for Bogdanov-Takens
    Δζ = Δ(Val(:der), J_at_xp, ζ, 0)
    𝐅.BT = BK.VI.inner(ζstar, Δζ)
    𝐅.CP = BK.getp(state.τ)

    return 𝐅.BT, 𝐅.CP
end
