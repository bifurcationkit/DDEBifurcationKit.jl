BK.jacobian(hopfpb::BK.FoldMAProblem{Tprob, BK.AutoDiff, Tu0, Tp, Tl, Tplot, Trecord}, x, p) where {Tprob <: BK.FoldProblemMinimallyAugmented{ <: ConstantDDEBifProblem}, Tu0, Tp, Tl <: Union{Lens, Nothing}, Tplot, Trecord} = JacobianCodim2DDE(hopfpb, ForwardDiff.jacobian(z -> hopfpb.prob(z, p), x), x, p)


function (eig::BK.FoldEigsolver)(Jdde::JacobianCodim2DDE, nev; kwargs...)
	xh = BK.getVec(Jdde.x)			# hopf point
	p1 = BK.getP(Jdde.x)			# first parameter
	newpar = set(Jdde.p, BK.getLens(Jdde.prob.prob), p1)

	J = BK.jacobian(Jdde.prob.prob.prob_vf, xh, newpar)
	eigenelts = eig.eigsolver(J, nev; kwargs...)

	return eigenelts
end
