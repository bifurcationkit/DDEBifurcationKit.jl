function BK.continuation(br::BK.AbstractResult{Tkind, Tprob},
					ind_bif::Int,
					_contParams::ContinuationPar,
					probPO::BK.AbstractPeriodicOrbitProblem ;
					alg = br.alg,
					δp = nothing,
					ampfactor = 1,
					usedeflation = false,
					nev = length(BK.eigenvalsfrombif(br, ind_bif)),
					kwargs...) where {Tkind, Tprob <: AbstractDDEBifurcationProblem}
	# compute the normal form of the branch point
	verbose = get(kwargs, :verbosity, 0) > 1 ? true : false
	verbose && (println("--> Considering bifurcation point:"); BK._show(stdout, br.specialpoint[ind_bif], ind_bif))

	cb = get(kwargs, :callbackN, BK.cbDefault)

	hopfpt = BK.hopfNormalForm(br.prob, br, ind_bif; nev = nev, verbose = verbose)
	@error "Careful here"
	@set! hopfpt.nf.a = 1.

	# compute predictor for point on new branch
	ds = isnothing(δp) ? _contParams.ds : δp
	Ty = typeof(ds)
	pred = predictor(hopfpt, ds; verbose = verbose, ampfactor = Ty(ampfactor))
	@error "Careful here"
	@set! pred.p = br.specialpoint[ind_bif].param + δp

	verbose && printstyled(color = :green, "#"^61*
			"\n┌─ Start branching from Hopf bif. point to periodic orbits.",
			"\n├─ Bifurcation type = ", hopfpt.type,
			"\n├─── Hopf param  p0 = ", br.specialpoint[ind_bif].param,
			"\n├─── new param    p = ", pred.p, ", p - p0 = ", pred.p - br.specialpoint[ind_bif].param,
			"\n├─── amplitude p.o. = ", pred.amp,
			"\n├─── period       T = ", pred.period,
			"\n├─ Method = \n", probPO, "\n")

	# we compute a phase so that the constraint equation
	# < u(0) − u_hopf, ψ > is satisfied, i.e. equal to zero.
	ζr = real.(hopfpt.ζ)
	ζi = imag.(hopfpt.ζ)
	# this phase is for POTrap problem constraint to be satisfied
	ϕ = atan(dot(ζr, ζr), dot(ζi, ζr))
	verbose && printstyled(color = :green, "├─── phase ϕ        = ", ϕ / pi, "⋅π\n")

	M = BK.getMeshSize(probPO)
	orbitguess_a = [pred.orbit(t - ϕ) for t in LinRange(0, 2pi, M + 1)[1:M]]
	# TODO THIS HAS BEEN ADDED FOR BETTER INITIAL GUESS
	orbitguess_a[M] .= orbitguess_a[1]

	# extract the vector field and use it possibly to affect the PO functional
	prob_vf = BK.reMake(br.prob, params = BK.setParam(br, pred.p))

	# build the variable to hold the functional for computing PO based on finite differences
	probPO, orbitguess = reMake(probPO, prob_vf, hopfpt, ζr, orbitguess_a, abs(2pi/pred.ω); orbit = pred.orbit)
	# ADDED
	probPO.xπ .= orbitguess[1:end-1]

	if _contParams.newtonOptions.linsolver isa GMRESIterativeSolvers
		_contParams = @set _contParams.newtonOptions.linsolver.N = length(orbitguess)
	end

	# perform continuation
	branch = continuation(
		probPO, orbitguess, alg,
		_contParams;
		kwargs...
	)

	return BK.Branch(branch, hopfpt)
end
