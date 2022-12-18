getP(u, hopfPb::HopfDDEProblem) = u.x[end]
getVec(u, hopfPb::HopfDDEProblem) = u.x[1]

# this function encodes the functional
function (hp::HopfDDEProblem)(x, vR, vI, p::T, ω::T, params) where T
	# These are the equations of the Hopf bifurcation point
	# input:
	# - x guess for the point at which the jacobian has a purely imaginary eigenvalue
	# - p guess for the parameter for which the jacobian has a purely imaginary eigenvalue
	a = hp.a
	b = hp.b
	# update parameter
	par = set(params, BK.getLens(hp), p)
	tmp = vR .+ 1im .* vI
	w = Δ(hp.prob_vf, x, par, tmp, Complex(zero(T), ω))
	ps = dot(b, tmp)
	return BK.residual(hp.prob_vf, x, par), real(w), imag(w), [real(ps) - 1, imag(ps)]
end

function (hopfpb::HopfDDEProblem)(x::ArrayPartition, params)
	res = hopfpb(x.x[1], x.x[2], x.x[3], x.x[4][1], x.x[4][2], params)
	ArrayPartition(res...)
end
###################################################################################################
# define a problem <: AbstractBifurcationProblem
# @inline hasAdjoint(hopfpb::HopfMAProblem) = hasAdjoint(hopfpb.prob)
# @inline isSymmetric(hopfpb::HopfMAProblem) = isSymmetric(hopfpb.prob)
# residual(hopfpb::HopfMAProblem, x, p) = hopfpb.prob(x, p)
# # jacobian(hopfpb::HopfMAProblem, x, p) = hopfpb.jacobian(x, p)
# jacobian(hopfpb::HopfMAProblem{Tprob, Nothing, Tu0, Tp, Tl, Tplot, Trecord}, x, p) where {Tprob, Tu0, Tp, Tl <: Union{Lens, Nothing}, Tplot, Trecord} = (x = x, params = p, hopfpb = hopfpb.prob)
#
struct JacobianHopfDDE{T1,T2,T3,T4}
	prob::T1
	J::T2
	x::T3
	p::T4
end

(l::BK.DefaultLS)(J::JacobianHopfDDE, args...; kw...) = l(J.J, args...; kw...)


BK.jacobian(hopfpb::BK.HopfMAProblem{Tprob, BK.AutoDiff, Tu0, Tp, Tl, Tplot, Trecord}, x, p) where {Tprob <: HopfDDEProblem, Tu0, Tp, Tl <: Union{Lens, Nothing}, Tplot, Trecord} = JacobianHopfDDE(hopfpb, ForwardDiff.jacobian(z -> hopfpb.prob(z, p), x), x, p)
#
# jacobian(hopfpb::HopfMAProblem{Tprob, FiniteDifferences, Tu0, Tp, Tl, Tplot, Trecord}, x, p) where {Tprob, Tu0, Tp, Tl <: Union{Lens, Nothing}, Tplot, Trecord} = dx -> (hopfpb.prob(x .+ 1e-8 .* dx, p) .- hopfpb.prob(x .- 1e-8 .* dx, p)) / (2e-8)
################################################################################################### Newton / Continuation functions
function BK.newtonHopf(prob::ConstantDDEBifProblem,
			hopfpointguess::ArrayPartition,
			par,
			eigenvec, eigenvec_ad,
			options::NewtonPar;
			normN = norm,
			bdlinsolver::BK.AbstractBorderedLinearSolver = MatrixBLS(),
			usehessian = true,
			kwargs...)
	# we first need to update d2F and d3F for them to accept complex arguments

	hopfproblem = HopfDDEProblem(
		prob,
		BK._copy(eigenvec_ad),	# this is pb.a ≈ null space of (J - iω I)^*
		BK._copy(eigenvec), 	# this is pb.b ≈ null space of  J - iω I
		options.linsolver,
		# do not change linear solver if user provides it
		@set bdlinsolver.solver = (isnothing(bdlinsolver.solver) ? options.linsolver : bdlinsolver.solver);
		usehessian = usehessian)

	prob_h = BK.HopfMAProblem(hopfproblem, BK.AutoDiff(), hopfpointguess, par, nothing, prob.plotSolution, prob.recordFromSolution)

	# options for the Newton Solver
	opt_hopf = @set options.linsolver = BK.DefaultLS()

	# return prob_h

	# solve the hopf equations
	return BK.newton(prob_h, opt_hopf, normN = normN, kwargs...)
end

function BK.newtonHopf(br::BK.AbstractBranchResult, ind_hopf::Int;
			prob::ConstantDDEBifProblem = br.prob,
			normN = norm,
			options = br.contparams.newtonOptions,
			verbose = true,
			nev = br.contparams.nev,
			startWithEigen = false,
			kwargs...)
	hopfpointguess = HopfPoint(br, ind_hopf)
	ω = hopfpointguess.p[2]
	bifpt = br.specialpoint[ind_hopf]
	options.verbose && println("--> Newton Hopf, the eigenvalue considered here is ", br.eig[bifpt.idx].eigenvals[bifpt.ind_ev])
	@assert bifpt.idx == bifpt.step + 1 "Error, the bifurcation index does not refer to the correct step"
	ζ = geteigenvector(options.eigsolver, br.eig[bifpt.idx].eigenvecs, bifpt.ind_ev)
	ζ ./= normN(ζ)
	ζad = LinearAlgebra.conj.(ζ)

	if startWithEigen
		# computation of adjoint eigenvalue. Recall that b should be a null vector of J-iω
		λ = Complex(0, ω)
		p = bifpt.param
		parbif = BK.setParam(br, p)

		# jacobian at bifurcation point
		L = BK.jacobian(prob, bifpt.x, parbif)

		# computation of adjoint eigenvector
		_Jt = ~BK.hasAdjoint(prob) ? adjoint(L) : jad(prob, bifpt.x, parbif)

		ζstar, λstar = BK.getAdjointBasis(_Jt, conj(λ), options.eigsolver; nev = nev, verbose = false)
		ζad .= ζstar ./ dot(ζstar, ζ)
	end

	# solve the hopf equations
	hopfpointguess = ArrayPartition(hopfpointguess.u, real(ζ), imag(ζ), hopfpointguess.p)
	return newtonHopf(prob, hopfpointguess, getParams(br), ζ, ζad, options; normN = normN, kwargs...)
end


function BK.continuationHopf(prob_vf::ConstantDDEBifProblem, alg::BK.AbstractContinuationAlgorithm,
				hopfpointguess::ArrayPartition, par,
				lens1::Lens, lens2::Lens,
				eigenvec, eigenvec_ad,
				options_cont::ContinuationPar ;
				updateMinAugEveryStep = 0,
				normC = norm,
				bdlinsolver::BK.AbstractBorderedLinearSolver = MatrixBLS(),
				jacobian_ma::Symbol = :autodiff,
				computeEigenElements = false,
				usehessian = true,
				massmatrix = LinearAlgebra.I,
				kwargs...) where {Tb, vectype}
	@assert lens1 != lens2 "Please choose 2 different parameters. You only passed $lens1"
	@assert lens1 == BK.getLens(prob_vf)

	# options for the Newton Solver inheritated from the ones the user provided
	options_newton = options_cont.newtonOptions
	threshBT = 100options_newton.tol

	hopfPb = HopfDDEProblem(
		prob_vf,
		BK._copy(eigenvec_ad),	# this is a ≈ null space of (J - iω I)^*
		BK._copy(eigenvec), 	# this is b ≈ null space of  J - iω I
		options_newton.linsolver,
		# do not change linear solver if user provides it
		@set bdlinsolver.solver = (isnothing(bdlinsolver.solver) ? options_newton.linsolver : bdlinsolver.solver);
		usehessian = usehessian,
		massmatrix = massmatrix)

	# Jacobian for the Hopf problem
	if jacobian_ma == :autodiff
		# hopfpointguess = vcat(hopfpointguess.u, hopfpointguess.p)
		prob_h = BK.HopfMAProblem(hopfPb, BK.AutoDiff(), hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
		opt_hopf_cont = @set options_cont.newtonOptions.linsolver = DefaultLS()
	elseif jacobian_ma == :finiteDifferencesMF
		hopfpointguess = vcat(hopfpointguess.u, hopfpointguess.p)
		prob_h = FoldMAProblem(hopfPb, FiniteDifferences(), hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
		opt_hopf_cont = @set options_cont.newtonOptions.linsolver = options_cont.newtonOptions.linsolver
	else
		prob_h = HopfMAProblem(hopfPb, nothing, hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
		opt_hopf_cont = @set options_cont.newtonOptions.linsolver = HopfLinearSolverMinAug()
	end

	# this functions allows to tackle the case where the two parameters have the same name
	lenses = BK.getLensSymbol(lens1, lens2)

	# current lyapunov coefficient
	eTb = eltype(hopfpointguess)
	hopfPb.l1 = Complex{eTb}(0, 0)
	hopfPb.BT = one(eTb)
	hopfPb.GH = one(eTb)

	# this function is used as a Finalizer
	# it is called to update the Minimally Augmented problem
	# by updating the vectors a, b
	function updateMinAugHopf(z, tau, step, contResult; kUP...)
		# we first check that the continuation step was successful
		# if not, we do not update the problem with bad information!
		success = get(kUP, :state, nothing).converged
		(~modCounter(step, updateMinAugEveryStep) || success == false) && return true
		x = getVec(z.u, hopfPb)	# hopf point
		p1, ω = getP(z.u, hopfPb)
		p2 = z.p		# second parameter
		newpar = set(par, lens1, p1)
		newpar = set(newpar, lens2, p2)

		# we stop continuation at Bogdanov-Takens points

		# CA NE DEVRAIT PAS ETRE ISSNOT?
		isbt = isnothing(contResult) ? true : isnothing(findfirst(x -> x.type in (:bt, :ghbt, :btgh), contResult.specialpoint))

		# if the frequency is null, this is not a Hopf point, we halt the process
		if abs(ω) < threshBT
			@warn "[Codim 2 Hopf - Finalizer] The Hopf curve seems to be close to a BT point: ω ≈ $ω. Stopping computations at ($p1, $p2). If the BT point is not detected, try lowering Newton tolerance or dsmax."
		end

		# call the user-passed finalizer
		finaliseUser = get(kwargs, :finaliseSolution, nothing)
		resFinal = isnothing(finaliseUser) ? true : finaliseUser(z, tau, step, contResult; prob = hopfPb, kUP...)

		return abs(ω) >= threshBT && isbt && resFinal
	end

	function testBT_GH(iter, state)
		z = getx(state)
		x = getVec(z, hopfPb)		# hopf point
		p1, ω = getP(z, hopfPb)		# first parameter
		p2 = getp(state)			# second parameter
		newpar = set(par, lens1, p1)
		newpar = set(newpar, lens2, p2)

		probhopf = iter.prob.prob

		# expression of the jacobian
		J_at_xp = BK.jacobian(probhopf.prob_vf, x, newpar)

		# compute new b
		T = typeof(p1)
		# ζ = probhopf.linbdsolver(J_at_xp, a, b, T(0), probhopf.zero, T(1); shift = Complex(0, -ω), Mass = hopfPb.massmatrix)[1]
		λ = Complex(0, ω)
		ζ = @. z.x[2] + im * z.x[3]
		ζ ./= normC(ζ)

		# compute new ζstar
		# JAd_at_xp = BK.hasAdjoint(probhopf.prob_vf) ? jad(probhopf.prob_vf, x, newpar) : transpose(J_at_xp)

		JAd_at_xp = BK.hasAdjoint(probhopf.prob_vf) ? jad(probhopf.prob_vf, x, newpar) : adjoint(J_at_xp)
		ζstar, λstar = BK.getAdjointBasis(JAd_at_xp, conj(λ), BK.getContParams(iter).newtonOptions.eigsolver.eigsolver)
		# ζstar = probhopf.linbdsolver(JAd_at_xp, b, a, T(0), hopfPb.zero, T(1); shift = Complex(0, ω), Mass = transpose(hopfPb.massmatrix))[1]
		# test function for Bogdanov-Takens
		probhopf.BT = ω
		# BT2 = real( dot(ζstar ./ normC(ζstar), ζ) )
		# ζstar ./= dot(ζ, ζstar)

		hp = BK.Hopf(x, p1, ω, newpar, lens1, ζ, ζstar, (a = Complex{T}(0,0), b = Complex{T}(0,0)), :hopf)
		BK.hopfNormalForm(probhopf.prob_vf, hp, options_newton.linsolver, verbose = false)

		# lyapunov coefficient
		probhopf.l1 = hp.nf.b
		# test for Bautin bifurcation.
		# If GH is too large, we take the previous value to avoid spurious detection
		# GH will be large close to BR points
		probhopf.GH = abs(real(hp.nf.b)) < 1e5 ? real(hp.nf.b) : state.eventValue[2][2]
		return probhopf.BT, probhopf.GH
	end

	# the following allows to append information specific to the codim 2 continuation to the user data
	_printsol = get(kwargs, :recordFromSolution, nothing)
	_printsol2 = isnothing(_printsol) ?
		(u, p; kw...) -> (; zip(lenses, (getP(u, hopfPb)[1], p))..., ω = getP(u, hopfPb)[2], l1 = hopfPb.l1, BT = hopfPb.BT, GH = hopfPb.GH, BK.namedprintsol(BK.recordFromSolution(prob_vf)(getVec(u, hopfPb), p; kw...))...) :
		(u, p; kw...) -> (; BK.namedprintsol(_printsol(getVec(u, hopfPb), p; kw...))..., zip(lenses, (getP(u, hopfPb)[1], p))..., ω = getP(u, hopfPb)[2], l1 = hopfPb.l1, BT = hopfPb.BT, GH = hopfPb.GH)

	prob_h = reMake(prob_h, recordFromSolution = _printsol2)

	# eigen solver
	eigsolver = HopfDDEEig(BK.getsolver(opt_hopf_cont.newtonOptions.eigsolver))

	# event for detecting codim 2 points
	if computeEigenElements
		event = BK.PairOfEvents(BK.ContinuousEvent(2, testBT_GH, computeEigenElements, ("bt", "gh"), threshBT), BifDetectEvent)
		# careful here, we need to adjust the tolerance for stability to avoid
		# spurious ZH or HH bifurcations
		@set! opt_hopf_cont.tolStability = 10opt_hopf_cont.newtonOptions.tol
	else
		event = BK.ContinuousEvent(2, testBT_GH, false, ("bt", "gh"), threshBT)
	end

	prob_h = reMake(prob_h, recordFromSolution = _printsol2)

	# solve the hopf equations
	br = continuation(
		prob_h, alg,
		(@set opt_hopf_cont.newtonOptions.eigsolver = eigsolver);
		kwargs...,
		kind = BK.HopfCont(),
		linearAlgo = BorderingBLS(solver = opt_hopf_cont.newtonOptions.linsolver, checkPrecision = false),
		normC = normC,
		finaliseSolution = updateMinAugEveryStep ==0 ? get(kwargs, :finaliseSolution, BK.finaliseDefault) : updateMinAugHopf,
		event = event
	)
	@assert ~isnothing(br) "Empty branch!"
	return correctBifurcation(br)
end

function BK.continuationHopf(prob::ConstantDDEBifProblem,
						br::BK.AbstractBranchResult, ind_hopf::Int64,
						lens2::Lens,
						options_cont::ContinuationPar = br.contparams;
						alg = br.alg,
						startWithEigen = false,
						normC = norm,
						kwargs...)
	hopfpointguess = HopfPoint(br, ind_hopf)
	ω = hopfpointguess.p[2]
	bifpt = br.specialpoint[ind_hopf]

	@assert ~isnothing(br.eig) "The branch contains no eigen elements. This is strange because a Hopf point was detected. Please open an issue on the website."

	@assert ~isnothing(br.eig[1].eigenvecs) "The branch contains no eigenvectors for the Hopf point. Please provide one."

	ζ = geteigenvector(br.contparams.newtonOptions.eigsolver, br.eig[bifpt.idx].eigenvecs, bifpt.ind_ev)
	ζ ./= normC(ζ)
	ζad = conj.(ζ)

	p = bifpt.param
	parbif = BK.setParam(br, p)

	hopfpointguess = ArrayPartition(hopfpointguess.u, real(ζ), imag(ζ), hopfpointguess.p)

	if startWithEigen
		# computation of adjoint eigenvalue
		λ = Complex(0, ω)
		# jacobian at bifurcation point
		L = BK.jacobian(prob, bifpt.x, parbif)

		# jacobian adjoint at bifurcation point
		_Jt = ~BK.hasAdjoint(prob) ? adjoint(L) : jad(prob, bifpt.x, parbif)

		ζstar, λstar = BK.getAdjointBasis(_Jt, conj(λ), br.contparams.newtonOptions.eigsolver; nev = br.contparams.nev, verbose = false)
		ζad .= ζstar ./ dot(ζstar, ζ)
	end

	return BK.continuationHopf(br.prob, alg,
					hopfpointguess, parbif,
					BK.getLens(br), lens2,
					ζ, ζad,
					options_cont ;
					normC = normC,
					kwargs...)
end

function BK.continuationHopf(prob_vf::ConstantDDEBifProblem, alg::BK.AbstractContinuationAlgorithm,
				hopfpointguess::ArrayPartition, par,
				lens1::Lens, lens2::Lens,
				eigenvec, eigenvec_ad,
				options_cont::ContinuationPar ;
				updateMinAugEveryStep = 0,
				normC = norm,
				bdlinsolver::BK.AbstractBorderedLinearSolver = MatrixBLS(),
				jacobian_ma::Symbol = :autodiff,
				computeEigenElements = false,
				usehessian = true,
				massmatrix = LinearAlgebra.I,
				kwargs...) where {Tb, vectype}
	@assert lens1 != lens2 "Please choose 2 different parameters. You only passed $lens1"
	@assert lens1 == BK.getLens(prob_vf)

	# options for the Newton Solver inheritated from the ones the user provided
	options_newton = options_cont.newtonOptions
	threshBT = 100options_newton.tol

	hopfPb = HopfDDEProblem(
		prob_vf,
		BK._copy(eigenvec_ad),	# this is a ≈ null space of (J - iω I)^*
		BK._copy(eigenvec), 	# this is b ≈ null space of  J - iω I
		options_newton.linsolver,
		# do not change linear solver if user provides it
		@set bdlinsolver.solver = (isnothing(bdlinsolver.solver) ? options_newton.linsolver : bdlinsolver.solver);
		usehessian = usehessian,
		massmatrix = massmatrix)

	# Jacobian for the Hopf problem
	if jacobian_ma == :autodiff
		# hopfpointguess = vcat(hopfpointguess.u, hopfpointguess.p)
		prob_h = BK.HopfMAProblem(hopfPb, BK.AutoDiff(), hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
		opt_hopf_cont = @set options_cont.newtonOptions.linsolver = DefaultLS()
	elseif jacobian_ma == :finiteDifferencesMF
		hopfpointguess = vcat(hopfpointguess.u, hopfpointguess.p)
		prob_h = FoldMAProblem(hopfPb, FiniteDifferences(), hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
		opt_hopf_cont = @set options_cont.newtonOptions.linsolver = options_cont.newtonOptions.linsolver
	else
		prob_h = HopfMAProblem(hopfPb, nothing, hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
		opt_hopf_cont = @set options_cont.newtonOptions.linsolver = HopfLinearSolverMinAug()
	end

	# this functions allows to tackle the case where the two parameters have the same name
	lenses = BK.getLensSymbol(lens1, lens2)

	# current lyapunov coefficient
	eTb = eltype(hopfpointguess)
	hopfPb.l1 = Complex{eTb}(0, 0)
	hopfPb.BT = one(eTb)
	hopfPb.GH = one(eTb)

	# this function is used as a Finalizer
	# it is called to update the Minimally Augmented problem
	# by updating the vectors a, b
	function updateMinAugHopf(z, tau, step, contResult; kUP...)
		# we first check that the continuation step was successful
		# if not, we do not update the problem with bad information!
		success = get(kUP, :state, nothing).converged
		(~modCounter(step, updateMinAugEveryStep) || success == false) && return true
		x = getVec(z.u, hopfPb)	# hopf point
		p1, ω = getP(z.u, hopfPb)
		p2 = z.p		# second parameter
		newpar = set(par, lens1, p1)
		newpar = set(newpar, lens2, p2)

		a = hopfPb.a
		b = hopfPb.b

		# expression of the jacobian
		J_at_xp = jacobian(hopfPb.prob_vf, x, newpar)

		# compute new b
		T = typeof(p1)
		newb = hopfPb.linbdsolver(J_at_xp, a, b, T(0), hopfPb.zero, T(1); shift = Complex(0, -ω), Mass = hopfPb.massmatrix)[1]

		# compute new a
		JAd_at_xp = hasAdjoint(hopfPb) ? jad(hopfPb.prob_vf, x, newpar) : adjoint(J_at_xp)
		newa = hopfPb.linbdsolver(JAd_at_xp, b, a, T(0), hopfPb.zero, T(1); shift = Complex(0, ω), Mass = adjoint(hopfPb.massmatrix))[1]

		hopfPb.a .= newa ./ normC(newa)
		# do not normalize with dot(newb, hopfPb.a), it prevents BT detection
		hopfPb.b .= newb ./ normC(newb)

		# we stop continuation at Bogdanov-Takens points

		# CA NE DEVRAIT PAS ETRE ISSNOT?
		isbt = isnothing(contResult) ? true : isnothing(findfirst(x -> x.type in (:bt, :ghbt, :btgh), contResult.specialpoint))

		# if the frequency is null, this is not a Hopf point, we halt the process
		if abs(ω) < threshBT
			@warn "[Codim 2 Hopf - Finalizer] The Hopf curve seems to be close to a BT point: ω ≈ $ω. Stopping computations at ($p1, $p2). If the BT point is not detected, try lowering Newton tolerance or dsmax."
		end

		# call the user-passed finalizer
		finaliseUser = get(kwargs, :finaliseSolution, nothing)
		resFinal = isnothing(finaliseUser) ? true : finaliseUser(z, tau, step, contResult; prob = hopfPb, kUP...)

		return abs(ω) >= threshBT && isbt && resFinal
	end

	function testBT_GH(iter, state)
		z = getx(state)
		x = getVec(z, hopfPb)		# hopf point
		p1, ω = getP(z, hopfPb)		# first parameter
		p2 = getp(state)			# second parameter
		newpar = set(par, lens1, p1)
		newpar = set(newpar, lens2, p2)

		probhopf = iter.prob.prob

		a = probhopf.a
		b = probhopf.b

		# expression of the jacobian
		J_at_xp = BK.jacobian(probhopf.prob_vf, x, newpar)

		# compute new b
		T = typeof(p1)
		# ζ = probhopf.linbdsolver(J_at_xp, a, b, T(0), probhopf.zero, T(1); shift = Complex(0, -ω), Mass = hopfPb.massmatrix)[1]
		ζ = @. z.x[2] + im * z.x[3]
		ζ ./= normC(ζ)

		# compute new a
		# JAd_at_xp = BK.hasAdjoint(probhopf) ? jad(probhopf.prob_vf, x, newpar) : transpose(J_at_xp)
		# ζstar = probhopf.linbdsolver(JAd_at_xp, b, a, T(0), hopfPb.zero, T(1); shift = Complex(0, ω), Mass = transpose(hopfPb.massmatrix))[1]
		# test function for Bogdanov-Takens
		# JE ME DEMANDE SI CA NE FOUT PAS LA MERDE AVEC LA BISSECTION. EN EFFET BT EST DE SIGNE CONSTANT ET DONC LA BISSSECTION ET LE CHANGEMENT DE SIGNES VONT MERDER. IL SUFFIT PE DE PRENDRE ω - ϵ Newton
		probhopf.BT = ω
		# BT2 = real( dot(ζstar ./ normC(ζstar), ζ) )
		# ζstar ./= dot(ζ, ζstar)

		# hp = Hopf(x, p1, ω, newpar, lens1, ζ, ζstar, (a = Complex{T}(0,0), b = Complex{T}(0,0)), :hopf)
		# hopfNormalForm(prob_vf, hp, options_newton.linsolver, verbose = false)

		# lyapunov coefficient
		# probhopf.l1 = hp.nf.b
		# test for Bautin bifurcation.
		# If GH is too large, we take the previous value to avoid spurious detection
		# GH will be large close to BR points
		# probhopf.GH = abs(real(hp.nf.b)) < 1e5 ? real(hp.nf.b) : state.eventValue[2][2]
		return probhopf.BT, probhopf.GH
	end

	# the following allows to append information specific to the codim 2 continuation to the user data
	_printsol = get(kwargs, :recordFromSolution, nothing)
	_printsol2 = isnothing(_printsol) ?
		(u, p; kw...) -> (; zip(lenses, (getP(u, hopfPb)[1], p))..., ω = getP(u, hopfPb)[2], l1 = hopfPb.l1, BT = hopfPb.BT, GH = hopfPb.GH, BK.namedprintsol(BK.recordFromSolution(prob_vf)(getVec(u, hopfPb), p; kw...))...) :
		(u, p; kw...) -> (; BK.namedprintsol(_printsol(getVec(u, hopfPb), p; kw...))..., zip(lenses, (getP(u, hopfPb)[1], p))..., ω = getP(u, hopfPb)[2], l1 = hopfPb.l1, BT = hopfPb.BT, GH = hopfPb.GH)

	prob_h = reMake(prob_h, recordFromSolution = _printsol2)

	# eigen solver
	eigsolver = HopfDDEEig(BK.getsolver(opt_hopf_cont.newtonOptions.eigsolver))

	# event for detecting codim 2 points
	if computeEigenElements
		event = PairOfEvents(ContinuousEvent(2, testBT_GH, computeEigenElements, ("bt", "gh"), threshBT), BifDetectEvent)
		# careful here, we need to adjust the tolerance for stability to avoid
		# spurious ZH or HH bifurcations
		@set! opt_hopf_cont.tolStability = 10opt_hopf_cont.newtonOptions.tol
	else
		event = ContinuousEvent(2, testBT_GH, false, ("bt", "gh"), threshBT)
	end

	prob_h = reMake(prob_h, recordFromSolution = _printsol2)

	# solve the hopf equations
	br = continuation(
		prob_h, alg,
		(@set opt_hopf_cont.newtonOptions.eigsolver = eigsolver);
		kwargs...,
		kind = BK.HopfCont(),
		linearAlgo = BorderingBLS(solver = opt_hopf_cont.newtonOptions.linsolver, checkPrecision = false),
		normC = normC,
		finaliseSolution = updateMinAugEveryStep ==0 ? get(kwargs, :finaliseSolution, BK.finaliseDefault) : updateMinAugHopf,
		event = event
	)
	@assert ~isnothing(br) "Empty branch!"
	return correctBifurcation(br)
end

function BK.continuationHopf(prob::ConstantDDEBifProblem,
						br::BK.AbstractBranchResult, ind_hopf::Int64,
						lens2::Lens,
						options_cont::ContinuationPar = br.contparams;
						alg = br.alg,
						startWithEigen = false,
						normC = norm,
						kwargs...)
	hopfpointguess = HopfPoint(br, ind_hopf)
	ω = hopfpointguess.p[2]
	bifpt = br.specialpoint[ind_hopf]

	@assert ~isnothing(br.eig) "The branch contains no eigen elements. This is strange because a Hopf point was detected. Please open an issue on the website."

	@assert ~isnothing(br.eig[1].eigenvecs) "The branch contains no eigenvectors for the Hopf point. Please provide one."

	ζ = geteigenvector(br.contparams.newtonOptions.eigsolver, br.eig[bifpt.idx].eigenvecs, bifpt.ind_ev)
	ζ ./= normC(ζ)
	ζad = conj.(ζ)

	p = bifpt.param
	parbif = BK.setParam(br, p)

	hopfpointguess = ArrayPartition(hopfpointguess.u, real(ζ), imag(ζ), hopfpointguess.p)

	if startWithEigen
		# computation of adjoint eigenvalue
		λ = Complex(0, ω)
		# jacobian at bifurcation point
		L = BK.jacobian(prob, bifpt.x, parbif)

		# jacobian adjoint at bifurcation point
		_Jt = ~BK.hasAdjoint(prob) ? adjoint(L) : jad(prob, bifpt.x, parbif)

		ζstar, λstar = BK.getAdjointBasis(_Jt, conj(λ), br.contparams.newtonOptions.eigsolver; nev = br.contparams.nev, verbose = false)
		ζad .= ζstar ./ dot(ζstar, ζ)
	end

	return BK.continuationHopf(br.prob, alg,
					hopfpointguess, parbif,
					BK.getLens(br), lens2,
					ζ, ζad,
					options_cont ;
					normC = normC,
					kwargs...)
end


# structure to compute the eigenvalues along the Hopf branch
struct HopfDDEEig{S} <: BK.AbstractCodim2EigenSolver
	eigsolver::S
end

function (eig::HopfDDEEig)(Jddehopf, nev; kwargs...)
	xh = Jddehopf.x.x[1]			# hopf point
	p1, ω = Jddehopf.x.x[4]			# first parameter

	newpar = set(Jddehopf.p, BK.getLens(Jddehopf.prob.prob), p1)

	J = BK.jacobian(Jddehopf.prob.prob.prob_vf, xh, newpar)

	eigenelts = eig.eigsolver(J, nev; kwargs...)

	@show nev
	# @infiltrate
	return eigenelts
end
