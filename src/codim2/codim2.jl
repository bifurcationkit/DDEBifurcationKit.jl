abstract type AbstractCodim2DDEEigenSolver <: BK.AbstractEigenSolver end

for op in (:HopfDDEProblem,)
	@eval begin
		"""
		$(TYPEDEF)

		Structure to encode Hopf functional based for a DDE problem with constant delays.

		# Fields

		$(FIELDS)
		"""
		mutable struct $op{Tprob <: BK.AbstractBifurcationProblem, vectype, T <: Real, S <: BK.AbstractLinearSolver, Sa <: BK.AbstractLinearSolver, Sbd <: BK.AbstractBorderedLinearSolver, Sbda <: BK.AbstractBorderedLinearSolver, Tmass}
			"Functional F(x, p) - vector field - with all derivatives"
			prob_vf::Tprob
			"close to null vector of Jᵗ"
			a::vectype
			"close to null vector of J"
			b::vectype
			"vector zero, to avoid allocating it many times"
			zero::vectype
			"Lyapunov coefficient"
			l1::Complex{T}
			"Cusp test value"
			CP::T
			"Bogdanov-Takens test value"
			BT::T
			"Bautin test values"
			GH::T
			"Zero-Hopf test values"
			ZH::Int
			"linear solver. Used to invert the jacobian of MA functional"
			linsolver::S
			"linear solver for the jacobian adjoint"
			linsolverAdjoint::Sa
			"bordered linear solver"
			linbdsolver::Sbd
			"linear bordered solver for the jacobian adjoint"
			linbdsolverAdjoint::Sbda
			"wether to use the hessian of prob_vf"
			usehessian::Bool
			"wether to use a mass matrix M for studying M∂tu = F(u), default = I"
			massmatrix::Tmass
		end

		@inline hasHessian(pb::$op) = hasHessian(pb.prob_vf)
		@inline BK.isSymmetric(pb::$op) = false
		@inline hasAdjoint(pb::$op) = hasAdjoint(pb.prob_vf)
		@inline hasAdjointMF(pb::$op) = hasAdjointMF(pb.prob_vf)
		@inline BK.isInplace(pb::$op) = BK.isInplace(pb.prob_vf)
		@inline BK.getLens(pb::$op) = BK.getLens(pb.prob_vf)
		jad(pb::$op, args...) = jad(pb.prob_vf, args...)

		# constructor
		function $op(prob::ConstantDDEBifProblem, a, b, linsolve::BK.AbstractLinearSolver, linbdsolver = BK.MatrixBLS(); usehessian = true, massmatrix = LinearAlgebra.I)
			# determine scalar type associated to vectors a and b
			α = norm(a) # this is valid, see https://jutho.github.io/KrylovKit.jl/stable/#Package-features-and-alternatives-1
			Ty = eltype(α)
			return $op(prob, a, b, 0*a,
						complex(zero(Ty)),   # l1
						real(one(Ty)),		# cp
						real(one(Ty)),		# bt
						real(one(Ty)),		# gh
						1,							# zh
						linsolve, linsolve, linbdsolver, linbdsolver, usehessian, massmatrix)
		end
	end
end
################################################################################
"""
$(SIGNATURES)

This function turns an initial guess for a Fold/Hopf point into a solution to the Fold/Hopf problem based on a Minimally Augmented formulation.

## Arguments
- `br` results returned after a call to [continuation](@ref Library-Continuation)
- `ind_bif` bifurcation index in `br`

# Optional arguments:
- `options::NewtonPar`, default value `br.contparams.newtonOptions`
- `normN = norm`
- `options` You can pass newton parameters different from the ones stored in `br` by using this argument `options`.
- `bdlinsolver` bordered linear solver for the constraint equation
- `startWithEigen = false` whether to start the Minimally Augmented problem with information from eigen elements.
- `kwargs` keywords arguments to be passed to the regular Newton-Krylov solver

!!! tip "ODE problems"
    For ODE problems, it is more efficient to use the Bordered Linear Solver using the option `bdlinsolver = MatrixBLS()`

!!! tip "startWithEigen"
    It is recommanded that you use the option `startWithEigen=true`
"""
function BK.newton(br::BK.AbstractResult{Tkind, Tprob}, ind_bif::Int64; normN = norm, options = br.contparams.newtonOptions, startWithEigen = false, lens2::Lens = (@lens _), kwargs...) where {Tkind, Tprob <: ConstantDDEBifProblem}
	@assert length(br.specialpoint) > 0 "The branch does not contain bifurcation points"
	if br.specialpoint[ind_bif].type == :hopf
		return newtonHopf(br, ind_bif; normN = normN, options = options, startWithEigen = startWithEigen, kwargs...)
	elseif br.specialpoint[ind_bif].type == :bt
		return newtonBT(br, ind_bif; lens2 = lens2, normN = normN, options = options, startWithEigen = startWithEigen, kwargs...)
	else
		return newtonFold(br, ind_bif; normN = normN, options = options, startWithEigen = startWithEigen, kwargs...)
	end
end
################################################################################
"""
$(SIGNATURES)

Codimension 2 continuation of Fold / Hopf points. This function turns an initial guess for a Fold/Hopf point into a curve of Fold/Hopf points based on a Minimally Augmented formulation. The arguments are as follows
- `br` results returned after a call to [continuation](@ref Library-Continuation)
- `ind_bif` bifurcation index in `br`
- `lens2` second parameter used for the continuation, the first one is the one used to compute `br`, e.g. `getLens(br)`
- `options_cont = br.contparams` arguments to be passed to the regular [continuation](@ref Library-Continuation)

# Optional arguments:
- `bdlinsolver` bordered linear solver for the constraint equation
- `updateMinAugEveryStep` update vectors `a,b` in Minimally Formulation every `updateMinAugEveryStep` steps
- `startWithEigen = false` whether to start the Minimally Augmented problem with information from eigen elements
- `detectCodim2Bifurcation ∈ {0,1,2}` whether to detect Bogdanov-Takens, Bautin and Cusp. If equals `1` non precise detection is used. If equals `2`, a bisection method is used to locate the bifurcations.
- `kwargs` keywords arguments to be passed to the regular [continuation](@ref Library-Continuation)

where the parameters are as above except that you have to pass the branch `br` from the result of a call to `continuation` with detection of bifurcations enabled and `index` is the index of Hopf point in `br` you want to refine.

!!! tip "ODE problems"
    For ODE problems, it is more efficient to pass the Bordered Linear Solver using the option `bdlinsolver = MatrixBLS()`

!!! tip "startWithEigen"
    It is recommanded that you use the option `startWithEigen = true`
"""
function BK.continuation(br::BK.AbstractResult{Tkind, Tprob},
					ind_bif::Int64,
					lens2::Lens,
					options_cont::ContinuationPar = br.contparams ;
					startWithEigen = false,
					detectCodim2Bifurcation::Int = 0,
					kwargs...) where {Tkind, Tprob <: ConstantDDEBifProblem}
	@assert length(br.specialpoint) > 0 "The branch does not contain bifurcation points"
	# options to detect codim2 bifurcations
	computeEigenElements = options_cont.detectBifurcation > 0
	_options_cont = BK.detectCodim2Parameters(detectCodim2Bifurcation, options_cont; kwargs...)

	if br.specialpoint[ind_bif].type == :hopf
		return continuationHopf(br.prob, br, ind_bif, lens2, _options_cont;
			startWithEigen = startWithEigen,
			computeEigenElements = computeEigenElements,
			kwargs...)
	else
		return continuationFold(br.prob, br, ind_bif, lens2, _options_cont;
			startWithEigen = startWithEigen,
			computeEigenElements = computeEigenElements,
			kwargs...)
	end
end
################################################################################
"""
$(SIGNATURES)

This function uses information in the branch to detect codim 2 bifurcations like BT, ZH and Cusp.
"""
function correctBifurcation(contres::ContResult)
	if contres.prob.prob isa HopfDDEProblem == false
		return contres
	end
	if contres.prob.prob isa FoldProblemMinimallyAugmented
		conversion = Dict(:bp => :bt, :hopf => :zh, :fold => :cusp, :nd => :nd)
	elseif contres.prob.prob isa HopfDDEProblem
		conversion = Dict(:bp => :zh, :hopf => :hh, :fold => :nd, :nd => :nd, :ghbt => :bt, :btgh => :bt)
	else
		throw("Error! this should not occur. Please open an issue on the website of BifurcationKit.jl")
	end
	for (ind, bp) in pairs(contres.specialpoint)
		if bp.type in keys(conversion)
			@set! contres.specialpoint[ind].type = conversion[bp.type]
		end
	end
	return contres
end