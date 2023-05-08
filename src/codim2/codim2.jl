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

		@inline BK.hasHessian(pb::$op) = BK.hasHessian(pb.prob_vf)
		@inline BK.isSymmetric(pb::$op) = false
		@inline BK.hasAdjoint(pb::$op) = BK.hasAdjoint(pb.prob_vf)
		@inline BK.hasAdjointMF(pb::$op) = BK.hasAdjointMF(pb.prob_vf)
		@inline BK.isInplace(pb::$op) = BK.isInplace(pb.prob_vf)
		@inline BK.getLens(pb::$op) = BK.getLens(pb.prob_vf)
		jad(pb::$op, args...) = jad(pb.prob_vf, args...)

		# constructor
		function $op(prob::AbstractDDEBifurcationProblem, a, b, linsolve::BK.AbstractLinearSolver, linbdsolver = BK.MatrixBLS(); usehessian = true, massmatrix = LinearAlgebra.I)
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

This function uses information in the branch to detect codim 2 bifurcations like BT, ZH and Cusp.
"""
function correctBifurcation(contres::ContResult)
	if contres.prob.prob isa HopfDDEProblem == false
		return contres
	end
	if contres.prob.prob isa HopfDDEProblem
		conversion = Dict(:bp => :zh, :hopf => :hh, :fold => :nd, :nd => :nd, :ghbt => :bt, :btgh => :bt, :ghbp => :zh)
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
