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

        @inline BK.has_hessian(pb::$op) = BK.hasHessian(pb.prob_vf)
        @inline BK.is_symmetric(pb::$op) = false
        @inline BK.has_adjoint(pb::$op) = BK.has_adjoint(pb.prob_vf)
        @inline BK.has_adjoint_MF(pb::$op) = BK.has_adjoint_MF(pb.prob_vf)
        @inline BK.isinplace(pb::$op) = BK.isinplace(pb.prob_vf)
        @inline BK.getlens(pb::$op) = BK.getlens(pb.prob_vf)
        jad(pb::$op, args...) = jad(pb.prob_vf, args...)
        @inline BK.getdelta(pb::$op) = BK.getdelta(pb.prob_vf)

        # constructor
        function $op(prob::AbstractDDEBifurcationProblem, a, b, linsolve::BK.AbstractLinearSolver, linbdsolver = BK.MatrixBLS(); usehessian = true, massmatrix = LinearAlgebra.I)
            # determine scalar type associated to vectors a and b
            α = norm(a) # this is valid, see https://jutho.github.io/KrylovKit.jl/stable/#Package-features-and-alternatives-1
            Ty = eltype(α)
            return $op(prob, a, b, 0*a,
                        complex(zero(Ty)),   # l1
                        real(one(Ty)),        # cp
                        real(one(Ty)),        # bt
                        real(one(Ty)),        # gh
                        1,                            # zh
                        linsolve, linsolve, linbdsolver, linbdsolver, usehessian, massmatrix)
        end
    end
end
################################################################################
#
# function BK.newton(br::BK.AbstractResult{Tkind, Tprob}, ind_bif::Int64; normN = norm, options = br.contparams.newtonOptions, startWithEigen = false, lens2::Lens = (@lens _), kwargs...) where {Tkind, Tprob <: ConstantDDEBifProblem}
#     @assert length(br.specialpoint) > 0 "The branch does not contain bifurcation points"
#     if br.specialpoint[ind_bif].type == :hopf
#         return newtonHopf(br, ind_bif; normN = normN, options = options, startWithEigen = startWithEigen, kwargs...)
#     elseif br.specialpoint[ind_bif].type == :bt
#         return newtonBT(br, ind_bif; lens2 = lens2, normN = normN, options = options, startWithEigen = startWithEigen, kwargs...)
#     else
#         return newtonFold(br, ind_bif; normN = normN, options = options, startWithEigen = startWithEigen, kwargs...)
#     end
# end
################################################################################
# function BK.continuation(br::BK.AbstractResult{Tkind, Tprob},
#                     ind_bif::Int64,
#                     lens2::Lens,
#                     options_cont::ContinuationPar = br.contparams ;
#                     startWithEigen = false,
#                     detectCodim2Bifurcation::Int = 0,
#                     kwargs...) where {Tkind, Tprob <: ConstantDDEBifProblem}
#     @assert length(br.specialpoint) > 0 "The branch does not contain bifurcation points"
#     # options to detect codim2 bifurcations
#     computeEigenElements = options_cont.detectBifurcation > 0
#     _options_cont = BK.detectCodim2Parameters(detectCodim2Bifurcation, options_cont; kwargs...)
#
#     if br.specialpoint[ind_bif].type == :hopf
#         return continuationHopf(br.prob, br, ind_bif, lens2, _options_cont;
#             startWithEigen = startWithEigen,
#             computeEigenElements = computeEigenElements,
#             kwargs...)
#     else
#         return continuationFold(br.prob, br, ind_bif, lens2, _options_cont;
#             startWithEigen = startWithEigen,
#             computeEigenElements = computeEigenElements,
#             kwargs...)
#     end
# end
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
            @reset contres.specialpoint[ind].type = conversion[bp.type]
        end
    end
    return contres
end
