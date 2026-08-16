abstract type AbstractCodim2DDEEigenSolver <: BK.AbstractEigenSolver end

for op in (:HopfDDEFormulation,)
    @eval begin
        """
        $(TYPEDEF)

        Structure to encode Hopf functional based for a DDE problem with constant delays.

        # Fields

        $(FIELDS)
        """
        mutable struct $op{Tprob <: BK.AbstractBifurcationProblem, vectype, T <: Real, S <: BK.AbstractLinearSolver, Sa <: BK.AbstractLinearSolver, Sbd <: BK.AbstractBorderedLinearSolver, Sbda <: BK.AbstractBorderedLinearSolver, Tmass, Tn}
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
            "norm to normalize vector in update or test"
            norm::Tn
            "Update the problem every such step"
            update_minaug_every_step::Int
        end

        @inline BK.has_hessian(pb::$op) = BK.hasHessian(pb.prob_vf)
        @inline BK.is_symmetric(pb::$op) = false
        @inline BK.has_adjoint(pb::$op) = BK.has_adjoint(pb.prob_vf)
        @inline BK.has_adjoint_MF(pb::$op) = BK.has_adjoint_MF(pb.prob_vf)
        @inline BK.isinplace(pb::$op) = BK.isinplace(pb.prob_vf)
        @inline BK.getlens(pb::$op) = BK.getlens(pb.prob_vf)
        jad(pb::$op, args...) = jad(pb.prob_vf, args...)
        @inline BK.getdelta(pb::$op) = BK.getdelta(pb.prob_vf)
        @inline BK.getparams(pb::$op) = BK.getparams(pb.prob_vf)

        # constructor
        function $op(prob::AbstractDDEBifurcationProblem, a, b,
                    linsolve::BK.AbstractLinearSolver,
                    linbdsolver = BK.MatrixBLS(); usehessian = true,
                    massmatrix = LA.I,
                    _norm = LA.norm,
                    update_minaug_every_step = 0)
            # determine scalar type associated to vectors a and b
            α = LA.norm(a) # this is valid, see https://jutho.github.io/KrylovKit.jl/stable/#Package-features-and-alternatives-1
            Ty = eltype(α)
            return $op(prob, a, b, 0*a,
                        complex(zero(Ty)),   # l1
                        real(one(Ty)),        # cp
                        real(one(Ty)),        # bt
                        real(one(Ty)),        # gh
                        1,                            # zh
                        linsolve, linsolve, linbdsolver, linbdsolver, usehessian, massmatrix, _norm, update_minaug_every_step)
        end
    end
end

@inline BK.getp(x, ::HopfDDEFormulation) = BK.get_par_bls(x, 2)
@inline BK.getvec(x, 𝐏𝐛::HopfDDEFormulation) = getVec(x, 𝐏𝐛)
@inline BK.get_parameter(x, 𝐏𝐛::HopfDDEFormulation) = BK.getp(x, 𝐏𝐛)[1]
@inline BK.get_frequency(x, 𝐏𝐛::HopfDDEFormulation) = BK.getp(x, 𝐏𝐛)[2]

BK.update!(𝐏𝐛::HopfDDEFormulation, iter::BK.ContIterable, state::BK.ContState) = true

function BK.save_solution(𝐏𝐛::HopfDDEFormulation, x, p2)
    p1 = BK.get_parameter(x, 𝐏𝐛)
    # TODO!! is it a copy or else?
    x_ma = BK.save_solution(𝐏𝐛.prob_vf, BK.getvec(x, 𝐏𝐛), p2)
    return BK.MASolutionFreq(x_ma, p1, BK.get_frequency(x, 𝐏𝐛))
end

function BK.re_make(𝐌𝐚::HopfDDEFormulation;
                 params = BK.getparams(𝐌𝐚))
    new_prob = BK.re_make(𝐌𝐚.prob_vf; params)
    return (@set 𝐌𝐚.prob_vf = new_prob)
end

################################################################################
struct JacobianCodim2DDE{T1,T2,T3,T4}
    prob::T1
    J::T2
    x::T3
    p::T4
end

(l::BK.DefaultLS)(J::JacobianCodim2DDE, args...; kw...) = l(J.J, args...; kw...)

BK.jacobian(hopfpb::BK.HopfMAProblem{Tprob, BK.AutoDiff, Tl, Tplot, Trecord}, x, p) where {Tprob <: HopfDDEFormulation, Tl <: Union{BK.AllOpticTypes, Nothing}, Tplot, Trecord} = JacobianCodim2DDE(hopfpb, ForwardDiff.jacobian(z -> hopfpb.prob(z, p), x), x, p)
################################################################################
"""
$(SIGNATURES)

This function uses information in the branch to detect codim 2 bifurcations like BT, ZH and Cusp.
"""
function correctBifurcation(contres::ContResult)
    if contres.prob.prob isa HopfDDEFormulation == false
        return contres
    end
    if contres.prob.prob isa HopfDDEFormulation
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
