abstract type AbstractDDEBifurcationProblem <: BK.AbstractBifurcationProblem end

"""
$(TYPEDEF)

Structure to hold the bifurcation problem. If don't have parameters, you can pass `nothing`.

## Internal fields

$(TYPEDFIELDS)

## Methods

- `getu0(pb)` calls `pb.u0`
- `getparams(pb)` calls `pb.params`
- `getlens(pb)` calls `pb.lens`
- `setparam(pb, p0)` calls `set(pb.params, pb.lens, p0)`
- `record_from_solution(pb)` calls `pb.record_from_solution`

## Constructors

- `ConstantDDEBifProblem(F, delays, u0, params, lens; J, Jᵗ, d2F, d3F, kwargs...)` and `kwargs` are the fields above.

"""
struct ConstantDDEBifProblem{Tvf, Tdf, Tu, Td, Tp, Tl <: Union{Nothing, BK.AllOpticTypes}, Tplot, Trec, Tgets, Tδ} <: AbstractDDEBifurcationProblem
    "Vector field, typically a [`BifFunction`](@ref). For more information, please look at the website https://bifurcationkit.github.io/DDEBifurcationKit.jl/dev/BifProblem."
    VF::Tvf
    "function delays. It takes the parameters and return the non-zero delays in an `AbstractVector` form. Example: `delays = par -> [1, 2]`."
    delays::Tdf
    "Initial guess"
    u0::Tu
    "[internal] initial delays (set internally by the constructor)."
    delays0::Td
    "parameters."
    params::Tp
    "Typically a `Accessors.@optic`. It specifies which parameter axis among `params` is used for continuation. For example, if `par = (α = 1.0, β = 1)`, we can perform continuation w.r.t. `α` by using `lens = (@optic _.α)`. If you have an array `par = [ 1.0, 2.0]` and want to perform continuation w.r.t. the first variable, you can use `lens = (@optic _[1])`. For more information, we refer to `Accessors.jl`."
    lens::Tl
    "user function to plot solutions during continuation. Signature: `plotSolution(x, p; kwargs...)`."
    plotSolution::Tplot
    "`record_from_solution = (x, p; k...) -> norm(x)` function used record a few indicators about the solution. It could be `norm` or `(x, p) -> x[1]`. This is also useful when saving several huge vectors is not possible for memory reasons (for example on GPU...). This function can return pretty much everything but you should keep it small. For example, you can do `(x, p) -> (x1 = x[1], x2 = x[2], nrm = norm(x))` or simply `(x, p) -> (sum(x), 1)`. This will be stored in `contres.branch` (see below). Finally, the first component is used to plot in the continuation curve."
    recordFromSolution::Trec
    "function to save the full solution on the branch. Some problem are mutable (like periodic orbit functional with adaptive mesh) and this function allows to save the state of the problem along with the solution itself. Signature: `save_solution(x, p)`."
    save_solution::Tgets
    "delta for Finite differences."
    δ::Tδ
end

BK.isinplace(::ConstantDDEBifProblem) = false
BK.is_symmetric(::ConstantDDEBifProblem) = false
BK._getvectortype(prob::ConstantDDEBifProblem{Tvf, Tdf, Tu, Td, Tp, Tl, Tplot, Trec}) where {Tvf, Tdf, Tu, Td, Tp, Tl, Tplot, Trec} = Tu
BK.getlens(prob::ConstantDDEBifProblem) = prob.lens
BK.has_adjoint(prob::ConstantDDEBifProblem) = true
BK.has_adjoint_MF(prob::ConstantDDEBifProblem) = false
BK.getdelta(prob::ConstantDDEBifProblem) = prob.δ
BK.save_solution(prob::ConstantDDEBifProblem, x, p) = prob.save_solution(x, p)
@inline delays(prob::ConstantDDEBifProblem, x, pars) = prob.delays(pars)

function Base.show(io::IO, prob::ConstantDDEBifProblem; prefix = "")
    print(io, prefix * "┌─ Constant Delays Bifurcation Problem with uType ")
    printstyled(io, BK._getvectortype(prob), color=:cyan, bold = true)
    print(io, prefix * "\n├─ # delays : ", length(prob.delays0))
    print(io, prefix * "\n├─ Inplace  : ")
    printstyled(io, BK.isinplace(prob), color=:cyan, bold = true)
    print(io, "\n" * prefix * "├─ Dimension: ")
    printstyled(io, length(BK.getu0(prob)), color = :cyan, bold = true)
    print(io, "\n" * prefix * "└─ Parameter: ")
    printstyled(io, BK.get_lens_symbol(getlens(prob)), color=:cyan, bold = true)
end

function ConstantDDEBifProblem(F, delayF, u0, parms, lens = (@optic _);
                F! = nothing,
                J! = nothing,
                jvp = nothing,
                vjp = nothing,
                J = nothing,
                Jᵗ = nothing,
                d2F = nothing,
                d2Fc = nothing,
                d3F = nothing,
                d3Fc = nothing,
                issymmetric::Bool = false,
                record_from_solution = BK.record_sol_default,
                plot_solution = BK.plot_default,
                save_solution = BK.save_solution_default,
                inplace = false,
                δ = convert(eltype(u0), 1e-8),
                kwargs_jet...
                )
    @assert lens isa Int || lens isa BK.AllOpticTypes

    VF = BifFunction(F, F!, jvp, vjp, J, Jᵗ, nothing, d2F, d3F, d2Fc, d3Fc, issymmetric, δ, inplace, BK.Jet(;kwargs_jet...)) ## TODO: it requires a specific DDEBifFunction
    return ConstantDDEBifProblem(VF,
                                 delayF,
                                 u0,
                                 delayF(parms),
                                 parms,
                                 lens,
                                 plot_solution,
                                 record_from_solution,
                                 save_solution,
                                 δ)
end

function F_voa(prob::ConstantDDEBifProblem, x, p)
    n_delays = length(delays(prob, x, p))
    F_voa(prob, VectorOfArray([copy(x) for _ in 1:(n_delays+1)]), p)
end

function F_voa(prob::ConstantDDEBifProblem, xd::VectorOfArray, p)
    prob.VF.F(xd.u[1], VectorOfArray(xd.u[2:end]), p)
end

function BK.dF(prob::ConstantDDEBifProblem{ <: BifFunction{Tf, TFinp, Nothing}}, x, p, dx) where {Tf, TFinp}
    𝒯 = BK.VI.scalartype(x)
    return ForwardDiff.derivative(t -> F_voa(prob, x .+ t .* dx, p), zero(𝒯))
end

BK.update!(prob::ConstantDDEBifProblem, args...; kwargs...) = BK.update_default(args...; kwargs...)

function BK.residual(prob::ConstantDDEBifProblem, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    prob.VF.F(x, xd, p)
end

function BK.residual!(prob::ConstantDDEBifProblem, o, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    o .= prob.VF.F(x, xd, p)
    o
end

function jacobian_forwarddiff(prob::ConstantDDEBifProblem, x, p)
    xd = VectorOfArray([copy(x) for _ in eachindex(prob.delays0)])
    J0 = ForwardDiff.jacobian(z -> prob.VF.F(z, xd, p), x)
    Jd = [ ForwardDiff.jacobian(z -> prob.VF.F(x, (@set xd.u[ii] = z), p), x) for ii in eachindex(prob.delays0)]
    return J0, Jd
end

BK.jacobian(prob::ConstantDDEBifProblem, x, pars) = jacobian_cst_dde(prob.VF, prob, x, pars)

function jacobian_cst_dde(VF::BifFunction, prob, x, pars)
    J0, Jd = VF.J(x, pars)
    return JacobianDDE(prob, J0, Jd, prob.delays(pars))
end

function jacobian_cst_dde(VF::BifFunction{Tf, TFinp, Tdf, Tdfad, Nothing}, prob, x, pars) where {Tf, TFinp, Tdf, Tdfad}
    J0, Jd = jacobian_forwarddiff(prob, x, pars)
    return JacobianDDE(prob, J0, Jd, prob.delays(pars))
end

function BK.jacobian_adjoint(prob::AbstractDDEBifurcationProblem, x, p)
    J = BK.jacobian(prob, x, p)
    if isnothing(J.Jall) == false
        J.Jall .= J.Jall'
    end
    J.J0 .= J.J0'
    for _J in J.Jd
        _J .= _J'
    end
    J
end

function jacobian(prob::ConstantDDEBifProblem, x, xd, p)
    J0 = ForwardDiff.jacobian(z -> prob.VF.F(z, xd, p), x)
    Jd = [ ForwardDiff.jacobian(z -> prob.VF.F(x, (@set xd.u[ii] = z), p), xd.u[ii]) for ii in eachindex(prob.delays0)]
    return JacobianDDE(prob, missing, J0, Jd, prob.delays(p))
end

function BK.d2F(prob::ConstantDDEBifProblem{ <: BifFunction{Tf, TFinp, Nothing}}, x, p, dx1, dx2) where {Tf, TFinp}
    𝒯 = BK.VI.scalartype(x)
    ForwardDiff.derivative(t -> BK.dF(prob, x .+ t .* dx2, p, dx1), zero(𝒯))
end

function BK.d3F(prob::ConstantDDEBifProblem{ <: BifFunction{Tf, TFinp, Nothing}}, x, p, dx1, dx2, dx3) where {Tf, TFinp}
    𝒯 = BK.VI.scalartype(x)
    ForwardDiff.derivative(t -> BK.d2F(prob, x .+ t .* dx3, p, dx1, dx2), zero(𝒯))
end
####################################################################################################
"""
$(TYPEDEF)

Structure to hold the bifurcation problem. If don't have parameters, you can pass `nothing`.

## Internal fields

$(TYPEDFIELDS)

## Methods

- `getu0(pb)` calls `pb.u0`
- `getparams(pb)` calls `pb.params`
- `getlens(pb)` calls `pb.lens`
- `setparam(pb, p0)` calls `set(pb.params, pb.lens, p0)`
- `record_from_solution(pb)` calls `pb.record_from_solution`

## Constructors

- `SDDDEBifProblem(F, delays, u0, params, lens; J, Jᵗ, d2F, d3F, kwargs...)` and `kwargs` are the fields above.

"""
struct SDDDEBifProblem{Tvf, Tdf, Tu, Td, Tp, Tl <: Union{Nothing, BK.AllOpticTypes}, Tplot, Trec, Tgets, Tδ} <: AbstractDDEBifurcationProblem
    "Vector field, typically a [`BifFunction`](@ref). For more information, please look at the website https://bifurcationkit.github.io/DDEBifurcationKit.jl/dev/BifProblem."
    VF::Tvf
    "function delays. It takes the state and the parameters and return the non-zero delays in an `AsbtractVector` form. Example: `delays = (u, pars) -> [1 + u[1]^2]`."
    delays::Tdf
    "Initial guess."
    u0::Tu
    "[internal] initial delays (set internally by the constructor)."
    delays0::Td
    "parameters"
    params::Tp
    "see ConstantDDEBifProblem for more information."
    lens::Tl
    "user function to plot solutions during continuation. Signature: `plotSolution(x, p; kwargs...)`."
    plotSolution::Tplot
    "`record_from_solution = (x, p; k...) -> norm(x)` function used record a few indicators about the solution. It could be `norm` or `(x, p) -> x[1]`. This is also useful when saving several huge vectors is not possible for memory reasons (for example on GPU...). This function can return pretty much everything but you should keep it small. For example, you can do `(x, p) -> (x1 = x[1], x2 = x[2], nrm = norm(x))` or simply `(x, p) -> (sum(x), 1)`. This will be stored in `contres.branch` (see below). Finally, the first component is used to plot in the continuation curve."
    recordFromSolution::Trec
    "function to save the full solution on the branch. Some problem are mutable (like periodic orbit functional with adaptive mesh) and this function allows to save the state of the problem along with the solution itself. Signature: `save_solution(x, p)`."
    save_solution::Tgets
    "delta for Finite differences."
    δ::Tδ
end

BK.isinplace(::SDDDEBifProblem) = false
BK.is_symmetric(::SDDDEBifProblem) = false
BK._getvectortype(prob::SDDDEBifProblem{Tvf, Tdf, Tu, Td, Tp, Tl, Tplot, Trec}) where {Tvf, Tdf, Tu, Td, Tp, Tl, Tplot, Trec} = Tu
BK.getlens(prob::SDDDEBifProblem) = prob.lens
BK.has_adjoint(prob::SDDDEBifProblem) = true
BK.getdelta(prob::SDDDEBifProblem) = prob.δ
BK.save_solution(prob::SDDDEBifProblem, x, p) = prob.save_solution(x, p)
@inline delays(prob::SDDDEBifProblem, x, pars) = prob.delays(x, pars)

function SDDDEBifProblem(F, delayF, u0, parms, lens = (@optic _);
                F! = nothing,
                jvp = nothing,
                vjp = nothing,
                J = nothing,
                Jᵗ = nothing,
                J! = nothing,
                d2F = nothing,
                d2Fc = nothing,
                d3F = nothing,
                d3Fc = nothing,
                issymmetric::Bool = false,
                record_from_solution = BK.record_sol_default,
                plot_solution = BK.plot_default,
                save_solution = BK.save_solution_default,
                inplace = false,
                δ = convert(eltype(u0), 1e-8),
                kwargs_jet...
                )
    @assert lens isa Int || lens isa BK.AllOpticTypes
    VF = BifFunction(F, F!, jvp, vjp, J, Jᵗ, nothing, d2F, d3F, d2Fc, d3Fc, issymmetric, δ, inplace, BK.Jet(;kwargs_jet...))
    return SDDDEBifProblem(VF,
                           delayF,
                           u0,
                           delayF(u0, parms),
                           parms,
                           lens,
                           plot_solution,
                           record_from_solution,
                           save_solution,
                           δ)
end

BK.update!(prob::SDDDEBifProblem, args...; kwargs...) = BK.update_default(args...; kwargs...)

function Base.show(io::IO, prob::SDDDEBifProblem; prefix = "")
    print(io, prefix * "┌─ State-dependent delays Bifurcation Problem with uType ")
    printstyled(io, BK._getvectortype(prob), color=:cyan, bold = true)
    print(io, prefix * "\n├─ Inplace:  ")
    printstyled(io, BK.isinplace(prob), color=:cyan, bold = true)
    print(io, "\n" * prefix * "├─ Dimension: ")
    printstyled(io, length(BK.getu0(prob)), color = :cyan, bold = true)
    print(io, "\n" * prefix * "└─ Parameter: ")
    printstyled(io, BK.get_lens_symbol(getlens(prob)), color=:cyan, bold = true)
end

function BK.residual(prob::SDDDEBifProblem, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    prob.VF.F(x, xd, p)
end

function BK.residual!(prob::SDDDEBifProblem, o, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    o .= prob.VF.F(x, xd, p)
    o
end

function BK.jacobian(prob::SDDDEBifProblem, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    J0 = ForwardDiff.jacobian(z -> prob.VF.F(z, xd, p), x)
    Jd = [ ForwardDiff.jacobian(z -> prob.VF.F(x, (@set xd.u[ii] = z), p), x) for ii in eachindex(prob.delays0)]
    return JacobianDDE(prob, J0 + sum(Jd), J0, Jd, prob.delays(x, p))
end

function jacobian(prob::SDDDEBifProblem, x, xd, p)
    J0 = ForwardDiff.jacobian(z -> prob.VF.F(z, xd, p), x)
    Jd = [ ForwardDiff.jacobian(z -> prob.VF.F(x, (@set xd.u[ii] = z), p), xd.u[ii]) for ii in eachindex(prob.delays0)]
    return JacobianDDE(prob, missing, J0, Jd, prob.delays(x, p))
end