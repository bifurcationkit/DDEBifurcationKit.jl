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
BK.d2F(prob::ConstantDDEBifProblem, x, p, dx1, dx2) = BK.d2F(prob.VF, x, p, dx1, dx2)
BK.d3F(prob::ConstantDDEBifProblem, x, p, dx1, dx2, dx3) = BK.d3F(prob.VF, x, p, dx1, dx2, dx3)
BK.save_solution(prob::ConstantDDEBifProblem, x, p) = prob.save_solution(x, p)
@inline delays(prob::ConstantDDEBifProblem, x, pars) = prob.delays(pars)

function Base.show(io::IO, prob::ConstantDDEBifProblem; prefix = "")
    print(io, prefix * "┌─ Constant Delays Bifurcation Problem with uType ")
    printstyled(io, BK._getvectortype(prob), color=:cyan, bold = true)
    print(io, prefix * "\n├─ Inplace:  ")
    printstyled(io, BK.isinplace(prob), color=:cyan, bold = true)
    print(io, "\n" * prefix * "└─ Parameter: ")
    printstyled(io, BK.get_lens_symbol(getlens(prob)), color=:cyan, bold = true)
end

function ConstantDDEBifProblem(F, delayF, u0, parms, lens = (@optic _);
                F! = nothing,
                J! = nothing,
                dF = nothing,
                dFad = nothing,
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
    F_voa = (xd, p) -> F(xd.u[1], VectorOfArray(xd.u[2:end]), p) # VectorOfArray version of F
    𝒯 = BK.VI.scalartype(u0)
    # J = isnothing(J) ? (x,p) -> ForwardDiff.jacobian(z -> F(z, p), x) : J
    dF = isnothing(dF) ? (x, p, dx) -> ForwardDiff.derivative(t -> F_voa(x .+ t .* dx, p), zero(𝒯)) : dF
    jvp(x, p, dx1) = ForwardDiff.derivative(t -> F_voa(x .+ t .* dx1, p), zero(𝒯))
    if isnothing(d2F)
        d2F = (x, p, dx1, dx2) -> ForwardDiff.derivative(t -> jvp(x .+ t .* dx2, p, dx1), zero(𝒯))
        d2Fc = (x, p, dx1, dx2) -> BilinearMap((_dx1, _dx2) -> d2F(x,p,_dx1,_dx2))(dx1,dx2)
    else
        d2Fc = d2F
    end

    if isnothing(d3F)
        d3F  = (x, p, dx1, dx2, dx3) -> ForwardDiff.derivative(t -> d2F(x .+ t .* dx3, p, dx1, dx2), zero(𝒯))
        d3Fc = (x, p, dx1, dx2, dx3) -> TrilinearMap((_dx1, _dx2, _dx3) -> d3F(x,p,_dx1,_dx2,_dx3))(dx1, dx2, dx3)
    else
        d3Fc = d3F
    end

    d3F = isnothing(d3F) ? (x,p,dx1,dx2,dx3) -> ForwardDiff.derivative(t -> d2F(x .+ t .* dx3, p, dx1, dx2), zero(𝒯)) : d3F
    # VF = BifFunction(F, dF, dFad, J, Jᵗ, d2F, d3F, d2Fc, d3Fc, issymmetric, 1e-8, inplace)
    VF = BifFunction(F, F!, dF, dFad, J, Jᵗ, nothing, d2F, d3F, d2Fc, d3Fc, issymmetric, δ, inplace, BK.Jet(;kwargs_jet...))
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

BK.dF(prob::ConstantDDEBifProblem, x,p, dx) = BK.dF(prob.VF, x, p, dx)
BK.update!(prob::ConstantDDEBifProblem, args...; kwargs...) = BK.update_default(args...; kwargs...)

struct JacobianDDE{Tp,T1,T2,T3,Td}
    prob::Tp
    Jall::T1
    J0::T2
    Jd::T3
    delays::Td
end

function BK.residual(prob::ConstantDDEBifProblem, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    prob.VF.F(x,xd,p)
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

function BK.jacobian(prob::ConstantDDEBifProblem, x, p)
    if isnothing(prob.VF.J)
        J0, Jd = jacobian_forwarddiff(prob, x, p)
    else
        J0, Jd = prob.VF.J(x, p)
    end
    return JacobianDDE(prob, J0 + sum(Jd), J0, Jd, prob.delays(p))
end

function BK.jacobian_adjoint(prob::ConstantDDEBifProblem, x, p)
    J = BK.jacobian(prob, x, p)
    J.Jall .= J.Jall'
    J.J0 .= J.J0'
    for _J in J.Jd
        _J .= _J'
    end
    J
end

"""
$(SIGNATURES)

Evaluate ∑ᵢ exp(-λᵢτᵢ)xᵢ
"""
function expθ(J::JacobianDDE, x, λ::T) where T
    buffer = [one(T) * x]
    for τ in J.delays
        push!(buffer, copy(x) * exp(λ * (-τ)))
    end
    VectorOfArray(buffer)
end

"""
$(SIGNATURES)

Evaluate Δ(λ)⋅v where
    Δ(λ) = λI - J₀ - exp(-λτ)J₁
"""
function Δ(prob::AbstractDDEBifurcationProblem, x, p, v, λ)
    J = BK.jacobian(prob, x, p)
    Δ(J, v, λ)
end

function Δ(J::JacobianDDE, v, λ)
    res = λ .* v
    LA.mul!(res, J.J0, v, -1, 1)
    for (ind, A) in pairs(J.Jd)
        LA.mul!(res, A, v, -exp(-λ * J.delays[ind]), 1)
    end
    res
end

function A(J::JacobianDDE, v, λ)
    res = (0λ) .* v
    LA.mul!(res, J.J0, v, 1, 1)
    for (ind, A) in pairs(J.Jd)
        LA.mul!(res, A, v, exp(-λ * J.delays[ind]), 1)
    end
    res
end

"""
$(SIGNATURES)

Evaluate Δ'(λ)⋅v
"""
function Δ(::Val{:der}, J::JacobianDDE, v, λ)
    res = Complex.(v)
    for (ind, A) in pairs(J.Jd)
        LA.mul!(res, A, v, J.delays[ind] * exp(-λ * J.delays[ind]), 1)
    end
    res
end

function Δ(J::JacobianDDE, λ)
    n = size(J.Jall, 1)
    res = λ .* LA.I(n) .- J.J0
    for (ind, A) in pairs(J.Jd)
        res .+= (-exp(-λ * J.delays[ind])) .* A
    end
    res
end

function (l::BK.DefaultLS)(J::JacobianDDE, args...; kwargs...)
    l(J.Jall, args...; kwargs...)
end

function (l::BK.MatrixBLS)(iter::BK.AbstractContinuationIterable, state::BK.AbstractContinuationState, J::JacobianDDE, args...; kwargs...)
    l(iter, state, J.Jall, args...; kwargs...)
end

function (l::BK.MatrixBLS)(J::JacobianDDE, args...; kwargs...)
    l(J.Jall, args...; kwargs...)
end

function (l::BK.MatrixBLS)(J::JacobianDDE, dR,
                        dzu, dzp::T, R::AbstractVecOrMat, n::T,
                        ξu::T = one(T), ξp::T = one(T) ; kwargs...) where {T <: Number}
    l(J.Jall, dR, dzu, dzp, R, n, ξu, ξp ; kwargs...)
end

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
                dF = nothing,
                dFad = nothing,
                J = nothing,
                Jᵗ = nothing,
                J! = nothing,
                d2F = nothing,
                d3F = nothing,
                issymmetric::Bool = false,
                record_from_solution = BifurcationKit.record_sol_default,
                plot_solution = BifurcationKit.plot_default,
                save_solution = BifurcationKit.save_solution_default,
                inplace = false,
                δ = convert(eltype(u0), 1e-8),
                kwargs_jet...
                )
    @assert lens isa Int || lens isa BK.AllOpticTypes
    J = isnothing(J) ? (x,p) -> ForwardDiff.jacobian(z -> F(z, p), x) : J
    dF = isnothing(dF) ? (x,p,dx) -> ForwardDiff.derivative(t -> F(x .+ t .* dx, p), 0.) : dF
    d1Fad(x,p,dx1) = ForwardDiff.derivative(t -> F(x .+ t .* dx1, p), 0.)
    if isnothing(d2F)
        d2F = (x,p,dx1,dx2) -> ForwardDiff.derivative(t -> d1Fad(x .+ t .* dx2, p, dx1), 0.)
        d2Fc = (x,p,dx1,dx2) -> BilinearMap((_dx1, _dx2) -> d2F(x,p,_dx1,_dx2))(dx1,dx2)
    else
        d2Fc = d2F
    end
    if isnothing(d3F)
        d3F  = (x,p,dx1,dx2,dx3) -> ForwardDiff.derivative(t -> d2F(x .+ t .* dx3, p, dx1, dx2), 0.)
        d3Fc = (x,p,dx1,dx2,dx3) -> TrilinearMap((_dx1, _dx2, _dx3) -> d3F(x,p,_dx1,_dx2,_dx3))(dx1,dx2,dx3)
    else
        d3Fc = d3F
    end

    d3F = isnothing(d3F) ? (x,p,dx1,dx2,dx3) -> ForwardDiff.derivative(t -> d2F(x .+ t .* dx3, p, dx1, dx2), 0.) : d3F
    VF = BifFunction(F, nothing, dF, dFad, J, Jᵗ, nothing, d2F, d3F, d2Fc, d3Fc, issymmetric, δ, inplace, BK.Jet(;kwargs_jet...))
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
    # printstyled(io, isSymmetric(prob), color=:cyan, bold = true)
    print(io, "\n" * prefix * "└─ Parameter: ")
    printstyled(io, BK.get_lens_symbol(getlens(prob)), color=:cyan, bold = true)
end

function BK.residual(prob::SDDDEBifProblem, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    prob.VF.F(x,xd,p)
end

function BK.residual!(prob::SDDDEBifProblem, o, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    o .= prob.VF.F(x,xd,p)
    o
end

function BK.jacobian(prob::SDDDEBifProblem, x, p)
    xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
    J0 = ForwardDiff.jacobian(z -> prob.VF.F(z, xd, p), x)
    Jd = [ ForwardDiff.jacobian(z -> prob.VF.F(x, (@set xd.u[ii] = z), p), x) for ii in eachindex(prob.delays0)]
    return JacobianDDE(prob, J0 + sum(Jd), J0, Jd, prob.delays(x, p))
end

function BK.jacobian_adjoint(prob::SDDDEBifProblem, x, p)
    J = BK.jacobian(prob, x, p)
    J.Jall .= J.Jall'
    J.J0 .= J.J0'
    for _J in J.Jd
        _J .= _J'
    end
    J
end
