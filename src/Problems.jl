abstract type AbstractDDEBifurcationProblem <: BK.AbstractBifurcationProblem end

"""
$(TYPEDEF)

Structure to hold the bifurcation problem. If don't have parameters, you can pass `nothing`.

## Fields

$(TYPEDFIELDS)

## Methods

- `getu0(pb)` calls `pb.u0`
- `getParams(pb)` calls `pb.params`
- `getLens(pb)` calls `pb.lens`
- `getParam(pb)` calls `get(pb.params, pb.lens)`
- `setParam(pb, p0)` calls `set(pb.params, pb.lens, p0)`
- `recordFromSolution(pb)` calls `pb.recordFromSolution`
- `plotSolution(pb)` calls `pb.plotSolution`
- `isSymmetric(pb)` calls `isSymmetric(pb.prob)`

## Constructors

- `ConstantDDEBifProblem(F, delays, u0, params, lens; J, Jᵗ, d2F, d3F, kwargs...)` and `kwargs` are the fields above.

"""
struct ConstantDDEBifProblem{Tvf, Tdf, Tu, Td, Tp, Tl <: Lens, Tplot, Trec, Tδ} <: AbstractDDEBifurcationProblem
	"Vector field, typically a [`BifFunction`](@ref)"
	VF::Tvf
	"function delays. It takes the parameters and return the non-zero delays in an `AsbtractVector` form. Example: `delays = par -> [1.]`"
	delays::Tdf
	"Initial guess"
	u0::Tu
	"initial delays (set internally by the constructor)"
	delays0::Td
	"parameters"
	params::Tp
	"Typically a `Setfield.Lens`. It specifies which parameter axis among `params` is used for continuation. For example, if `par = (α = 1.0, β = 1)`, we can perform continuation w.r.t. `α` by using `lens = (@lens _.α)`. If you have an array `par = [ 1.0, 2.0]` and want to perform continuation w.r.t. the first variable, you can use `lens = (@lens _[1])`. For more information, we refer to `SetField.jl`."
	lens::Tl
	"user function to plot solutions during continuation. Signature: `plotSolution(x, p; kwargs...)`"
	plotSolution::Tplot
	"`recordFromSolution = (x, p) -> norm(x)` function used record a few indicators about the solution. It could be `norm` or `(x, p) -> x[1]`. This is also useful when saving several huge vectors is not possible for memory reasons (for example on GPU...). This function can return pretty much everything but you should keep it small. For example, you can do `(x, p) -> (x1 = x[1], x2 = x[2], nrm = norm(x))` or simply `(x, p) -> (sum(x), 1)`. This will be stored in `contres.branch` (see below). Finally, the first component is used to plot in the continuation curve."
	recordFromSolution::Trec
	"delta for Finite differences"
	δ::Tδ
end

BK.isInplace(::ConstantDDEBifProblem) = false
BK.isSymmetric(::ConstantDDEBifProblem) = false
BK.getVectorType(prob::ConstantDDEBifProblem{Tvf, Tdf, Tu, Td, Tp, Tl, Tplot, Trec}) where {Tvf, Tdf, Tu, Td, Tp, Tl <: Lens, Tplot, Trec} = Tu
BK.getLens(prob::ConstantDDEBifProblem) = prob.lens
BK.hasAdjoint(prob::ConstantDDEBifProblem) = true
BK.getDelta(prob::ConstantDDEBifProblem) = prob.δ
BK.d2F(prob::ConstantDDEBifProblem, x, p, dx1, dx2) = BK.d2F(prob.VF, x, p, dx1, dx2)
BK.d3F(prob::ConstantDDEBifProblem, x, p, dx1, dx2, dx3) = BK.d3F(prob.VF, x, p, dx1, dx2, dx3)

function Base.show(io::IO, prob::ConstantDDEBifProblem; prefix = "")
	print(io, prefix * "┌─ Constant Delays Bifurcation Problem with uType ")
	printstyled(io, BK.getVectorType(prob), color=:cyan, bold = true)
	print(io, prefix * "\n├─ Inplace:  ")
	printstyled(io, BK.isInplace(prob), color=:cyan, bold = true)
	print(io, "\n" * prefix * "└─ Parameter: ")
	printstyled(io, BK.getLensSymbol(getLens(prob)), color=:cyan, bold = true)
end

function ConstantDDEBifProblem(F, delayF, u0, parms, lens = (@lens _);
				dF = nothing,
				dFad = nothing,
				J = nothing,
				Jᵗ = nothing,
				d2F = nothing,
				d3F = nothing,
				issymmetric::Bool = false,
				recordFromSolution = BifurcationKit.recordSolDefault,
				plotSolution = BifurcationKit.plotDefault,
				inplace = false,
				δ = convert(eltype(u0), 1e-8)
				)

	Fc = (xd, p) -> F(xd[1], xd[2:end], p)
	# J = isnothing(J) ? (x,p) -> ForwardDiff.jacobian(z -> F(z, p), x) : J
	dF = isnothing(dF) ? (x,p,dx) -> ForwardDiff.derivative(t -> Fc(x .+ t .* dx, p), 0.) : dF
	d1Fad(x,p,dx1) = ForwardDiff.derivative(t -> Fc(x .+ t .* dx1, p), 0.)
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
	VF = BifFunction(F, dF, dFad, J, Jᵗ, d2F, d3F, d2Fc, d3Fc, issymmetric, 1e-8, inplace)
	return ConstantDDEBifProblem(VF, delayF, u0, delayF(parms), parms, lens, plotSolution, recordFromSolution, δ)
end

struct JacobianConstantDDE{Tp,T1,T2,T3,Td}
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

function jacobianFD(prob::ConstantDDEBifProblem, x, p)
	xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
	J0 = ForwardDiff.jacobian(z -> prob.VF.F(z, xd, p), x)

	Jd = [ ForwardDiff.jacobian(z -> prob.VF.F(x, (@set xd[ii] = z), p), x) for ii in eachindex(prob.delays0)]
	return J0, Jd
end

function BK.jacobian(prob::ConstantDDEBifProblem, x, p)
	if isnothing(prob.VF.J)
		J0, Jd = jacobianFD(prob, x, p)
	else
		J0, Jd = prob.VF.J(x, p)
	end
	return JacobianConstantDDE(prob, J0 + sum(Jd), J0, Jd, prob.delays(p))
end

function BK.jad(prob::ConstantDDEBifProblem, x, p)
	J = BK.jacobian(prob, x, p)
	J.Jall .= J.Jall'
	J.J0 .= J.J0'
	for _J in J.Jd
		_J .= _J'
	end
	J
end

function expθ(J::JacobianConstantDDE, x, λ::T) where T
	buffer = [one(T)*x]
	for τ in J.delays
		push!(buffer, copy(x) * exp(λ*(-τ)))
	end
	VectorOfArray(buffer)
end

function Δ(prob::ConstantDDEBifProblem, x, p, v, λ)
	J = BK.jacobian(prob, x, p)
	Δ(J, v, λ)
end

function Δ(J::JacobianConstantDDE, v, λ)
	res = λ .* v
	mul!(res, J.J0, v, -1, 1)
	for (ind, A) in pairs(J.Jd)
		mul!(res, A, v, -exp(-λ * J.delays[ind]), 1)
	end
	res
end

function Δ(::Val{:der}, J::JacobianConstantDDE, v, λ)
	res = Complex.(v)
	for (ind, A) in pairs(J.Jd)
		mul!(res, A, v, J.delays[ind] * exp(-λ * J.delays[ind]), 1)
	end
	res
end

function Δ(J::JacobianConstantDDE, λ)
	n = size(J.Jall, 1)
	res = λ .* I(n) .- J.J0
	for (ind, A) in pairs(J.Jd)
		res .+= (-exp(-λ * J.delays[ind])) .* A
	end
	res
end

function (l::BK.DefaultLS)(J::JacobianConstantDDE, args...; kwargs...)
	l(J.Jall, args...; kwargs...)
end

function (l::BK.MatrixBLS)(iter::BK.AbstractContinuationIterable, state::BK.AbstractContinuationState, J::JacobianConstantDDE, args...; kwargs...)
	l(iter, state, J.Jall, args...; kwargs...)
end

function (l::BK.MatrixBLS)(J::JacobianConstantDDE, args...; kwargs...)
	l(J.Jall, args...; kwargs...)
end

function (l::BK.MatrixBLS)(J::JacobianConstantDDE, dR,
						dzu, dzp::T, R::AbstractVecOrMat, n::T,
						ξu::T = T(1), ξp::T = T(1) ; kwargs...) where {T <: Number, Tξ}
	l(J.Jall, dR, dzu, dzp, R, n, ξu, ξp ; kwargs...)
end

"""
$(TYPEDEF)

Structure to hold the bifurcation problem. If don't have parameters, you can pass `nothing`.

## Fields

$(TYPEDFIELDS)

## Methods

- `getu0(pb)` calls `pb.u0`
- `getParams(pb)` calls `pb.params`
- `getLens(pb)` calls `pb.lens`
- `getParam(pb)` calls `get(pb.params, pb.lens)`
- `setParam(pb, p0)` calls `set(pb.params, pb.lens, p0)`
- `recordFromSolution(pb)` calls `pb.recordFromSolution`
- `plotSolution(pb)` calls `pb.plotSolution`
- `isSymmetric(pb)` calls `isSymmetric(pb.prob)`

## Constructors

- `SDDDEBifProblem(F, delays, u0, params, lens; J, Jᵗ, d2F, d3F, kwargs...)` and `kwargs` are the fields above.

"""
struct SDDDEBifProblem{Tvf, Tdf, Tu, Td, Tp, Tl <: Lens, Tplot, Trec, Tδ} <: AbstractDDEBifurcationProblem
	"Vector field, typically a [`BifFunction`](@ref)"
	VF::Tvf
	"function delays. It takes the parameters and the state and return the non-zero delays in an `AsbtractVector` form. Example: `delays = (par, u) -> [1. + u[1]^2]`"
	delays::Tdf
	"Initial guess"
	u0::Tu
	"initial delays (set internally by the constructor)"
	delays0::Td
	"parameters"
	params::Tp
	"Typically a `Setfield.Lens`. It specifies which parameter axis among `params` is used for continuation. For example, if `par = (α = 1.0, β = 1)`, we can perform continuation w.r.t. `α` by using `lens = (@lens _.α)`. If you have an array `par = [ 1.0, 2.0]` and want to perform continuation w.r.t. the first variable, you can use `lens = (@lens _[1])`. For more information, we refer to `SetField.jl`."
	lens::Tl
	"user function to plot solutions during continuation. Signature: `plotSolution(x, p; kwargs...)`"
	plotSolution::Tplot
	"`recordFromSolution = (x, p) -> norm(x)` function used record a few indicators about the solution. It could be `norm` or `(x, p) -> x[1]`. This is also useful when saving several huge vectors is not possible for memory reasons (for example on GPU...). This function can return pretty much everything but you should keep it small. For example, you can do `(x, p) -> (x1 = x[1], x2 = x[2], nrm = norm(x))` or simply `(x, p) -> (sum(x), 1)`. This will be stored in `contres.branch` (see below). Finally, the first component is used to plot in the continuation curve."
	recordFromSolution::Trec
	"delta for Finite differences"
	δ::Tδ
end

BK.isInplace(::SDDDEBifProblem) = false
BK.isSymmetric(::SDDDEBifProblem) = false
BK.getVectorType(prob::SDDDEBifProblem{Tvf, Tdf, Tu, Td, Tp, Tl, Tplot, Trec}) where {Tvf, Tdf, Tu, Td, Tp, Tl <: Lens, Tplot, Trec} = Tu
BK.getLens(prob::SDDDEBifProblem) = prob.lens
BK.hasAdjoint(prob::SDDDEBifProblem) = true

function SDDDEBifProblem(F, delayF, u0, parms, lens = (@lens _);
				dF = nothing,
				dFad = nothing,
				J = nothing,
				Jᵗ = nothing,
				d2F = nothing,
				d3F = nothing,
				issymmetric::Bool = false,
				recordFromSolution = BifurcationKit.recordSolDefault,
				plotSolution = BifurcationKit.plotDefault,
				inplace = false,
				δ = convert(eltype(u0), 1e-8)
				)
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
	VF = BifFunction(F, dF, dFad, J, Jᵗ, d2F, d3F, d2Fc, d3Fc, issymmetric, 1e-8, inplace)
	return SDDDEBifProblem(VF, delayF, u0, delayF(u0, parms), parms, lens, plotSolution, recordFromSolution, δ)
end

function Base.show(io::IO, prob::SDDDEBifProblem; prefix = "")
	print(io, prefix * "┌─ State-dependent delays Bifurcation Problem with uType ")
	printstyled(io, BK.getVectorType(prob), color=:cyan, bold = true)
	print(io, prefix * "\n├─ Inplace:  ")
	printstyled(io, BK.isInplace(prob), color=:cyan, bold = true)
	# printstyled(io, isSymmetric(prob), color=:cyan, bold = true)
	print(io, "\n" * prefix * "└─ Parameter: ")
	printstyled(io, BK.getLensSymbol(getLens(prob)), color=:cyan, bold = true)
end

function BK.residual(prob::SDDDEBifProblem, x, p)
	xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
	prob.VF.F(x,xd,p)
end

function BK.jacobian(prob::SDDDEBifProblem, x, p)
	xd = VectorOfArray([x for _ in eachindex(prob.delays0)])
	J0 = ForwardDiff.jacobian(z -> prob.VF.F(z, xd, p), x)

	Jd = [ ForwardDiff.jacobian(z -> prob.VF.F(x, (@set xd[ii] = z), p), x) for ii in eachindex(prob.delays0)]
	return JacobianConstantDDE(prob, J0 + sum(Jd), J0, Jd, prob.delays(x, p))
end
