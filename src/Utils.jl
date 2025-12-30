struct JacobianDDE{Tp,T1,T2,T3,Td}
    prob::Tp
    Jall::T1
    J0::T2
    Jd::T3
    delays::Td
end

# for matrix-free operators, we do not sum the arrays
JacobianDDE(prob, J0, Jd, delays) = JacobianDDE(prob, nothing, J0, Jd, delays)
JacobianDDE(prob, J0::AbstractArray, Jd::AbstractVector{ <: AbstractArray}, delays) = JacobianDDE(prob, J0 + sum(Jd), J0, Jd, delays)

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

function Δ(J::JacobianDDE, λ)
    n = size(J.Jall, 1)
    res = λ .* LA.I(n) .- J.J0
    for (ind, A) in pairs(J.Jd)
        res .+= (-exp(-λ * J.delays[ind])) .* A
    end
    res
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
####################################################################################################
# we need to add this dispatch as AbstractVectorOfArray is not a AbstractArray
(b::BK.TrilinearMap)(dx1::T, dx2::T, dx3::T) where {T <: RecursiveArrayTools.AbstractVectorOfArray{<: Real}} = b.tl(dx1, dx2, dx3)

(b::BK.BilinearMap)(dx1::T, dx2::T) where {T <: RecursiveArrayTools.AbstractVectorOfArray{<: Real}} = b.bl(dx1, dx2)
