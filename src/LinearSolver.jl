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

function BK.solve_bls_block(lbs::MatrixBLS,
                           J::JacobianDDE,
                           a::Tuple,
                           b::Tuple,
                           c::AbstractMatrix,
                           rhst,
                           rhsb)
    BK.solve_bls_block(lbs, J.Jall, a, b, c, rhst, rhsb)
end
