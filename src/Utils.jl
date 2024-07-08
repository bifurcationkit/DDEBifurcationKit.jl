# we need to add this dispatch as AbstractVectorOfArray is not a AbstractArray
(b::BK.TrilinearMap)(dx1::T, dx2::T, dx3::T) where {T <: RecursiveArrayTools.AbstractVectorOfArray{<: Real}} = b.tl(dx1, dx2, dx3)

(b::BK.BilinearMap)(dx1::T, dx2::T) where {T <: RecursiveArrayTools.AbstractVectorOfArray{<: Real}} = b.bl(dx1, dx2)