using DDEBifurcationKit
using Test

@testset "DDEBifurcationKit.jl" begin
    include("normalform.jl")
    include("codim2.jl")
    include("pocoll.jl")
    include("pocoll-sdde.jl")
end
