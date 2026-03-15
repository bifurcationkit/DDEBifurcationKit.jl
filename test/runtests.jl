using DDEBifurcationKit
using Test

@testset "DDEBifurcationKit.jl" begin
    @testset "Normal forms" begin
        include("normalform.jl")
    end
    @testset "Codim 2" begin
        include("codim2.jl")
    end
    @testset "Collocation DDE" begin
        include("pocoll.jl")
    end
    @testset "Collocation SD-DDE" begin
        include("pocoll-sdde.jl")
    end
    @testset "Floquet" begin
        include("testfloquet.jl")
    end
end
