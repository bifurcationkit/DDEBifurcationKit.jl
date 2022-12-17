module DDEBifurcationKit

	using BifurcationKit, Setfield, DocStringExtensions, RecursiveArrayTools
	using ForwardDiff, Parameters, LinearAlgebra
	using NonlinearEigenproblems
	const BK = BifurcationKit

	abstract type AbstractDDEBifurcationProblem <: BK.AbstractBifurcationProblem end

	include("problems.jl")
	include("diffeq.jl")
	include("EigSolver.jl")
	include("codim2/codim2.jl")
	include("codim2/Hopf.jl")



end
