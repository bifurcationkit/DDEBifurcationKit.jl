module DDEBifurcationKit

	using BifurcationKit, Setfield, DocStringExtensions, RecursiveArrayTools
	using ForwardDiff, Parameters, LinearAlgebra
	using NonlinearEigenproblems
	const BK = BifurcationKit

	abstract type AbstractDDEBifurcationProblem <: BK.AbstractBifurcationProblem end

	include("Problems.jl")
	include("NormalForms.jl")
	include("EigSolver.jl")
	include("codim2/codim2.jl")
	include("codim2/Hopf.jl")
	include("codim2/Fold.jl")

	include("periodicorbit/PeriodicOrbits.jl")
	include("periodicorbit/PeriodicOrbitCollocation.jl")

	export ConstantDDEBifProblem, SDDDEBifProblem
	export DDE_DefaultEig



end
