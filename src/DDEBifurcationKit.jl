module DDEBifurcationKit
    using BifurcationKit, DocStringExtensions, RecursiveArrayTools
    using ForwardDiff, Parameters
    import LinearAlgebra as LA
    using NonlinearEigenproblems
    const BK = BifurcationKit


    include("Problems.jl")
    include("Utils.jl")
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
