module DDEBifurcationKit
    using BifurcationKit, DocStringExtensions, RecursiveArrayTools
    using ForwardDiff, Parameters
    import LinearAlgebra as LA
    import NonlinearEigenproblems as NLE
    const BK = BifurcationKit


    include("Problems.jl")
    include("Utils.jl")
    include("LinearSolver.jl")
    include("NormalForms.jl")
    include("EigSolver.jl")
    include("codim2/codim2.jl")
    include("codim2/Hopf.jl")
    include("codim2/Fold.jl")

    include("periodicorbit/PeriodicOrbits.jl")
    include("periodicorbit/PeriodicOrbitCollocation.jl")

    export ConstantDDEBifProblem, SDDDEBifProblem
    export DDE_DefaultEig

    # re-export methods, macros from BK
    export continuation, @optic, NewtonPar, ContinuationPar, PALC, MoorePenrose, Natural, Polynomial, Secant, Bordered
    export PeriodicOrbitTrapProblem, PeriodicOrbitOCollProblem, ShootingProblem, FloquetCollGEV
    export @reset, norminf, get_periodic_orbit, getperiod, get_normal_form, setparam
end
