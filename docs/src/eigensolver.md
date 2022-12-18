# Eigen solvers (Eig)

See also [Eigen solvers](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/eigensolver/) for more information, for example to implement your own.



## List of implemented eigen solvers
1. Default [`DDE_DefaultEig`](@ref) eigensolver for DDE. You can create it via `eig = DDE_DefaultEig()`. It is based on the package [NonlinearEigenproblems.jl](https://github.com/nep-pack/NonlinearEigenproblems.jl).