# Periodic orbits based on orthogonal collocation

We compute `Ntst` time slices of a periodic orbit using orthogonal collocation. This is implemented in the structure `BifurcationKit.PeriodicOrbitOCollProblem`.

!!! warning "Large scale"
    The current implementation is not yet optimized for large scale problems. This will be improved in the future.    

The general method is explained in [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbitCollocation/).
