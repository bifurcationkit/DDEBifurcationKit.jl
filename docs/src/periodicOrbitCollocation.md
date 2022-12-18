# Periodic orbits based on orthogonal collocation

We compute `Ntst` time slices of a periodic orbit using orthogonal collocation. This is implemented in the structure `BifurcationKit.PeriodicOrbitOCollProblem`.

!!! warning "Large scale"
    The current implementation is not yet optimized for large scale problems. This will be improved in the future.    

The general method is explained in [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbitCollocation/).

## Phase condition

To ensure uniqueness of the solution to the functional, we add the following phase condition

$$\frac{1}{T} \int_{0}^{T}\left\langle x(s), \dot x_0(s)\right\rangle d s \approx  \sum_{j=1}^{N_{tst}}\sum_{i=1}^{m}\omega_i\left\langle x_{i,j}, \phi_{i,j}\right\rangle=0$$

> During continuation at step $k$, we use $\frac{1}{T} \int_{0}^{T}\left\langle x(s), \dot x_{k-1}(s)\right\rangle d s$


