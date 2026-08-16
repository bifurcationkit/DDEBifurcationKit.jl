# [Branch switching](@id Branch-switching-page)

The precise definition of the methods are given in [Branch switching (branch point)](@ref) and [Branch switching (Hopf point)](@ref).

```@contents
Pages = ["branchswitching.md"]
Depth = 3
```

## Branch switching from simple branch point to equilibria

You can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar; kwargs...)
```

where `br` is a branch computed after a call to `continuation` with detection of bifurcation points enabled. This call computes the branch bifurcating from the `ind_bif `th bifurcation point in `br`. An example of use is provided in [2d generalized Bratu–Gelfand problem](@ref).

> See [Branch switching (branch point)](@ref) precise method definition

## Branch switching from non simple branch point to equilibria

We provide an automatic branch switching method in this case. The method is to first compute the reduced equation (see [Non-simple branch point](@ref)) and use it to compute the nearby solutions. These solutions are seeded as initial guess for [`continuation`](@ref). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar; kwargs...)
```

An example of use is provided in [2d generalized Bratu–Gelfand problem](@ref).

> See [Branch switching (branch point)](@ref) for the precise method definition

## Branch switching from Hopf point to periodic orbits

In order to compute the bifurcated branch of periodic solutions at a Hopf bifurcation point, you need to choose a method to compute periodic orbits among:

- [Periodic orbits based on orthogonal collocation](@ref)

Once you have decided which method to use, you use the following call:

```julia
continuation(br::ContResult, ind_HOPF::Int, _contParams::ContinuationPar,
	prob::AbstractPeriodicOrbitProblem ;
	δp = nothing, ampfactor = 1, kwargs...)
```

We refer to [`continuation`](@ref) for more information about the arguments. Here, we just say a few words about how we can specify `prob::AbstractPeriodicOrbitProblem`.

- For [Periodic orbits based on orthogonal collocation](@ref), you can pass `PeriodicOrbitOCollProblem(M, m)` where `M` is the number of times slices in the periodic orbit and `m` is the degree of the collocation polynomials.

## From Hopf-Hopf (HH) point to Fold / Hopf curve

We provide an automatic branch switching method in this case (see for example [Neuron model](@ref lorenz)). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_HH::Int,
	options_cont::ContinuationPar = br.contparams;
	detailed = Val(false), # MUST BE PASSED
	alg = getalg(br),
	δp = nothing, ampfactor::Real = 1,
	nev = options_cont.nev,
	detect_codim2_bifurcation::Int = 0,
	start_with_eigen = false,
	scaleζ = norm,
	kwargs...)
```

where `ind_HH` is the index of the HH point in `br`. Note that the HH has been detected during Hopf continuation. Calling the above method thus switches from Hopf continuation to another Hopf branch automatically with the same parameter axis.

> Check the docs of [Fold / Hopf Continuation](@ref) and particularly [Setting the jacobian](@ref jac-fold) for improving the speed of computation for large scale systems.
	
