# [Branch switching](@id Branch-switching-page)

The precise definition of the methods are given in [Branch switching (branch point)](@ref) and [Branch switching (Hopf point)](@ref).

```@contents
Pages = ["branchswitching.md"]
Depth = 3
```

!!! tip "BifurcationKit.jl"
    We only highlight what is specific to DDE. For more information go to [BifurcationKit](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/intro-abs/).


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
	
