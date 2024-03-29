# Library

```@contents
Pages = ["library.md"]
Depth = 3
```

## Parameters

```@docs
BifurcationKit.NewtonPar
```

```@docs
BifurcationKit.ContinuationPar
```

## Problems

```@docs
ConstantDDEBifProblem
```

```@docs
SDDDEBifProblem
```

## Eigen solvers

```@docs
DDEBifurcationKit.DDE_DefaultEig
```

## Branch switching (branch point)

```@docs
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar ; kwargs...)
```

## Branch switching (Hopf point)
```@docs
continuation(br::BifurcationKit.AbstractBranchResult, ind_bif::Int, _contParams::ContinuationPar, prob::BifurcationKit.AbstractPeriodicOrbitProblem ; kwargs...)
```

## Utils for periodic orbits

```@docs
getperiod
```

```@docs
getamplitude
```

```@docs
getmaximum
```

## Misc.
```@docs
guess_from_hopf(br, ind_hopf, eigsolver::AbstractEigenSolver, M, amplitude; phase = 0)
```

```@docs
get_normal_form
```