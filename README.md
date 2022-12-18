# DDEBifurcationKit

`DDEBifurcationKit.jl` is a component package in the `BifurcationKit` ecosystem. It holds the delay differential equation (DDE) utilities. While completely independent
and usable on its own, users interested in using this
functionality should check out [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl).

## Installation

Assuming that you already have Julia correctly installed, it suffices to import
`DDEBifurcationKit.jl` in the standard way:

```julia
import Pkg; Pkg.add("DDEBifurcationKit")
```

## Support and citation
If you use this package for your work, we ask that you cite the following paper. Open source development as part of academic research strongly depends on this. Please also consider starring this repository if you like our work, this will help us to secure funding in the future. It is referenced on HAL-Inria as follows:

```
@misc{veltz:hal-02902346,
  TITLE = {{BifurcationKit.jl}},
  AUTHOR = {Veltz, Romain},
  URL = {https://hal.archives-ouvertes.fr/hal-02902346},
  INSTITUTION = {{Inria Sophia-Antipolis}},
  YEAR = {2020},
  MONTH = Jul,
  KEYWORDS = {pseudo-arclength-continuation ; periodic-orbits ; floquet ; gpu ; bifurcation-diagram ; deflation ; newton-krylov},
  PDF = {https://hal.archives-ouvertes.fr/hal-02902346/file/354c9fb0d148262405609eed2cb7927818706f1f.tar.gz},
  HAL_ID = {hal-02902346},
  HAL_VERSION = {v1},
}
```