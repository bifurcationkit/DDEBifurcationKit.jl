# DDEBifurcationKit

| **Documentation** | **Build Status** | **Downloads** |
|:-----------------:|:----------------:|:-------------:|
| [![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] |   |  |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev


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

## Main features

Type of delay: Constant (D), state-dependent (S), nested (N)

|Features| delay type | Matrix Free|Custom state| [Tutorial](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/) | GPU |
|---|---|---|---|---|---|
| (Deflated) Krylov-Newton| C/SD |  Yes | Yes| | |
| Continuation PALC (Natural, Secant, Tangent, Polynomial) | C/SD| | | | |
| Bifurcation / Fold / Hopf point detection | C/SD | Yes|   |  | |
| Fold Point continuation |C/SD | Yes|  |  |
| Hopf Point continuation | C/SD | Y | `AbstractArray` | | |
| ~~Bogdanov-Takens Point newton~~ | C/SD | Y | `AbstractArray` | | |
| Branch point / Fold / Hopf normal form | C/SD | Yes|  | |  | |
| Branch switching at Branch / Hopf points | C/SD | Y | `AbstractArray` |  |  
| <span style="color:red">**Automatic bifurcation diagram computation of equilibria**</span> | C/SD| Yes| `AbstractArray` |  | |
| ~~Periodic Orbit (Trapezoid) Newton / continuation~~ | | | `AbstractVector` |  | |
| Periodic Orbit (Collocation) Newton / continuation | C |  | `AbstractVector` |  | |
| ~~Periodic Orbit (Parallel Poincaré / Standard Shooting) Newton / continuation~~ | | | `AbstractArray` |   | |
| ~~Fold, Neimark-Sacker, Period doubling detection~~ | | | `AbstractVector` |   | |
| ~~Continuation of Fold of periodic orbits~~ | | | `AbstractVector` |  |  |
| Bogdanov-Takens / Bautin / Cusp / Zero-Hopf / Hopf-Hopf point detection | C/SD| Yes|  |  |
|~~Bogdanov-Takens / Bautin / Cusp / Zero-Hopf / Hopf-Hopf normal forms~~ | | Yes|  |  |
| ~~Branching from Bogdanov-Takens points to Fold / Hopf curve~~ |  | |  `AbstractVector` | |  |
