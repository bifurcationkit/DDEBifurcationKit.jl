# DDEBifurcationKit.jl

This Julia package aims at performing **automatic bifurcation analysis** of possibly large dimensional equations delay differential equations (DDE) by taking advantage of iterative methods, dense / sparse formulation and specific hardwares (*e.g.* GPU).

It builds upon [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl) with version > 0.2 to perform continuation and numerical bifurcation analysis.

## 📦 Installation

Assuming that you already have Julia correctly installed, it suffices to import DDEBifurcationKit.jl in the standard way:

`import Pkg; Pkg.add("https://github.com/bifurcationkit/DDEBifurcationKit.jl")`

## 📚 Citing this work
If you use this package for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib).

## 🧑‍💻 Other softwares

There are several good softwares already available.

- For continuation in small dimension, most softwares are listed on [DSWeb](https://ddebiftool.sourceforge.net). One can mention the widely used [DDE-BIFTOOL](http://www.math.pitt.edu/~bard/xpp/xpp.html), [Knut](https://rs1909.github.io/knut/). All these are very reliable and some address high codimension bifurcations.

- For large scale problems, there is very little.

In Julia, the present package seems to be the only one.

## A word on performance

The examples which follow have not **all** been written with the goal of performance but rather simplicity (for now).

## Main features

- Newton-Krylov solver with generic linear / eigen *preconditioned* solver. Idem for the arc-length continuation.
- Newton-Krylov solver with nonlinear deflation and preconditioner. It can be used for branch switching for example.
- Continuation written as an [iterator](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/iterator/)
- Monitoring user functions along curves computed by continuation, see [events](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/EventCallback/).
- Continuation methods: PALC, Moore-Penrose, etc. See [methods](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/IntroContinuation/).
- Bifurcation points are located using a bisection algorithm
- detection of Branch, Fold, Hopf bifurcation point of stationary solutions and computation of their normal form.
- <s>Automatic branch switching at branch points (whatever the dimension of the kernel)</s>
- Automatic branch switching at simple Hopf points to periodic orbits
- **Automatic bifurcation diagram computation of equilibria**
- Fold / Hopf continuation.
- detection all codim 2 bifurcations of equilibria and <s>computation of the normal forms of Bogdanov-Takens, Bautin and Cusp</s>
- <s>Branching from Bogdanov-Takens points to Fold / Hopf curve</s>
- Periodic orbit computation and continuation using <s>Shooting, Finite Differences or </s>Orthogonal Collocation.
- <s>detection of Branch, Fold, Neimark-Sacker, Period Doubling bifurcation point of periodic orbits.</s>
- <s>Continuation of Fold of periodic orbits</s>

Custom state means, we can use something else than `AbstractArray`, for example your own `struct`.

Type of delay: Constant (C), state-dependent (SD), nested (N)

|Features| delay type | Matrix Free|Custom state| [Tutorial](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/) | GPU |
|---|---|---|---|---|---|
| (Deflated) Krylov-Newton| C/SD |  Yes | Yes| | |
| Continuation PALC (Natural, Secant, Tangent, Polynomial) | C/SD| | | | |
| Bifurcation / Fold / Hopf point detection | C/SD | Y|   |  | |
| Fold Point continuation |C/SD | Y |  |  |
| Hopf Point continuation | C/SD |  | `AbstractArray` | | |
| ~~Bogdanov-Takens Point newton~~ | C/SD | Y | `AbstractArray` | | |
| Branch point / Fold / Hopf normal form | C/SD | Y|  | |  | |
| Branch switching at Branch / Hopf points | C/SD | Y | `AbstractArray` |  |  
| <span style="color:red">**Automatic bifurcation diagram computation of equilibria**</span> | C/SD| Y| `AbstractArray` |  | |
| ~~Periodic Orbit (Trapezoid) Newton / continuation~~ | | | `AbstractVector` |  | |
| Periodic Orbit (Collocation) Newton / continuation | C/SD |  | `AbstractVector` |  | |
| ~~Periodic Orbit (Parallel Poincaré / Standard Shooting) Newton / continuation~~ | | | `AbstractArray` |   | |
| ~~Fold, Neimark-Sacker, Period doubling detection~~ | | | `AbstractVector` |   | |
| ~~Continuation of Fold of periodic orbits~~ | | | `AbstractVector` |  |  |
| Bogdanov-Takens / Bautin / Cusp / Zero-Hopf / Hopf-Hopf point detection | C/SD| Y|  |  |
|~~Bogdanov-Takens / Bautin / Cusp / Zero-Hopf / Hopf-Hopf normal forms~~ | | Y|  |  |
| ~~Branching from Bogdanov-Takens points to Fold / Hopf curve~~ |  | |  `AbstractVector` | |  |

## Requested methods for Custom State
Needless to say, if you use regular arrays, you don't need to worry about what follows.

We make the same requirements as `KrylovKit.jl`. Hence, we refer to its [docs](https://jutho.github.io/KrylovKit.jl/stable/#Package-features-and-alternatives-1) for more information. We additionally require the following methods to be available:

- `Base.length(x)`: it is used in the constraint equation of the pseudo arclength continuation method (see [`continuation`](@ref) for more details). If `length` is not available for your "vector", define it `length(x) = 1` and adjust tuning the parameter `theta` in `ContinuationPar`.
- `Base.copyto!(dest, in)` this is required to reduce the allocations by avoiding too many copies
- `Base.eltype` must be extended to your vector type. It is mainly used for branching.

## Citations
The papers citing this work are collected on [google scholar](https://scholar.google.fr/scholar?hl=fr&as_sdt=2005&cites=159498619004863176%2C8662907770106865595&scipsc=&as_ylo=&as_yhi=).