# DDEBifurcationKit.jl

This Julia package aims at performing **automatic bifurcation analysis** of possibly large dimensional equations delay differential equations (DDE) by taking advantage of iterative methods, dense / sparse formulation and specific hardwares (*e.g.* GPU).

It builds upon [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl) with version > 0.2 to perform continuation and numerical bifurcation analysis.

## 📦 Installation

Assuming that you already have Julia correctly installed, it suffices to add `DDEBifurcationKit.jl` in the standard way:

`] add DDEBifurcationKit.jl`

## 📚 Citing this work
If you use this package for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib).

## 🧑‍💻 Other softwares

There are several good softwares already available.

- For continuation in small dimension, most softwares are listed on [DSWeb](https://dsweb.siam.org/Software/PID/1956/SearchID/1964/cfs/True/cblcid_96_84/84). One can mention the widely used [DDE-BIFTOOL](https://github.com/DDE-BifTool/DDE-Biftool), [Knut](https://rs1909.github.io/knut/). All these are very reliable and some address high codimension bifurcations.

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

### Capabilities related to equilibria
- detection of Branch, Fold, Hopf bifurcation point of stationary solutions and computation of their normal form.
```@raw html
<ul> 
<li><del>Automatic branch switching at branch points (whatever the dimension of the kernel)</del></li></ul>
```
- **Automatic bifurcation diagram computation of equilibria**
- Fold / Hopf continuation.
```@raw html
<ul> 
<li>detection all codim 2 bifurcations of equilibria and <del>computation of the normal forms of Bogdanov-Takens, Bautin and Cusp</del></li>
<li><del>Branching from Bogdanov-Takens points to Fold / Hopf curve</del></li>
<li>Periodic orbit computation and continuation using <del>Shooting, Finite Differences or </del>Orthogonal Collocation.</li>
<li><del>detection of Branch, Fold, Neimark-Sacker, Period Doubling bifurcation point of periodic orbits.</del></li>
<li><del>Continuation of Fold of periodic orbits</del></li>
</ul>
```

Custom state means, we can use something else than `AbstractArray`, for example your own `struct`.

Type of delay: Constant (C), state-dependent (SD), nested (N)

|Features| delay type | Matrix Free|Custom state| [Tutorials](https://bifurcationkit.github.io/DDEBifurcationKit.jl/dev/tutorials/tutorials/) | GPU |
|---|---|---|---|---|---|
| (Deflated) Krylov-Newton| C/SD |  Yes | Yes| | |
| Continuation PALC (Natural, Secant, Tangent, Polynomial) | C/SD| | | | |
| Bifurcation / Fold / Hopf point detection | C/SD | Y|   |  | |
| Fold Point continuation |C/SD | Y |  |  |
| Hopf Point continuation | C/SD |  | `AbstractArray` | | |
| Branch point / Fold form | C/SD | |  | |  | |
| Hopf normal form | C | |  | |  | |
| Branch switching at Branch / Hopf points | C/SD |  | `AbstractArray` |  |  |
| Automatic bifurcation diagram computation of equilibria | C/SD| Y| `AbstractArray` |  | |
| Periodic Orbit (Collocation) Newton / continuation | C/SD |  | `AbstractVector` |  | |
| Bogdanov-Takens / Bautin / Cusp / Zero-Hopf / Hopf-Hopf point detection | NA | |  |  |

### Capabilities related to Periodic orbits (PO)

- PO computation and continuation Orthogonal Collocation (mesh adaptive).
- Computation of Floquet exponents
- Automatic branch switching at simple Hopf points (DDE, SD-DDE) to periodic orbits
- Detection of Branch, Fold, Neimark-Sacker (NS), Period Doubling (PD) bifurcation points of PO.
- Assisted branch switching from simple Period-Doubling points to PO
- Assisted branch switching from simple Branch points to PO

|Features| delay type | Matrix Free|Custom state| [Tutorials](https://bifurcationkit.github.io/DDEBifurcationKit.jl/dev/tutorials/tutorials/) | GPU |
|---|---|---|---|---|---|
| Periodic Orbit (Collocation) Newton / continuation | C/SD |  | `AbstractVector` |  | |
| Floquet exponents |C| | `AbstractVector` |  | |
| Fold, Neimark-Sacker, Period doubling detection |C| | `AbstractVector` |  | |


## Requested methods for Custom State

If you use standard arrays, you can skip this section.

We make the same requirements as `KrylovKit.jl`. Hence, we refer to its [docs](https://jutho.github.io/KrylovKit.jl/stable/#Package-features-and-alternatives-1) and the package [VectorInterface.jl](https://github.com/Jutho/VectorInterface.jl?tab=readme-ov-file) for more information. We additionally require the following methods to be available:

- `Base.length(x)`: it is used in the constraint equation of the pseudo arclength continuation method (see [`continuation`](@ref) for more details). If `length` is not available for your "vector" type, define `length(x) = 1` and adjust the parameter `θ` in `PALC`.
- `Base.copyto!(dest, in)` this is used to reduce the allocations.
- `real` this is used to compute normal forms.
- `conj` this is used to compute normal forms.

## Citations
Papers citing this work are collected on [Zotero](https://www.zotero.org/groups/6097154/citations_of_bifurcationkit/library).

These citations are aggregated from [Google Scholar (search)](https://scholar.google.com/scholar?q=bifurcationkit&hl=en&as_sdt=0,5) and [Google Scholar (citations)](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=159498619004863176,12573642401780006854,8662907770106865595). Note that each link may reference different subsets of papers.