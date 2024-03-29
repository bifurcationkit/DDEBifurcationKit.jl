# Simple Hopf point


At a Hopf branch point $(x_0,p_0)$ for the problem $F(x,p)=0$, the spectrum of the linear operator $dF(x_0,p_0)$ contains two purely imaginary $\pm i\omega,\ \omega > 0$ which are simple. At such point, we can compute the **normal form** to transform the DDE problem

$$\dot x = \mathbf{F}(x_t,p)$$

in large dimensions to a **complex** polynomial vector field ($\delta p\equiv p-p_0$):

$$\dot z = z\left(a \cdot\delta p + i\omega + l_1|z|^2\right)\quad\text{(E)}$$

whose solutions give access to the solutions of the Cauchy problem in a neighborhood of $(x,p)$.

!!! tip "Coefficient $l_1$"
    The coefficient $l_1$ above is called the **Lyapunov** coefficient

## Normal form computation

The normal form (E) is automatically computed as follows

```julia
get_normal_form(br::ContResult, ind_bif::Int ;
	verbose = false, ζs = nothing, lens = br.param_lens)
```

where `prob` is a bifurcation problem. `br` is a branch computed after a call to `continuation` with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters, we refer to [`get_normal_form`](@ref). The above call returns a point with information needed to compute the bifurcated branch.

```julia
mutable struct Hopf{Tv, T, Tω, Tevr, Tevl, Tnf} <: BifurcationPoint
	"Hopf point"
	x0::Tv

	"Parameter value at the Hopf point"
	p::T

	"Frequency of the Hopf point"
	ω::Tω

	"Right eigenvector"
	ζ::Tevr

	"Left eigenvector"
	ζstar::Tevl

	"Normal form coefficient (a = 0., b = 1 + 1im)"
	nf::Tnf

	"Type of Hopf bifurcation"
	type::Symbol
end
```

!!! info "Note"
    You should not need to call `get_normal_form ` except if you need the full information about the branch point.

## References

[^Haragus]: > Haragus, Mariana, and Gérard Iooss. Local Bifurcations, Center Manifolds, and Normal Forms in Infinite-Dimensional Dynamical Systems. London: Springer London, 2011. https://doi.org/10.1007/978-0-85729-112-7.
