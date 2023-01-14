# Fold / Hopf Continuation

In this page, we explain how to perform continuation of Fold / Hopf points and detect the associated bifurcations.

### List of detected codim 2 bifurcation points
|Bifurcation|symbol used|
|---|---|
| Bogdanov-Takens | bt |
| Bautin | gh |
| Cusp | cusp |
| Zero-Hopf | zh |
| Hopf-Hopf | hh |

In a nutshell, all you have to do (see below) is to call `continuation(br, ind_bif, lens2)` to continue the bifurcation point stored in `br.specialpoint[ind_bif]` and set proper options.

## Fold continuation

The continuation of Fold bifurcation points is based on a **Minimally Augmented**[^Govaerts] formulation which is an efficient way to detect singularities. See [docs in BifurcationKit](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/codim2Continuation/) for more information.

### Detection of codim 2 bifurcation points

You can detect the following codim 2 bifurcation points by using the option `detectCodim2Bifurcation` in the method `continuation`.

- the detection of Cusp (Cusp) is done by the detection of Fold bifurcation points along the curve of Folds by monitoring the parameter component of the tangent.
- the detection of Bogdanov-Takens (BT) is performed using the test function[^Bindel] $\psi_{BT}(p) = \langle w(p),v(p)\rangle$
- the detection of Zero-Hopf (ZH) is performed by monitoring the number of eigenvalues $\lambda$ such that $\Re\lambda > \min\limits_{\nu\in\Sigma(dF)}|\Re\nu|$ and $\Im\lambda > \epsilon$ where $\epsilon$ is the Newton tolerance.

## Hopf continuation

The continuation of Fold bifurcation points is based on solving the extended system for $(u^*, v, \omega)$

$$\begin{aligned}
& 0=\mathbf F\left(u^*, u^*; p\right) \\
& 0=\Delta\left(u^*, p, \mathrm{i} \omega\right) v \\
& 0=v^{\mathrm{H}} v-1
\end{aligned}$$

where $\Delta(\lambda)\cdot v := \lambda v - d_1\mathbf F(u^*,u^*; p)\cdot v-d_2\mathbf F(u^*, u^*; p)\cdot(e^{\lambda\cdot}v)$.

### Detection of codim 2 bifurcation points

You can detect the following codim 2 bifurcation points by using the option `detectCodim2Bifurcation` in the method `continuation`.

- the detection of Bogdanov-Takens (BT) is performed using the test function $\psi_{BT}(p) = 	\omega$
- the detection of Bautin (GH) is based on the test function $\psi_{GH}(p) = \Re(l_1(p))$ where $l_1$ is the Lyapunov coefficient defined in [Simple Hopf point](@ref).
- the detection of Zero-Hopf (ZH) is performed by monitoring the eigenvalues.
- the detection of Hopf-Hopf (HH) is performed by monitoring the eigenvalues.

> The continuation of Hopf points is stopped at BT and when $\omega<100\epsilon$ where $\epsilon$ is the newton tolerance.


## Codim 2 continuation

To compute the codim 2 curve of Fold / Hopf points, one can call `continuation` with the following options

```@docs
 continuation(br::BifurcationKit.AbstractBranchResult, ind_bif::Int64,
				lens2::Lens, options_cont::ContinuationPar = br.contparams ;
				kwargs...)
```

where the options are as above except with have an additional parameter axis `lens2` which is used to locate the bifurcation points.
