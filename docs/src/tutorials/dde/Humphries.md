# Humphries model (codim 2)

```@contents
Pages = ["Humphries.md"]
Depth = 3
```
Consider the model [^Hum] as an example of state-dependent delays

$$x^{\prime}(t)=-\gamma x(t)-\kappa_1 x\left(t-a_1-c x(t)\right)-\kappa_2 x\left(t-a_2-c x(t)\right)$$


## Continuation and codim 1 bifurcations

We first instantiate the model

```@example TUTHumphries
using Revise, DDEBifurcationKit, Parameters, Setfield, Plots
using BifurcationKit
const BK = BifurcationKit

function humpriesVF(x, xd, p)
   @unpack κ1,κ2,γ,a1,a2,c = p
   [
      -γ * x[1] - κ1 * xd[1][1] - κ2 * xd[2][1]
   ]
end

function delaysF(x, par)
   [
      par.a1 + par.c * x[1],
      par.a2 + par.c * x[1],
   ]
end


pars = (κ1=0.,κ2=2.3,a1=1.3,a2=6,γ=4.75,c=1.)
x0 = zeros(1)

prob = SDDDEBifProblem(humpriesVF, delaysF, x0, pars, (@lens _.κ1))
```

We then compute the branch

```@example TUTHumphries
optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 13., pMin = 0., newtonOptions = optn, ds = -0.01, detectBifurcation = 3, nev = 3, )
br = continuation(prob, PALC(), opts; verbosity = 0, bothside = true)
```

and plot it

```@example TUTHumphries
scene = plot(br)
```

## Continuation of Hopf point

We follow the Hopf points in the parameter plane $(\kappa_1,\kappa_2)$.
We tell the solver to consider br.specialpoint[2] and continue it.

```@example TUTHumphries
brhopf = continuation(br, 2, (@lens _.κ2),
         setproperties(br.contparams, detectBifurcation = 2, dsmax = 0.04, maxSteps = 230, pMax = 5., pMin = -1.,ds = -0.02);
         verbosity = 0, plot = false,
         # we disable detection of Bautin bifurcation as the
         # Hopf normal form is not implemented for SD-DDE
         detectCodim2Bifurcation = 0,
         bothside = true,
         startWithEigen = true)

scene = plot(brhopf, vars = (:κ1, :κ2))
```

## References
[^Hum]: > Humphries et al. (2012), Dynamics of a delay differential equation with multiple state-dependent delays, Discrete and Continuous Dynamical Systems 32(8) pp. 2701-2727 http://dx.doi.org/10.3934/dcds.2012.32.2701)