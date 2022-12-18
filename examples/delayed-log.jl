cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."
using Revise, DDEBifurcationKit, Parameters, Setfield, RecursiveArrayTools, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

norminf(x) = norm(x, Inf)

function delayedlogVF(x, xd, p)
   @unpack λ = p
   y = xd[1][1]
   [
      (λ - y) * x[1]
   ]
end

function delaysF(par)
   [
      1.0,
   ]
end

pars = (λ=1.1,b=0.)
x0 = [1.]

prob = ConstantDDEBifProblem(delayedlogVF, delaysF, x0, pars, (@lens _.λ), recordFromSolution=(x,p)-> (x=x[1], _x=1))

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 2., pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 6, nInversion = 6 )
br = BK.continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = false)
plot(br)
