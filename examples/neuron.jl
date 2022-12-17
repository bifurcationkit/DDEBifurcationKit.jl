using Revise, DDEBifurcationKit, Parameters, Setfield, RecursiveArrayTools
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

function neuronVF(x,xd,p)
   @unpack κ, β, a12, a21, τs, τ1, τ2 = p
   [
      -κ * x[1] + β * tanh(xd[3][1]) + a12 * tanh(xd[2][2]),
      -κ * x[2] + β * tanh(xd[3][2]) + a21 * tanh(xd[1][1])
   ]
end

function delaysF(dest, par)
   dest[1] = par.τ1
   dest[2] = par.τ2
   dest[3] = par.τs
   dest
end


pars = (κ = 0.5, β=-1, a12=1, a21=2.34, τ1=0.2, τ2=0.2, τs=1.5)
x0 = zeros(2)
delays0 = [0.2,0.2,1.5]
optn = NewtonPar(verbose = true, eigsolver = DDEBK.DDE_NLEVEigSolver(100))
opts = ContinuationPar(pMax = 3., pMin = 0., newtonOptions = optn, ds = -0.01, detectBifurcation = 3, nev = 3, )
br = BK.continuation(prob, PALC(), opts; verbosity = 1, plot = true)
