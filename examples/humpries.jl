cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."

# https://ddebiftool.sourceforge.net/demos/neuron/html/demo1_stst.html
using Revise, DDEBifurcationKit, Parameters, Setfield, RecursiveArrayTools
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

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

prob = DDEBK.SDDDEBifProblem(humpriesVF, delaysF, x0, pars, (@lens _.κ1))

optn = NewtonPar(verbose = true, eigsolver = DDEBK.DDE_NLEVEigSolver(maxit=100))
opts = ContinuationPar(pMax = 13., pMin = 0., newtonOptions = optn, ds = -0.01, detectBifurcation = 3, nev = 3, )
br = BK.continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = true)
plot(br)
