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

function ikedaVF(x, xd, p)
   @unpack Λ = p
   y = xd[1][1]
   [
      -pi/2 + Λ/2 * y^2;
   ]
end

function delaysF(par)
   [
      1.0,
   ]
end

pars = (Λ=0.1,b=0.)
x0 = [-sqrt(pi)]

prob = DDEBK.ConstantDDEBifProblem(ikedaVF, delaysF, x0, pars, (@lens _.Λ), recordFromSolution=(x,p)-> (x=x[1], _x=1))

optn = NewtonPar(verbose = false, eigsolver = DDEBK.DDE_NLEVEigSolver(maxit=100))
opts = ContinuationPar(pMax = 2., pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 4, nInversion = 12 )
br = BK.continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = false)
plot(br)

BK.getNormalForm(br, 1) # l1= -0.0591623057, b = 0.09293196762669392
