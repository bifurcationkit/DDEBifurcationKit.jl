cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."
using Revise, DDEBifurcationKit, Parameters, Setfield, RecursiveArrayTools
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

g(z) = (tanh(z − 1) + tanh(1))*cosh(1)^2
function neuron2VF(x, xd, p)
   @unpack a,b,c,d = p
   [
      -x[1] - a * g(b*xd[1][1]) + c * g(d*xd[2][2]),
      -x[2] - a * g(b*xd[1][2]) + c * g(d*xd[2][1])
   ]
end

function delaysF(par)
   [par.τ1, par.τ2]
end

pars = (a = 0.25, b = 2., c = 15/29, d = 1.2, τ1 = 12.7, τ2 = 20.2)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuron2VF, delaysF, x0, pars, (@lens _.a))

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig(maxit=100))
opts = ContinuationPar(pMax = 1., pMin = 0., newtonOptions = optn, ds = 0.01, detectBifurcation = 3, nev = 9, dsmax = 0.2, nInversion = 4)
br = continuation(prob, PALC(tangent=Bordered()), opts; verbosity = 0, plot = true, bothside = false)

plot(br)

hpnf = BK.getNormalForm(br, 2)
################################################################################
brhopf = continuation(br, 1, (@lens _.c),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.01, maxSteps = 100, pMax = 1.1, pMin = -0.1,ds = 0.01, nInversion = 2);
         verbosity = 0, plot = true,
         detectCodim2Bifurcation = 2,
         bothside = true,
         startWithEigen = true)

brhopf2 = continuation(br, 2, (@lens _.c),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.01, maxSteps = 100, pMax = 1.1, pMin = -0.1,ds = -0.01);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 2,
         bothside = true,
         startWithEigen = true)

plot(brhopf, vars = (:a, :c), xlims = (0,0.7), ylims = (0,1))
   plot!(brhopf2, vars = (:a, :c), xlims = (-0,0.7), ylims = (-0.1,1))

################################################################################
prob2 = ConstantDDEBifProblem(neuron2VF, delaysF, x0, (@set pars.a = 0.12), (@lens _.c))

br2 = continuation(prob2, PALC(), setproperties(opts, pMax = 1.22); verbosity = 1, plot = true, bothside = false)

plot(br2)
