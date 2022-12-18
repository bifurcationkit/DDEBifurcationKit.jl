using Revise, DDEBifurcationKit, Parameters, Setfield, RecursiveArrayTools
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

using Plots

function neuronVF(x, xd, p)
   @unpack κ, β, a12, a21, τs, τ1, τ2 = p
   [
      -κ * x[1] + β * tanh(xd[3][1]) + a12 * tanh(xd[2][2]),
      -κ * x[2] + β * tanh(xd[3][2]) + a21 * tanh(xd[1][1])
   ]
end

delaysF(par) = [par.τ1, par.τ2, par.τs]

pars = (κ = 0.5, β = -1, a12 = 1, a21 = 0.5, τ1 = 0.2, τ2 = 0.2, τs = 1.5)
x0 = [0.01, 0.001]

prob = ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@lens _.τs))

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(pMax = 13., pMin = 0., newtonOptions = optn, ds = -0.01, detectBifurcation = 3, nev = 5, dsmax = 0.2, nInversion = 4)
br = BK.continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = true, normC = norminf)

plot(br)
################################################################################
prob2 = DDEBK.ConstantDDEBifProblem(neuronVF, delaysF, x0, pars, (@lens _.a21))
br2 = BK.continuation(prob2, PALC(), opts; verbosity = 1, plot = true, bothside = true)
BK.getNormalForm(br2, 3)
#L1 ≈ −0.0601.
################################################################################
brhopf = continuation(br, 3, (@lens _.a21),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.04, maxSteps = 230, pMax = 15., pMin = -1.,ds = -0.02);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 2,
         # bothside = true,
         startWithEigen = true)

plot(brhopf, vars = (:a21, :τs))
plot(brhopf, vars = (:τs, :ω))

brhopf2 = continuation(br, 2, (@lens _.a21),
         setproperties(br.contparams, detectBifurcation = 1, dsmax = 0.1, maxSteps = 56, pMax = 15., pMin = -1.,ds = -0.01, nInversion = 4);
         verbosity = 2, plot = true,
         detectCodim2Bifurcation = 2,
         startWithEigen = true,
         bothside=true)

plot(brhopf, brhopf2, legend = :top)
