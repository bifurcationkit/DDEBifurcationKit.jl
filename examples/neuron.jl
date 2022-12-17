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

function delaysF(par)
   [par.τ1, par.τ2, par.τs]
end


pars = (κ = 0.5, β=-1, a12=1, a21=0.5, τ1=0.2, τ2=0.2, τs=1.5)
x0 = zeros(2)
delays0 = [0.2,0.2,1.5]

# xd0 = VectorOfArray([copy(x0) for _ in eachindex(delays0)])
# neuronVF(x0,xd0,pars)

prob = DDEBK.ConstantDDEBifProblem(neuronVF, delaysF, x0, delays0, pars, (@lens _.τs))

# BK.residual(prob, x0, pars)

# BK.jacobian(prob, x0, pars)

optn = NewtonPar(verbose = true, eigsolver = DDEBK.DDE_NLEVEigSolver(maxit=100))
opts = ContinuationPar(pMax = 13., pMin = 0., newtonOptions = optn, ds = -0.01, detectBifurcation = 3, nev = 3, )
br = BK.continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = true)
plot(br)
################################################################################
pthopf = newton(br, 3,; lens2 = (@lens _.τs),
         startWithEigen = true)

brhopf = continuation(br, 3, (@lens _.a21),
         setproperties(br.contparams, detectBifurcation = 0, dsmax = 0.04, maxSteps = 230, pMax = 15., pMin = -1.,ds = 0.01);
         verbosity = 3,
         detectCodim2Bifurcation = 2,
         startWithEigen = true)

plot(brhopf, vars = (:a21, :τs))
plot(brhopf, vars = (:τs, :ω))

brhopf2 = continuation(br, 2, (@lens _.a21),
         setproperties(br.contparams, detectBifurcation = 0, dsmax = 0.04, maxSteps = 230, pMax = 15., pMin = -1.,ds = -0.01);
         verbosity = 3,
         detectCodim2Bifurcation = 2,
         startWithEigen = true,
         bothside=true)

plot(brhopf, brhopf2)
