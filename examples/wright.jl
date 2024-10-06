cd(@__DIR__)
cd("..")
# using Pkg, LinearAlgebra, Test
# pkg"activate ."

using Revise, DDEBifurcationKit
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

function wrightVF(x, xd, p)
   (;a) = p
   y = xd[1][1]
   [
      -a * y * (1 + x[1])
   ]
end

delaysF(par) = [1.0]

pars = (a=0.1,b=0.)
x0 = [0.]

prob = ConstantDDEBifProblem(wrightVF, delaysF, x0, pars, (@optic _.a), record_from_solution=(x,p;k...)-> (x=x[1], _x=1))

optn = NewtonPar(verbose = false, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 9., p_min = 0., newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 4, n_inversion = 12 )
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = false)
plot(br)

BK.get_normal_form(br, 1)

function NFcoeff(N)
   # Faria 2006
   ωN = pi/2 + N*pi
   aN = ωN*(-1)^N
   K1 = aN / (1+aN^2)
   K2 = ωN / (5(1+aN^2))*((-1)^N-3ωN)
   return (a=K1, b=K2)
end

h1 = BK.get_normal_form(br, 1)
@test isapprox(real(h1.nf.a), NFcoeff(0).a; rtol = 1e-5)
@test isapprox(real(h1.nf.b), NFcoeff(0).b; rtol = 1e-5)

h2 = BK.get_normal_form(br, 2)
@test isapprox(real(h2.nf.a), NFcoeff(2).a; rtol = 1e-5)
@test isapprox(real(h2.nf.b), NFcoeff(2).b; rtol = 1e-5)
