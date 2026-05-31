cd(@__DIR__)
using Pkg, LinearAlgebra, Test
pkg"activate ."
using Revise, BifurcationKit, DDEBifurcationKit, Plots
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

g(z) = (tanh(z − 1) + tanh(1)) * cosh(1)^2

function neuron2VF(x, xd, p)
   (; a,b,c,d) = p
   [
      -x[1] - a * g(b * xd.u[1][1]) + c * g(d * xd.u[2][2]),
      -x[2] - a * g(b * xd.u[1][2]) + c * g(d * xd.u[2][1])
   ]
end

delaysF(par) = [par.τ1, par.τ2]

pars = (a = 0.069, b = 2., c = 0.6, d = 1.2, τ1 = 11.6, τ2 = 20.3)
prob = ConstantDDEBifProblem(neuron2VF, delaysF, zeros(2), pars, (@optic _.c))

optn = NewtonPar(eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 1.0, p_min = 0.25, newton_options = optn, ds = 0.01, detect_bifurcation = 3, nev = 9, dsmax = 0.05, n_inversion = 8)
br = continuation(prob, PALC(), opts; verbosity = 0, bothside = true)

plot(br)

hpnf = BK.get_normal_form(br, 2)
################################################################################
# computation periodic orbit
@views function plot_solution(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, nothing)
    plot!(xtt.t, xtt[1,:]; label = "V1", marker = :d, markersize = 2, k...)
    plot!(xtt.t, xtt[2,:]; label = "V2", k...)
    plot!(br; subplot = 1, putspecialptlegend = false)
end

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.02, ds = -0.0001, p_max = 1.2, p_min=0.100661, max_steps = 200, nev = 10, tol_stability = 1e-4, plot_every_step = 10)

# arguments for periodic orbits
args_po = (  ;
    plot_solution,
    normC = norminf)

probpo = Collocation(20, 5, 
            # jacobian = BK.AutoDiffDense(),
            # jacobian = BK.DenseAnalytical(),
            jacobian = BK.DenseAnalyticalInplace(),
            # meshadapt = true,
            )
br_pocoll = @time continuation(
            br, 2, ContinuationPar(opts_po_cont; max_steps = 200, detect_bifurcation = 2),
            probpo;
            # verbosity = 2,
            plot = true,
            args_po...,
            # eigsolver = BK.FloquetGEV(DDE_DefaultEig(maxit=100, tol = 1e-12, σ = 1e-3)),
            )


- #  1, ns at c ≈ +0.69504545 ∈ (+0.69504545, +0.71231214), |δp|=2e-02
- #  2, pd at c ≈ +0.64290758 ∈ (+0.64290758, +0.66022045), |δp|=2e-02
- #  3, bp at c ≈ +0.46195230 ∈ (+0.46195230, +0.46196651), |δp|=1e-05
- #  4, pd at c ≈ +0.46585936 ∈ (+0.46461857, +0.46585936), |δp|=1e-03
- #  5, pd at c ≈ +0.59688847 ∈ (+0.59454913, +0.59688847), |δp|=2e-03
- #  6, ns at c ≈ +0.61515883 ∈ (+0.61515883, +0.61515932), |δp|=5e-07,
- #  7, bp at c ≈ +0.61512835 ∈ (+0.61512835, +0.61513186), |δp|=4e-06,
- #  8, pd at c ≈ +0.52150492 ∈ (+0.52150492, +0.52212080), |δp|=6e-04
- #  9, pd at c ≈ +0.52237601 ∈ (+0.52166120, +0.52237601), |δp|=7e-04

_color = [:blue, :red, :maroon, :violet, :black]
plot(br_pocoll, vars=(:param, :amplitude), linecolor = [_color[min(5,u+1)] for u in br_pocoll.n_unstable])

plot(br_pocoll.param, pars.τ2 ./ br_pocoll.period, ylabel = "max(tau) / Period")

ind_po = 1
br_pocoll.sol[ind_po].p
period = br_pocoll.sol[ind_po].x[end]
_pars = BK.setparam(br_pocoll,br_pocoll.sol[ind_po].p)
_po = br_pocoll.sol[ind_po].x
_sol = BK.get_periodic_orbit(br_pocoll, ind_po)

_J = @time BK.jacobian(br_pocoll.prob, BK.AutoDiffDense(), br_pocoll.sol[ind_po].x, BK.setparam(br,br_pocoll.sol[ind_po].p));
heatmap(iszero.(_J), yflip=true, color = :viridis)

_J2 = @time DDEBifurcationKit.analytical_jacobian_dde_cst(br_pocoll.prob.prob, br_pocoll.sol[ind_po].x, BK.setparam(br,br_pocoll.sol[ind_po].p));
heatmap(abs.(_J-_J2) .>1e-7, yflip=true, color = :viridis)
(_J-_J2)[1:end-1,1:end-1] |> norminf


_J2 = @time DDEBifurcationKit.analytical_jacobian_dde_cst_floquetcoll(br_pocoll.prob.prob, br_pocoll.sol[ind_po].x, BK.setparam(br,br_pocoll.sol[ind_po].p));

heatmap(abs.(_J2.J0) .>0, yflip=true, color = :viridis)
heatmap(abs.(_J2.Jd) .>0, yflip=true, color = :viridis)

(_J2.J0 + _J2.Jd -_J)[1:end-1,1:end-1] |> norminf

begin
    mu = @time DDEBK.__floquet_coll_gev(BK.FloquetGEV(DDE_DefaultEig(maxit=100, tol = 1e-12, σ = 1e-3, logger = 0), length(probpo), 0), br_pocoll.prob, _po, _pars, 10)[1]
    mu = @time DDEBK.__floquet_coll(BK.FloquetColl(), br_pocoll.prob.prob, _po, _pars, 24)[1];
    plot(layout = 2)
    scatter!(exp.(mu), subplot = 1)
    plot!(cos.(LinRange(0,2pi,100)), sin.(LinRange(0,2pi,100)), subplot = 1, label = "")
    scatter!(mu0, marker=:cross, color=:black, subplot = 1)

    scatter!((mu), subplot = 2); scatter!(log.(mu0), subplot = 2, marker=:cross, color=:black,)
    plot!(0*cos.(LinRange(0,2pi,100)), LinRange(-2pi,2pi,100), label = "", subplot = 2, xlims = (-1,.2))
end



mu = @time DDEBK.__floquet_coll(BK.FloquetColl(), br_pocoll.prob.prob, _po, _pars, 20)[1]

DDEBK.__floquet_coll_gev(BK.FloquetGEV(DDE_DefaultEig(maxit=100, tol = 1e-12, σ = 1e-3, logger = 0)), br_pocoll.prob, _po, _pars, 10)[1]

scatter(vals)


_sol = BK.get_periodic_orbit(br_pocoll, 11)
plot(_sol.t, _sol[1,:])
_sol[1,1]

mu0=[1.002119820872172 + 0.000000000000000im,
  1.000000000016908 + 0.000000000000000im,
 -0.939571779727754 + 0.163504369938072im,
 -0.939571779727754 - 0.163504369938072im,
 -0.804768036913029 + 0.242759690179555im,
 -0.804768036913029 - 0.242759690179555im,
  0.795147997061884 + 0.000000000000000im,
 -0.756988550643102 + 0.022360710023499im,
 -0.756988550643102 - 0.022360710023499im,
  0.708964798330404 + 0.229587565692066im,
  0.708964798330404 - 0.229587565692066im,
  0.658004630535091 + 0.153781117001026im,
  0.658004630535091 - 0.153781117001026im,
  0.482246232875615 + 0.394706666785414im,
  0.482246232875615 - 0.394706666785414im,
 -0.475251580963130 + 0.339692214571253im,
 -0.475251580963130 - 0.339692214571253im,
 -0.485777418170274 + 0.271236751955848im,
 -0.485777418170274 - 0.271236751955848im,
  0.201977571266859 + 0.437443780222598im]


mu0 = [     1.207563057234381 + 0.000000000000000im,
 -1.087640083550283 + 0.000000000000000im,
 -1.014564264466880 + 0.000000000000000im,
  0.999999878103086 + 0.000000000000000im,
 -0.845487685322784 + 0.176000280545985im,
 -0.845487685322784 - 0.176000280545985im,
  0.779185386517539 + 0.156341576349155im,
  0.779185386517539 - 0.156341576349155im,
  0.544783176802123 + 0.336560592534749im,
  0.544783176802123 - 0.336560592534749im,
 -0.562255480173753 + 0.295083355371657im,
 -0.562255480173753 - 0.295083355371657im,
  0.593726906939041 + 0.000000000000000im,
 -0.570000386390988 + 0.023199862309099im,
 -0.570000386390988 - 0.023199862309099im,
  0.490328812338157 + 0.109495381458279im,
  0.490328812338157 - 0.109495381458279im,
  0.290839677165292 + 0.405339617963459im,
  0.290839677165292 - 0.405339617963459im,
 -0.271987089605255 + 0.388379581749829im,]


lmu0 = [    
    0.1886043264932222 + 0.0im,
    0.08401028814511996 + 3.141592653589793im,
   0.014459224228058705 + 3.141592653589793im,
 -1.2189692143653188e-7 + 0.0im,
   -0.14663179355919825 + 2.9363593530623406im,
   -0.14663179355919825 - 2.9363593530623406im,
   -0.22977122529405136 + 0.19801804942500548im,
   -0.22977122529405136 - 0.19801804942500548im,
    -0.4457237699766191 + 0.5533964508569221im,
    -0.4457237699766191 - 0.5533964508569221im,
    -0.4541545703405431 + 2.658286164206433im,
    -0.4541545703405431 - 2.658286164206433im,
    -0.5213358179660129 + 0.0im,
     -0.561290620157923 + 3.100913621552018im,
     -0.561290620157923 - 3.100913621552018im,
    -0.6883471294192669 + 0.21970540063974678im,
    -0.6883471294192669 - 0.21970540063974678im,
    -0.6953762947251434 + 0.94840785108968im,
    -0.6953762947251434 - 0.94840785108968im,
    -0.7462372141878477 + 2.1817320157109563im,]


########
plotH(x) = heatmap(x, yflip=true, color = :viridis)
_coll = br_pocoll.prob.prob

ext_sol = DDEBK.extended_sol(_coll, _po, _pars)
ext_sol(-0.2)

plot(LinRange(ext_sol.interval[1], ext_sol.interval[2], 100), [ext_sol(t)[1] for t in LinRange(ext_sol.interval[1], ext_sol.interval[2], 100)])

_poc = BK.get_time_slices(_coll, _po)
_outc = zero(_poc) .+ 1

DDEBK._residual_for_extended_meshv0!(_coll,
                                    _outc,
                                    _poc,
                                    _po[end],
                                    BK.get_Ls(_coll.mesh_cache),
                                    _pars,
                                    _po)
plot(_outc[1,1:end-1])

DDEBK.jacobian_extended_mesh(_coll, _po, _pars)|>plot

Je = DDEBK.jacobian_extended_mesh(_coll, _po, _pars)
plotH(Je)

_Je, _A, _B = DDEBK.jacobian_extended_mesh(_coll, _po, _pars)
_A1, _B1 = DDEBK.__floquet_coll(BK.FloquetColl(), br_pocoll.prob.prob, _po, _pars, 15)

plotH(_A)
plotH(_B)
plotH(iszero.(_Je))
plotH(iszero.(_A))
plotH(iszero.(_B))

_A1

_A[1:5,1:5]

_B[1:5,1:5]

_J2.Jd

plotH(_J2.Jd)

_J2.Jd[:,end-26:end]