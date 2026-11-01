
using Revise, DDEBifurcationKit
using BifurcationKit
const BK = BifurcationKit
const DDEBK = DDEBifurcationKit

using Plots

function MKG(x, xd, p)
   (;a, b) = p
   y = xd.u[1][1]
   [
      a * x[1] + b * y / (1 + y^10)
   ]
end

delaysF(par) = [2.0]

pars = (a=−1., b=1.5)
x0 = [((1 - 1.5)/(-1))^(1/10)]

prob = ConstantDDEBifProblem(MKG, delaysF, x0, pars, (@optic _.b), record_from_solution=(x,p;k...)-> (x=x[1], _x=1), R01 = BK.AutoDiff(), R11 = BK.AutoDiff(), R02 = BK.AutoDiff())

optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit = 150, tol=1e-8))
opts = ContinuationPar(p_max = 1.6, p_min = 1.2, newton_options = optn, ds = -0.01, dsmax = 0.1, detect_bifurcation = 3, nev = 5, n_inversion = 6, max_steps = 20 )
br = continuation(prob, PALC(), opts)
plot(br)
##########################################################################################
brhopf = continuation(br, 1, (@optic _.a), ContinuationPar(opts, p_min = -2., ds = 0.001, max_steps = 50); 
         bothside = Val(true))
plot(brhopf)
##########################################################################################
BK.get_normal_form(br, 1)
# arguments for periodic orbits
args_po = (    record_from_solution = (x, p; iter, k...) -> begin
        xtt = BK.get_periodic_orbit(BK.getprob(iter), x, nothing)
        mi, ma = extrema(xtt[1,:])
        return (max = ma,
                min = mi,
                amp = ma-mi,
                period = getperiod(BK.getprob(iter), x, nothing))
    end,
    plot_solution = (x, p; iter, k...) -> begin
        xtt = BK.get_periodic_orbit(BK.getprob(iter), x, nothing)
        plot!(xtt.t, xtt[1,:]; label = "V1", k...)
        plot!(br; subplot = 1, putspecialptlegend = false)
        end,
    normC = norminf)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.01, ds= 0.001, dsmin = 1e-4, p_max = 1.9, nev = 20, tol_stability = 1e-5, detect_bifurcation = 3, newton_options = NewtonPar(tol = 1e-12), n_inversion = 4)

probpo = Collocation(20, 4; N = 1, jacobian = BK.DenseAnalyticalInplace(), meshadapt = false)
br_pocoll = @time continuation(
            br, 1, ContinuationPar(opts_po_cont;),
            probpo;
            verbosity = 3,
            plot = true,
            args_po...,
            normC = norminf,
            )

plot(br, br_pocoll)

# witching to the periodic solution branch at the HOPF point (TF).
#     T = 5.490064E+00, arg(Z) = 3.642945E-01 * 2Pi, DST=5.409710E-06

# LABEL TP  NORM           b              T             U  IT
#     1     8.969028E-01   1.336864E+00   5.490064E+00  0  0 +2
#    10     8.972851E-01   1.341726E+00   5.492598E+00  0  1 +0
#    20     8.988527E-01   1.361686E+00   5.501855E+00  0  1 +0
#    30     9.008289E-01   1.386866E+00   5.510990E+00  0  1 +0o|
#    40     9.029960E-01   1.414447E+00   5.517881E+00  0  1 +0
#    50     9.052803E-01   1.443448E+00   5.521765E+00  0  1 +0
#    60     9.076338E-01   1.473227E+00   5.522297E+00  0  1 +0o|
#    70     9.100142E-01   1.503249E+00   5.519381E+00  0  1 +0
#    80     9.123834E-01   1.533058E+00   5.513112E+00  0  1 +0
#    88 PD  9.142482E-01   1.556504E+00   5.505846E+00  1  1 +0
#    90     9.147083E-01   1.562291E+00   5.503740E+00  1  1 +0o|
#   100     9.169638E-01   1.590710E+00   5.491626E+00  1  1 +0
#   110     9.191337E-01   1.618206E+00   5.477211E+00  1  1 +0
#   120 EP  9.212125E-01   1.644807E+00   5.460977E+00  1  1 +0o|

#############################################################################
# branching from period doubling point
BK.get_normal_form(br_pocoll, 1; detailed = Val(false))

bo_pd = continuation(deepcopy(br_pocoll), 1;
               verbosity = 3,
               plot = true,
               args_po...,
               δp = 0.001, ampfactor = 0.1,
               normal_form_kwargs = (detailed = Val(false),),
               normC = norminf,
               )

plot(bo_pd, br_pocoll)

bo_pd2 = continuation(deepcopy(bo_pd), 1;
               verbosity = 3,
               plot = true,
               args_po...,
               δp = 0.001, ampfactor = 0.1,
               normal_form_kwargs = (detailed = Val(false),),
               normC = norminf,
               )

plot(bo_pd, bo_pd2, br_pocoll)