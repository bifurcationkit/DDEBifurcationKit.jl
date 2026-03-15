getP(u, hopfPb::HopfDDEProblem) = u.x[end]
getVec(u, hopfPb::HopfDDEProblem) = u.x[1]

# this function encodes the functional
function (hp::HopfDDEProblem)(x, vR, vI, p::T, ω::T, params) where T
    # These are the equations of the Hopf bifurcation point
    # input:
    # - x guess for the point at which the jacobian has a purely imaginary eigenvalue
    # - p guess for the parameter for which the jacobian has a purely imaginary eigenvalue
    b = hp.b
    # update parameter
    par = set(params, BK.getlens(hp), p)
    tmp = vR .+ 1im .* vI
    w = Δ(hp.prob_vf, x, par, tmp, Complex(zero(T), ω))
    ps = LA.dot(b, tmp)
    return BK.residual(hp.prob_vf, x, par), real(w), imag(w), [real(ps) - 1, imag(ps)]
end

function (hopfpb::HopfDDEProblem)(x::ArrayPartition, params)
    res = hopfpb(x.x[1], x.x[2], x.x[3], x.x[4][1], x.x[4][2], params)
    ArrayPartition(res...)
end
###################################################################################################
struct JacobianCodim2DDE{T1,T2,T3,T4}
    prob::T1
    J::T2
    x::T3
    p::T4
end

(l::BK.DefaultLS)(J::JacobianCodim2DDE, args...; kw...) = l(J.J, args...; kw...)

BK.jacobian(hopfpb::BK.HopfMAProblem{Tprob, BK.AutoDiff, Tu0, Tp, Tl, Tplot, Trecord}, x, p) where {Tprob <: HopfDDEProblem, Tu0, Tp, Tl <: Union{BK.AllOpticTypes, Nothing}, Tplot, Trecord} = JacobianCodim2DDE(hopfpb, ForwardDiff.jacobian(z -> hopfpb.prob(z, p), x), x, p)
################################################################################################### Newton / Continuation functions
function BK.newton_hopf(prob::AbstractDDEBifurcationProblem,
                        hopfpointguess::ArrayPartition,
                        par,
                        eigenvec, eigenvec_ad,
                        options::NewtonPar;
                        normN = LA.norm,
                        bdlinsolver::BK.AbstractBorderedLinearSolver = MatrixBLS(),
                        usehessian = true,
                        kwargs...)
    # we first need to update d2F and d3F for them to accept complex arguments

    hopfproblem = HopfDDEProblem(
        prob,
        BK._copy(eigenvec_ad),    # this is pb.a ≈ null space of (J - iω I)^*
        BK._copy(eigenvec),     # this is pb.b ≈ null space of  J - iω I
        options.linsolver,
        # do not change linear solver if user provides it
        @set bdlinsolver.solver = (isnothing(bdlinsolver.solver) ? options.linsolver : bdlinsolver.solver);
        usehessian = usehessian)

    prob_h = BK.HopfMAProblem(hopfproblem, BK.AutoDiff(), hopfpointguess, par, nothing, prob.plotSolution, prob.recordFromSolution)

    # options for the Newton Solver
    opt_hopf = @set options.linsolver = BK.DefaultLS()

    # solve the hopf equations
    return BK.solve(prob_h, Newton(), opt_hopf; normN, kwargs...)
end

function BK.newton_hopf(br::BK.AbstractResult{Tk, Tp}, 
                        ind_hopf::Int;
                        prob::AbstractDDEBifurcationProblem = br.prob,
                        normN = LA.norm,
                        options = br.contparams.newton_options,
                        verbose = true,
                        nev = br.contparams.nev,
                        start_with_eigen = false,
                        kwargs...) where {Tk, Tp <: AbstractDDEBifurcationProblem}
    hopfpointguess = BK.hopf_point(br, ind_hopf)
    ω = hopfpointguess.p[2]
    bifpt = br.specialpoint[ind_hopf]
    options.verbose && println("--> Newton Hopf, the eigenvalue considered here is ", br.eig[bifpt.idx].eigenvals[bifpt.ind_ev])
    @assert bifpt.idx == bifpt.step + 1 "Error, the bifurcation index does not refer to the correct step"
    ζ = geteigenvector(options.eigsolver, br.eig[bifpt.idx].eigenvecs, bifpt.ind_ev)
    ζ ./= normN(ζ)
    ζad = LA.conj.(ζ)

    if start_with_eigen
        # computation of adjoint eigenvalue. Recall that b should be a null vector of J-iω
        λ = Complex(0, ω)
        p = bifpt.param
        parbif = BK.setparam(br, p)

        # jacobian at bifurcation point
        L = BK.jacobian(prob, bifpt.x, parbif)

        # computation of adjoint eigenvector
        _Jt = ~BK.hasAdjoint(prob) ? adjoint(L) : jad(prob, bifpt.x, parbif)

        ζstar, λstar = BK.getAdjointBasis(_Jt, conj(λ), options.eigsolver; nev = nev, verbose = false)
        ζad .= ζstar ./ LA.dot(ζstar, ζ)
    end

    # solve the hopf equations
    hopfpointguess = ArrayPartition(hopfpointguess.u, real(ζ), imag(ζ), hopfpointguess.p)
    return newton_hopf(prob, hopfpointguess, getparams(br), ζ, ζad, options; normN, kwargs...)
end

function BK.continuation_hopf(prob_vf::AbstractDDEBifurcationProblem, alg::BK.AbstractContinuationAlgorithm,
                            hopfpointguess::ArrayPartition, par,
                            lens1::BK.AllOpticTypes, lens2::BK.AllOpticTypes,
                            eigenvec, eigenvec_ad,
                            options_cont::ContinuationPar ;
                            update_minaug_every_step = 0,
                            normC = norm,
                            bdlinsolver::BK.AbstractBorderedLinearSolver = MatrixBLS(),
                            jacobian_ma::Symbol = :autodiff,
                            compute_eigen_elements = false,
                            usehessian = true,
                            massmatrix = LA.I,
                            kwargs...)
    @assert lens1 != lens2 "Please choose 2 different parameters. You only passed $lens1"
    @assert lens1 == BK.getlens(prob_vf)

    # options for the Newton Solver inheritated from the ones the user provided
    options_newton = options_cont.newton_options
    threshBT = 100options_newton.tol

    𝐇 = HopfDDEProblem(
            prob_vf,
            BK._copy(eigenvec_ad),    # this is a ≈ null space of (J - iω I)^*
            BK._copy(eigenvec),       # this is b ≈ null space of  J - iω I
            options_newton.linsolver,
            # do not change linear solver if user provides it
            @set bdlinsolver.solver = (isnothing(bdlinsolver.solver) ? options_newton.linsolver : bdlinsolver.solver);
            usehessian,
            massmatrix,
            _norm = normC,
            update_minaug_every_step)

    # Jacobian for the Hopf problem
    if jacobian_ma == :autodiff
        # hopfpointguess = vcat(hopfpointguess.u, hopfpointguess.p)
        prob_h = BK.HopfMAProblem(𝐇, BK.AutoDiff(), hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
        opt_hopf_cont = @set options_cont.newton_options.linsolver = DefaultLS()
    elseif jacobian_ma == :finiteDifferencesMF
        hopfpointguess = vcat(hopfpointguess.u, hopfpointguess.p)
        prob_h = FoldMAProblem(𝐇, FiniteDifferences(), hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
        opt_hopf_cont = @set options_cont.newton_options.linsolver = options_cont.newton_options.linsolver
    else
        prob_h = BK.HopfMAProblem(𝐇, nothing, hopfpointguess, par, lens2, prob_vf.plotSolution, prob_vf.recordFromSolution)
        opt_hopf_cont = @set options_cont.newton_options.linsolver = HopfLinearSolverMinAug()
    end

    # this functions allows to tackle the case where the two parameters have the same name
    lenses = BK.get_lens_symbol(lens1, lens2)

    # current lyapunov coefficient
    eTb = eltype(hopfpointguess)
    𝐇.l1 = Complex{eTb}(0, 0)
    𝐇.BT = one(eTb)
    𝐇.GH = one(eTb)

    # this function is used as a Finalizer
    # it is called to update the Minimally Augmented problem
    # by updating the vectors a, b
    function updateMinAugHopf(z, tau, step, contResult; kUP...)
        # we first check that the continuation step was successful
        # if not, we do not update the problem with bad information!
        success = get(kUP, :state, nothing).converged
        (~BK.mod_counter(step, update_minaug_every_step) || success == false) && return true
        x = getVec(z.u, 𝐇)    # hopf point
        p1, ω = getP(z.u, 𝐇)
        p2 = z.p        # second parameter
        newpar = set(par, lens1, p1)
        newpar = set(newpar, lens2, p2)

        # we stop continuation at Bogdanov-Takens points

        # CA NE DEVRAIT PAS ETRE ISSNOT?
        isbt = isnothing(contResult) ? true : isnothing(findfirst(x -> x.type in (:bt, :ghbt, :btgh), contResult.specialpoint))

        # if the frequency is null, this is not a Hopf point, we halt the process
        if abs(ω) < threshBT
            @warn "[Codim 2 Hopf - Finalizer] The Hopf curve seems to be close to a BT point: ω ≈ $ω. Stopping computations at ($p1, $p2). If the BT point is not detected, try lowering Newton tolerance or dsmax."
        end

        # call the user-passed finalizer
        finaliseUser = get(kwargs, :finaliseSolution, nothing)
        resFinal = isnothing(finaliseUser) ? true : finaliseUser(z, tau, step, contResult; prob = 𝐇, kUP...)

        return abs(ω) >= threshBT && isbt && resFinal
    end

    # the following allows to append information specific to the codim 2 continuation to the user data
    _printsol = get(kwargs, :record_from_solution, nothing)
    _printsol2 = isnothing(_printsol) ?
        (u, p; kw...) -> (; zip(lenses, (getP(u, 𝐇)[1], p))..., ω = getP(u, 𝐇)[2], l1 = 𝐇.l1, BT = 𝐇.BT, GH = 𝐇.GH, BK._namedrecordfromsol(BK.record_from_solution(prob_vf)(getVec(u, 𝐇), p; kw...))...) :
        (u, p; kw...) -> (; BK._namedrecordfromsol(_printsol(getVec(u, 𝐇), p; kw...))..., zip(lenses, (getP(u, 𝐇)[1], p))..., ω = getP(u, 𝐇)[2], l1 = 𝐇.l1, BT = 𝐇.BT, GH = 𝐇.GH)

    prob_h = re_make(prob_h, record_from_solution = _printsol2)

    # eigen solver
    eigsolver = HopfDDEEig(BK.getsolver(opt_hopf_cont.newton_options.eigsolver))

    # event for detecting codim 2 points
    if compute_eigen_elements
        event = BK.PairOfEvents(BK.ContinuousEvent(2, testBT_GH, compute_eigen_elements, ("bt", "gh"), threshBT), BifDetectEvent)
        # careful here, we need to adjust the tolerance for stability to avoid
        # spurious ZH or HH bifurcations
        @reset opt_hopf_cont.tol_stability = 10opt_hopf_cont.newton_options.tol
    else
        event = BK.ContinuousEvent(2, testBT_GH, false, ("bt", "gh"), threshBT)
    end

    prob_h = re_make(prob_h, record_from_solution = _printsol2)

    # solve the hopf equations
    br = continuation(
        prob_h, alg,
        (@set opt_hopf_cont.newton_options.eigsolver = eigsolver);
        kwargs...,
        kind = BK.HopfCont(),
        linear_algo = BorderingBLS(solver = opt_hopf_cont.newton_options.linsolver, check_precision = false),
        normC = normC,
        finalise_solution = update_minaug_every_step == 0 ? get(kwargs, :finalise_solution, BK.finalise_default) : updateMinAugHopf,
        event = event
    )
    @assert ~isnothing(br) "Empty branch!"
    return correctBifurcation(br)
end

function BK.continuation_hopf(prob::AbstractDDEBifurcationProblem,
                        br::BK.AbstractBranchResult, ind_hopf::Int64,
                        lens2::BK.AllOpticTypes,
                        options_cont::ContinuationPar = br.contparams;
                        alg = br.alg,
                        start_with_eigen = false,
                        normC = LA.norm,
                        kwargs...)
    hopfpointguess = BK.hopf_point(br, ind_hopf)
    ω = hopfpointguess.p[2]
    bifpt = br.specialpoint[ind_hopf]

    @assert ~isnothing(br.eig) "The branch contains no eigen elements. This is strange because a Hopf point was detected. Please open an issue on the website."

    @assert ~isnothing(br.eig[1].eigenvecs) "The branch contains no eigenvectors for the Hopf point. Please provide one."

    ζ = geteigenvector(br.contparams.newton_options.eigsolver, br.eig[bifpt.idx].eigenvecs, bifpt.ind_ev)
    ζ ./= normC(ζ)
    ζad = conj.(ζ)

    p = bifpt.param
    parbif = BK.setparam(br, p)

    hopfpointguess = ArrayPartition(hopfpointguess.u, real(ζ), imag(ζ), hopfpointguess.p)

    if start_with_eigen
        # computation of adjoint eigenvalue
        λ = Complex(0, ω)
        # jacobian at bifurcation point
        L = BK.jacobian(prob, bifpt.x, parbif)

        # jacobian adjoint at bifurcation point
        _Jt = ~BK.has_adjoint(prob) ? adjoint(L) : BK.jacobian_adjoint(prob, bifpt.x, parbif)

        ζstar, λstar = BK.get_adjoint_basis(_Jt, conj(λ), br.contparams.newton_options.eigsolver; nev = br.contparams.nev, verbose = false)
        ζad .= ζstar ./ LA.dot(ζstar, ζ)
    end

    return BK.continuation_hopf(BK.getprob(br), alg,
                    hopfpointguess, parbif,
                    BK.getlens(br), lens2,
                    ζ, ζad,
                    options_cont ;
                    normC = normC,
                    kwargs...)
end

function testBT_GH(iter, state)
    probma = BK.getprob(iter)
    lens1, lens2 = BK.get_lenses(probma)

    𝐇 = probma.prob
    𝒯 = eltype(𝐇) 
    z = getx(state)
    x = getVec(z, 𝐇)          # hopf point
    p1, ω = getP(z, 𝐇)        # first parameter
    p2 = getp(state)           # second parameter
    par = getparams(probma)
    newpar = set(par, lens1, p1)
    newpar = set(newpar, lens2, p2)

    # expression of the jacobian
    J_at_xp = BK.jacobian(𝐇.prob_vf, x, newpar)

    # compute new b
    T = typeof(p1)
    # ζ = 𝐇.linbdsolver(J_at_xp, a, b, T(0), 𝐇.zero, T(1); shift = Complex(0, -ω), Mass = hopfPb.massmatrix)[1]
    λ = Complex(0, ω)
    ζ = @. z.x[2] + im * z.x[3]
    ζ ./= 𝐇.norm(ζ)

    # compute new ζstar
    # JAd_at_xp = BK.hasAdjoint(𝐇.prob_vf) ? jad(𝐇.prob_vf, x, newpar) : transpose(J_at_xp)

    JAd_at_xp = BK.has_adjoint(𝐇.prob_vf) ? BK.jacobian_adjoint(𝐇.prob_vf, x, newpar) : adjoint(J_at_xp)
    ζ★, _ = BK.get_adjoint_basis(JAd_at_xp, conj(λ), BK.getcontparams(iter).newton_options.eigsolver.eigsolver)
    # ζ★ = 𝐇.linbdsolver(JAd_at_xp, b, a, T(0), hopfPb.zero, T(1); shift = Complex(0, ω), Mass = transpose(hopfPb.massmatrix))[1]

    # test function for Bogdanov-Takens
    𝐇.BT = ω
    # BT2 = real( LA.dot(ζ★ ./ 𝐇.norm(ζ★), ζ) )
    # ζ★ ./= LA.dot(ζ, ζ★)

    hp0 = BK.Hopf(x, nothing, p1, ω, newpar, lens1, ζ, ζ★, (a = Complex{T}(0, 0), b = Complex{T}(0, 0)), :hopf)
    hp = BK.hopf_normal_form(𝐇.prob_vf, hp0, 𝐇.linsolver; verbose = false) #increase nev?

    # lyapunov coefficient
    𝐇.l1 = hp.nf.b
    # test for Bautin bifurcation.
    # If GH is too large, we take the previous value to avoid spurious detection
    # GH will be large close to BR points
    𝐇.GH = abs(real(hp.nf.b)) < 1e5 ? real(hp.nf.b) : state.eventValue[2][2]
    return 𝐇.BT, 𝐇.GH
    end

# structure to compute the eigenvalues along the Hopf branch
struct HopfDDEEig{S} <: BK.AbstractCodim2EigenSolver
    eigsolver::S
end

function (eig::HopfDDEEig)(Jddehopf, nev; kwargs...)
    xh = Jddehopf.x.x[1]            # hopf point
    p1, ω = Jddehopf.x.x[4]         # first parameter
    newpar = set(Jddehopf.p, BK.getlens(Jddehopf.prob.prob), p1)
    J = BK.jacobian(Jddehopf.prob.prob.prob_vf, xh, newpar)
    eigenelts = eig.eigsolver(J, nev; kwargs...)
    return eigenelts
end
