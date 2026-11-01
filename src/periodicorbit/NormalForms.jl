function BK.period_doubling_normal_form_iooss(pbwrap::BK.PeriodicOrbitFunctionalColl{ <: BK.Collocation{Tprob}},
                                            pd::BK.PeriodDoubling;
                                            nev::Int = 3,
                                            verbose::Bool = false,
                                            lens = BK.getlens(pbwrap),
                                            detailed::Val{detailed_type} = Val(true),
                                            kwargs_nf...) where {Tprob <: AbstractDDEBifurcationProblem, detailed_type}
    if detailed_type
        error("This is not yet implemented.")
    end

    coll = BK.get_discretization(pbwrap)
    x0 = BK.saved_solution(pd.x0)
    par = pd.params
    period = abs(BK.getperiod(coll, x0, par))

    ζ, μ = _pd_floquet_eigenvector(coll, x0, par)
    verbose && println("├─ [PD-DDE] Floquet multiplier μ = ", μ)

    # normalise as the generic Iooss normal form (L² norm over one period)
    n, m, Ntst = size(coll)
    ζm = reshape(ζ, n, m * Ntst + 1)
    nrm = sqrt(BK.∫(coll, ζm, ζm))
    nrm > 0 && (ζ ./= nrm)

    nf = @set pd.nf = (a = 0, b3 = 0)
    return BK.PeriodDoublingPO(pd.x0, period, real.(ζ), nothing, nf, coll, false)
end
