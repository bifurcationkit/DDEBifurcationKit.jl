@views function BK.residual!(coll::PeriodicOrbitOCollProblem{Tprob},
                            result,
                            u::AbstractVector, 
                            pars) where {Tprob <: AbstractDDEBifurcationProblem}
    uc = BK.get_time_slices(coll, u)
    period = BK.getperiod(coll, u, nothing)
    resultc = BK.get_time_slices(coll, result)
    functional_coll!(coll, resultc, uc, period, BK.get_Ls(coll.mesh_cache), pars, u)
    # add the phase condition
    result[end] = BK.phase_condition(coll, uc, BK.get_Ls(coll.mesh_cache), period)
    return result
end

function __po_coll_bc!(coll::PeriodicOrbitOCollProblem, dest, ∂u, u, ud, par, h, tmp)
    tmp .= coll.prob_vf.VF.F(u, ud, par)
    @. dest = ∂u - h * tmp
end

# function for collocation problem
@views function functional_coll!(coll::PeriodicOrbitOCollProblem{Tprob},
                                 out,
                                 u,
                                 period,
                                 (L, ∂L), 
                                 pars, 
                                 result) where {Tprob <: AbstractDDEBifurcationProblem}
    𝒯 = eltype(u)
    n, ntimes = size(u)
    m = coll.mesh_cache.degree
    Ntst = coll.mesh_cache.Ntst
    # we want slices at fixed times, hence gj[:, j] is the fastest
    # temporaries to reduce allocations
    gj  = BK.get_tmp(coll.cache.gj, u)  # zeros(𝒯, n, m)
    ∂gj = BK.get_tmp(coll.cache.∂gj, u) # zeros(𝒯, n, m)
    uj  = zeros(𝒯, n, m+1)

    # get P.O. interpolation which allows to get result(t)
    interp = BK.POSolution(coll, result)
    VF = coll.prob_vf
    _delays = delays(VF, gj[:, 1], pars)

    # get the mesh of the collocation problem
    mesh = BK.getmesh(coll)
    σs = LinRange{𝒯}(0, 2, m) # TODO: better to rely on BK.get_mesh_coll(coll)
    udj = VectorOfArray([copy(uj[:, 1]) for _ in _delays])

    # TODO: there is an issue here. If we use `σs = BK.get_mesh_coll(coll)` which is equivalent to
    # choosing `σs = LinRange(0, 2, m+1)`, do we take `σs[l+1]` below which runs for l in 1:m ?
    # if we do, newton does not converge which indicates an issue with `interp`

    # range for locating time slices
    rg = UnitRange(1, m+1)
    for j in 1:Ntst
        uj .= u[:, rg]
        LA.mul!(gj, uj, L)
        LA.mul!(∂gj, uj, ∂L)

        # get the delayed states
        τj = mesh[j]
        dτj = (mesh[j+1] - mesh[j]) / 2

        # compute the collocation residual
        for l in 1:m
            τ = τj + dτj * σs[l]
            if VF isa SDDDEBifProblem
                _delays = delays(VF, gj[:, l], pars)
            end
            # udj = VectorOfArray([interp(mod(τ * period - d, period)) for d in _delays])
            for (ind, d) in enumerate(_delays)
                udj.u[ind] .= interp(τ * period - d)
            end
            __po_coll_bc!(coll, out[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dτj, out[:, end])

        end
        rg = rg .+ m
    end
    # add the periodicity condition
    @. out[:, end] = u[:, end] - u[:, 1]
end

# analytical jacobian for constant DDE
# there is a function specific to the computation of Floquet coefficients which is slightly different
# and relies on the same code
for (fname, floquet) in ((:analytical_jacobian_dde_cst, false), 
                         (:analytical_jacobian_dde_cst_floquetgev, true),
                         (:analytical_jacobian_dde_cst_floquetcoll, true),
                         )
    @eval begin
    function $fname(wrap::BK.WrapPOColl{ <: PeriodicOrbitOCollProblem{Tprob}}, 
                            u::AbstractVector{𝒯}, 
                            pars;
                            ρD = one(𝒯),
                            ρF = one(𝒯),
                            ρI = zero(𝒯)) where {Tprob <: ConstantDDEBifProblem, 𝒯 }
        coll = wrap.prob # TODO: use getdisc
        J = zeros(𝒯, length(coll)+1, length(coll)+1)

        n, m, Ntst = size(coll)
        nJ = length(coll) + 1
        L, ∂L = BK.get_Ls(coll.mesh_cache) # L is of size (m+1, m)
        mesh = BK.getmesh(coll)
        σs = LinRange{𝒯}(0, 2, m) # TODO: better to rely on BK.get_mesh_coll(coll)
        ω = coll.mesh_cache.gauss_weight
        period = BK.getperiod(coll, u, nothing)
        phase = zero(𝒯)
        uc = BK.get_time_slices(coll, u)
        pj = BK.get_tmp(coll.cache.gi, u) # zeros(𝒯, n, m)
        In = coll.cache.In # this helps greatly the for loop for J0 below

        # vector field
        interp = BK.POSolution(coll, u)
        VF = coll.prob_vf

        # put boundary condition
        J[nJ-n:nJ-1, nJ-n:nJ-1] .= In
        J[nJ-n:nJ-1, 1:n] .= (-1) .* In

        # loop over the mesh intervals
        rg = UnitRange(1, m+1)
        rgNx = UnitRange(1, n)
        rgNy = UnitRange(1, n)

        delays = VF.delays(pars)
        udj = VectorOfArray([zeros(𝒯, n) for d in delays])
        J0 = zeros(𝒯, n, n)

        if $(fname == :analytical_jacobian_dde_cst_floquetgev)
            # arrays to store the jacobian of the delayed terms
            Jd = [zeros(𝒯, length(coll)+1, length(coll)+1) for _ in 1:length(delays)]
        elseif $(fname == :analytical_jacobian_dde_cst_floquetcoll)
            Jd = zeros(𝒯, length(coll)+1, length(coll)+1)
        end

        # !!!!!  SD-DDE  luzyanina_computing_2002

        for j in 1:Ntst
            LA.mul!(pj, uc[:, rg], L) # pj ≈ (L * uj')'
            τj = mesh[j]
            dτj = (mesh[j+1] - mesh[j]) / 2
            α = period * dτj
            for l in 1:m
                _rgX = rgNx .+ (l-1)*n
                τ = τj + dτj * σs[l]
                # udj = VectorOfArray([interp(mod(τ * period - d, period)) for d in delays])
                for (ind, d) in enumerate(delays)
                    udj.u[ind] .= interp(mod(τ * period - d, period))
                end
                JacDDE = BK.jacobian(VF, pj[:, l], udj, pars)
                J0 .= JacDDE.J0
                for l2 in 1:m+1
                    J[_rgX, rgNy .+ (l2-1)*n ] .+= @. (-α * L[l2, l] * ρF) * J0 +
                                                    (ρD * ∂L[l2, l] - α * L[l2, l] * ρI) * In
                    for (idelay, d) in enumerate(delays)
                        # find interval where t-τ/period belongs
                        t0 = τ * period - d
                        τd = mod(t0, period) / period
                        index_t = searchsortedfirst(mesh, τd) - 1
                        # index_t = max(1, min(index_t, Ntst))
                        @assert 1 <= index_t <= Ntst
                        rgNy_delay = UnitRange(1, n) .+ ((m * n) * (index_t-1))
                        σ = BK.σj(τd, mesh, index_t)
                        β = BK.lagrange(l2, σ, BK.get_mesh_coll(coll))
                        if $(fname == :analytical_jacobian_dde_cst_floquetgev)
                            Jd[idelay][_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[idelay] .* β
                        elseif ($(fname == :analytical_jacobian_dde_cst_floquetcoll) && t0 < 0)
                            Jd[_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[idelay] .* β
                        else
                            J[_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[idelay] .* β
                        end
                    end
                end
            end
            rg = rg .+ m
            rgNx = rgNx .+ (m * n)
            rgNy = rgNy .+ (m * n)
        end
        if $floquet
            return JacobianDDE(missing, missing, J, Jd, delays)
        else
            return J
        end
    end # function end

    end # begin
end # for-loop end
