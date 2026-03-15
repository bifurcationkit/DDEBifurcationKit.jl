# TODO use getter from BK
_get_gauss_nodes(coll) = coll.mesh_cache.gauss_nodes

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

function (sol::BK.POSolution{ <: PeriodicOrbitOCollProblem})(::Val{:der}, t0)
    ForwardDiff.derivative(sol, t0)
end

# function for collocation problem
@views function functional_coll!(coll::PeriodicOrbitOCollProblem{Tprob},
                                 outc::AbstractMatrix{𝒯},
                                 uc::AbstractMatrix{𝒯},
                                 period,
                                 (L, ∂L), 
                                 pars, 
                                 u, # uc is a view of u[1:end-1] 
                                 ) where {Tprob <: AbstractDDEBifurcationProblem, 𝒯}
    n, m, Ntst = size(coll)
    # we want slices at fixed times, hence gj[:, j] is the fastest
    # temporaries to reduce allocations
    gj  = BK.get_tmp(coll.cache.gj, uc)  # zeros(𝒯, n, m)
    ∂gj = BK.get_tmp(coll.cache.∂gj, uc) # zeros(𝒯, n, m)
    uj  = zeros(𝒯, n, m+1)

    # get P.O. interpolation which allows to get interp(t)
    interp = BK.POSolution(coll, u)
    VF = coll.prob_vf
    _delays = delays(VF, gj[:, 1], pars)

    # get the mesh of the collocation problem
    mesh = BK.getmesh(coll)
    σs = _get_gauss_nodes(coll)
    udj = VectorOfArray([copy(uj[:, 1]) for _ in _delays])

    # range for locating time slices
    rg = UnitRange(1, m+1)
    for j in 1:Ntst
        uj .= uc[:, rg]
        LA.mul!(gj, uj, L)
        LA.mul!(∂gj, uj, ∂L)

        # get the delayed states
        dτj = (mesh[j+1] - mesh[j]) / 2

        # compute the collocation residual
        for l in 1:m
            τ = BK.τj(σs[l], mesh, j)
            if VF isa SDDDEBifProblem
                _delays = delays(VF, gj[:, l], pars)
            end
            # udj = VectorOfArray([interp(mod(τ * period - d, period)) for d in _delays])
            for (ind, d) in enumerate(_delays)
                udj.u[ind] .= interp(τ * period - d)
            end
            __po_coll_bc!(coll, outc[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dτj, outc[:, end])
        end
        rg = rg .+ m
    end
    # add the periodicity condition
    @. outc[:, end] = uc[:, end] - uc[:, 1]
end

# analytical jacobian for constant DDE
for (fname, floquet) in ((:analytical_jacobian_dde_cst, false), 
                         (:analytical_jacobian_dde_cst_floquetgev, true),
                         (:analytical_jacobian_dde_cst_floquetcoll, true),
                         )
    @eval begin
    @views function $fname(coll::PeriodicOrbitOCollProblem{Tprob}, 
                            u::AbstractVector{𝒯}, 
                            pars;
                            ρD = one(𝒯),
                            ρF = one(𝒯),
                            ρI = zero(𝒯)) where {Tprob <: AbstractDDEBifurcationProblem, 𝒯 }
        n, m, Ntst = size(coll)
        nJ = length(coll) + 1
        L, ∂L = BK.get_Ls(coll.mesh_cache) # L is of size (m+1, m)
        mesh = BK.getmesh(coll)            # coarse mesh of size Ntst + 1
        σs = _get_gauss_nodes(coll)
        ω = coll.mesh_cache.gauss_weight
        period = BK.getperiod(coll, u, nothing)
        phase = zero(𝒯)
        uc = BK.get_time_slices(coll, u)
        pj = BK.get_tmp(coll.cache.gi, u) # zeros(𝒯, n, m)
        In = coll.cache.In # this helps greatly the for loop for J0 below

        # vector field
        interp = BK.POSolution(coll, u)
        VF = coll.prob_vf

        # loop over the mesh intervals
        rg = UnitRange(1, m+1)
        rgNx = UnitRange(1, n)
        rgNy = UnitRange(1, n)

        delays_v = delays(VF, u[1:n], pars) # vector of delays
        udj = VectorOfArray([zeros(𝒯, n) for d in delays_v])
        J = zeros(𝒯, length(coll) + 1, length(coll) + 1)
        J0 = zeros(𝒯, n, n)

        # put boundary condition
        J[nJ-n:nJ-1, nJ-n:nJ-1] .= In
        J[nJ-n:nJ-1, 1:n] .= (-1) .* In

        if $(fname == :analytical_jacobian_dde_cst_floquetgev)
            # arrays to store the jacobian of the delayed terms
            Jd = [zeros(𝒯, length(coll)+1, length(coll)+1) for _ in 1:length(delays_v)]
        elseif $(fname == :analytical_jacobian_dde_cst_floquetcoll)
            # this part contains the times t - d/period which are negative
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
                τ = BK.τj(σs[l], mesh, j) # collocation nodes
                # udj = VectorOfArray([interp(mod(τ * period - d, period)) for d in delays_v])
                for (ind, d) in enumerate(delays_v)
                    udj.u[ind] .= interp(mod(τ * period - d, period))
                end
                JacDDE = jacobian(VF, pj[:, l], udj, pars)
                J0 .= JacDDE.J0
                for l2 in 1:m+1
                    J[_rgX, rgNy .+ (l2-1)*n ] .+= @. (-α * L[l2, l] * ρF) * J0 +
                                                    (ρD * ∂L[l2, l] - α * L[l2, l] * ρI) * In
                    for (idelay, d) in enumerate(delays_v)
                        # find interval where t-τ/period belongs
                        t0 = τ * period - d
                        τd = mod(t0, period) / period
                        index_t = searchsortedfirst(mesh, τd) - 1
                        @assert 1 <= index_t <= Ntst "We have index_t = $index_t, which is out of bounds for mesh of size $(length(mesh)) and τd = $τd. Please open an issue on the website of BifurcationKit.jl"

                        rgNy_delay = UnitRange(1, n) .+ ((m * n) * (index_t - 1))
                        σ = BK.σj(τd, mesh, index_t)
                        β = BK.lagrange(l2, σ, BK.get_mesh_coll(coll)) * ρF

                        if $(fname == :analytical_jacobian_dde_cst_floquetgev)
                            Jd[idelay][_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[idelay] .* β
                        elseif ($(fname == :analytical_jacobian_dde_cst_floquetcoll) && t0 < 0)
                            rgNy_delay = UnitRange(1, n) .+ ((m * n) * (index_t - 1))
                            Jd[_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[idelay] .* β
                        else # case analytical_jacobian_dde_cst
                            J[_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[idelay] .* β
                        end
                    end
                end
                # ================================
                # add derivative w.r.t. the period
                J[_rgX, nJ] .= coll.prob_vf.VF.F(pj[:, l], udj, pars) .* (-dτj)
                for (idelay, d) in enumerate(delays_v)
                    J[_rgX, nJ] .+= -(α * d/period ) .* (JacDDE.Jd[idelay] * interp(Val(:der), τ * period - d))
                end
                # ================================
                phase += LA.dot(pj[:, l], coll.∂ϕ[:, (j-1)*m + l]) * ω[l]
            end
            rg = rg .+ m
            rgNx = rgNx .+ (m * n)
            rgNy = rgNy .+ (m * n)
        end
        if $floquet
            return JacobianDDE(missing, missing, J, Jd, delays_v)
        else
            J[end, begin:end-1] .= coll.cache.∇phase ./ period
            J[nJ, nJ] = -phase / period^2
            return J
        end
    end # function end

    end # begin
end # for-loop end
