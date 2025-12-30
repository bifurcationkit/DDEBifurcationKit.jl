@views function BK.residual!(coll::PeriodicOrbitOCollProblem{Tprob},
                            result,
                            u::AbstractVector, 
                            pars) where {Tprob <: AbstractDDEBifurcationProblem}
    uc = BK.get_time_slices(coll, u)
    T = BK.getperiod(coll, u, nothing)
    resultc = BK.get_time_slices(coll, result)
    functional_coll!(coll, resultc, uc, T, BK.get_Ls(coll.mesh_cache), pars, u)
    # add the phase condition
    result[end] = BK.phase_condition(coll, uc, BK.get_Ls(coll.mesh_cache), T)
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
            τ = τj + dτj * (σs[l])
            if VF isa SDDDEBifProblem
                _delays = delays(VF, gj[:, l], pars)
            end
            udj = VectorOfArray([interp(mod(τ * period - d, period)) for d in _delays])
            # for (ind, d) in enumerate(_delays)
                # udj.u[ind] .= interp(τ * period - d)
            # end
            __po_coll_bc!(coll, out[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dτj, out[:, end])

        end
        rg = rg .+ m
    end
    # add the periodicity condition
    @. out[:, end] = u[:, end] - u[:, 1]
end

