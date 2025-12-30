@views function BK.residual!(prob::PeriodicOrbitOCollProblem{Tprob},
                            result,
                            u::AbstractVector, 
                            pars) where {Tprob <: AbstractDDEBifurcationProblem}
    uc = BK.get_time_slices(prob, u)
    T = BK.getperiod(prob, u, nothing)
    resultc = BK.get_time_slices(prob, result)
    functional_coll!(prob, resultc, uc, T, BK.get_Ls(prob.mesh_cache), pars, u)
    # add the phase condition
    result[end] = BK.phase_condition(prob, uc, BK.get_Ls(prob.mesh_cache), T)
    return result
end

function _po_coll_bc!(pb::PeriodicOrbitOCollProblem, dest, ∂u, u, ud, par, h, tmp)
    tmp .= pb.prob_vf.VF.F(u, ud, par)
    dest .= @. ∂u - h * tmp
end

# function for collocation problem
@views function functional_coll!(pb::PeriodicOrbitOCollProblem{Tprob},
                                 out,
                                 u,
                                 period,
                                 (L, ∂L), 
                                 pars, 
                                 result) where {Tprob <: ConstantDDEBifProblem}
    Ty = eltype(u)
    n, ntimes = size(u)
    m = pb.mesh_cache.degree
    Ntst = pb.mesh_cache.Ntst
    # we want slices at fixed  times, hence gj[:, j] is the fastest
    # temporaries to reduce allocations
    gj  = BK.get_tmp(pb.cache.gj, u)  #zeros(𝒯, n, m)
    ∂gj = BK.get_tmp(pb.cache.∂gj, u) #zeros(𝒯, n, m)
    uj  = zeros(Ty, n, m+1)

    # get P.O. interpolation which allows to get result(t)
    interp = BK.POSolution(pb, result)
    delays = pb.prob_vf.delays(pars)

    # get the mesh of the OCollProblem
    mesh = BK.getmesh(pb)
    σ = LinRange(0, 2, m)
    # udj = [copy(uj[:,1]) for _ in delays]

    # range for locating time slices
    rg = UnitRange(1, m+1)
    for j in 1:Ntst
        uj .= u[:, rg]
        LA.mul!(gj, uj, L)
        LA.mul!(∂gj, uj, ∂L)

        # get the delayed states
        tj = mesh[j]
        dtj = (mesh[j+1] - mesh[j]) / 2

        # compute the collocation residual
        for l in 1:m
            tσ = tj + dtj * σ[l]
            udj = [interp(mod(tσ * period - d, period)) for d in delays]
            # for (ind, d) in enumerate(delays)
                # udj[ind] .= interp(mod(tσ*period - d, period))
            # end
            _po_coll_bc!(pb, out[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dtj, out[:, end])

        end
        # carefull here https://discourse.julialang.org/t/is-this-a-bug-scalar-ranges-with-the-parser/70670/4"
        rg = rg .+ m
    end
    # add the periodicity condition
    @. out[:, end] = u[:, end] - u[:, 1]
end


# function for collocation problem
@views function functional_coll!(pb::PeriodicOrbitOCollProblem{Tprob},
                                 out,
                                 u,
                                 period,
                                 (L, ∂L),
                                 pars,
                                 result) where {Tprob <: SDDDEBifProblem}
    𝒯 = eltype(u)
    n, ntimes = size(u)
    m = pb.mesh_cache.degree
    Ntst = pb.mesh_cache.Ntst
    # we want slices at fixed  times, hence gj[:, j] is the fastest
    # temporaries to reduce allocations
    gj  = zeros(𝒯, n, m)
    ∂gj = zeros(𝒯, n, m)
    uj  = zeros(𝒯, n, m+1)

    # get interpolation which allows to get result(t)
    interp = BK.POSolution(pb, result)

    if period <= 0
        out .= Inf
        return out
    end

    # get the mesh of the OCollProblem
    mesh = BK.getmesh(pb)
    σ = LinRange(0, 2, m)

    # range for locating time slices
    rg = UnitRange(1, m+1)
    for j in 1:Ntst
        uj .= u[:, rg]
        LA.mul!(gj, uj, L)
        LA.mul!(∂gj, uj, ∂L)

        # get the delayed states
        tj = mesh[j]
        dtj = (mesh[j+1]-mesh[j]) / 2

        # compute the collocation residual
        for l in 1:m
            tσ = tj + dtj * σ[l]
            delays = pb.prob_vf.delays(gj[:, l], pars)
            udj = [interp(mod(tσ*period - d, period)) for d in delays]
            # out[:, end] can serve as buffer for now in the following function
            _po_coll_bc!(pb, out[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dtj, out[:, end])

        end
        # carefull here https://discourse.julialang.org/t/is-this-a-bug-scalar-ranges-with-the-parser/70670/4"
        rg = rg .+ m
    end
    # add the periodicity condition
    @. out[:, end] = u[:, end] - u[:, 1]
end
