# TODO the bottleneck is computhing the jacobian.
# J: 0.011025 seconds (58.15 k allocations: 68.364 MiB)
# L: 0.001172 seconds (16 allocations: 3.752 MiB)

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
                # udj.u[ind] .= interp(τ * period - d)
                udj.u[ind] .= BK.__interpolate_posolution(coll, τ - d/period, u, 1)
            end
            __po_coll_bc!(coll, outc[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dτj, outc[:, end])
        end
        rg = rg .+ m
    end
    # add the periodicity condition
    @. outc[:, end] = uc[:, end] - uc[:, 1]
end

function BK.jacobian(coll::PeriodicOrbitOCollProblem{Tprob}, 
                ::BK.DenseAnalytical,
                x, 
                p) where {Tprob <: ConstantDDEBifProblem}
    return analytical_jacobian_dde_cst(coll, x, p)
end

function BK.jacobian(coll::PeriodicOrbitOCollProblem{Tprob}, 
                ::BK.FullSparse,
                x, 
                p) where {Tprob <: ConstantDDEBifProblem}
    return analytical_jacobian_dde_cst(coll, x, p)
end

"""
using DifferentiationInterface to automatically derive the sparse jacobian.
"""
struct AutoSparseDI <: BK.AbstractJacobianSparseMatrix end

function BK._generate_jacobian(coll::PeriodicOrbitOCollProblem{Tprob}, ::AutoSparseDI, orbitguess, pars; k...) where {Tprob <: AbstractDDEBifurcationProblem}
    error("You need to import `DifferentiationInterface, SparseConnectivityTracer, SparseMatrixColorings` in order to use this jacobian")
end
########################################################################################
# analytical jacobians for constant DDE
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
        if coll.jacobian isa BK.FullSparse
            J = SA.spzeros(𝒯, length(coll) + 1, length(coll) + 1)
            J0 = SA.spzeros(𝒯, n, n)
            In = SA.sparse(In)
        else
            J = zeros(𝒯, length(coll) + 1, length(coll) + 1)
            J0 = zeros(𝒯, n, n)
        end

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
                        t0 = τ - d/period
                        τd = mod(t0, 1) / 1
                        index_t = searchsortedfirst(mesh, τd) - 1
                        @assert 1 <= index_t <= Ntst "We have index_t = $index_t, which is out of bounds for mesh of size $(length(mesh)) and τd = $τd. Please open an issue on the website of BifurcationKit.jl"

                        rgNy_delay = UnitRange(1, n) .+ ((m * n) * (index_t - 1))
                        σ = BK.σj(τd, mesh, index_t)
                        β = BK.lagrange(l2, σ, BK.get_mesh_coll(coll)) * ρF

                        if $(fname == :analytical_jacobian_dde_cst_floquetgev)
                            Jd[idelay][_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[idelay] .* β
                        elseif ($(fname == :analytical_jacobian_dde_cst_floquetcoll) && t0 < 0)
                            rgNy_delay = UnitRange(1, n) .+ ((m * n) * (index_t - 1))
                            # fullmesh = coll.mesh_cache.full_mesh
                            # index_tau = searchsortedlast(fullmesh, τd)- 2
                            # rgNy_delay = UnitRange(1, n) .+ index_tau * n
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
########################################################################################
struct ExtendedSolution{Tpb, Tx, T}
    extended_coll::Tpb 
    xc::Tx # AbstractMatrix
    interval::Tuple{T, T} # (-τ/D, 1])
    initial_n::Int # number of unknowns initial
end

function (sol::ExtendedSolution)(t)
    if ~(sol.interval[1] <= t <= sol.interval[2])
        error("You passed t=$t and \n $(sol.interval)")
    end
    # extended problem
    extended_coll = sol.extended_coll
    n, m, Ntst = size(extended_coll)
    xc = sol.xc
    mesh = BK.getmesh(extended_coll)
    index_t = searchsortedfirst(mesh, t) - 1
    if index_t <= 0
        return xc[:, 1]
    elseif index_t > Ntst
        return xc[:, end]
    end
    @assert mesh[index_t] <= t <= mesh[index_t+1] "Please open an issue on the website of BifurcationKit.jl"
    σ = BK.σj(t, mesh, index_t)
    # @assert -1 <= σ <= 1 "Strange value of $σ"
    σs = BK.get_mesh_coll(extended_coll)
    out = zeros(typeof(t), n)
    rg = (1:m+1) .+ (index_t - 1) * m
    for l in 1:m+1
        out .+= xc[:, rg[l]] .* BK.lagrange(l, σ, σs)
    end
    out
end

@views function extended_sol(coll::BK.PeriodicOrbitOCollProblem,
                            periodic_sol::AbstractVector, 
                            pars)
    result = copy(periodic_sol)
    periodic_solc = BK.get_time_slices(coll, periodic_sol)
    period = BK.getperiod(coll, periodic_sol, nothing)
    ratio = maximum(delays(coll.prob_vf, nothing, pars)) / period

    # from the mesh 0 = τ₁ < ... < τₙₜₛₜ₊₁ = 1, we build the new one:
    # τ₋ₙ < ratio < ... < τₙₜₛₜ-1 < τ₁ < ... < τₙₜₛₜ₊₁ = 1
    mesh = BK.getmesh(coll)      # τᵢ
    times = BK.get_times(coll)   # tᵢ

    extended_mesh = copy(mesh)
    extended_times = copy(times)
    extended_solution = copy(periodic_solc)

    for i = reverse(eachindex(mesh))
        new_t = mesh[i] - 1
        if new_t < 0
            pushfirst!(extended_mesh, new_t)
        end
        if new_t < -ratio
            break
        end
    end

    for i = reverse(eachindex(times))
        new_t = times[i] - 1
        if new_t < 0
            pushfirst!(extended_times, new_t)
            extended_solution = hcat(periodic_solc[:, i], extended_solution) 
        end
        if new_t < -ratio
            break
        end
    end
    # @error "" ratio
    # return extended_mesh

    extended_coll = deepcopy(coll)
    @reset extended_coll.mesh_cache.τs = extended_mesh
    @reset extended_coll.mesh_cache.Ntst = length(extended_mesh) - 1
    @reset extended_coll.mesh_cache.full_mesh = extended_times

    return ExtendedSolution(extended_coll, extended_solution, (-ratio, 1.0), length(coll))
end

# Continuation and Bifurcation Analysis of Delay Differential Equations page 10

@views function _residual_for_extended_meshv0!(coll::PeriodicOrbitOCollProblem{Tprob},
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
    gj  = BK.get_tmp(coll.cache.gj, u)  # zeros(𝒯, n, m)
    ∂gj = BK.get_tmp(coll.cache.∂gj, u) # zeros(𝒯, n, m)
    uj  = zeros(𝒯, n, m+1)

    # get P.O. interpolation which allows to get result(t)
    # interp = BK.POSolution(coll, u, pars)
    interp = extended_sol(coll, u, pars)
    VF = coll.prob_vf
    _delays = delays(VF, gj[:, 1], pars)

    # get the mesh of the collocation problem
    mesh = BK.getmesh(coll)
    σs = _get_gauss_nodes(coll)
    udj = VectorOfArray([copy(uj[:, 1]) for _ in _delays])

    # range for locating time slices
    rg = UnitRange(1, m+1)
    eq = 1
    for j in 1:Ntst
        uj .= uc[:, rg]
        LA.mul!(gj, uj, L)
        LA.mul!(∂gj, uj, ∂L)

        # get the delayed states
        dτj = (mesh[j+1] - mesh[j]) / 2

        # compute the collocation residual
        if mesh[j]>=0
            for l in 1:m
                τ = BK.τj(σs[l], mesh, j)
                # udj = VectorOfArray([interp(mod(τ * period - d, period)) for d in _delays])
                for (ind, d) in enumerate(_delays)
                    udj.u[ind] .= interp(τ - d/period)
                    # udj.u[ind] .= BK.__interpolate_posolution(coll, τ - d/period, u, 1)
                end
                __po_coll_bc!(coll, outc[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dτj, outc[:, end])
            end
        end
        rg = rg .+ m
        eq += 1
    end
    # add the periodicity condition
    index = interp.initial_n
    # @. outc[:, end] = uc[:, end] - uc[:, end-index+1]
    return outc
end

function _residual_for_extended_mesh(coll_ext::PeriodicOrbitOCollProblem,
                                    interp,
                                    u,
                                    pars,
                                    )
    uc = BK.get_time_slices(coll_ext, u)
    period = BK.getperiod(coll_ext, u, nothing)
    outc = 0*(uc) .+ 0
    out = vec(outc)
    _residual_for_extended_mesh!(coll_ext, interp, outc, uc, period, BK.get_Ls(coll_ext.mesh_cache), pars, u)
    return out
end

@views function _residual_for_extended_mesh!(coll::PeriodicOrbitOCollProblem,
                                    interp,
                                    outc::AbstractMatrix{𝒯},
                                    uc::AbstractMatrix{𝒯},
                                    period,
                                    (L, ∂L), 
                                    pars, 
                                    u, # uc is a view of u[1:end-1] 
                                    ) where {𝒯}
    n, m, Ntst = size(coll)
    # we want slices at fixed times, hence gj[:, j] is the fastest
    # temporaries to reduce allocations
    gj  = BK.get_tmp(coll.cache.gj, u)  # zeros(𝒯, n, m)
    ∂gj = BK.get_tmp(coll.cache.∂gj, u) # zeros(𝒯, n, m)
    uj  = zeros(𝒯, n, m+1)

    # get P.O. interpolation which allows to get result(t)
    # interp = BK.POSolution(coll, u, pars)
    # interp = extended_sol(coll, u, pars)
    VF = coll.prob_vf
    _delays = delays(VF, gj[:, 1], pars)

    # get the mesh of the collocation problem
    mesh = BK.getmesh(coll)
    σs = _get_gauss_nodes(coll)
    udj = VectorOfArray([copy(uj[:, 1]) for _ in _delays])

    # range for locating time slices
    rg = UnitRange(1, m+1)
    eq = 1
    for j in 1:Ntst
        uj .= uc[:, rg]
        LA.mul!(gj, uj, L)
        LA.mul!(∂gj, uj, ∂L)

        # get the delayed states
        dτj = (mesh[j+1] - mesh[j]) / 2

        # compute the collocation residual
        if mesh[j]>=0
            for l in 1:m
                τ = BK.τj(σs[l], mesh, j)
                # udj = VectorOfArray([interp(mod(τ * period - d, period)) for d in _delays])
                for (ind, d) in enumerate(_delays)
                    udj.u[ind] .= interp(τ - d/period)
                    # udj.u[ind] .= BK.__interpolate_posolution(coll, τ - d/period, u, 1)
                end
                __po_coll_bc!(coll, outc[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dτj, outc[:, end])
            end
        end
        rg = rg .+ m
        eq += 1
    end
    # add the periodicity condition
    index = interp.initial_n
    # @. outc[:, end] = uc[:, end] - uc[:, end-index+1]
end

function jacobian_extended_mesh(coll::PeriodicOrbitOCollProblem,
                            periodic_sol::AbstractVector, 
                            pars)
    period = BK.getperiod(coll, periodic_sol, pars)
    interp = extended_sol(coll, periodic_sol, pars)
    coll_ext = interp.extended_coll
    extended_uc = interp.xc
    extended_u = vcat(vec(extended_uc), period)
    extended_outc = zero(extended_uc) .+ 0
    index = interp.initial_n
    @error "" index size(extended_uc)

    # return _residual_for_extended_mesh(coll_ext, interp, extended_u, pars)

    # _residual_for_extended_mesh!(coll_ext, interp, outc, uc, period, BK.get_Ls(coll.mesh_cache), pars, u)
    # return outc

    n, m, Ntst = size(coll)

    function residual(extended_u0)
        interp0 = ExtendedSolution(coll_ext, extended_u0[1:end-1]', interp.interval, interp.initial_n)
        _residual_for_extended_mesh(coll_ext, interp0, extended_u0, pars)[1:end-1-n][end-index+1:end]
    end

    J = ForwardDiff.jacobian(residual, extended_u)
    ncol = size(J,2)-size(J,1)
    B = J[end-index+3:end, end-index:end-1]
    A = J[3:index, 1:ncol+1]
    
    @error "" residual(extended_u) size(A) size(B)
    return J, A, B


    Mₜ = -B\A
    Nₜ, N = size(Mₜ)
    @error "" size(Mₜ)
    if N <= Nₜ
        @error "1" N <= Nₜ
        It = Itilde(N, Nₜ)
        M = It * Mₜ # same as Mₜ[end-N+1:end, :]
    else
        @error "2" N <= Nₜ
        It = Itilde(N - Nₜ, N)
        M = vcat(It, Mₜ)
    end

    vals = LA.eigvals(M)
    logvals = log.(complex.(vals))
    I = sortperm(logvals, by = real, rev = true)
    # floquet exponents
    σ = logvals[I] #.* 40.10727283620028
end