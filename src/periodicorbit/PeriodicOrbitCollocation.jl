_get_gauss_nodes(coll) = coll.mesh_cache.gauss_nodes

@views function BK.po_residual!(coll::Collocation{Tprob},
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

function __po_coll_bc!(coll::Collocation, dest, ∂u, u, ud, par, h, tmp)
    tmp .= coll.prob_vf.VF.F(u, ud, par)
    @. dest = ∂u - h * tmp
end

function (sol::BK.POInterpolation{ <: Collocation})(::Val{:der}, t0)
    ForwardDiff.derivative(sol, t0)
end

# function for collocation problem
@views function functional_coll!(coll::Collocation{Tprob},
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
    interp = BK.POInterpolation(coll, u)
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

function BK._generate_jacobian(coll::BK.Collocation{ <: ConstantDDEBifProblem}, 
                                J::BK.DenseAnalyticalInplace,
                                orbitguess,
                                pars; 
                                Jcoll_matrix = nothing,
                                k...)
    _Jcoll_matrix = isnothing(Jcoll_matrix) ? analytical_jacobian_dde_cst(coll, orbitguess, pars) : Jcoll_matrix
    return (BK.DenseAnalyticalInplace(), _Jcoll_matrix)
end

function BK._generate_jacobian(coll::BK.Collocation{ <: SDDDEBifProblem}, 
                                J::BK.DenseAnalyticalInplace,
                                orbitguess,
                                pars; 
                                Jcoll_matrix = nothing,
                                k...)
    _Jcoll_matrix = isnothing(Jcoll_matrix) ? analytical_jacobian_dde(coll, orbitguess, pars) : Jcoll_matrix
    return (BK.DenseAnalyticalInplace(), _Jcoll_matrix)
end

function BK._jacobian_po(wrap::BK.PeriodicOrbitFunctionalColl{ <: BK.Collocation{ <: ConstantDDEBifProblem}}, 
                ::BK.DenseAnalytical,
                x, 
                p)
    coll = BK.get_discretization(wrap)
    return analytical_jacobian_dde_cst(coll, x, p)
end

function BK._jacobian_po(wrap::BK.PeriodicOrbitFunctionalColl{ <: BK.Collocation{ <: SDDDEBifProblem}}, 
                ::BK.DenseAnalytical,
                x, 
                p)
    coll = BK.get_discretization(wrap)
    return analytical_jacobian_dde(coll, x, p)
end

function BK._jacobian_po(wrap::BK.PeriodicOrbitFunctionalColl{ <: BK.Collocation{ <: ConstantDDEBifProblem}}, 
                J::Tuple{BK.DenseAnalyticalInplace, Tj},
                x, 
                p) where {Tj}
    _Jcoll_matrix = J[2]
    coll = BK.get_discretization(wrap)
    return analytical_jacobian_dde_cst(coll, x, p; Jcoll = _Jcoll_matrix)
end

function BK._jacobian_po(wrap::BK.PeriodicOrbitFunctionalColl{ <: BK.Collocation{ <: SDDDEBifProblem}}, 
                J::Tuple{BK.DenseAnalyticalInplace, Tj},
                x, 
                p) where {Tj}
    _Jcoll_matrix = J[2]
    coll = BK.get_discretization(wrap)
    return analytical_jacobian_dde(coll, x, p; Jcoll = _Jcoll_matrix)
end

function BK._jacobian_po(wrap::BK.PeriodicOrbitFunctionalColl{ <: BK.Collocation{ <: ConstantDDEBifProblem}}, 
                ::BK.FullSparse,
                x, 
                p)
    coll = BK.get_discretization(wrap)
    return analytical_jacobian_dde_cst(coll, x, p)
end

function BK._jacobian_po(wrap::BK.PeriodicOrbitFunctionalColl{ <: BK.Collocation{ <: SDDDEBifProblem}}, 
                ::BK.FullSparse,
                x, 
                p)
    coll = BK.get_discretization(wrap)
    return analytical_jacobian_dde(coll, x, p)
end

"""
using DifferentiationInterface to automatically derive the sparse jacobian.
"""
struct AutoSparseDI <: BK.AbstractJacobianSparseMatrix end

function BK._generate_jacobian(coll::Collocation{Tprob}, ::AutoSparseDI, orbitguess, pars; k...) where {Tprob <: AbstractDDEBifurcationProblem}
    error("You need to import `DifferentiationInterface, SparseConnectivityTracer, SparseMatrixColorings` in order to use this jacobian")
end
########################################################################################
# analytical jacobians for constant DDE
for (fname, floquet) in ((:analytical_jacobian_dde_cst, false), 
                         (:analytical_jacobian_dde_cst_floquetgev, true),
                         (:analytical_jacobian_dde_floquet, true),
                         )
    @eval begin
    @views function $fname(coll::Collocation{Tprob}, 
                            u::AbstractVector{𝒯}, 
                            pars;
                            ρD = one(𝒯),
                            ρF = one(𝒯),
                            ρI = zero(𝒯),
                            Jcoll = nothing) where {Tprob <: AbstractDDEBifurcationProblem, 𝒯}
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
        interp = BK.POInterpolation(coll, u)
        VF = coll.prob_vf

        # loop over the mesh intervals
        rg = UnitRange(1, m+1)
        rgNx = UnitRange(1, n)
        rgNy = UnitRange(1, n)

        delays_v = delays(VF, u[1:n], pars) # vector of delays
        udj = VectorOfArray([zeros(𝒯, n) for d in delays_v])
        if coll.jacobian isa BK.FullSparse
            J  = SA.spzeros(𝒯, length(coll) + 1, length(coll) + 1)
            J0 = SA.spzeros(𝒯, n, n)
            In = SA.sparse(In)
        else
            # careful, the sparsity pattern changes so we need to zero the output array Jcoll
            J = isnothing(Jcoll) ? zeros(𝒯, length(coll) + 1, length(coll) + 1) : Jcoll
            if isnothing(Jcoll) == false 
                J .= 0
            end
            J0 = zeros(𝒯, n, n)
        end

        # put boundary condition
        J[nJ-n:nJ-1, nJ-n:nJ-1] .= In
        J[nJ-n:nJ-1, 1:n] .= (-1) .* In

        if $(fname == :analytical_jacobian_dde_cst_floquetgev)
            # arrays to store the jacobian of the delayed terms
            Jd = [zeros(𝒯, length(coll)+1, length(coll)+1) for _ in 1:length(delays_v)]
        elseif $(fname == :analytical_jacobian_dde_floquet)
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
                        # the delayed time may fall exactly on t = 0 (τd = 0) which is
                        # the left boundary of the first interval
                        index_t = max(searchsortedfirst(mesh, τd) - 1, 1)
                        @assert index_t <= Ntst "We have index_t = $index_t, which is out of bounds for mesh of size $(length(mesh)) and τd = $τd. Please open an issue on the website of BifurcationKit.jl"

                        rgNy_delay = UnitRange(1, n) .+ ((m * n) * (index_t - 1))
                        σ = BK.σj(τd, mesh, index_t)
                        β = BK.lagrange(l2, σ, BK.get_mesh_coll(coll)) * ρF

                        if $(fname == :analytical_jacobian_dde_cst_floquetgev)
                            Jd[idelay][_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[idelay] .* β
                        elseif ($(fname == :analytical_jacobian_dde_floquet) && t0 < 0)
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
########################################################################################
# Floquet collocation Jacobian for state-dependent delays.
# The linearization of an SD-DDE around a T-periodic orbit x* reads
#   ẏ(t) = [A₀(t) - Σⱼ Aⱼ(t)·ẋ*(t-τⱼ*(t))·cⱼ(t)]·y(t) + Σⱼ Aⱼ(t)·y(t-τⱼ*(t))
# with τⱼ*(t) = τⱼ(x*(t), p) and cⱼ(t) = ∇τⱼ(x*(t)). Here the delays are evaluated at
# each collocation point and the delay-derivative (state-dependence) correction
# -Σⱼ (Aⱼ·ẋ*(t-τⱼ))·cⱼ is added to the undelayed coefficient. The returned JacobianDDE
# has the same J0 / Jd split as analytical_jacobian_dde_floquet for constant delays
# (Jd holds the delayed contributions with t0 < 0), so the Verheyden-Lust monodromy applies.
@views function analytical_jacobian_dde_floquet(coll::Collocation{Tprob},
                                                u::AbstractVector{𝒯},
                                                pars;
                                                Jcoll = nothing) where {Tprob <: SDDDEBifProblem, 𝒯}
    n, m, Ntst = size(coll)
    nJ = length(coll) + 1
    L, ∂L = BK.get_Ls(coll.mesh_cache)
    mesh = BK.getmesh(coll)
    σs = _get_gauss_nodes(coll)
    ω = coll.mesh_cache.gauss_weight
    period = BK.getperiod(coll, u, nothing)
    phase = zero(𝒯)
    uc = BK.get_time_slices(coll, u)
    pj = BK.get_tmp(coll.cache.gi, u)
    In = coll.cache.In
    interp = BK.POInterpolation(coll, u)
    VF = coll.prob_vf

    delays_v = delays(VF, u[1:n], pars) # reference delays (at t = 0)

    J  = zeros(𝒯, nJ, nJ)
    J0 = zeros(𝒯, n, n)
    # history part (delayed terms reaching into the previous period), single matrix
    Jd = zeros(𝒯, nJ, nJ)

    # periodic boundary condition
    J[nJ-n:nJ-1, nJ-n:nJ-1] .= In
    J[nJ-n:nJ-1, 1:n] .= (-1) .* In

    rg = UnitRange(1, m+1)
    rgNx = UnitRange(1, n)
    rgNy = UnitRange(1, n)
    for j in 1:Ntst
        LA.mul!(pj, uc[:, rg], L)
        dτj = (mesh[j+1] - mesh[j]) / 2
        α = period * dτj
        for l in 1:m
            _rgX = rgNx .+ (l-1)*n
            τ = BK.τj(σs[l], mesh, j)
            # state-dependent delays evaluated at the collocation point
            delays_l = delays(VF, pj[:, l], pars)
            udj = VectorOfArray([interp(mod(τ * period - d, period)) for d in delays_l])
            JacDDE = jacobian(VF, pj[:, l], udj, pars)
            J0 .= JacDDE.J0
            # delay-derivative correction: B₀ = J₀ - Σⱼ (Aⱼ·ẋ*(t-τⱼ))·cⱼ
            for (ind, d) in enumerate(delays_l)
                cj = ForwardDiff.gradient(z -> delays(VF, z, pars)[ind], pj[:, l])
                ẋd = interp(Val(:der), τ * period - d)
                J0 .-= (JacDDE.Jd[ind] * ẋd) * cj'
            end
            for l2 in 1:m+1
                J[_rgX, rgNy .+ (l2-1)*n] .+= @. (-α * L[l2, l]) * J0 + (∂L[l2, l]) * In
                for (ind, d) in enumerate(delays_l)
                    # find interval where t-τ/period belongs
                    t0 = τ * period - d
                    τd = mod(t0, period) / period
                    index_t = clamp(searchsortedfirst(mesh, τd) - 1, 1, Ntst)
                    rgNy_delay = UnitRange(1, n) .+ ((m * n) * (index_t - 1))
                    σ = BK.σj(τd, mesh, index_t)
                    β = BK.lagrange(l2, σ, BK.get_mesh_coll(coll))
                    if t0 < 0
                        Jd[_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[ind] .* β
                    else
                        J[_rgX, rgNy_delay .+ (l2-1)*n] .+= -α .* JacDDE.Jd[ind] .* β
                    end
                end
            end
            # derivative w.r.t. the period
            J[_rgX, nJ] .= VF.VF.F(pj[:, l], udj, pars) .* (-dτj)
            for (ind, d) in enumerate(delays_l)
                J[_rgX, nJ] .+= -(α * d / period) .* (JacDDE.Jd[ind] * interp(Val(:der), τ * period - d))
            end
            phase += LA.dot(pj[:, l], coll.∂ϕ[:, (j-1)*m + l]) * ω[l]
        end
        rg = rg .+ m
        rgNx = rgNx .+ (m * n)
        rgNy = rgNy .+ (m * n)
    end
    return JacobianDDE(missing, missing, J, Jd, delays_v)
end

# full analytical jacobian for the DDE collocation (constant or state-dependent delays).
# It is the Floquet jacobian (J0 + Jd) completed with the phase condition row, so that
# it can be used as a Newton jacobian. The phase condition code stays in the generated
# constant-delay jacobian (analytical_jacobian_dde_cst); it is added here for the Floquet
# version (whose monodromy extraction ignores the last n+1 rows anyway).
@views function analytical_jacobian_dde(coll::Collocation{ <: AbstractDDEBifurcationProblem},
                                        u::AbstractVector,
                                        pars;
                                        Jcoll = nothing)
    Jdde = analytical_jacobian_dde_floquet(coll, u, pars; Jcoll)
    M = Jdde.J0 + Jdde.Jd
    nJ = size(M, 1)
    period = BK.getperiod(coll, u, nothing)
    Ls = BK.get_Ls(coll.mesh_cache)
    phase = BK.phase_condition(coll, BK.get_time_slices(coll, u), Ls, period)
    M[end, begin:end-1] .= coll.cache.∇phase ./ period
    M[nJ, nJ] = -phase / period^2
    return M
end
########################################################################################
# BK.jacobian on a DDE collocation: use the analytical jacobian when the collocation is
# configured with DenseAnalytical, otherwise fall back to automatic differentiation.
function BK.jacobian(coll::BK.Collocation{<:AbstractDDEBifurcationProblem}, u, par)
    return _jacobian_dde_coll(coll, coll.jacobian, u, par)
end

_jacobian_dde_coll(coll, ::BK.DenseAnalytical, u, par) = analytical_jacobian_dde(coll, u, par)
_jacobian_dde_coll(coll, ::BK.AutoDiffDense, u, par) = ForwardDiff.jacobian(z -> BK.po_residual(coll, z, par), u)
