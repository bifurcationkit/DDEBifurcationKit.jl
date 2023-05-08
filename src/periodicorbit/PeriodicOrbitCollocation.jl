@views function (prob::PeriodicOrbitOCollProblem{Tprob})(u::AbstractVector, pars) where {Tprob <: AbstractDDEBifurcationProblem}
	uc = BK.getTimeSlices(prob, u)
	T = BK.getPeriod(prob, u, nothing)
	result = zero(u)
	resultc = BK.getTimeSlices(prob, result)
	functionalColl!(prob, resultc, uc, T, BK.getLs(prob.mesh_cache), pars, u)
	# add the phase condition
	result[end] = BK.phaseCondition(prob, (u, uc), BK.getLs(prob.mesh_cache))
	return result
end


function _POOCollScheme!(pb::PeriodicOrbitOCollProblem, dest, ∂u, u, ud, par, h, tmp)
	tmp2 = pb.prob_vf.VF.F(u, ud, par)
	dest .= @. ∂u - h * tmp2
end

# function for collocation problem
@views function functionalColl!(pb::PeriodicOrbitOCollProblem{Tprob}, out, u, period, (L, ∂L), pars, result) where {Tprob <: ConstantDDEBifProblem}
	Ty = eltype(u)
	n, ntimes = size(u)
	m = pb.mesh_cache.degree
	Ntst = pb.mesh_cache.Ntst
	# we want slices at fixed  times, hence gj[:, j] is the fastest
	# temporaries to reduce allocations
	gj  = zeros(Ty, n, m)
	∂gj = zeros(Ty, n, m)
	uj  = zeros(Ty, n, m+1)

	# get interpolator which allows to get result(t)
	interp = BK.POSolution(pb, result)
	delays = pb.prob_vf.delays(pars)

	# get the mesh of the OCollProblem
	mesh = BK.getMesh(pb)
	σ = LinRange(0,2,m)

	# range for locating time slices
	rg = UnitRange(1, m+1)
	for j in 1:Ntst
		uj .= u[:, rg]
		mul!(gj, uj, L')
		mul!(∂gj, uj, ∂L')

		# get the delayed states
		tj = mesh[j]
		dtj = (mesh[j+1]-mesh[j]) / 2

		# compute the collocation residual
		for l in 1:m
			tσ = tj + dtj * σ[l]
			udj = [interp(mod(tσ*period - d, period)) for d in delays]
			# out[:, end] can serve as buffer for now in the following function
			_POOCollScheme!(pb, out[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dtj, out[:, end])

		end
		# carefull here https://discourse.julialang.org/t/is-this-a-bug-scalar-ranges-with-the-parser/70670/4"
		rg = rg .+ m
	end
	# add the periodicity condition
	out[:, end] .= u[:, end] .- u[:, 1]
end


# function for collocation problem
@views function functionalColl!(pb::PeriodicOrbitOCollProblem{Tprob}, out, u, period, (L, ∂L), pars, result) where {Tprob <: SDDDEBifProblem}
	Ty = eltype(u)
	n, ntimes = size(u)
	m = pb.mesh_cache.degree
	Ntst = pb.mesh_cache.Ntst
	# we want slices at fixed  times, hence gj[:, j] is the fastest
	# temporaries to reduce allocations
	gj  = zeros(Ty, n, m)
	∂gj = zeros(Ty, n, m)
	uj  = zeros(Ty, n, m+1)

	# get interpolator which allows to get result(t)
	interp = BK.POSolution(pb, result)

	if period <= 0
		out .= 1e9
		return out
	end

	# get the mesh of the OCollProblem
	mesh = BK.getMesh(pb)
	σ = LinRange(0,2,m)

	# range for locating time slices
	rg = UnitRange(1, m+1)
	for j in 1:Ntst
		uj .= u[:, rg]
		mul!(gj, uj, L')
		mul!(∂gj, uj, ∂L')

		# get the delayed states
		tj = mesh[j]
		dtj = (mesh[j+1]-mesh[j]) / 2

		# compute the collocation residual
		for l in 1:m
			tσ = tj + dtj * σ[l]
			delays = pb.prob_vf.delays(gj[:, l], pars)
			udj = [interp(mod(tσ*period - d, period)) for d in delays]
			# out[:, end] can serve as buffer for now in the following function
			_POOCollScheme!(pb, out[:, rg[l]], ∂gj[:, l], gj[:, l], udj, pars, period * dtj, out[:, end])

		end
		# carefull here https://discourse.julialang.org/t/is-this-a-bug-scalar-ranges-with-the-parser/70670/4"
		rg = rg .+ m
	end
	# add the periodicity condition
	out[:, end] .= u[:, end] .- u[:, 1]
end
