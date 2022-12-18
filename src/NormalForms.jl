function BK.hopfNormalForm(prob::ConstantDDEBifProblem, pt::BK.Hopf, ls; verbose::Bool = false)
	x0 = pt.x0
	p = pt.p
	lens = pt.lens
	parbif = set(pt.params, lens, p)
	ω = pt.ω
	λ0 = Complex(0, ω)
	δ = BK.getDelta(prob)

	# jacobian at the bifurcation point
	# c'est recalcule ici!!!! 2x
	L = BK.jacobian(prob, x0, parbif)
	Δ0 = Δ(L, 0λ0)
	Δ2ω = Δ(L, 2λ0)

	ζ = pt.ζ
	cζ = conj.(pt.ζ)
	ζstar = copy(pt.ζstar)
	ζstar ./= conj(dot(ζstar, Δ(Val(:der), L, ζ, λ0)))
	if ~isapprox(dot(ζstar, Δ(Val(:der), L, ζ, λ0)), 1; rtol = 1e-2)
		@warn "We found instead $(dot(ζstar, Δ(Val(:der), L, ζ, λ0)))"
	end

	x0c = VectorOfArray([copy(x0) for _ in 1:length(prob.delays0)+1])

	ζθ = expθ(L, ζ, λ0)
	ζθc = conj.(ζθ)

	# we use BilinearMap to be able to call on complex valued arrays
	R2 = BK.BilinearMap( (dx1, dx2)      -> BK.d2F(prob, x0c, parbif, dx1, dx2) ./2)
	R3 = BK.TrilinearMap((dx1, dx2, dx3) -> BK.d3F(prob, x0c, parbif, dx1, dx2, dx3) ./6 )

	# −LΨ001 = R01
	R01 = (BK.residual(prob, x0, set(parbif, lens, p + δ)) .- BK.residual(prob, x0, set(parbif, lens, p - δ))) ./ (2δ)
	Ψ001 = Complex.(expθ(L, ls(Δ0, -R01)[1], 0))

	# (2iω−L)Ψ200 = R20(ζ,ζ)
	R20 = R2(ζθ, ζθ)
	Ψ200 = expθ(L, ls(Δ2ω, R20)[1], 2λ0)
	# @assert Ψ200 ≈ (Complex(0, 2ω)*I - L) \ R20

	# −LΨ110 = 2R20(ζ,cζ).
	R20 = 2 .* R2(ζθ, ζθc)
	Ψ110 = Complex.(expθ(L, ls(Δ0, -R20)[1], 0))

	# a = ⟨R11(ζ) + 2R20(ζ,Ψ001),ζ∗⟩
	# av = (apply(jacobian(prob, x0, set(parbif, lens, p + δ)), ζ) .- apply(jacobian(prob, x0, set(parbif, lens, p - δ)), ζ)) ./ (2δ)
	av = zero(ζ)
	av .+= 2 .* R2(ζθ, Ψ001)
	a = dot(ζstar, av)
	a = 0

	# b = ⟨2R20(ζ,Ψ110) + 2R20(cζ,Ψ200) + 3R30(ζ,ζ,cζ), ζ∗⟩)
	bv = 2 .* R2(ζθ, Ψ110) .+ 2 .* R2(ζθc, Ψ200) .+ 3 .* R3(ζθ, ζθ, ζθc)
	b = dot(ζstar, bv)

	# @info "info" b real(b)/ω/2 parbif δ #Ψ110 Ψ200 2λ0

	# return coefficients of the normal form
	verbose && println((a = a, b = b))
	pt.nf = (a = a, b = b)
	if real(a) * real(b) < 0
		pt.type = :SuperCritical
	elseif real(a) * real(b) > 0
		pt.type = :SubCritical
	else
		pt.type = :Singular
	end
	verbose && printstyled(color = :red,"--> Hopf bifurcation point is: ", pt.type, "\n")
	return pt
end

function BK.hopfNormalForm(prob::ConstantDDEBifProblem,
					br::BK.AbstractBranchResult, ind_hopf::Int;
					nev = length(BK.eigenvalsfrombif(br, id_bif)),
					verbose::Bool = false,
					lens = BK.getLens(br),
					Teigvec = BK.getvectortype(br),
					scaleζ = norm)
	@assert br.specialpoint[ind_hopf].type == :hopf "The provided index does not refer to a Hopf Point"
	verbose && println("#"^53*"\n--> Hopf Normal form computation")

	options = br.contparams.newtonOptions

	# bifurcation point
	bifpt = br.specialpoint[ind_hopf]
	eigRes = br.eig

	# eigenvalue
	λ = eigRes[bifpt.idx].eigenvals[bifpt.ind_ev]
	ω = imag(λ)

	# parameter for vector field
	p = bifpt.param
	parbif = set(getParams(br), lens, p)
	L = BK.jacobian(prob, convert(Teigvec, bifpt.x), parbif)

	# right eigenvector
	if BK.haseigenvector(br) == false
		# we recompute the eigen-elements if there were not saved during the computation of the branch
		_λ, _ev, _ = options.eigsolver(L, bifpt.ind_ev + 2)
		@assert _λ[bifpt.ind_ev] ≈ λ "We did not find the correct eigenvalue $λ. We found $(_λ)"
		ζ = geteigenvector(options.eigsolver, _ev, bifpt.ind_ev)
	else
		ζ = copy(geteigenvector(options.eigsolver ,br.eig[bifpt.idx].eigenvecs, bifpt.ind_ev))
	end
	ζ ./= scaleζ(ζ)

	# left eigen-elements
	_Jt = BK.hasAdjoint(prob) ? jad(prob, convert(Teigvec, bifpt.x), parbif) : adjoint(L)
	ζstar, λstar = BK.getAdjointBasis(_Jt, conj(λ), options.eigsolver; nev = nev, verbose = verbose)

	# check that λstar ≈ conj(λ)
	abs(λ + λstar) > 1e-2 && @warn "We did not find the left eigenvalue for the Hopf point to be very close to the imaginary part:\nλ ≈ $λ,\nλstar ≈ $λstar?\n You can perhaps increase the number of computed eigenvalues, the number is nev = $nev"

	# normalise left eigenvector
	ζstar ./= dot(ζ, ζstar)
	@assert dot(ζ, ζstar) ≈ 1

	hopfpt = BK.Hopf(bifpt.x, bifpt.param,
		ω,
		parbif, lens,
		ζ, ζstar,
		(a = zero(Complex{eltype(bifpt.x)}), b = zero(Complex{eltype(bifpt.x)}) ),
		:SuperCritical
	)
	return BK.hopfNormalForm(prob, hopfpt, options.linsolver ; verbose = verbose)
end