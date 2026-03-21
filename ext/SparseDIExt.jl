module SparseDIExt

using SparseConnectivityTracer, SparseMatrixColorings
import DifferentiationInterface as DI
import DDEBifurcationKit as DDEBK
import BifurcationKit as BK

coll_residual_for_di!(result, u, coll, pars) = BK.residual!(coll, result, u, pars)

function DDEBK.BK._generate_jacobian(coll::BK.PeriodicOrbitOCollProblem{Tprob},
                    ::DDEBK.AutoSparseDI, 
                    orbitguess, 
                    pars; 
                    k...) where {Tprob <: DDEBK.AbstractDDEBifurcationProblem}
    backend = DI.AutoForwardDiff()
    sparse_forward_backend = DI.AutoSparse(
        backend;
        # we use the following instead of TracerSparsityDetector
        # because of searchsortedfirst in POSolution
        # sparsity_detector = TracerLocalSparsityDetector(),
        sparsity_detector = DI.DenseSparsityDetector(backend, atol = 1e-8),
        coloring_algorithm = GreedyColoringAlgorithm(),
    )
    out = copy(orbitguess)

    jac_prep_sparse_nonallocating = DI.prepare_jacobian(coll_residual_for_di!,
                                    out,
                                    sparse_forward_backend,
                                    orbitguess,
                                    DI.Constant(coll),
                                    DI.Constant(pars),
                                    # strict = Val(false)
                                    )

    jac_buffer = similar(sparsity_pattern(jac_prep_sparse_nonallocating), eltype(ones(length(out))))
    L1 = copy(jac_buffer)
    return (DDEBK.AutoSparseDI(),
            L1,
            out,
            jac_prep_sparse_nonallocating,
            sparse_forward_backend)
end

function BK.jacobian(coll::BK.PeriodicOrbitOCollProblem, 
            J::Tuple{DDEBK.AutoSparseDI, T1, T2, T3, T4}, 
            x, 
            p) where {T1, T2, T3, T4}
    L1 = J[2]
    out = J[3]
    jac_prep_sparse_nonallocating = J[4]
    sparse_forward_backend = J[5]
    DI.jacobian!(coll_residual_for_di!,
                out,
                L1,
                jac_prep_sparse_nonallocating,
                sparse_forward_backend,
                x,
                DI.Constant(p))
    return L1
end

end # module