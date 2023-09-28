#include "solver/solver_ode.h"

Sparse QDACC::Numerics::ODESolver::calculate_propagator_single( System &s, size_t tensor_dim, double t0, double t_step, int i, int j, std::vector<QDACC::SaveState> &output, const Sparse &one ) {
    Log::L3( "[PathIntegral] Calculating Single Propagator at t = {} to t+dt = {} for i = {}, j = {}\n", t0, t0+s.parameters.t_step_pathint, i, j );
    Sparse projector = Sparse( tensor_dim, tensor_dim );
    projector.coeffRef( i, j ) = 1;
    Sparse M = iterate( projector, s, t0, t_step, output );//.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );

    Sparse map;
    if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
        map = QDACC::Matrix::sparse_projector( M );
    }
    for ( double tau = t_step; tau < s.parameters.t_step_pathint; tau += t_step ) {
        M = iterate( M, s, t0 + tau, t_step, output );//.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
    }
    if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
        return M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
    }
    M.makeCompressed();
    return M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
}