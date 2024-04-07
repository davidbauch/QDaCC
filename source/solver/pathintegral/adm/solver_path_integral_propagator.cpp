#include "solver/solver_ode.h"

MatrixMain QDACC::Numerics::ODESolver::calculate_propagator_single( System &s, size_t tensor_dim, double t0, double t_step, int i, int j, std::vector<QDACC::SaveState> &output ) {
    Log::L3( "[PathIntegral] Calculating Single Propagator at t = {} to t+dt = {} for i = {}, j = {}\n", t0, t0 + s.parameters.t_step_pathint, i, j );
    MatrixMain projector = MatrixMain( tensor_dim, tensor_dim );
    projector.setZero();
    projector.coeffRef( i, j ) = 1;
    MatrixMain M = iterate( projector, s, t0, t_step, output );

    MatrixMain map( tensor_dim, tensor_dim);
    if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
        #ifdef USE_SPARSE_MATRIX
        map = QDACC::Matrix::sparse_projector( M );
        #else
        map = QDACC::Matrix::dense_projector( M );
        #endif
    }
    for ( double tau = t_step; tau < s.parameters.t_step_pathint; tau += t_step ) {
        // Make sure to exactly hit the end time
        if ( tau > s.parameters.t_step_pathint ) {
            t_step = s.parameters.t_step_pathint - ( tau - t_step );
        }
        M = iterate( M, s, t0 + tau, t_step, output ); //;
    }
    if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
        return M.cwiseProduct( map );
    }
    return M;
}