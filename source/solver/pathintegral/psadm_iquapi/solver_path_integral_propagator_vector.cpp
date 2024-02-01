#include "solver/solver_ode.h"

// TODO: maybe we can only evaluate the oupper triangle matrix, then write a getter function that checks if j>i, return mat[i,j].dagger(). Test with samples if (i,j) = (j,i).dagger()!!
std::vector<std::vector<MatrixMain>> &QDACC::Numerics::ODESolver::calculate_propagator_vector( System &s, size_t tensor_dim, double t0, double t_step, std::vector<QDACC::SaveState> &output ) {
    if ( pathint_propagator.contains( t0 ) ) {
        return pathint_propagator[t0];
    }
    MatrixMain one( tensor_dim, tensor_dim );
    one.setIdentity();
    std::vector<std::vector<MatrixMain>> ret( tensor_dim, { tensor_dim, MatrixMain( tensor_dim, tensor_dim ) } );
    // Calculate first by hand to ensure the Hamilton gets calculated correclty
    ret[0][0] = calculate_propagator_single( s, tensor_dim, t0, t_step, 0, 0, output, one ); // pathint_propagator[-1][0][0] );
// Calculate the remaining propagators
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_secondary_threads )
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) { // <= i
            if ( i == 0 && j == 0 ) continue;
            ret[i][j] = calculate_propagator_single( s, tensor_dim, t0, t_step, i, j, output, one ); // pathint_propagator[-1][i][j] );
        }
    }
    // Only save propagator vector if correlation functions are calculated.
    if ( cache.empty() ) {
        pathint_propagator[-1] = ret;
    } else {
        Log::L3( "[PathIntegral] Caching propagator vector for t = {}\n", t0 );
#pragma omp critical
        pathint_propagator[t0] = ret;
        return pathint_propagator[t0];
    }
    return pathint_propagator[-1];
}