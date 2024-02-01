#include "solver/solver_ode.h"

/**
 * REQUIRES: tpm_purpose in output cache
*/

bool QDACC::Numerics::ODESolver::calculate_timebin_generator_expectation_values( System &s, const std::string &purpose ) {
    Log::L2( "[TimeBin] (Z)XZ Generator Expectation values for {}\n", purpose );

    // Progress
    auto progressbar = ProgressBar();

    // This is most certainly never the case, but we may as well make sure.
    if ( not to_output_m["TwoPMat"].contains( purpose ) ) {
        Log::Warning( "[TimeBin] Required TPM is not in cache!\n" );
        return false;
    }

    // Cache Progressbar Size
    auto cache_size = to_output_m["TwoPMat"][purpose].size();
    auto dim = to_output_m["TwoPMat"][purpose].front().rows();
    auto order = ( dim == 4 ) ? 2 : 3; // = log_2(dim)
    int pbsize = 2 * cache_size;
    Timer &timer_c = Timers::create( "Generator Expectation Values (" + purpose + ")" ).start();

    // Maximum Time / Dimension of the cache matrices
    auto maximum_time = std::min<size_t>( cache_size, savedStates.size() );

    // Generate Referenced Shortcuts
    auto &output = add_to_output( "GeneratorExpectationG", purpose, std::vector<Scalar>( maximum_time, 0 ), to_output );
    // Get the TPM References
    auto &twophotonmatrix = to_output_m["TwoPMat"][purpose];

#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t k = 0; k < maximum_time; k++ ) {
        Scalar upper, lower;
        // Neumann-Iterated TPDM
        const auto rho = twophotonmatrix[k];
        // Iterate upper and lower depending on the order (G2 or G3 TPM)
        if ( order == 2 ) {
            upper = 2.0 * ( rho( 1, 0 ) - rho( 3, 1 ) );
            lower = rho( 2, 1 ) + rho( 1, 2 ) + rho( 0, 0 ) + rho( 3, 3 );
        } else {
            upper = 2.0 * ( rho( 0, 2 ) - rho( 1, 6 ) - rho( 3, 4 ) + rho( 5, 7 ) );
            lower = rho( 0, 0 ) + rho( 7, 7 ) + rho( 1, 3 ) + rho( 3, 1 ) + rho( 2, 2 ) + rho( 5, 5 ) + rho( 6, 4 ) + rho( 4, 6 );
        }
        // Output
        output[k] = { std::real( std::real( upper ) / lower ), std::real( std::imag( upper ) / lower ) };
        // Progress
        timer_c.iterate();
        Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Generator Expectation Values (" + purpose + "): " );
    }

    Log::L1( "[TimeBin] (Z)XZ Generator Expectation value: {}.\n", std::real( output.back() ) );
    return true;
}