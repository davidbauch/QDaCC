#include "solver/solver_ode.h"

void QDACC::Numerics::ODESolver::calculate_g2( System &s, const std::string &s_op_i, const std::vector<std::string> &s_op_j, const std::vector<std::string> &s_op_k, const std::string &s_op_l, const std::vector<std::string> &purposes ) {
    Log::L2( "[CorrelationFunction] Preparing to calculate Correlation function\n" );
    Log::L2( "[CorrelationFunction] Generating Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_l = get_operators_matrix( s, s_op_l );
    Log::L2( "[CorrelationFunction] Iterated Operators are op_i = {} and op_l = {}\n", s_op_i, s_op_l );

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Generator super-purpose
    std::string super_purpose = std::accumulate( std::next( purposes.begin() ), purposes.end(), purposes.front(), []( const std::string &a, const std::string &b ) { return a + " and " + b; } );

    std::vector<std::pair<MatrixMain, std::string>> eval_operators;

    // Preconstruct
    for ( auto current = 0; current < s_op_k.size(); current++ ) {
        const auto &purpose = purposes[current];
        // Cancel if Purpose already exists
        if ( cache.contains( purpose ) ) {
            Log::L2( "[CorrelationFunction] Matrix for {} already exists! Skipping!\n", purpose );
            continue;
        }
        const auto &op_j = get_operators_matrix( s, s_op_j[current] );
        const auto &op_k = get_operators_matrix( s, s_op_k[current] );
        // Construct Evaluation Operators
        eval_operators.emplace_back( op_j * op_k, purpose );

        Log::L2( "[CorrelationFunction] Preparing Cache Matrices for {}. Using Eval Operators op_j = {}, op_k = {}\n", purpose, s_op_j[current], s_op_k[current] );
        cache[purpose] = MultidimensionalCacheMatrix( { matdim, matdim }, purpose );
        auto &mat = cache[purpose];
    }

    // Return if no Functions need to be evaluated
    if ( eval_operators.empty() )
        return;

    // Create Timer and Progresbar
    Timer &timer = Timers::create( "Correlation-Loop (" + super_purpose + ")" ).start();
    auto progressbar = ProgressBar();

    // Calculate G2 Function
    Log::L2( "[CorrelationFunction] Calculating G(tau)... purpose: {}, saving to matrix of size {}x{},  iterating over {} saved states...\n", super_purpose, matdim, matdim, std::min<size_t>( matdim, savedStates.size() ) );
    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t t = 0; t < std::min<size_t>( matdim, savedStates.size() ); t++ ) {
        const auto& current_state = savedStates.at( t );
        // Create and reserve past rho's vector
        std::vector<QDACC::SaveState> savedRhos; 
        // Get Time from saved State
        double t_t = current_state.t; 
        // Calculate New Modified Density Matrix
        MatrixMain rho_tau = s.dgl_calc_rhotau( current_state.mat, op_l, op_i, t_t );
        // Calculate Runge Kutta or PI
        if ( s.parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
            calculate_path_integral_correlation( rho_tau, t_t, s.parameters.t_end, timer, s, savedRhos, op_l, op_i, 1 /*ADM Cores*/ );
        } else {
            calculate_runge_kutta( rho_tau, t_t, s.parameters.t_end, timer, progressbar, super_purpose, s, savedRhos, false /*Output*/ );
        }
        // Interpolate saved states to equidistant timestep
        savedRhos = Numerics::interpolate_curve( savedRhos, t_t, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        for ( const auto &[eval, purpose] : eval_operators ) {
            auto &gmat = cache[purpose];
            for ( size_t tau = 0; tau < savedRhos.size(); tau++ ) {
                const double t_tau = savedRhos.at( tau ).t;
                gmat.get( t, tau ) = s.dgl_expectationvalue<MatrixMain>( savedRhos.at( tau ).mat, eval, t_tau );
            }
        }
        Timers::outputProgress( timer, progressbar, t, savedStates.size(), super_purpose );
    }

    timer.end();
    Timers::outputProgress( timer, progressbar, savedStates.size(), savedStates.size(), super_purpose, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[CorrelationFunction] G ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", super_purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );

    // Manually Apply the detector function
    for ( const auto &[eval, purpose] : eval_operators ) {
        auto &gmat = cache[purpose];
        if ( gmat.hasBeenFourierTransformed() ) {
            Log::L2( "[CorrelationFunction] Detector Function already applied to {}. Skipping!\n", purpose );
            continue;
        }
        apply_detector_function( s, gmat, purpose );
    }
}