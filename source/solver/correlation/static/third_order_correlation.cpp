#include "solver/solver_ode.h"

// G3 will NOT be cached like G1 and G2, because it is too memory intensive. Additionally,
// it is not required, because no other function depends on it / needs G3s (yet...).
void QDACC::Numerics::ODESolver::calculate_g3( System &s, const std::string &s_op_i, const std::vector<std::string> &s_op_j, const std::vector<std::string> &s_op_k, const std::vector<std::string> &s_op_l, const std::vector<std::string> &s_op_m, const std::string &s_op_n, const std::vector<std::string> &purposes ) {
    Log::L2( "[CorrelationFunction] Preparing to calculate Correlation function\n" );
    Log::L2( "[CorrelationFunction] Generating Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_n = get_operators_matrix( s, s_op_n );
    //Log::L2( "[CorrelationFunction] Iterated Operators are op_i = {} and op_l = {}\n", s_op_i, s_op_l );

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Generator super-purpose
    std::string super_purpose = std::accumulate( std::next( purposes.begin() ), purposes.end(), purposes.front(), []( const std::string &a, const std::string &b ) { return a + " and " + b; } );

    std::vector<std::tuple<MatrixMain, std::vector<MatrixMain>, std::vector<std::string>>> eval_operators;

    // Preconstruct
    // We dont calculate a time matrix here, because we store the time as a complex number,
    // which we cannot do for nD with n > 2
    // In the future, we should change the time tracking to the following:
    // save t_i_min and t_i_max for each i, and then calculate the time on the fly
    // I will try this in the output of G3 and then move it to G1 and G2 too.

    // Looping works as follows:
    // We pass fixed operators like in G1/2, which have the index i and n.
    // We pass the outermost operators, which have the index j and m.
    // We pass the innermost operators, which have the index k and l.
    // We then loop over the outermost operators, and for each outermost operator,
    // we loop over the innermost operators.
    // This way, we can calculate the G3 function for each combination of outermost
    // and innermost operators.
    // The purpose vector is REQUIRED to pass purposes in this order too.
    // This means, purposes for the input i, [j1,j2], [k1,k2], [l1,l2], [m1,m2], n should be as follows:
    // P(i,j1,k1,l1,m1,n), P(i,j1,k2,l2,m1), P(i,j2,k1,l1,m2), P(i,j2,k2,l2,m2)
    for ( auto outer = 0; outer < s_op_j.size(); outer++ ) {
        const auto &op_j = get_operators_matrix( s, s_op_j[outer] );
        const auto &op_m = get_operators_matrix( s, s_op_m[outer] );
        std::vector<MatrixMain> inners;
        std::vector<std::string> inner_purposes;
        for ( auto inner = 0; inner < s_op_k.size(); inner++ ) {
            auto current = outer * s_op_j.size() + inner;
            const auto &purpose = purposes[current];
            // Cancel if Purpose already exists. This should never be true for the current implmentation.
            if ( cache.contains( purpose ) ) {
                Log::L2( "[CorrelationFunction] Matrix for {} already exists! Skipping!\n", purpose );
                continue;
            }
            Log::L2( "[CorrelationFunction] Preparing Cache Matrices for {}.", purpose );
            const auto &op_k = get_operators_matrix( s, s_op_k[current] );
            const auto &op_l = get_operators_matrix( s, s_op_l[current] );
            inners.emplace_back( op_k * op_l );
            inner_purposes.emplace_back( purpose );
            // Construct Evaluation Operators
            cache[purpose] = MultidimensionalCacheMatrix( { matdim, matdim, matdim }, purpose );
            // Preconstruct t,tau time. 
            auto& mat = cache[purpose];
        }
        eval_operators.emplace_back( op_j * op_m, inners, inner_purposes );
    }

    // Create Timer and Progresbar
    Timer &timer = Timers::create( "Correlation-Loop (" + super_purpose + ")" ).start();
    auto progressbar = ProgressBar();

    // Calculate G3 Function
    Log::L2( "[CorrelationFunction] Calculating G3(tau)... purpose: {}, saving to matrix of size {}x{}x{}, iterating over {} saved states...\n", super_purpose, matdim, matdim, matdim, std::min<size_t>( matdim, savedStates.size() ) );
    // Main G3 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t i = 0; i < std::min<size_t>( matdim, savedStates.size() ); i++ ) {
        // Create and reserve past rho's vector.
        std::vector<QDACC::SaveState> savedRhos_tau1, savedRhos_tau2;
        // Get state corresponding to i
        const auto& current_state = savedStates.at( i );
        // Get Time from saved State
        double t_t = current_state.t;
        // Calculate New Modified Density Matrix
        MatrixMain rho_tau = s.dgl_calc_rhotau( current_state.mat, op_n, op_i, t_t );
        // Calculate Runge Kutta or PI
        if ( s.parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
            calculate_path_integral_correlation( rho_tau, t_t, s.parameters.t_end, timer, s, savedRhos_tau1, op_n, op_i, 1 /*ADM Cores*/ );
        } else {
            calculate_runge_kutta( rho_tau, t_t, s.parameters.t_end, timer, progressbar, super_purpose, s, savedRhos_tau1, false /*Output*/ );
        } 
        // Interpolate saved states to equidistant timestep
        savedRhos_tau1 = Numerics::interpolate_curve( savedRhos_tau1, t_t, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        // Iterate through all outermost operators
        for ( const auto &[outer, inner_vec, purpose_vec] : eval_operators ) {
            auto &op_j_times_op_m = outer;
            for ( size_t tau1 = 0; tau1 < savedRhos_tau1.size(); tau1++ ) {
                const double t_tau_1 = savedRhos_tau1.at( tau1 ).t;
                MatrixMain rho_tau_2 = savedRhos_tau1.at( tau1 ).mat * op_j_times_op_m;
                // In a G1 or G2 loop, we would start calculating expectation values. Here, for G3, we instead iterate the resulting matrix again.
                if ( s.parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
                    MatrixMain unity_matrix( matdim, matdim );
                    unity_matrix.setIdentity();
                    calculate_path_integral_correlation( rho_tau_2, t_tau_1, s.parameters.t_end, timer, s, savedRhos_tau2, unity_matrix, unity_matrix, 1 /*ADM Cores*/ );
                } else {
                    calculate_runge_kutta( rho_tau, t_tau_1, s.parameters.t_end, timer, progressbar, super_purpose, s, savedRhos_tau2, false /*Output*/ );
                }
                // Interpolate saved states to equidistant timestep
                savedRhos_tau2 = Numerics::interpolate_curve( savedRhos_tau2, t_tau_1, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
                // Iterate through all innermost operators
                for ( int inner_index = 0; inner_index < inner_vec.size(); inner_index++ ) {
                    auto &inner_purpose = purpose_vec[inner_index];
                    auto &op_k_times_op_l = inner_vec[inner_index];
                    auto &gmat = cache[inner_purpose];
                    for ( size_t tau2 = 0; tau2 < std::min<int>( matdim, savedRhos_tau2.size() ); tau2++ ) {
                        const double t_tau_2 = savedRhos_tau2.at( tau2 ).t;
                        gmat.get( { i, tau1, tau2 } ) = s.dgl_expectationvalue<MatrixMain>( savedRhos_tau2.at( tau2 ).mat, op_k_times_op_l, t_tau_2 );
                    }
                }
            }
        }
        Timers::outputProgress( timer, progressbar, i, savedStates.size(), super_purpose );
    }

    timer.end();
    Timers::outputProgress( timer, progressbar, savedStates.size(), savedStates.size(), super_purpose, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[CorrelationFunction] G^(3) ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", super_purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );

    // Don't Apply any detector function for now.
}