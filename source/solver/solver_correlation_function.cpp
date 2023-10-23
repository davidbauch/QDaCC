#include "solver/solver_ode.h"

double QDACC::Numerics::get_tdelta( const Dense &gmat_time, size_t fixed_index, size_t var_index ) {
    return var_index == 0 ? std::real( gmat_time( var_index + 1, fixed_index ) - gmat_time( var_index, fixed_index ) ) : std::real( gmat_time( var_index, fixed_index ) - gmat_time( var_index - 1, fixed_index ) );
}
double QDACC::Numerics::get_taudelta( const Dense &gmat_time, size_t fixed_index, size_t var_index ) {
    return var_index == 0 ? std::imag( gmat_time( fixed_index, var_index + 1 ) - gmat_time( fixed_index, var_index ) ) : std::imag( gmat_time( fixed_index, var_index ) - gmat_time( fixed_index, var_index - 1 ) );
}
double QDACC::Numerics::get_tdelta( const std::vector<SaveState> &savedStates, size_t var_index ) {
    return var_index == 0 ? savedStates[var_index + 1].t - savedStates[var_index].t : savedStates[var_index].t - savedStates[var_index - 1].t;
}

/**
 * These functions calculate the n-Dimensional G^(n) correlation function.
 *
 * The n-th order correlation function is defined as:
 * `G^(n)(t1,...,tn) = < a_1^d(t1) * a_2^d(t2) * ... * a_n^d(tn) * a_n+1(tn) * a_n+2(tn-1) * ... * a_2n(t) >`
 * `= Tr( rho * a_1^d(t1) * a_2^d(t2) * ... * a_n^d(tn) * a_n+1(tn) * a_n+2(tn-1) * ... * a_2n(t) )`
 *
 * The most used correlation functions are G^(1), G^(2) and G^(3). Using the unitary
 * time transformation `U(t0,t1) = exp(-i int_t0^t1 H(t) dt)` and the Quantum
 * Regression Theorem, we can rewrite the correlation functions as:
 *
 * `G^(1)(t,tau) = < a_1^d(t) * a_2(t+tau) >`
 * `             = Tr (rho * a_1^d(t) * a_2(t+tau) )`
 * `             = Tr( [a_2 * rho(t)]->tau * a_1^d )`
 * `G^(2)(t,tau) = < a_1^d(t) * a_2^d(t+tau) * a_3(t+tau) * a_4(t) >`
 * `             = Tr( rho * a_1^d(t) * a_2^d(t+tau) * a_3(t+tau) * a_4(t))`
 * `             = Tr( [a_4 * rho(t) * a_1^d]->tau * a_2^d * a_3 )`
 * `G^(3)(t,tau1,tau2) = < a_1^d(t) * a_2^d(t+tau1) * a_3^d(t+tau2) * a_4(t+tau2) * a_5(t+tau1) * a_6(t) >`
 * `                   = Tr( rho * a_1^d(t) * a_2^d(t+tau1) * a_3^d(t+tau2) * a_4(t+tau2) * a_5(t+tau1) * a_6(t) )`
 * `                   = Tr ( [ [ a_6 * rho(t) * a_1^d ]->tau1 a_2^d * a_5 ]->tau2 * a_3^d * a_4 ) )`
 * ...
 * `G^(n)(t,tau1,...,taun-1) = < a_1^d(t) * a_2^d(t+tau1) * ... * a_n^d(t+tau_n-1) * a_n+1(t+tau_n-1) * ... * a_2n(t) >`
 * `                         = Tr( rho * a_1^d(t) * a_2^d(t+tau1) * ... * a_n^d(t+tau_n-1) * a_n+1(t+tau_n-1) * ... * a_2n(t) )`
 * `                         = Tr( [ [ [ [ [ a_2n * rho(t) * a_1^d ]->tau1 * a_2^d * a_2n-1 ]->tau2 ] * a_3^d * a_2n-2 ]->tau3 ... * a_n-1^d a_n+2 ]->taun-1 * a_n^d * a_n+1 )`
 *
 * Effectively, we are chaining regular G(2) correlation calculations together.
 * With the MultidimensionalCacheMatrix, we can, in theory, calculate any G^(n)
 * as long as we have enough memory to save the resulting N^n matrices.
 */

void QDACC::Numerics::ODESolver::calculate_g1( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[CorrelationFunction] G1(tau) for {} already exists.\n", purpose );
    }
    calculate_g1( s, std::vector<std::string>{ s_op_i }, s_op_j, { purpose } );
}

void QDACC::Numerics::ODESolver::calculate_g1( System &s, const std::vector<std::string> &s_op_i, const std::string &s_op_j, const std::vector<std::string> &purposes ) {
    const auto size = s_op_i.size();
    calculate_g2( s, "internal_identitymatrix", std::vector<std::string>{ size, "internal_identitymatrix" }, s_op_i, s_op_j, purposes );
}

void QDACC::Numerics::ODESolver::calculate_g2( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[CorrelationFunction] G2(tau) for {} already exists.\n", purpose );
    }
    calculate_g2( s, s_op_i, std::vector<std::string>{ s_op_j }, { s_op_k }, s_op_l, { purpose } );
}

void QDACC::Numerics::ODESolver::calculate_g3( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &s_op_m, const std::string &s_op_n, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[CorrelationFunction] G2(tau) for {} already exists.\n", purpose );
    }
    calculate_g3( s, s_op_i, std::vector<std::string>{ s_op_j }, { s_op_k }, {s_op_l}, {s_op_m}, s_op_n, { purpose } );
}

void QDACC::Numerics::ODESolver::calculate_g2( System &s, const std::string &s_op_i, const std::vector<std::string> &s_op_j, const std::vector<std::string> &s_op_k, const std::string &s_op_l, const std::vector<std::string> &purposes ) {
    Log::L2( "[CorrelationFunction] Preparing to calculate Correlation function\n" );
    Log::L2( "[CorrelationFunction] Generating Sparse Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_l = get_operators_matrix( s, s_op_l );
    Log::L2( "[CorrelationFunction] Iterated Operators are op_i = {} and op_l = {}\n", s_op_i, s_op_l );

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Generator super-purpose
    std::string super_purpose = std::accumulate( std::next( purposes.begin() ), purposes.end(), purposes.front(), []( const std::string &a, const std::string &b ) { return a + " and " + b; } );

    std::vector<std::pair<Sparse, std::string>> eval_operators;

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
        cache[purpose] = MultidimensionalCacheMatrix( { (int)matdim, (int)matdim }, purpose );
        auto &mat = cache[purpose];
        // Fill Time Matrix
#pragma omp parallel for collapse( 2 ) schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
        for ( size_t i = 0; i < matdim; i++ ) {
            for ( size_t j = 0; j < matdim; j++ ) {
                double tau = i + j < matdim ? s.parameters.grid_values[i + j] : s.parameters.grid_values.back() + s.parameters.grid_steps.back() * ( i + j - matdim + 1 );
                mat.get_time( i, j ) = s.parameters.grid_values[i] + 1.0i * tau;
            }
        }
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
    for ( size_t i = 0; i < std::min<size_t>( matdim, savedStates.size() ); i++ ) {
        const auto& current_state = savedStates.at( i );
        // Create and reserve past rho's vector
        std::vector<QDACC::SaveState> savedRhos; 
        // Get Time from saved State
        double t_t = current_state.t; 
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau( current_state.mat, op_l, op_i, t_t );
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
            for ( size_t j = 0; j < savedRhos.size(); j++ ) {
                const double t_tau = savedRhos.at( j ).t;
                gmat.get( i, j ) = s.dgl_expectationvalue<Sparse, Scalar>( savedRhos.at( j ).mat, eval, t_tau );
            }
        }
        Timers::outputProgress( timer, progressbar, i, savedStates.size(), super_purpose );
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

// G3 will NOT be cached like G1 and G2, because it is too memory intensive. Additionally,
// it is not required, because no other function depends on it / needs G3s (yet...).
void QDACC::Numerics::ODESolver::calculate_g3( System &s, const std::string &s_op_i, const std::vector<std::string> &s_op_j, const std::vector<std::string> &s_op_k, const std::vector<std::string> &s_op_l, const std::vector<std::string> &s_op_m, const std::string &s_op_n, const std::vector<std::string> &purposes ) {
    Log::L2( "[CorrelationFunction] Preparing to calculate Correlation function\n" );
    Log::L2( "[CorrelationFunction] Generating Sparse Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_n = get_operators_matrix( s, s_op_n );
    //Log::L2( "[CorrelationFunction] Iterated Operators are op_i = {} and op_l = {}\n", s_op_i, s_op_l );

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Generator super-purpose
    std::string super_purpose = std::accumulate( std::next( purposes.begin() ), purposes.end(), purposes.front(), []( const std::string &a, const std::string &b ) { return a + " and " + b; } );

    std::vector<std::tuple<Sparse, std::vector<Sparse>, std::vector<std::string>>> eval_operators;

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
        std::vector<Sparse> inners;
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
            cache[purpose] = MultidimensionalCacheMatrix( { (int)matdim, (int)matdim, (int)matdim }, purpose );
            // Preconstruct t,tau time. 
            auto& mat = cache[purpose];
            // Fill Time Matrix
#pragma omp parallel for collapse( 2 ) schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
            for ( size_t i = 0; i < matdim; i++ ) {
                for ( size_t j = 0; j < matdim; j++ ) {
                    double tau = i + j < matdim ? s.parameters.grid_values[i + j] : s.parameters.grid_values.back() + s.parameters.grid_steps.back() * ( i + j - matdim + 1 );
                    mat.get_time( i, j ) = s.parameters.grid_values[i] + 1.0i * tau;
                }
            }
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
    for ( int i = 0; i < std::min<int>( matdim, savedStates.size() ); i++ ) {
        // Create and reserve past rho's vector.
        std::vector<QDACC::SaveState> savedRhos_tau1, savedRhos_tau2;
        // Get state corresponding to i
        const auto& current_state = savedStates.at( i );
        // Get Time from saved State
        double t_t = current_state.t;
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau( current_state.mat, op_n, op_i, t_t );
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
            for ( int tau1 = 0; tau1 < savedRhos_tau1.size(); tau1++ ) {
                const double t_tau_1 = savedRhos_tau1.at( tau1 ).t;
                Sparse rho_tau_2 = savedRhos_tau1.at( tau1 ).mat * op_j_times_op_m;
                // In a G1 or G2 loop, we would start calculating expectation values. Here, for G3, we instead iterate the resulting matrix again.
                if ( s.parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
                    Sparse unity_matrix = Dense::Identity( matdim, matdim ).sparseView();
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
                    for ( int tau2 = 0; tau2 < std::min<int>( matdim, savedRhos_tau2.size() ); tau2++ ) {
                        const double t_tau_2 = savedRhos_tau2.at( tau2 ).t;
                        gmat.get( { i, tau1, tau2 } ) = s.dgl_expectationvalue<Sparse, Scalar>( savedRhos_tau2.at( tau2 ).mat, op_k_times_op_l, t_tau_2 );
                    }
                }
            }
        }
        Timers::outputProgress( timer, progressbar, i, savedStates.size() * eval_operators.size(), super_purpose );
    }

    timer.end();
    Timers::outputProgress( timer, progressbar, savedStates.size(), savedStates.size(), super_purpose, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[CorrelationFunction] G^(3) ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", super_purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );

    // Don't Apply any detector function for now.
}