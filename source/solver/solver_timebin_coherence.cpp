#include "solver/solver_ode.h"

/**
 * @brief Calculates the G2 time bin coherence of a given operator set.
 * The operator set is supposed to represent the early and late photons of two 
 * Qubits. The G2 function is then calculated as follows:
 * 
 * Qubit_coherence_matrix(t1,t2) = <a^d(t1) a^d(t2) a(t2 + T) a(t1 + T)>
 * which with t1=t and t2 = 2T+tau turns out to be
 * Tr( a_k [ [ a_l [ rho(t) a_i^d ]_{t->t+T} ]_{t+T ->2T+tau} a_j^d ]_{2T+tau -> 3T+tau} )
*/

// EEEL
void QDACC::Numerics::ODESolver::calculate_timebin_coherence( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &purpose, const double t_start, const double time_bin_length ) {
    Log::L2( "[CorrelationFunction] Preparing to calculate Correlation function for time bin coherence\n" );
    Log::L2( "[CorrelationFunction] Generating Sparse Operator Matrices from String input...\n" );
    
    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_j = get_operators_matrix( s, s_op_j );
    const auto op_k = get_operators_matrix( s, s_op_k );
    const auto op_l = get_operators_matrix( s, s_op_l );
    const Sparse identity_matrix = Dense::Identity( op_i.rows(), op_i.cols() ).sparseView();
    Log::L2( "[CorrelationFunction] Iterated Operators are op_i = {}, op_j = {}, op_k = {} and op_l = {}\n", s_op_i, s_op_j, s_op_k, s_op_l );

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );


    // Create Timer and Progresbar
    Timer &timer = Timers::create( "Correlation-Loop (" + purpose + ")" ).start();
    auto progressbar = ProgressBar();

    // Time bin configuration
    auto end_of_first_bin = t_start + time_bin_length;
    auto end_of_second_bin = end_of_first_bin + time_bin_length;
    auto end_of_third_bin = end_of_second_bin + time_bin_length;
    auto end_of_fourth_bin = end_of_third_bin + time_bin_length;
    auto index_t_start = rho_index_map.lower_bound(t_start)->second;
    auto index_t_end_of_first_bin = rho_index_map.upper_bound(end_of_first_bin)->second;
    auto index_range = index_t_end_of_first_bin - index_t_start;

    cache[purpose] = MultidimensionalCacheMatrix( { (int)index_range, (int)index_range }, purpose );
    auto &gmat = cache[purpose];
    
    Log::L2( "[CorrelationFunction] Calculating G^2(tau)... purpose: {}, saving to matrix of size {}x{},  iterating over {} saved states...\n", purpose, matdim, matdim, index_range );

    Log::L2( "[CorrelationFunction] Time Bin configuration is: start = {}, bin length = {}\n", t_start, time_bin_length);
    Log::L2( "[CorrelationFunction] Which results in: starting_index = {}, ending_index = {}, index_range = {}\n", index_t_start, index_t_end_of_first_bin, index_range);
    Log::L2( "[CorrelationFunction] The resulting bin times are: T1 = {}, T2 = {}, T3 = {}\n", end_of_first_bin, end_of_second_bin, end_of_third_bin);

    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t index_t1 = index_t_start; index_t1 < index_t_end_of_first_bin; index_t1++ ) {
        // Rho(t)
        const auto& current_state = savedStates.at( index_t1 );
        // Create and reserve past rho's vector
        std::vector<QDACC::SaveState> savedRhos; 
        // Get Time from saved State
        double t_t = current_state.t; 
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau( current_state.mat, identity_matrix, op_i, t_t );
        // Calculate Runge Kutta. For simplicity, we dont even offer the option to use the PI here.

        // First, iterate rho(t) a_i^d from t_start to t+time_bin_length
        calculate_runge_kutta( rho_tau, t_t, t_t + time_bin_length, timer, progressbar, purpose, s, savedRhos, false  );
        auto time_rho_tau_transformed = savedRhos.back().t;
        Sparse rho_times_a_i_d = savedRhos.back().mat;
        // Erase cached rhos, we dont need them anymore.
        savedRhos.clear();

        // Calculate a_k * rho_times_a_i_d. savedRhos will then contain the tau-dependent results, which we will then again iteratre towars the final time
        Sparse a_l_times_rho_times_a_i_d = s.dgl_calc_rhotau( rho_times_a_i_d, op_l, identity_matrix, time_rho_tau_transformed );
        calculate_runge_kutta( a_l_times_rho_times_a_i_d, time_rho_tau_transformed, t_start+2.0*time_bin_length, timer, progressbar, purpose, s, savedRhos, false  );
        a_l_times_rho_times_a_i_d = savedRhos.back().mat;
        time_rho_tau_transformed = savedRhos.back().t;
        savedRhos.clear();
        // "Tau" direction. Delta tau is in between 0 and T.
        calculate_runge_kutta( a_l_times_rho_times_a_i_d, time_rho_tau_transformed, time_rho_tau_transformed + 1.0*time_bin_length, timer, progressbar, purpose, s, savedRhos, false  );

        // Calculate (a_k * rho_times_a_i_d)(t = 2T+tau) * a_j^d and iterate to 3T+tau
        savedRhos = Numerics::interpolate_curve( savedRhos, time_rho_tau_transformed, time_rho_tau_transformed + 1.0*time_bin_length, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );

        Sparse eval = op_j*op_k;
        //std::cout << "Adding contributions to gmat.\n";
        //std::cout << "For t1 = " << t_t << " -> index = " << index_t1 << std::endl;
        //std::cout << "    t2 in [" << savedRhos.front().t << "," << savedRhos.back().t << "]" << " -> " << savedRhos.size() << " elements." << std::endl;
        for ( size_t tau_index = 0; tau_index < savedRhos.size(); tau_index++ ) {
            const double t_tau = savedRhos.at( tau_index ).t ;
            const auto& rho_tau = savedRhos.at( tau_index ).mat;
            gmat.get( index_t1 - index_t_start, tau_index ) = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, eval, t_tau );
        }
        
        Timers::outputProgress( timer, progressbar, index_t1-index_t_start, index_range, purpose );
    }

    timer.end();
    Timers::outputProgress( timer, progressbar, savedStates.size(), savedStates.size(), purpose, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[CorrelationFunction] G ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );

}

/*
void QDACC::Numerics::ODESolver::calculate_timebin_coherence( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &purpose, const double t_start, const double time_bin_length ) {
    Log::L2( "[CorrelationFunction] Preparing to calculate Correlation function for time bin coherence\n" );
    Log::L2( "[CorrelationFunction] Generating Sparse Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_l = get_operators_matrix( s, s_op_l );
    const auto &op_j = get_operators_matrix( s, s_op_j );
    const auto &op_k = get_operators_matrix( s, s_op_k );
    const Sparse identity_matrix = Dense::Identity( op_i.rows(), op_i.cols() ).sparseView();
    Log::L2( "[CorrelationFunction] Iterated Operators are op_i = {} and op_l = {}\n", s_op_i, s_op_l );

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );


    // Create Timer and Progresbar
    Timer &timer = Timers::create( "Correlation-Loop (" + purpose + ")" ).start();
    auto progressbar = ProgressBar();

    // Time bin configuration
    auto end_of_first_bin = t_start + time_bin_length;
    auto end_of_second_bin = end_of_first_bin + time_bin_length;
    auto end_of_third_bin = end_of_second_bin + time_bin_length;
    auto end_of_fourth_bin = end_of_third_bin + time_bin_length;
    auto index_t_start = rho_index_map.lower_bound(t_start)->second;
    auto index_t_end_of_first_bin = rho_index_map.upper_bound(end_of_first_bin)->second;
    auto index_range = index_t_end_of_first_bin - index_t_start;

    cache[purpose] = MultidimensionalCacheMatrix( { (int)index_range, (int)index_range }, purpose );
    auto &gmat = cache[purpose];
    
    Log::L2( "[CorrelationFunction] Calculating G^2(tau)... purpose: {}, saving to matrix of size {}x{},  iterating over {} saved states...\n", purpose, matdim, matdim, index_range );

    Log::L2( "[CorrelationFunction] Time Bin configuration is: start = {}, bin length = {}\n", t_start, time_bin_length);
    Log::L2( "[CorrelationFunction] Which results in: starting_index = {}, ending_index = {}, index_range = {}\n", index_t_start, index_t_end_of_first_bin, index_range);
    Log::L2( "[CorrelationFunction] The resulting bin times are: T1 = {}, T2 = {}, T3 = {}\n", end_of_first_bin, end_of_second_bin, end_of_third_bin);

    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t index_t1 = index_t_start; index_t1 < index_t_end_of_first_bin; index_t1++ ) {
        // Rho(t)
        const auto& current_state = savedStates.at( index_t1 );
        // Create and reserve past rho's vector
        std::vector<QDACC::SaveState> savedRhos; 
        // Get Time from saved State
        double t_t = current_state.t; 
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau( current_state.mat, identity_matrix, op_i, t_t );
        // Calculate Runge Kutta. For simplicity, we dont even offer the option to use the PI here.

        // First, iterate rho(t) a_i^d from t_start to t+time_bin_length
        //std::cout << "1. - First iteration, we start at t1 = t_t = " << t_t << " and we iterate a single matrix to t1 + T = " << t_t + time_bin_length << std::endl;
        //std::cout << "     The maximum time we can reach this way is t_start + 2T = " << t_start + 2.0*time_bin_length;
        calculate_runge_kutta( rho_tau, t_t, t_t + time_bin_length, timer, progressbar, purpose, s, savedRhos, false );
        //std::cout << ", the resulting time is " << savedRhos.back().t << ". We now have a_k *[rho(t1)*a_i^d]." << std::endl;
        auto time_rho_tau_transformed = savedRhos.back().t;
        Sparse rho_times_a_i_d = savedRhos.back().mat;
        // Erase cached rhos, we dont need them anymore.
        savedRhos.clear();

        // Calculate a_k * rho_times_a_i_d. savedRhos will then contain the tau-dependent results, which we will then again iteratre towars the final time
        Sparse a_k_times_rho_times_a_i_d = s.dgl_calc_rhotau( rho_times_a_i_d, op_k, identity_matrix, time_rho_tau_transformed );
        //std::cout << "2. - Second iteration, we start at t_t + T = " << time_rho_tau_transformed << ", which should be smaller than t_start + 2T = " << t_start + 2.0*time_bin_length << std::endl;
        //std::cout << "     We integrate the matrix to t_start + 2T = " << t_start + 2.0*time_bin_length << std::endl;
        calculate_runge_kutta( a_k_times_rho_times_a_i_d, time_rho_tau_transformed, t_start+2.0*time_bin_length, timer, progressbar, purpose, s, savedRhos, false );
        a_k_times_rho_times_a_i_d = savedRhos.back().mat;
        time_rho_tau_transformed = savedRhos.back().t;
        savedRhos.clear();
        //std::cout << "3. - Third iteration, we should now be at t_start + 2T = " << t_start + 2.0*time_bin_length <<". Our real matrix is at " << time_rho_tau_transformed << "." << std::endl;
        //std::cout << "     We now evaluate the tau integral. This means, we iterate the current matrix by another deltaTau in (0,T). We save the taus in a vector." << std::endl;
        //std::cout << "     The maximum time we should be able to reach this way is t_start + 3T " << t_start + 3.0*time_bin_length << std::endl;
        // "Tau" direction. Delta tau is in between 0 and T.
        calculate_runge_kutta( a_k_times_rho_times_a_i_d, time_rho_tau_transformed, time_rho_tau_transformed + time_bin_length, timer, progressbar, purpose, s, savedRhos, false );
        //std::cout << "   - The resulting time minimum is " << savedRhos.front().t << " and the maximum is " << savedRhos.back().t << ", which should lie between " << t_start+2.0*time_bin_length << " and " << t_start+3.0*time_bin_length << std::endl;

        // Calculate (a_k * rho_times_a_i_d)(t = 2T+tau) * a_j^d and iterate to 3T+tau
        //std::cout << "4. - Fourth iteration, we now iterate another time bin for each tau. This way, we should only reach a maximum time of t_start + 4T = " << t_start + 4.0*time_bin_length << std::endl; 
        //std::cout << "     We still have " << savedRhos.size() << " elements. Since we want a NxN matrix, we interpolate to ";
        savedRhos = Numerics::interpolate_curve( savedRhos, time_rho_tau_transformed, time_rho_tau_transformed + time_bin_length, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        //std::cout << savedStates.size() << " elements." << std::endl;
        std::ranges::for_each( savedRhos, [&]( auto &state ) {
            state.mat = (state.mat*op_j).eval(); // Prevent possible aliasing by using .eval()
            std::vector<QDACC::SaveState> temp;
            //std::cout << "   - We start at t = " << state.t << " and iterate to t = " << state.t + time_bin_length << ", which should be smaller than " << t_start+4.0*time_bin_length << std::endl;
            calculate_runge_kutta( state.mat, state.t, state.t + time_bin_length, timer, progressbar, purpose, s, temp, false );
            state.mat = temp.back().mat;
            state.t = temp.back().t; // These times should now all be between 3T and 4T.
            //std::cout << "    We now have a matrix at t = " << state.t << std::endl;
        } );
        
        for ( size_t tau_index = 0; tau_index < savedRhos.size(); tau_index++ ) {
            const double t_tau = savedRhos.at( tau_index ).t ;
            gmat.get( index_t1 - index_t_start, tau_index ) = s.dgl_expectationvalue<Sparse, Scalar>( savedRhos.at( tau_index ).mat, op_l, t_tau );
        }
        
        Timers::outputProgress( timer, progressbar, index_t1-index_t_start, index_range, purpose );
    }

    timer.end();
    Timers::outputProgress( timer, progressbar, savedStates.size(), savedStates.size(), purpose, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[CorrelationFunction] G ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );

}
*/