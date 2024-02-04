#include "solver/solver_ode.h"

/**
 * @brief Calculates the G3 time bin coherence of a given operator set.
 * The operator set is supposed to represent the early and late photons of two
 * Qubits. The G3 function is then calculated as follows:
 *
 * Qubit_coherence_matrix(t1,t2,t3) = <a^d(t1) a^d(t2) a^d(t3) a(t3) a(t2) a(t1)>
 * where each operator can be either early (E = a(t)) or late (L = a(t+T)).
 * Other than in the static case, we do not calculate individual correlation functions but instead calculate the
 * complete photon matrix instead.
 *
 * The base of the matrix is build using: EEE EEL ELE ELL LLE LEE LEL LLL
 * The matrix then reads
 * | EEEEEE EEEEEL EEEELE EEEELL EEELLE EEELEE EEELEL EEELLL |
 * | EELEEE EELEEL EELELE EELELL EELLLE EELLEE EELLEL EELLLL |
 * | ELEEEE ELEEEL ELEELE ELEELL ELELLE ELELEE ELELEL ELELLL |
 * | ELLEEE ELLEEL ELLELE ELLELL ELLLLE ELLLEE ELLLEL ELLLLL |
 * | LLEEEE LLEEEL LLEELE LLEELL LLELLE LLELEE LLELEL LLELLL |
 * | LEEEEE LEEEEL LEEELE LEEELL LEELLE LEELEE LEELEL LEELLL |
 * | LELEEE LELEEL LELELE LELELL LELLLE LELLEE LELLEL LELLLL |
 * | LLLEEE LLLEEL LLLELE LLLELL LLLLLE LLLLEE LLLLEL LLLLLL |
 *
 *
 * We required 64 correlation functions in total, which are defined by the inner product of the matrix base
 * The general forumla for any of the G3 time bin functions for 0 < t1 < t1+T1/6 < 2T < t2 < t2+T2/5 < 4T < t3 < t3+T3/4 < 6T is:
 *
 *
 * Note, that in this case, all of the entries are required since this is inherintly a 4 photon process.
 *
 * We group the functions into 4 main groups, each with 3 subgroups to reduce calculation efforts as follows:
 * 1: ELEELE ELLELE ELELLE ELLLLE EELLEE EELEEE EEEEEE EEELEE
 *    ELEEEE ELLEEE ELLLEE ELELEE
 *    EEELLE EELELE EELLLE EEEELE
 * 2: LELEEL LEEEEL LEELEL LLEELL LLELLL LLLLLL LLLELL LELLEL
 *    LLEEEL LLELEL LLLLEL LLLEEL
 *    LEEELL LELLLL LELELL LEELLL
 * 3: ELLLLL EEEEEL EELEEL EEELEL EEELLL ELELLL ELLELL EELLEL EEEELL ELEELL
 *    ELEEEL ELELEL ELLLEL ELLEEL
 *    EELLLL EELELL
 * 4: LLEELE LLELLE LEEEEE LLLEEE LLLELE LEELEE LELLEE LELEEE LLLLLE LLLLEE
 *    LEEELE LEELLE LELELE LELLLE
 *    LLELEE LLEEEE
 * 
 * I decided to not use conditional statements to reduce the amount of if statements in the code, but instead write out
 * the full matrix. This is not the most efficient way, but it is the most readable and maintainable way. This way we can
 * directly compare any of the functions to the theoretical formula and ensure the implementation is correct.
 * Maybe, after everything is implemented correctly, ill change this to a more efficient way.
 * For example, we can create cases for the different branches of the matrix, and ensure the evaluated functions are also
 * cached. Then, instead of blindly recalculating everything, we create a sort of look-up-table for the intermediate results.
 */

void QDACC::Numerics::ODESolver::calculate_timebin_g3_correlations( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &s_op_m,
                                                                    const std::string &s_op_n, const std::string &purpose, double t_start, double time_bin_length, int deph ) {
    Log::L2( "[TimeBinCorrelation] Preparing to calculate Correlation function for time bin coherence\n" );
    Log::L2( "[TimeBinCorrelation] Generating Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_j = get_operators_matrix( s, s_op_j );
    const auto op_k = get_operators_matrix( s, s_op_k );
    const auto op_l = get_operators_matrix( s, s_op_l );
    const auto op_m = get_operators_matrix( s, s_op_m );
    const auto op_n = get_operators_matrix( s, s_op_n );

    Log::L2( "[TimeBinCorrelation] Iterated Operators are op_i = {}, op_j = {}, op_k = {}, op_l = {}, op_m = {}, op_n = {}\n", s_op_i, s_op_j, s_op_k, s_op_l, s_op_m, s_op_n );

    std::vector modes = {
        "EEEEEE", "EEEEEL", "EEEELE", "EEEELL", "EEELLE", "EEELEE", "EEELEL", "EEELLL", "EELEEE", "EELEEL", "EELELE", "EELELL", "EELLLE", "EELLEE", "EELLEL", "EELLLL", "ELEEEE", "ELEEEL", "ELEELE", "ELEELL", "ELELLE", "ELELEE",
        "ELELEL", "ELELLL", "ELLEEE", "ELLEEL", "ELLELE", "ELLELL", "ELLLLE", "ELLLEE", "ELLLEL", "ELLLLL", "LLEEEE", "LLEEEL", "LLEELE", "LLEELL", "LLELLE", "LLELEE", "LLELEL", "LLELLL", "LEEEEE", "LEEEEL", "LEEELE", "LEEELL",
        "LEELLE", "LEELEE", "LEELEL", "LEELLL", "LELEEE", "LELEEL", "LELELE", "LELELL", "LELLLE", "LELLEE", "LELLEL", "LELLLL", "LLLEEE", "LLLEEL", "LLLELE", "LLLELL", "LLLLLE", "LLLLEE", "LLLLEL", "LLLLLL",
    };

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Create Timer and Progresbar
    Timer &timer = Timers::create( "TimeBin TPM-Loop (" + purpose + ")" ).start();
    auto progressbar = ProgressBar();

    auto dt = s.parameters.grid_values[1] - s.parameters.grid_values[0];
    // Correct t_start and the time bin length so they are on the grid
    t_start = s.parameters.grid_values[std::max<size_t>( 0, std::min<size_t>( s.parameters.grid_values.size() - 1, std::round( t_start / dt ) ) )] - s.parameters.grid_values[0];
    time_bin_length = s.parameters.grid_values[std::max<size_t>( 0, std::min<size_t>( s.parameters.grid_values.size() - 1, std::round( time_bin_length / dt ) ) )] - s.parameters.grid_values[0];

    // 6 Time bin configuration
    auto end_of_first_bin = t_start + time_bin_length;        // end t1
    auto end_of_second_bin = t_start + 2.0 * time_bin_length; // end t1+T
    auto end_of_third_bin = t_start + 3.0 * time_bin_length;  // end t2
    auto end_of_fourth_bin = t_start + 4.0 * time_bin_length; // end t2+T
    auto end_of_fifth_bin = t_start + 5.0 * time_bin_length;  // end t3
    auto end_of_sixth_bin = t_start + 6.0 * time_bin_length;  // end t3+T
    auto index_t_start_of_first_bin = rho_index_map.lower_bound( t_start )->second;
    auto index_t_end_of_first_bin = rho_index_map.upper_bound( end_of_first_bin )->second;
    auto index_t_start_of_second_bin = index_t_end_of_first_bin;
    auto index_t_end_of_second_bin = rho_index_map.upper_bound( end_of_second_bin )->second;
    auto index_t_start_of_third_bin = index_t_end_of_second_bin;
    auto index_t_end_of_third_bin = rho_index_map.upper_bound( end_of_third_bin )->second;
    auto index_t_start_of_fourth_bin = index_t_end_of_third_bin;
    auto index_t_end_of_fourth_bin = rho_index_map.upper_bound( end_of_fourth_bin )->second;
    auto index_t_start_of_fith_bin = index_t_end_of_fourth_bin;
    auto index_t_end_of_fifth_bin = rho_index_map.upper_bound( end_of_fifth_bin )->second;
    auto index_t_start_of_sixth_bin = index_t_end_of_fifth_bin;
    auto index_t_end_of_sixth_bin = rho_index_map.upper_bound( end_of_sixth_bin )->second;

    size_t index_range = index_t_end_of_first_bin - index_t_start_of_first_bin - 1;

    for ( size_t i = 0; i < s.parameters.grid_values.size(); i++ ) {
        if ( i % 10 == 0 ) {
        }
    }

    std::map<std::string, std::string> mode_lookup;

    Log::L2( "[TimeBinCorrelation] Reserving memory for {0} matrices of size {1}x{1}x{1}.\n", modes.size(), index_range + 1 );
    for ( const auto &mode : modes ) {
        const auto name = purpose + "_" + mode;
        mode_lookup[mode] = name;
        cache[name] = MultidimensionalCacheMatrix( { index_range + 1, index_range + 1, index_range + 1 }, name );
    }

    auto &EEEEEE_cache = cache[mode_lookup["EEEEEE"]];
    auto &EEEEEL_cache = cache[mode_lookup["EEEEEL"]];
    auto &EEEELE_cache = cache[mode_lookup["EEEELE"]];
    auto &EEEELL_cache = cache[mode_lookup["EEEELL"]];
    auto &EEELLE_cache = cache[mode_lookup["EEELLE"]];
    auto &EEELEE_cache = cache[mode_lookup["EEELEE"]];
    auto &EEELEL_cache = cache[mode_lookup["EEELEL"]];
    auto &EEELLL_cache = cache[mode_lookup["EEELLL"]];
    auto &EELEEE_cache = cache[mode_lookup["EELEEE"]];
    auto &EELEEL_cache = cache[mode_lookup["EELEEL"]];
    auto &EELELE_cache = cache[mode_lookup["EELELE"]];
    auto &EELELL_cache = cache[mode_lookup["EELELL"]];
    auto &EELLLE_cache = cache[mode_lookup["EELLLE"]];
    auto &EELLEE_cache = cache[mode_lookup["EELLEE"]];
    auto &EELLEL_cache = cache[mode_lookup["EELLEL"]];
    auto &EELLLL_cache = cache[mode_lookup["EELLLL"]];
    auto &ELEEEE_cache = cache[mode_lookup["ELEEEE"]];
    auto &ELEEEL_cache = cache[mode_lookup["ELEEEL"]];
    auto &ELEELE_cache = cache[mode_lookup["ELEELE"]];
    auto &ELEELL_cache = cache[mode_lookup["ELEELL"]];
    auto &ELELLE_cache = cache[mode_lookup["ELELLE"]];
    auto &ELELEE_cache = cache[mode_lookup["ELELEE"]];
    auto &ELELEL_cache = cache[mode_lookup["ELELEL"]];
    auto &ELELLL_cache = cache[mode_lookup["ELELLL"]];
    auto &ELLEEE_cache = cache[mode_lookup["ELLEEE"]];
    auto &ELLEEL_cache = cache[mode_lookup["ELLEEL"]];
    auto &ELLELE_cache = cache[mode_lookup["ELLELE"]];
    auto &ELLELL_cache = cache[mode_lookup["ELLELL"]];
    auto &ELLLLE_cache = cache[mode_lookup["ELLLLE"]];
    auto &ELLLEE_cache = cache[mode_lookup["ELLLEE"]];
    auto &ELLLEL_cache = cache[mode_lookup["ELLLEL"]];
    auto &ELLLLL_cache = cache[mode_lookup["ELLLLL"]];
    auto &LLEEEE_cache = cache[mode_lookup["LLEEEE"]];
    auto &LLEEEL_cache = cache[mode_lookup["LLEEEL"]];
    auto &LLEELE_cache = cache[mode_lookup["LLEELE"]];
    auto &LLEELL_cache = cache[mode_lookup["LLEELL"]];
    auto &LLELLE_cache = cache[mode_lookup["LLELLE"]];
    auto &LLELEE_cache = cache[mode_lookup["LLELEE"]];
    auto &LLELEL_cache = cache[mode_lookup["LLELEL"]];
    auto &LLELLL_cache = cache[mode_lookup["LLELLL"]];
    auto &LEEEEE_cache = cache[mode_lookup["LEEEEE"]];
    auto &LEEEEL_cache = cache[mode_lookup["LEEEEL"]];
    auto &LEEELE_cache = cache[mode_lookup["LEEELE"]];
    auto &LEEELL_cache = cache[mode_lookup["LEEELL"]];
    auto &LEELLE_cache = cache[mode_lookup["LEELLE"]];
    auto &LEELEE_cache = cache[mode_lookup["LEELEE"]];
    auto &LEELEL_cache = cache[mode_lookup["LEELEL"]];
    auto &LEELLL_cache = cache[mode_lookup["LEELLL"]];
    auto &LELEEE_cache = cache[mode_lookup["LELEEE"]];
    auto &LELEEL_cache = cache[mode_lookup["LELEEL"]];
    auto &LELELE_cache = cache[mode_lookup["LELELE"]];
    auto &LELELL_cache = cache[mode_lookup["LELELL"]];
    auto &LELLLE_cache = cache[mode_lookup["LELLLE"]];
    auto &LELLEE_cache = cache[mode_lookup["LELLEE"]];
    auto &LELLEL_cache = cache[mode_lookup["LELLEL"]];
    auto &LELLLL_cache = cache[mode_lookup["LELLLL"]];
    auto &LLLEEE_cache = cache[mode_lookup["LLLEEE"]];
    auto &LLLEEL_cache = cache[mode_lookup["LLLEEL"]];
    auto &LLLELE_cache = cache[mode_lookup["LLLELE"]];
    auto &LLLELL_cache = cache[mode_lookup["LLLELL"]];
    auto &LLLLLE_cache = cache[mode_lookup["LLLLLE"]];
    auto &LLLLEE_cache = cache[mode_lookup["LLLLEE"]];
    auto &LLLLEL_cache = cache[mode_lookup["LLLLEL"]];
    auto &LLLLLL_cache = cache[mode_lookup["LLLLLL"]];

    Log::L2(
        "[TimeBinCorrelation] Calculating G^3(tau)... purpose: {0}, saving to matrix of size {1}x{1}x{1}, iterating over {2} saved "
        "states...\n",
        purpose, index_range, savedStates.size() );

    Log::L2( "[TimeBinCorrelation] Time Bin configuration is: start = {}, bin length = {}\n", t_start, time_bin_length );
    Log::L2( "[TimeBinCorrelation] The resulting bin times are: T1 = {}, T2 = {}, T3 = {}, T4 = {}, T5 = {}, T6 = {}\n", t_start, end_of_first_bin, end_of_second_bin, end_of_third_bin, end_of_fourth_bin, end_of_fifth_bin );

    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t index_t1 = index_t_start_of_first_bin; index_t1 < index_t_end_of_first_bin; index_t1++ ) {
        auto mat_t_index = index_t1 - index_t_start_of_first_bin;
        // Rho(t)
        const auto &current_state = savedStates.at( index_t1 );
        const auto &current_state_shifted = savedStates.at( index_t1 + index_range );
        // Create and reserve past rho's vector
        std::vector<QDACC::SaveState> saved_rho, saved_rho_inner, saved_rho_temp;
        // Get Time from saved State
        double t1 = current_state.t;
        double t1_p_T = current_state_shifted.t;

        // Calculate New Modified Density Matrices
        MatrixMain a_rho_t1 = s.dgl_calc_rhotau( current_state.mat, op_n, s.operatorMatrices.identity, t1 );
        MatrixMain rho_t1_ad = s.dgl_calc_rhotau( current_state.mat, s.operatorMatrices.identity, op_i, t1 );
        MatrixMain a_rho_t1_ad = s.dgl_calc_rhotau( current_state.mat, op_n, op_i, t1 );
        MatrixMain a_rho_t1pT_ad = s.dgl_calc_rhotau( current_state_shifted.mat, op_n, op_i, t1_p_T );

        // Calculate Runge Kutta. For simplicity, we dont even offer the option to use the PI here.

        MatrixMain eval = op_k * op_l;

        // Define temporary variables outside of the loops. Preallocate a maximum size for matrices
        MatrixMain temp( op_i.rows(), op_i.cols() );
        MatrixMain outer( op_i.rows(), op_i.cols() );

        // Reserve memory for the past rho vectors
        saved_rho.reserve( 2 * index_range );
        saved_rho_inner.reserve( 2 * index_range );
        saved_rho_temp.reserve( 2 * index_range );

        // Precalculate the t_index,t2,t3 dense index of the cache matrix
        std::vector<std::vector<size_t>> cache_index( index_range + 1 );
        for ( size_t t2 = 0; t2 < index_range; t2++ ) {
            cache_index[t2].reserve( index_range + 1 );
            cache_index[t2] = std::vector<size_t>( index_range + 1 );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                cache_index[t2][t3] = EELEEE_cache.get_index( { mat_t_index, t2, t3 } );
            }
        }

        // ELEELE ELLELE ELELLE ELLLLE EELLEE EELEEE EEEEEE EEELEE & ELEEEE ELLEEE ELLLEE ELELEE & EEELLE EELELE EELLLE EEEELE
        // ===================================================================================================================
        calculate_runge_kutta( a_rho_t1_ad, t1, end_of_fourth_bin + dt, timer, progressbar, purpose, s, saved_rho_inner, false );
        saved_rho_inner = Numerics::interpolate_curve( saved_rho_inner, end_of_second_bin, end_of_fourth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                       s.parameters.numerics_interpolate_method_tau );
        auto rho_size = saved_rho_inner.size() - 1;
        for ( size_t t2 = 0; t2 < index_range; t2++ ) {
            // Calculate a __ a^+ for t2
            MatrixMain outer = op_m * saved_rho_inner.at( t2 ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // EEEEEE -> Calculate Exp. Value for eval
                EEEEEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // EELLEE -> Calculate Exp. Value for eval
                EELLEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // EELEEE -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                EELEEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // EEELEE -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                EEELEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }

            saved_rho.clear();

            // Calculate a __ a^+ for t2+T
            auto t2pT = rho_size - t2;
            outer = op_m * saved_rho_inner.at( t2pT ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // ELEELE -> Calculate Exp. Value for eval
                ELEELE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // ELLLLE -> Calculate Exp. Value for eval
                ELLLLE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // ELLELE -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                ELLELE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // ELELLE -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                ELELLE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
            // Calculate a __ for t2+T
            outer = op_m * saved_rho_inner.at( t2 ).mat;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, saved_rho_inner.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho, false );
            outer = saved_rho_inner.back().mat * op_j;
            saved_rho.clear();
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // ELEEEE
                ELEEEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // ELLLEE
                ELLLEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // ELLEEE
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                ELLEEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // ELELEE
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                ELELEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
            // Calculate __ a^+ for t2+T
            outer = saved_rho_inner.at( t2 ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, saved_rho_inner.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho, false );
            outer = op_m * saved_rho.back().mat;
            saved_rho.clear();
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // EEEELE
                EEEELE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // EELLLE
                EELLLE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // EEELLE
                MatrixMain temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                EEELLE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // EELELE
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                EELELE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
        }
        saved_rho_inner.clear();
        Timers::outputProgress( timer, progressbar, 4 * mat_t_index, 4 * index_range, purpose );

        // LLEELL LLELLL LLLLLL LLLELL LELEEL LEEEEL LEELEL LELLEL & LLEEEL LLELEL LLLLEL LLLEEL & LEEELL LELLLL LELELL LEELLL
        // ===================================================================================================================
        calculate_runge_kutta( a_rho_t1pT_ad, t1_p_T, end_of_fourth_bin + dt, timer, progressbar, purpose, s, saved_rho_inner, false );
        saved_rho_inner =
            Numerics::interpolate_curve( saved_rho_inner, t1_p_T, end_of_fourth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        rho_size = saved_rho_inner.size() - 1;
        for ( size_t t2 = 0; t2 < index_range; t2++ ) {
            // Calculate a __ a^+ for t2
            MatrixMain outer = op_m * saved_rho_inner.at( t2 ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );

                // LEEEEL -> Calculate Exp. Value for eval
                LEEEEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // LELLEL -> Calculate Exp. Value for eval
                LELLEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // LELEEL -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LELEEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // LEELEL -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LEELEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
            // Calculate a __ a^+ for t2+T
            auto t2pT = rho_size - t2;
            outer = op_m * saved_rho_inner.at( t2pT ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // LLEELL -> Calculate Exp. Value for eval
                LLEELL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // LLLLLL -> Calculate Exp. Value for eval
                LLLLLL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // LLLELL -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                LLLELL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // LLELLL -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                LLELLL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
            // Calculate a __ for t2+T
            outer = op_m * saved_rho_inner.at( t2 ).mat;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, saved_rho_inner.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho, false );
            outer = saved_rho_inner.back().mat * op_j;
            saved_rho.clear();
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // LLEEEL
                LLEEEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // LLLLEL
                LLLLEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // LLLEEL
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                LLLEEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // LLELEL
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                LLELEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
            // Calculate __ a^+ for t2+T
            outer = saved_rho_inner.at( t2 ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, saved_rho_inner.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho, false );
            outer = op_m * saved_rho.back().mat;
            saved_rho.clear();
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // LEEELL
                LEEELL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // LELLLL
                LELLLL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // LEELLL
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                LEELLL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // LELELL
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, end_of_sixth_bin, timer, progressbar, purpose, s, saved_rho_temp, false );
                LELELL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
        }
        saved_rho_inner.clear();
        Timers::outputProgress( timer, progressbar, 4 * mat_t_index + 1, 4 * index_range, purpose );

        // ELLLLL EEEEEL EELEEL EEELEL EEELLL ELELLL ELLELL EELLEL EEEELL ELEELL & ELEEEL ELELEL ELLLEL ELLEEL & EELLLL EELELL
        // ===================================================================================================================
        calculate_runge_kutta( rho_t1_ad, t1, t1_p_T, timer, progressbar, purpose, s, saved_rho_inner, false );
        MatrixMain inner = op_n * saved_rho_inner.back().mat;
        calculate_runge_kutta( inner, saved_rho_inner.back().t, end_of_fourth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
        saved_rho_inner = Numerics::interpolate_curve( saved_rho, saved_rho_inner.back().t, end_of_fourth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                       s.parameters.numerics_interpolate_method_tau );
        rho_size = saved_rho_inner.size() - 1;
        for ( size_t t2 = 0; t2 < index_range; t2++ ) {
            // Calculate a __ a^+ for t2
            MatrixMain outer = op_m * saved_rho_inner.at( t2 ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // EEEEEEL -> Calculate Exp. Value for eval
                EEEEEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // EELLEL
                EELLEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // EELEEL -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                EELEEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // EEELEL -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                EEELEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
            // Calculate a __ a^+ for t2+T
            auto t2pT = rho_size - t2;
            outer = op_m * saved_rho_inner.at( t2pT ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );

                // ELEELL -> Calculate Exp. Value for eval
                ELEELL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // ELLLLL
                ELLLLL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // ELLELL -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                ELLELL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // ELELLL -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                ELELLL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();

            // Calculate ( a __ ) a^+ for ( t2+T ) t_3(+T)
            outer = op_m * saved_rho_inner.at( t2 ).mat;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, saved_rho_inner.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho, false );
            outer = saved_rho.back().mat * op_j;
            saved_rho.clear();
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // ELEEEL -> Calculate Exp. Value for eval
                ELEEEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // ELLLEL
                ELLLEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // ELLEEL -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                ELLEEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // ELELEL -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                ELELEL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();

            // Calculate a ( __ a^+ ) for (t2+T) t_3(+T)
            outer = saved_rho_inner.at( t2 ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, saved_rho_inner.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho, false );
            outer = op_m * saved_rho.back().mat;
            saved_rho.clear();
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // EEEELL -> Calculate Exp. Value for eval
                EEEELL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // EELLLL
                EELLLL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // EELELL -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                EELELL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // EEELLL -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                EEELLL_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
        }
        saved_rho_inner.clear();
        Timers::outputProgress( timer, progressbar, 4 * mat_t_index + 2, 4 * index_range, purpose );

        // LLEELE LLELLE LLLELE LLLLLE LEELEE LELLEE LELEEE LEEEEE LLELEE LLEEEE LLLEEE LEEELE LEELLE LELELE LELLLE LLLLEE
        // ===================================================================================================================
        calculate_runge_kutta( rho_t1_ad, t1, t1_p_T, timer, progressbar, purpose, s, saved_rho_inner, false );
        inner = op_n * saved_rho_inner.back().mat;
        calculate_runge_kutta( inner, saved_rho_inner.back().t, end_of_fourth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
        saved_rho_inner = Numerics::interpolate_curve( saved_rho, saved_rho_inner.back().t, end_of_fourth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                       s.parameters.numerics_interpolate_method_tau );
        rho_size = saved_rho_inner.size() - 1;
        for ( size_t t2 = 0; t2 < index_range; t2++ ) {
            // Calculate a __ a^+ for t2
            MatrixMain outer = op_m * saved_rho_inner.at( t2 ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // LEEEEE -> Calculate Exp. Value for eval
                LEEEEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // LELLEE
                LELLEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // LELEEE -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LELEEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // LEELEE -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LEELEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();

            // Calculate a __ a^+ for t2+T
            auto t2pT = rho_size - t2;
            outer = op_m * saved_rho_inner.at( t2pT ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );

                // LLEELE -> Calculate Exp. Value for eval
                LLEELE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // LLLLLE
                LLLLLE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // LLLELE -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LLLELE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // LLELLE -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LLELLE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();

            // Calculate ( a __ ) a^+ for ( t2+T ) t_3(+T)
            outer = op_m * saved_rho_inner.at( t2 ).mat;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, saved_rho_inner.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho, false );
            outer = saved_rho.back().mat * op_j;
            saved_rho.clear();
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // LEEELE -> Calculate Exp. Value for eval
                LEEELE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // LELLLE
                LELLLE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // LELELE -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LELELE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // LEELLE -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LEELLE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
            // Calculate a ( __ a^+ ) for (t2+T) t_3(+T)
            outer = saved_rho_inner.at( t2 ).mat * op_j;
            calculate_runge_kutta( outer, saved_rho_inner.at( t2 ).t, saved_rho_inner.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho, false );
            outer = op_m * saved_rho.back().mat;
            saved_rho.clear();
            calculate_runge_kutta( outer, saved_rho_inner.at( t2pT ).t, end_of_sixth_bin + dt, timer, progressbar, purpose, s, saved_rho, false );
            saved_rho = Numerics::interpolate_curve( saved_rho, end_of_fourth_bin, end_of_sixth_bin + dt, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                     s.parameters.numerics_interpolate_method_tau );
            for ( size_t t3 = 0; t3 < index_range; t3++ ) {
                auto t3pT = t3 + index_range;
                auto &state_at_t3 = saved_rho.at( t3 );
                auto &state_at_t3pT = saved_rho.at( t3pT );
                // LLEEEE -> Calculate Exp. Value for eval
                LLEEEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3.mat, eval, state_at_t3.t );
                // LLLLEE
                LLLLEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( state_at_t3pT.mat, eval, state_at_t3pT.t );
                if ( deph <= 1 ) continue;
                // LLLEEE -> Iterate a*rho to t3+T, then calculate Exp Value for a^+
                temp = op_l * state_at_t3.mat;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LLLEEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_k, saved_rho_temp.back().t );
                saved_rho_temp.clear();
                // LLELEE -> Iterate rho*a^+ to t3+T, then calculate Exp Value for a
                temp = state_at_t3.mat * op_k;
                calculate_runge_kutta( temp, state_at_t3.t, state_at_t3.t + time_bin_length, timer, progressbar, purpose, s, saved_rho_temp, false );
                LLELEE_cache.get( cache_index[t2][t3] ) = s.dgl_expectationvalue<MatrixMain>( saved_rho_temp.back().mat, op_l, saved_rho_temp.back().t );
                saved_rho_temp.clear();
            }
            saved_rho.clear();
        }
        saved_rho_inner.clear();

        Timers::outputProgress( timer, progressbar, 4 * mat_t_index + 3, 4 * index_range, purpose );
    }
    timer.end();
    Timers::outputProgress( timer, progressbar, index_range, index_range, purpose, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[TimeBinCorrelation] G ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );
}
