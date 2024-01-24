#include "solver/solver_ode.h"

/**
 * @brief Calculates the G2 time bin coherence of a given operator set.
 * The operator set is supposed to represent the early and late photons of two
 * Qubits. The G2 function is then calculated as follows:
 *
 * Qubit_coherence_matrix(t1,t2) = <a^d(t1) a^d(t2) a(t2) a(t1)>
 * where each operator can be either early (E = a(t)) or late (L = a(t+T)).
 * Other than in the static case, we do not calculate individual correlation functions but instead calculate the
 * complete photon matrix instead. The Two Photon Matrix is actually a two qubit matrix and is equal to the static case
 * (H==E and V==L) with
 *
 * | EEEE EEEL EELE EELL |
 * | ELEE ELEL ELLE ELLL |
 * | LEEE LEEL LELE LELL |
 * | LLEE LLEL LLLE LLLL |
 *
 * We required 16 correlation functions in total, which are defined by the inner product of EE,EL,LE,LL
 * The general forumla for any of the G2 time bin functions for 0 < t1 < t1+T1/4 < 2T < t2 < t2+T2/3 < 4T is:
 *
 *                                                       / T2 < T3: U^+_{t2+T2->t2+T3} \
 * Tr( rho(t1) U^+_{t1->t1+T1} a^+ U^+_{t1+T1->t2+T2} * | T2 = T3: 1                    | * a U_{t1+T4->t2+T3} a
 * U_{t1->t1+T4} ) \ T2 > T3: U_{t2+T3->t2+T2}   /
 *
 * Note, that in this case, all of the entries are required since this is inherintly a 4 photon process.
 *
 * We group the functions to reduce calculation efforts as follows:
 *  EELL ELLL EEEL ELEL -> [rho(t1) a^+]_{t1,t1+T}
 *  LLEE LLLE LEEE LELE -> [a rho(t1)]_{t1,t1+T}
 *  EEEE ELEE EELE ELLE -> [a rho(t1) a^+]_{t1,t2+T}
 *  LLLL LELL LLEL LEEL -> [a rho(t1 + T) a^+]_{t1,t2+T}
 */

void QDACC::Numerics::ODESolver::calculate_timebin_g2_correlations( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &purpose, double t_start,
                                                                    double time_bin_length ) {
    Log::L2( "[TimeBinCorrelation] Preparing to calculate Correlation function for time bin coherence\n" );
    Log::L2( "[TimeBinCorrelation] Generating Sparse Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_i = get_operators_matrix( s, s_op_i );
    const auto op_j = get_operators_matrix( s, s_op_j );
    const auto op_k = get_operators_matrix( s, s_op_k );
    const auto op_l = get_operators_matrix( s, s_op_l );
    const Sparse identity_matrix = Dense::Identity( op_i.rows(), op_i.cols() ).sparseView();
    Log::L2( "[TimeBinCorrelation] Iterated Operators are op_i = {}, op_j = {}, op_k = {} and op_l = {}\n", s_op_i, s_op_j, s_op_k, s_op_l );

    // We always evaluate the full matrix.
    std::vector modes = {
        "EEEE", "EEEL", "EELE", "EELL", "ELEE", "ELEL", "ELLE", "ELLL", "LEEE", "LEEL", "LELE", "LELL", "LLEE", "LLEL", "LLLE", "LLLL",
    };
    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Create Timer and Progresbar
    Timer &timer = Timers::create( "TimeBin TPM-Loop (" + purpose + ")" ).start();
    auto progressbar = ProgressBar();

    // 4 Time bin configuration
    auto end_of_first_bin = t_start + time_bin_length;
    auto end_of_second_bin = end_of_first_bin + time_bin_length;
    auto end_of_third_bin = end_of_second_bin + time_bin_length;
    auto end_of_fourth_bin = end_of_third_bin + time_bin_length;
    auto index_t_start_of_first_bin = rho_index_map.lower_bound( t_start )->second;
    auto index_t_end_of_first_bin = rho_index_map.upper_bound( end_of_first_bin )->second;
    auto index_t_start_of_second_bin = index_t_end_of_first_bin;
    auto index_t_end_of_second_bin = rho_index_map.upper_bound( end_of_second_bin )->second;
    auto index_t_start_of_third_bin = index_t_end_of_second_bin;
    auto index_t_end_of_third_bin = rho_index_map.upper_bound( end_of_third_bin )->second;
    auto index_t_start_of_fourth_bin = index_t_end_of_third_bin;
    auto index_t_end_of_fourth_bin = rho_index_map.upper_bound( end_of_fourth_bin )->second;

    size_t index_range = index_t_end_of_first_bin - index_t_start_of_first_bin;

    std::map<std::string, std::string> mode_lookup;

    Log::L2( "[TimeBinCorrelation] Reserving memory for {} matrices of size {}x{}.\n", modes.size(), index_range, index_range );
    for ( const auto &mode : modes ) {
        const auto name = purpose + "_" + mode;
        mode_lookup[mode] = name;
        cache[name] = MultidimensionalCacheMatrix( { index_range, index_range }, name );
    }

    auto &LLLE_cache = cache[mode_lookup["LLLE"]];
    auto &LEEE_cache = cache[mode_lookup["LEEE"]];
    auto &LELE_cache = cache[mode_lookup["LELE"]];
    auto &LLEE_cache = cache[mode_lookup["LLEE"]];

    auto &EELL_cache = cache[mode_lookup["EELL"]];
    auto &ELLL_cache = cache[mode_lookup["ELLL"]];
    auto &EEEL_cache = cache[mode_lookup["EEEL"]];
    auto &ELEL_cache = cache[mode_lookup["ELEL"]];

    auto &EEEE_cache = cache[mode_lookup["EEEE"]];
    auto &ELEE_cache = cache[mode_lookup["ELEE"]];
    auto &EELE_cache = cache[mode_lookup["EELE"]];
    auto &ELLE_cache = cache[mode_lookup["ELLE"]];
    auto &LLLL_cache = cache[mode_lookup["LLLL"]];
    auto &LELL_cache = cache[mode_lookup["LELL"]];
    auto &LLEL_cache = cache[mode_lookup["LLEL"]];
    auto &LEEL_cache = cache[mode_lookup["LEEL"]];

    Log::L2(
        "[TimeBinCorrelation] Calculating G^2(tau)... purpose: {}, saving to matrix of size {}x{}, iterating over {} saved "
        "states...\n",
        purpose, index_range, index_range, index_range );

    Log::L2( "[TimeBinCorrelation] Time Bin configuration is: start = {}, bin length = {}\n", t_start, time_bin_length );
    Log::L2( "[TimeBinCorrelation] The resulting bin times are: T1 = {}, T2 = {}, T3 = {}, T4 = {}\n", t_start, end_of_first_bin, end_of_second_bin, end_of_third_bin );

    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t index_t1 = index_t_start_of_first_bin; index_t1 < index_t_end_of_first_bin; index_t1++ ) {
        auto mat_t_index = index_t1 - index_t_start_of_first_bin;
        // Rho(t)
        const auto &current_state = savedStates.at( index_t1 );
        const auto &current_state_shifted = savedStates.at( index_t1 + index_range );
        // Create and reserve past rho's vector
        std::vector<QDACC::SaveState> saved_rho, saved_rho_inner;
        // Get Time from saved State
        double t1 = current_state.t;
        double t1_p_T = current_state_shifted.t;

        // Calculate New Modified Density Matrices
        Sparse a_rho_t1 = s.dgl_calc_rhotau( current_state.mat, op_l, identity_matrix, t1 );
        Sparse rho_t1_ad = s.dgl_calc_rhotau( current_state.mat, identity_matrix, op_i, t1 );
        Sparse a_rho_t1_ad = s.dgl_calc_rhotau( current_state.mat, op_l, op_i, t1 );
        Sparse a_rho_t1pT_ad = s.dgl_calc_rhotau( current_state_shifted.mat, op_l, op_i, t1_p_T );

        // Calculate Runge Kutta. For simplicity, we dont even offer the option to use the PI here.

        // LLLE,LLEE,LEEEE,LELE
        // =============================================================================================
        calculate_runge_kutta( a_rho_t1, t1, t1 + time_bin_length, timer, progressbar, purpose, s, saved_rho,
                               false ); // L..E. We only need the last matrix here
        Sparse a_rho_t1_ad_t1pT = saved_rho.back().mat * op_i;
        saved_rho.clear();
        // LEEE and LLLE
        calculate_runge_kutta( a_rho_t1_ad_t1pT, t1 + time_bin_length, end_of_fourth_bin, timer, progressbar, purpose, s, saved_rho, false );
        saved_rho = Numerics::interpolate_curve( saved_rho, t1 + time_bin_length, end_of_fourth_bin, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                 s.parameters.numerics_interpolate_method_tau );
        auto rho_size = saved_rho.size() - 1;
        Sparse eval = op_j * op_k;
        for ( size_t tau = 0; tau < index_range; tau++ ) {
            auto t2 = rho_size - index_range - tau;
            auto t2pT = rho_size - tau;
            LLLE_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho.at( t2pT ).mat, eval, saved_rho.at( t2pT ).t );
            LEEE_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho.at( t2 ).mat, eval, saved_rho.at( t2 ).t );
            // LLEE and LELE; Iterate further
            Sparse a_a_rho_t1_t1pT_ad_t2 = op_k * saved_rho.at( t2 ).mat;
            Sparse a_rho_t1_t1pT_ad_t2_ad = saved_rho.at( t2 ).mat * op_j;
            // Advance further to t2+T
            calculate_runge_kutta( a_a_rho_t1_t1pT_ad_t2, saved_rho.at( t2 ).t, saved_rho.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho_inner, false );
            LLEE_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho_inner.back().mat, op_j, saved_rho_inner.back().t );
            calculate_runge_kutta( a_rho_t1_t1pT_ad_t2_ad, saved_rho.at( t2 ).t, saved_rho.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho_inner, false );
            LELE_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho_inner.back().mat, op_k, saved_rho_inner.back().t );
            saved_rho_inner.clear();
        }
        saved_rho.clear();

        // EELL ELLL EEEL ELEL
        // =============================================================================================
        calculate_runge_kutta( rho_t1_ad, t1, t1 + time_bin_length, timer, progressbar, purpose, s, saved_rho,
                               false ); // E..L. We only need the last matrix here
        a_rho_t1_ad_t1pT = op_l * saved_rho.back().mat;
        saved_rho.clear();
        // ELLL and EEEL
        calculate_runge_kutta( a_rho_t1_ad_t1pT, t1 + time_bin_length, end_of_fourth_bin, timer, progressbar, purpose, s, saved_rho, false );
        saved_rho = Numerics::interpolate_curve( saved_rho, t1 + time_bin_length, end_of_fourth_bin, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false,
                                                 s.parameters.numerics_interpolate_method_tau );
        rho_size = saved_rho.size() - 1;
        for ( size_t tau = 0; tau < index_range; tau++ ) {
            auto t2 = rho_size - index_range - tau;
            auto t2pT = rho_size - tau;
            ELLL_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho.at( t2pT ).mat, eval, saved_rho.at( t2pT ).t );
            EEEL_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho.at( t2 ).mat, eval, saved_rho.at( t2 ).t );
            // EELL and ELEL; iterate further
            Sparse a_rho_t1_t1pT_ad_t2_ad = saved_rho.at( t2 ).mat * op_j;
            Sparse a_a_rho_t1_t1pT_ad_t2 = op_k * saved_rho.at( t2 ).mat;
            // Advance further to t2+T
            calculate_runge_kutta( a_rho_t1_t1pT_ad_t2_ad, saved_rho.at( t2 ).t, saved_rho.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho_inner, false );
            EELL_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho_inner.back().mat, op_k, saved_rho_inner.back().t );
            saved_rho_inner.clear();
            calculate_runge_kutta( a_a_rho_t1_t1pT_ad_t2, saved_rho.at( t2 ).t, saved_rho.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho_inner, false );
            ELEL_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho_inner.back().mat, op_j, saved_rho_inner.back().t );
            saved_rho_inner.clear();
        }
        saved_rho.clear();

        // EEEE ELEE EELE ELLE
        // =============================================================================================
        calculate_runge_kutta( a_rho_t1_ad, t1, end_of_fourth_bin, timer, progressbar, purpose, s, saved_rho, false );
        saved_rho = Numerics::interpolate_curve( saved_rho, t1, end_of_fourth_bin, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        rho_size = saved_rho.size() - 1;
        for ( size_t tau = 0; tau < index_range; tau++ ) {
            auto t2 = rho_size - index_range - tau;
            auto t2pT = rho_size - tau;
            EEEE_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho.at( t2 ).mat, eval, saved_rho.at( t2 ).t );
            ELLE_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho.at( t2pT ).mat, eval, saved_rho.at( t2pT ).t );
            // ELLE and EELE require additional time transformation
            Sparse a_a_rhot1_ad_t2 = op_k * saved_rho.at( t2 ).mat;
            Sparse a_rho_t1_ad_t2_ad = saved_rho.at( t2 ).mat * op_j;
            // Advance further to t2+T
            calculate_runge_kutta( a_a_rhot1_ad_t2, saved_rho.at( t2 ).t, saved_rho.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho_inner, false );
            EELE_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho_inner.back().mat, op_j, saved_rho_inner.back().t );
            saved_rho_inner.clear();
            calculate_runge_kutta( a_rho_t1_ad_t2_ad, saved_rho.at( t2 ).t, saved_rho.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho_inner, false );
            ELLE_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho_inner.back().mat, op_k, saved_rho_inner.back().t );
            saved_rho_inner.clear();
        }
        saved_rho.clear();

        // LLLL LELL LLEL LEEL
        // =============================================================================================
        calculate_runge_kutta( a_rho_t1pT_ad, t1_p_T, end_of_fourth_bin, timer, progressbar, purpose, s, saved_rho, false );
        saved_rho = Numerics::interpolate_curve( saved_rho, t1_p_T, end_of_fourth_bin, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        rho_size = saved_rho.size() - 1;
        for ( size_t tau = 0; tau < index_range; tau++ ) {
            auto t2 = rho_size - index_range - tau;
            auto t2pT = rho_size - tau;
            LLLL_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho.at( t2pT ).mat, eval, saved_rho.at( t2pT ).t );
            LLEL_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho.at( t2 ).mat, eval, saved_rho.at( t2 ).t );
            // LLEL and LEEL require additional time transformation
            Sparse a_rho_t1pT_ad_t2_ad = saved_rho.at( t2 ).mat * op_j;
            Sparse a_a_rhot_ad_t1pT_t2 = op_k * saved_rho.at( t2 ).mat;
            // Advance further to t2+T
            calculate_runge_kutta( a_rho_t1pT_ad_t2_ad, saved_rho.at( t2 ).t, saved_rho.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho_inner, false );
            LLEL_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho_inner.back().mat, op_k, saved_rho_inner.back().t );
            saved_rho_inner.clear();
            calculate_runge_kutta( a_a_rhot_ad_t1pT_t2, saved_rho.at( t2 ).t, saved_rho.at( t2pT ).t, timer, progressbar, purpose, s, saved_rho_inner, false );
            LEEL_cache.get( mat_t_index, index_range - 1 - tau ) = s.dgl_expectationvalue<Sparse, Scalar>( saved_rho_inner.back().mat, op_j, saved_rho_inner.back().t );
            saved_rho_inner.clear();
        }
        saved_rho.clear();

        Timers::outputProgress( timer, progressbar, index_t1 - index_t_start_of_first_bin, index_range, purpose );
    }
    timer.end();
    Timers::outputProgress( timer, progressbar, index_range, index_range, purpose, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[TimeBinCorrelation] G ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );
}