#include "solver/solver_ode.h"

double QDLC::Numerics::get_tdelta( const Dense &gmat_time, size_t fixed_index, size_t var_index ) {
    return var_index == 0 ? std::real( gmat_time( var_index + 1, fixed_index ) - gmat_time( var_index, fixed_index ) ) : std::real( gmat_time( var_index, fixed_index ) - gmat_time( var_index - 1, fixed_index ) );
}
double QDLC::Numerics::get_taudelta( const Dense &gmat_time, size_t fixed_index, size_t var_index ) {
    return var_index == 0 ? std::imag( gmat_time( fixed_index, var_index + 1 ) - gmat_time( fixed_index, var_index ) ) : std::imag( gmat_time( fixed_index, var_index ) - gmat_time( fixed_index, var_index - 1 ) );
}
double QDLC::Numerics::get_tdelta( const std::vector<SaveState> &savedStates, size_t var_index ) {
    return var_index == 0 ? savedStates[var_index + 1].t - savedStates[var_index].t : savedStates[var_index].t - savedStates[var_index - 1].t;
}

void QDLC::Numerics::ODESolver::calculate_g1( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[G1Correlation] G1(tau) for {} already exists.\n", purpose );
    }
    calculate_g1( s, std::vector<std::string>{ s_op_creator }, s_op_annihilator, { purpose } ); // String in braced initializer list = string, not vector!
}

// TODO: code repetition vermeiden, generische G_x funktion? ODer wenigstens nur den uniquen teil dublizieren
void QDLC::Numerics::ODESolver::calculate_g1( System &s, const std::vector<std::string> &s_op_creator, const std::string &s_op_annihilator, const std::vector<std::string> &purposes ) {
    // Find Operator Matrices
    Log::L2( "[G1Correlation] Preparing to calculate G1 Correlation function\n" );
    Log::L2( "[G1Correlation] Generating Sparse Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_annihilator = get_operators_matrix( s, s_op_annihilator );

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Generator super-purpose
    std::string super_purpose = std::accumulate( std::next( purposes.begin() ), purposes.end(), purposes.front(), []( const std::string &a, const std::string &b ) { return a + " and " + b; } );

    std::vector<std::pair<Sparse, std::string>> eval_operators;

    // Preconstruct
    for ( auto current = 0; current < s_op_creator.size(); current++ ) {
        const auto &op_creator = get_operators_matrix( s, s_op_creator[current] );
        const auto &purpose = purposes[current];
        if ( cache.contains( purpose ) ) {
            Log::L2( "[G1Correlation] Matrix for {} already exists! Skipping!\n", purpose );
            continue;
        }
        // Construct Evaluation Operators
        eval_operators.emplace_back( op_creator, purpose );

        Log::L2( "[G1Correlation] Preparing Cache Matrices for {}...\n", purpose );
        cache[purpose] = Dense::Zero( matdim, matdim );
        cache[purpose + "_time"] = Dense::Zero( matdim, matdim );
        auto &gmat_time = cache[purpose + "_time"];
// Fill Time Matrix
#pragma omp parallel for collapse( 2 ) schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
        for ( size_t i = 0; i < matdim; i++ ) {
            for ( size_t j = 0; j < matdim; j++ ) {
                double tau = i + j < matdim ? s.parameters.grid_values[i + j] : s.parameters.grid_values.back() + s.parameters.grid_steps.back() * ( i + j - matdim + 1 );
                gmat_time( i, j ) = s.parameters.grid_values[i] + 1.0i * tau;
            }
        }
    }
    // Create Timer and Progressbar
    std::string progressstring = "G1(" + super_purpose + "): ";
    Timer &timer = Timers::create( "RungeKutta-G1-Loop (" + super_purpose + ")" );
    ProgressBar progressbar = ProgressBar();
    timer.start();

    if ( eval_operators.empty() )
        return;
    // Calculate G1 Function
    Log::L2( "[G1Correlation] Calculating G1(tau)... purpose: {}, saving to matrix of size {}x{}, iterating over {} saved states...\n", super_purpose, matdim, matdim, savedStates.size() );
    // Main G1 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t i = 0; i < std::min<size_t>( matdim, savedStates.size() ); i++ ) {
        std::vector<QDLC::SaveState> savedRhos;
        // Get Time from saved State
        double t_t = get_time_at( i );
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau( get_rho_at( i ), op_annihilator, t_t );
        // Calculate Runge Kutta
        calculate_runge_kutta( rho_tau, t_t, s.parameters.t_end, timer, progressbar, progressstring, s, savedRhos, false );
        // Interpolate saved states to equidistant timestep
        savedRhos = Numerics::interpolate_curve( savedRhos, t_t, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        for ( const auto &[eval, purpose] : eval_operators ) {
            auto &gmat = cache[purpose];
            for ( size_t j = 0; j < savedRhos.size() and i + j < std::min<size_t>( matdim, savedRhos.size() ); j++ ) {
                const double t_tau = savedRhos.at( j ).t;
                gmat( i, j ) = s.dgl_expectationvalue<Sparse, Scalar>( savedRhos.at( j ).mat, eval, t_tau );
            }
        }
        Timers::outputProgress( timer, progressbar, i, savedStates.size(), progressstring );
    }

    timer.end();
    Timers::outputProgress( timer, progressbar, savedStates.size(), savedStates.size(), progressstring, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[G1Correlation] G1 ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", super_purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );

    // Manually Apply the detector function
    for ( const auto &[eval, purpose] : eval_operators ) {
        auto &gmat = cache[purpose];
        const auto &gmat_time = cache[purpose + "_time"];
        apply_detector_function( s, gmat, gmat_time, purpose );
    }
}

void QDLC::Numerics::ODESolver::calculate_g2( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[G2Correlation] G2(tau) for {} already exists.\n", purpose );
    }
    calculate_g2( s, s_op_creator_1, std::vector<std::string>{ s_op_annihilator_1 }, { s_op_creator_2 }, s_op_annihilator_2, { purpose } );
}

void QDLC::Numerics::ODESolver::calculate_g2( System &s, const std::string &s_op_creator_1, const std::vector<std::string> &s_op_annihilator_1, const std::vector<std::string> &s_op_creator_2, const std::string &s_op_annihilator_2, const std::vector<std::string> &purposes ) {
    Log::L2( "[G2Correlation] Preparing to calculate G2 Correlation function\n" );
    Log::L2( "[G2Correlation] Generating Sparse Operator Matrices from String input...\n" );

    // Find Operator Matrices
    const auto op_creator_1 = get_operators_matrix( s, s_op_creator_1 );
    const auto op_annihilator_2 = get_operators_matrix( s, s_op_annihilator_2 );

    // Matrix Dimension
    const size_t matdim = s.parameters.grid_values.size(); // int( savedStates.size() / s.parameters.iterations_t_skip );

    // Generator super-purpose
    std::string super_purpose = std::accumulate( std::next( purposes.begin() ), purposes.end(), purposes.front(), []( const std::string &a, const std::string &b ) { return a + " and " + b; } );

    std::vector<std::pair<Sparse, std::string>> eval_operators;

    // Preconstruct
    for ( auto current = 0; current < s_op_annihilator_1.size(); current++ ) {
        const auto &op_creator_2 = get_operators_matrix( s, s_op_creator_2[current] );
        const auto &op_annihilator_1 = get_operators_matrix( s, s_op_annihilator_1[current] );
        const auto &purpose = purposes[current];
        // Cancel if Purpose already exists
        if ( cache.contains( purpose ) ) {
            Log::L2( "[G2Correlation] Matrix for {} already exists! Skipping!\n", purpose );
            continue;
        }
        // Construct Evaluation Operators
        eval_operators.emplace_back( op_creator_2 * op_annihilator_1, purpose );

        Log::L2( "[G2Correlation] Preparing Cache Matrices for {}...\n", purpose );
        cache[purpose] = Dense::Zero( matdim, matdim );
        cache[purpose + "_time"] = Dense::Zero( matdim, matdim );
        auto &gmat_time = cache[purpose + "_time"];
        // Fill Time Matrix
#pragma omp parallel for collapse( 2 ) schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
        for ( size_t i = 0; i < matdim; i++ ) {
            for ( size_t j = 0; j < matdim; j++ ) {
                double tau = i + j < matdim ? s.parameters.grid_values[i + j] : s.parameters.grid_values.back() + s.parameters.grid_steps.back() * ( i + j - matdim + 1 );
                gmat_time( i, j ) = s.parameters.grid_values[i] + 1.0i * tau;
            }
        }
    }
    if ( eval_operators.empty() )
        return;
    // Create Timer and Progresbar
    std::string progressstring = "G2(" + super_purpose + "): ";
    Timer &timer = Timers::create( "RungeKutta-G2-Loop (" + super_purpose + ")" );
    auto progressbar = ProgressBar();
    timer.start();

    // Calculate G2 Function
    Log::L2( "[G2Correlation] Calculating G2(tau)... purpose: {}, saving to matrix of size {}x{},  iterating over {} saved states...\n", super_purpose, matdim, matdim, std::min<size_t>( matdim, savedStates.size() ) );
    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t i = 0; i < std::min<size_t>( matdim, savedStates.size() ); i++ ) {
        // Create and reserve past rho's vector
        std::vector<QDLC::SaveState> savedRhos;
        // Get Time from saved State
        double t_t = get_time_at( i );
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau_2( get_rho_at( i ), op_annihilator_2, op_creator_1, t_t );
        // Calculate Runge Kutta
        calculate_runge_kutta( rho_tau, t_t, s.parameters.t_end, timer, progressbar, progressstring, s, savedRhos, false );
        // Interpolate saved states to equidistant timestep
        savedRhos = Numerics::interpolate_curve( savedRhos, t_t, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
        for ( const auto &[eval, purpose] : eval_operators ) {
            auto &gmat = cache[purpose];
            for ( size_t j = 0; j < savedRhos.size() and i + j < std::min<size_t>( matdim, savedRhos.size() ); j++ ) {
                const double t_tau = savedRhos.at( j ).t;
                gmat( i, j ) = s.dgl_expectationvalue<Sparse, Scalar>( savedRhos.at( j ).mat, eval, t_tau );
            }
        }
        Timers::outputProgress( timer, progressbar, i, savedStates.size(), progressstring );
    }

    timer.end();
    Timers::outputProgress( timer, progressbar, savedStates.size(), savedStates.size(), progressstring, Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[G2Correlation] G2 ({}) Hamilton Statistics: Attempts w/r: {}, Write: {}, Read: {}, Calc: {}.\n", super_purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );

    // Manually Apply the detector function
    for ( const auto &[eval, purpose] : eval_operators ) {
        auto &gmat = cache[purpose];
        const auto &gmat_time = cache[purpose + "_time"];
        apply_detector_function( s, gmat, gmat_time, purpose );
    }
}