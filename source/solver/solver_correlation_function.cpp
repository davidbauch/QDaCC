#include "solver/solver_ode.h"

bool ODESolver::scale_grid( System &s, Dense &cache, std::vector<std::vector<SaveScalar>> &cache_noneq ) {
    Log::L2( "Scaling grid... " );
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
    for ( long unsigned int i = 0; i < cache_noneq.size(); i++ ) {
        // Interpolant for real and imag of cache.at(i)
        // fill cache(i,j) with evaluated interpolant
        std::vector<double> real, imag, times;
        Interpolant interpolant_real, interpolant_imag;
        real.reserve( cache_noneq.size() );
        imag.reserve( cache_noneq.size() );
        times.reserve( cache_noneq.size() );
        for ( long unsigned int j = 0; j < cache_noneq.size(); j++ ) {
            if ( i < cache_noneq.at( j ).size() ) {
                real.emplace_back( std::real( cache_noneq.at( j ).at( i ).scalar ) );
                imag.emplace_back( std::imag( cache_noneq.at( j ).at( i ).scalar ) );
                times.emplace_back( cache_noneq.at( j ).at( i ).t );
            } else {
                break;
            }
        }

        if ( real.size() > 1 ) {
            interpolant_real = Interpolant( times, real, "linear" );
            interpolant_imag = Interpolant( times, imag, "linear" );
        }

        Scalar scalar;
        for ( long unsigned int i = 0; i < dim; i++ ) {
            for ( long unsigned int j = 0; j < dim; j++ ) {
                if ( i + j >= dim ) {
                    break;
                }
                double t = ( (double)j ) * s.parameters.t_step;
                if ( real.size() > 1 )
                    scalar = Scalar( interpolant_real.evaluate( t ), interpolant_imag.evaluate( t ) );
                else
                    scalar = Scalar( real.front(), imag.front() );
                cache( j, i ) = scalar;
            }
        }
    }
    Log::L2( "Done!\n" );
    return true;
}

// Description: Calculates the G1(tau) function. Uses akf_mat temporary variable to save the tau-direction expectation values. Calculates <b^+(t) * b(t+tau)> via quantum regression theorem. Logs and outputs progress.
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator: [&Sparse] Creator operator (adjunct of annihilator)
// @param op_annihilator: [&Sparse] Annihilator operator
// @return: [bool] True if calculations were sucessfull, else false
std::tuple<Sparse, Sparse> ODESolver::calculate_g1( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator, std::string purpose ) {
    // Find Operator Matrices
    Log::L2( " : Preparing to calculate G1 Correlation function\n : Generating Sparse Operator Matrices from String input...\n" );
    Sparse op_creator = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_annihilator = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    for ( auto &split_s_op_creator : splitline( s_op_creator, '+' ) )
        op_creator += s.operatorMatrices.el_transitions.count( split_s_op_creator ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_creator].hilbert : s.operatorMatrices.ph_transitions[split_s_op_creator].hilbert;
    for ( auto &split_s_op_annihilator : splitline( s_op_annihilator, '+' ) )
        op_annihilator += s.operatorMatrices.el_transitions.count( split_s_op_annihilator ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_annihilator].hilbert : s.operatorMatrices.ph_transitions[split_s_op_annihilator].hilbert;

    if ( cache.count( purpose ) != 0 ) {
        Log::L2( " : G1(tau) for {} already exists.\n", purpose );
        return {op_creator, op_annihilator};
    }

    // Generate Cache Matrices
    Log::L2( " : Preparing Cache Matrices...\n" );
    cache[purpose] = Dense::Zero( dim, dim );
    Dense &gmat = cache[purpose];

    Timer &timer = Timers::create( "RungeKutta-G1-Loop (" + purpose + ")" );
    int totalIterations = getIterationNumberTau( s );
    ProgressBar progressbar = ProgressBar( totalIterations );
    timer.start();
    Log::L2( " : Calculating G1(tau)... purpose: {}, saving to matrix of size {}x{}, maximum iterations: {}...\n", purpose, dim, dim, totalIterations );
    std::string progressstring = "G1(" + purpose + "): ";
    // Reserve Cache Matrix/Vector
    Log::L2( " : Preparing cache vector...\n" );
    std::vector<std::vector<SaveScalar>> cache_noneq;
    cache_noneq.reserve( savedStates.size() );
    for ( long unsigned int j = 0; j < savedStates.size(); j += s.parameters.iterations_t_skip ) {
        cache_noneq.emplace_back( std::vector<SaveScalar>() );
    }
    Log::L2( " : Done! Calculating g1(t,tau)...\n" );
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( long unsigned int i = 0; i < savedStates.size(); i += s.parameters.iterations_t_skip ) {
        std::vector<SaveState> savedRhos;
        // Get Time from saved State
        double t_t = getTimeAt( i );
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau( getRhoAt( i ), op_annihilator, t_t );
        // Calculate Runge Kutta
        calculate_runge_kutta( rho_tau, t_t, s.parameters.t_end, s.parameters.t_step, timer, progressbar, progressstring, s, savedRhos, false );
        for ( long unsigned int j = 0; j < savedRhos.size(); j++ ) {
            double t_tau = savedRhos.at( j ).t;
            cache_noneq.at( i / s.parameters.iterations_t_skip ).emplace_back( SaveScalar( s.dgl_expectationvalue<Sparse, Scalar>( savedRhos.at( j ).mat, op_creator, t_tau ), t_t, t_tau ) );
        }
        Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring );
    }
    // Because G1,G2 are better to be calculated on an equidistant grid, we interpolate the values to fit the grid.
    // The grid is determined by the initial timestep chosen (+ skip step)
    if ( s.parameters.numerics_stretch_correlation_grid && s.parameters.iterations_t_skip > 1 ) {
        scale_grid( s, gmat, cache_noneq );
    } else {
        for ( long unsigned int i = 0; i < cache_noneq.size(); i++ ) {
            for ( long unsigned int j = 0; j < cache_noneq.at( i ).size(); j++ ) {
                gmat( i * s.parameters.iterations_t_skip, j ) = cache_noneq.at( i ).at( j ).scalar;
            }
        }
    }

    Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring, PROGRESS_FORCE_OUTPUT );
    timer.end();
    Log::L2( ": Done! G1 ({}): Attempts w/r: {}, Write: {}, Read: {}, Calc: {}. Done!\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );
    return {op_creator, op_annihilator};
}

std::tuple<Sparse, Sparse, Sparse, Sparse> ODESolver::calculate_g2( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2, std::string purpose ) {
    // Find Operator Matrices
    Log::L2( " : Preparing to calculate G1 Correlation function\n : Generating Sparse Operator Matrices from String input...\n" );
    Sparse op_creator_1 = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_creator_2 = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_annihilator_1 = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_annihilator_2 = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    for ( auto &split_s_op_creator_1 : splitline( s_op_creator_1, '+' ) )
        op_creator_1 += s.operatorMatrices.el_transitions.count( split_s_op_creator_1 ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_creator_1].hilbert : s.operatorMatrices.ph_transitions[split_s_op_creator_1].hilbert;
    for ( auto &split_s_op_creator_2 : splitline( s_op_creator_2, '+' ) )
        op_creator_2 += s.operatorMatrices.el_transitions.count( split_s_op_creator_2 ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_creator_2].hilbert : s.operatorMatrices.ph_transitions[split_s_op_creator_2].hilbert;
    for ( auto &split_s_op_annihilator_1 : splitline( s_op_annihilator_1, '+' ) )
        op_annihilator_1 += s.operatorMatrices.el_transitions.count( split_s_op_annihilator_1 ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_annihilator_1].hilbert : s.operatorMatrices.ph_transitions[split_s_op_annihilator_1].hilbert;
    for ( auto &split_s_op_annihilator_2 : splitline( s_op_annihilator_2, '+' ) )
        op_annihilator_2 += s.operatorMatrices.el_transitions.count( split_s_op_annihilator_2 ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_annihilator_2].hilbert : s.operatorMatrices.ph_transitions[split_s_op_annihilator_2].hilbert;

    if ( cache.count( purpose ) != 0 ) {
        Log::L2( " : G2(tau) for {} already exists.\n", purpose );
        return {op_creator_1, op_annihilator_1, op_creator_2, op_annihilator_2};
    }

    Log::L2( " : Preparing Cache Matrices...\n" );
    cache[purpose] = Dense::Zero( dim, dim );
    Dense &gmat = cache[purpose];

    // Create Timer and Progresbar
    Timer &timer = Timers::create( "RungeKutta-G2-Loop (" + purpose + ")" );
    int totalIterations = getIterationNumberTau( s );
    ProgressBar progressbar = ProgressBar( totalIterations );
    timer.start();
    Log::L2( ": Calculating G2(tau)... purpose: {}, saving to matrix of size {}x{}...\n", purpose, dim, dim ); ////cache.size(), cache.at(0).size()
    Sparse evalOperator = op_creator_2 * op_annihilator_1;
    std::string progressstring = "G2(" + purpose + "): ";
    // Reserve Cache Matrix/Vector
    std::vector<std::vector<SaveScalar>> cache_noneq;
    cache_noneq.reserve( savedStates.size() );
    for ( long unsigned int j = 0; j < savedStates.size(); j += s.parameters.iterations_t_skip ) {
        cache_noneq.emplace_back( std::vector<SaveScalar>() );
    }
    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( long unsigned int i = 0; i < savedStates.size(); i += s.parameters.iterations_t_skip ) {
        // Create and reserve past rho's vector
        std::vector<SaveState> savedRhos;
        // Get Time from saved State
        double t_t = getTimeAt( i );
        // Calculate New Modified Density Matrix
        Sparse rho_tau = s.dgl_calc_rhotau_2( getRhoAt( i ), op_annihilator_2, op_creator_1, t_t );
        // Calculate Runge Kutta
        calculate_runge_kutta( rho_tau, t_t, s.parameters.t_end, s.parameters.t_step, timer, progressbar, progressstring, s, savedRhos, false );
        for ( long unsigned int j = 0; j < savedRhos.size(); j++ ) {
            double t_tau = savedRhos.at( j ).t;
            cache_noneq.at( i / s.parameters.iterations_t_skip ).emplace_back( s.dgl_expectationvalue<Sparse, Scalar>( savedRhos.at( j ).mat, evalOperator, t_tau ), t_t, t_tau );
        }
        Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring );
    }
    // Because G1,G2 are better to be calculated on an equidistant grid, we interpolate the values to fit the grid.
    // The grid is determined by the initial timestep chosen (+ skip step)
    if ( s.parameters.numerics_stretch_correlation_grid && s.parameters.iterations_t_skip > 1 ) {
        scale_grid( s, gmat, cache_noneq );
    } else {
        for ( long unsigned int i = 0; i < cache_noneq.size(); i++ ) {
            for ( long unsigned int j = 0; j < cache_noneq.at( i ).size(); j++ ) {
                gmat( i * s.parameters.iterations_t_skip, j ) = cache_noneq.at( i ).at( j ).scalar;
            }
        }
    }

    Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring, PROGRESS_FORCE_OUTPUT );
    timer.end();
    Log::L2( ": G2 ({}): Attempts w/r: {}, Write: {}, Read: {}, Calc: {}. Done!\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );
    return {op_creator_1, op_annihilator_1, op_creator_2, op_annihilator_2};
}

// Description: Calculates the G2(tau=0) function. Calculates <b^+(t) * b^+(t) * b(t) * b(t)> / <b^+(t) * b(t)>^2 . Logs and outputs progress. Saves resulting function.
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator: [&Sparse] Creator operator (adjunct of annihilator)
// @param op_annihilator: [&Sparse] Annihilator operator
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false
/*bool ODESolver::calculate_g2_0( System &s, const Sparse &op_creator, const Sparse &op_annihilator, std::string fileOutputName = "g2(0).txt" ) {
    if ( !( (int)savedStates.size() > 0 ) ) {
        Log::L1( "Need to calculate t-direction first!\n" );
        return false;
    }

    Timer &timer = Timers::create( "G2-0-Loop" );
    int totalIterations = (int)savedStates.size() / s.parameters.iterations_t_skip;
    ProgressBar progressbar = ProgressBar( totalIterations);
    timer.start();
    Log::L2( "Calculating G2(0)... " );

    std::vector<Scalar> g2Values;
    g2Values.reserve( totalIterations );
    for ( int i = 0; i < (int)savedStates.size(); i += s.parameters.iterations_t_skip ) {
        double t_t = getTimeAt( i );
        Sparse rho = getRhoAt( i );
        Sparse M1 = op_creator * op_creator * op_annihilator * op_annihilator;
        Sparse M2 = op_creator * op_annihilator;
        g2Values.emplace_back( s.dgl_expectationvalue<Sparse,Scalar>( rho, M1, t_t ) / std::pow( s.dgl_expectationvalue<Sparse,Scalar>( rho, M2, t_t ), 2 ) );
        timer.iterate();
        Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "G2: " );
    }
    timer.end();
    Log::L2( "Done, saving to {}... ", fileOutputName );
    std::string filepath = s.parameters.subfolder + fileOutputName;
    FILE *g2file = std::fopen( filepath.c_str(), "w" );
    if ( !g2file ) {
        Log::L2( "Failed to open outputfile for g2(0)!\n" );
        return false;
    }
    for ( int i = 0; i < totalIterations; i++ ) {
        fmt::print( g2file, "{0:.15e}\t{1:.15e}\n", getTimeAt( i ), std::real( g2Values.at( i ) ) );
    }
    std::fclose( g2file );
    Log::L2( "Done!\n" );
    return true;
}*/