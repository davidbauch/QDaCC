#include "solver.h"

// Description: Calculates the G1(tau) function. Uses akf_mat temporary variable to save the tau-direction expectation values. Calculates <b^+(t) * b(t+tau)> via quantum regression theorem. Logs and outputs progress.
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator: [&MatType] Creator operator (adjunct of annihilator)
// @param op_annihilator: [&MatType] Annihilator operator
// @return: [bool] True if calculations were sucessfull, else false
bool ODESolver::calculate_g1( System &s, const MatType &op_creator, const MatType &op_annihilator, Dense &cache, std::string purpose ) {
    if ( !( (int)savedStates.size() > 0 ) ) {
        logs( "Need to calculate t-direction first!\n" );
        return false;
    }
    reset( s );
    Timer &timer = createTimer( "RungeKutta-G1-Loop (" + purpose + ")" );
    int totalIterations = getIterationNumberTau( s );
    ProgressBar progressbar = ProgressBar( totalIterations );
    timer.start();
    logs.level2( "Calculating G1(tau)... purpose: {}, saving to matrix of size {}x{}... ", purpose, cache.rows(), cache.cols() );
    std::string progressstring = "G1(" + purpose + "): ";
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < (int)savedStates.size(); i += s.parameters.iterations_t_skip ) {
        std::vector<SaveState> past_rhos;
        past_rhos.reserve( (int)( ( savedStates.size() - i ) / s.parameters.iterations_t_skip ) );
        int k = i / s.parameters.iterations_t_skip;
        double t_t = getTimeAt( i );
        MatType rho_tau = s.dgl_calc_rhotau( getRhoAt( i ), op_annihilator, t_t );
        saveState( rho_tau, t_t, past_rhos );

        cache( k, 0 ) = s.dgl_expectationvalue<MatType, Scalar>( rho_tau, op_creator, t_t );
        int j = 1;
        int curIt_tau = 1;
        for ( double t_tau = t_t + s.parameters.t_step; t_tau < s.parameters.t_end; t_tau += s.parameters.t_step ) { // t + +s.parameters.t_step
            rho_tau = iterate( rho_tau, s, t_tau, past_rhos, DIR_TAU );
            saveState( rho_tau, t_tau, past_rhos );
            timer.iterate();
            if ( queueNow( s, curIt_tau ) ) {
                cache( k, j ) = s.dgl_expectationvalue<MatType, Scalar>( rho_tau, op_creator, t_tau );
                j++; // equivalent to s.parameters.akf_vecIndex
            }
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring );
        }
    }
    outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring, PROGRESS_FORCE_OUTPUT );
    timer.end();
    logs.level2( "G1 ({}): Attempts w/r: {}, Write: {}, Read: {}, Calc: {}. Done!\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );
    return true;
}

bool ODESolver::calculate_g2( System &s, const MatType &op_creator_1, const MatType &op_annihilator_1, const MatType &op_creator_2, const MatType &op_annihilator_2, Dense &cache, std::string purpose ) {
    // Ensuring T-Direction was calculated first
    if ( !( (int)savedStates.size() > 0 ) ) {
        logs( "Need to calculate t-direction first!\n" );
        return false;
    }
    reset( s );
    // Create Timer and Progresbar
    Timer &timer = createTimer( "RungeKutta-G2-Loop (" + purpose + ")" );
    int totalIterations = getIterationNumberTau( s );
    ProgressBar progressbar = ProgressBar( totalIterations );
    timer.start();
    logs.level2( "Calculating G2(tau)... purpose: {}, saving to matrix of size {}x{}... ", purpose, cache.rows(), cache.cols() );
    MatType evalOperator = op_creator_2 * op_annihilator_1;
    std::string progressstring = "G2(" + purpose + "): ";
    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < (int)savedStates.size(); i += s.parameters.iterations_t_skip ) {
        // Create and reserve past rho's vector
        std::vector<SaveState> past_rhos;
        past_rhos.reserve( (int)( ( savedStates.size() - i ) / s.parameters.iterations_t_skip ) );
        // Get index incrementing by 1 from global index
        int k = i / s.parameters.iterations_t_skip;
        // Get Time from saved State
        double t_t = getTimeAt( i );
        // Calculate rho_tau
        MatType rho_tau = s.dgl_calc_rhotau_2( getRhoAt( i ), op_annihilator_2, op_creator_1, t_t );
        saveState( rho_tau, t_t, past_rhos );

        cache( k, 0 ) = s.dgl_expectationvalue<MatType, Scalar>( rho_tau, evalOperator, t_t );
        int j = 1;
        int curIt_tau = 1;
        for ( double t_tau = t_t + s.parameters.t_step; t_tau < s.parameters.t_end; t_tau += s.parameters.t_step ) { // t + +s.parameters.t_step
            rho_tau = iterate( rho_tau, s, t_tau, past_rhos, DIR_TAU );
            saveState( rho_tau, t_tau, past_rhos );
            timer.iterate();
            if ( queueNow( s, curIt_tau ) ) {
                cache( k, j ) = s.dgl_expectationvalue<MatType, Scalar>( rho_tau, evalOperator, t_tau );
                j++;
            }
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring );
        }
    }
    outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring, PROGRESS_FORCE_OUTPUT );
    timer.end();
    logs.level2( "G2 ({}): Attempts w/r: {}, Write: {}, Read: {}, Calc: {}. Done!\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );
    return true;
}

// Description: Calculates the G2(tau=0) function. Calculates <b^+(t) * b^+(t) * b(t) * b(t)> / <b^+(t) * b(t)>^2 . Logs and outputs progress. Saves resulting function.
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator: [&MatType] Creator operator (adjunct of annihilator)
// @param op_annihilator: [&MatType] Annihilator operator
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false
/*bool ODESolver::calculate_g2_0( System &s, const MatType &op_creator, const MatType &op_annihilator, std::string fileOutputName = "g2(0).txt" ) {
    if ( !( (int)savedStates.size() > 0 ) ) {
        logs( "Need to calculate t-direction first!\n" );
        return false;
    }

    Timer &timer = createTimer( "G2-0-Loop" );
    int totalIterations = (int)savedStates.size() / s.parameters.iterations_t_skip;
    ProgressBar progressbar = ProgressBar( totalIterations);
    timer.start();
    logs.level2( "Calculating G2(0)... " );

    std::vector<Scalar> g2Values;
    g2Values.reserve( totalIterations );
    for ( int i = 0; i < (int)savedStates.size(); i += s.parameters.iterations_t_skip ) {
        double t_t = getTimeAt( i );
        MatType rho = getRhoAt( i );
        MatType M1 = op_creator * op_creator * op_annihilator * op_annihilator;
        MatType M2 = op_creator * op_annihilator;
        g2Values.emplace_back( s.dgl_expectationvalue<MatType,Scalar>( rho, M1, t_t ) / std::pow( s.dgl_expectationvalue<MatType,Scalar>( rho, M2, t_t ), 2 ) );
        timer.iterate();
        outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "G2: " );
    }
    timer.end();
    logs.level2( "Done, saving to {}... ", fileOutputName );
    std::string filepath = s.parameters.subfolder + fileOutputName;
    FILE *g2file = std::fopen( filepath.c_str(), "w" );
    if ( !g2file ) {
        logs.level2( "Failed to open outputfile for g2(0)!\n" );
        return false;
    }
    for ( int i = 0; i < totalIterations; i++ ) {
        fmt::print( g2file, "{0:.15e}\t{1:.15e}\n", getTimeAt( i ), std::real( g2Values.at( i ) ) );
    }
    std::fclose( g2file );
    logs.level2( "Done!\n" );
    return true;
}*/