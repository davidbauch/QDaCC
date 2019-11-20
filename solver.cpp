#pragma once
#include "global.h"
#include "solver.h"

// Description: ODESolver class provides both Runge-Kutta functions of different orders and functions for different numerical operations
// Type: ODESolver Class Constructor
// @param s: [&System] Class providing set of system functions
ODESolver::ODESolver( System_Parent &s ) {
    logs.level2( "Creating ODESolver Class... " );
    reset( s );
    savedStates.clear();
    savedStates.reserve( dim );
    logs.level2( "Done!\n" );
}

// Description: Gatheres a Hamiltonoperator using a systems get_hamilton() function. This workaround enables the saving of already calculated matrices for dublicate uses. Seems to have almost no influence on runtime.
// Type: ODESolver private function
// @param s: [&System] Class providing set of system functions
// @param t: [double] Time to return Hamiltonoperator at
// @return: [Eigent::MatrixXcd] hamilton matrix of type Eigen::MatrixXcd
MatType ODESolver::getHamilton( System_Parent &s, const double t, bool use_saved_hamiltons = false ) {
    if ( s.parameters.numerics_use_saved_hamiltons || use_saved_hamiltons ) {
        if ( savedHamiltons.size() == 0 || t > savedHamiltons.back().t ) {
            track_gethamilton_calcattempt++;
            // Calculate new Hamilton for t, save it, return it
            saveHamilton( s.dgl_getHamilton( t ), t );
            track_gethamilton_write++;
            return savedHamiltons.back().mat;
        } else {
            int approx = std::floor( t / s.getTimeStep() / 2.0 );
            // Look for corresponding Hamilton. t-direction has to be calculated first if used with multiple threads
            while ( t > savedHamiltons.at( approx ).t )
                approx++;
            while ( t < savedHamiltons.at( approx ).t )
                approx--;
            if ( approx < (int)savedHamiltons.size() ) {
                track_gethamilton_read++;
                return savedHamiltons.at( approx ).mat;
            }
        }
    }
    track_gethamilton_calc++;
    return s.dgl_getHamilton( t );
}

// Description: Iterates Runge-Kutta of order 4 at time t onto rho using the systems hamilton operator.
// Type: ODESolver private function
// @param rho: [&Eigen::MatrixXcd] Input (density-) matrix
// @param s: [&System] Class providing set of system functions
// @param t: [double] Time to iterate at
// @return: [Eigen::MatrixXcd] rho at time t+t_step
MatType ODESolver::iterateRungeKutta4( const MatType &rho, System_Parent &s, const double t, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-4 ausrechnen
    MatType H_calc_k1 = getHamilton( s, t );
    MatType H_calc_k23 = getHamilton( s, t + s.getTimeStep() * 0.5 );
    MatType H_calc_k4 = getHamilton( s, t + s.getTimeStep() );
    // k1-4 ausrechnen
    MatType rk1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    MatType rk2 = s.dgl_rungeFunction( rho + s.getTimeStep() * 0.5 * rk1, H_calc_k23, t + s.getTimeStep() * 0.5, savedStates );
    MatType rk3 = s.dgl_rungeFunction( rho + s.getTimeStep() * 0.5 * rk2, H_calc_k23, t + s.getTimeStep() * 0.5, savedStates );
    MatType rk4 = s.dgl_rungeFunction( rho + s.getTimeStep() * rk3, H_calc_k4, t + s.getTimeStep(), savedStates );
    // Dichtematrix
    return rho + s.getTimeStep() / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
}

// Description: Iterates Runge-Kutta of order 5 at time t onto rho using the systems hamilton operator.
// Type: ODESolver private function
// @param rho: [&Eigen::MatrixXcd] Input (density-) matrix
// @param s: [&System] Class providing set of system functions
// @param t: [double] Time to iterate at
// @return: [Eigen::MatrixXcd] rho at time t+t_step
MatType ODESolver::iterateRungeKutta5( const MatType &rho, System_Parent &s, const double t, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-6 ausrechnen
    MatType H_calc_k1 = getHamilton( s, t );
    MatType H_calc_k2 = getHamilton( s, t + a2 * s.getTimeStep() );
    MatType H_calc_k3 = getHamilton( s, t + a3 * s.getTimeStep() );
    MatType H_calc_k4 = getHamilton( s, t + a4 * s.getTimeStep() );
    MatType H_calc_k5 = getHamilton( s, t + a5 * s.getTimeStep() );
    MatType H_calc_k6 = getHamilton( s, t + a6 * s.getTimeStep() );
    // k1-6 ausrechnen
    MatType k1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    MatType k2 = s.dgl_rungeFunction( rho + s.getTimeStep() * b11 * k1, H_calc_k2, t + a2 * s.getTimeStep(), savedStates );
    MatType k3 = s.dgl_rungeFunction( rho + s.getTimeStep() * ( b21 * k1 + b22 * k2 ), H_calc_k3, t + a3 * s.getTimeStep(), savedStates );
    MatType k4 = s.dgl_rungeFunction( rho + s.getTimeStep() * ( b31 * k1 + b32 * k2 + b33 * k3 ), H_calc_k4, t + a4 * s.getTimeStep(), savedStates );
    MatType k5 = s.dgl_rungeFunction( rho + s.getTimeStep() * ( b41 * k1 + b42 * k2 + b43 * k3 + b44 * k4 ), H_calc_k5, t + a5 * s.getTimeStep(), savedStates );
    MatType k6 = s.dgl_rungeFunction( rho + s.getTimeStep() * ( b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4 + b55 * k5 ), H_calc_k6, t + a6 * s.getTimeStep(), savedStates );
    // Dichtematrix
    return rho + s.getTimeStep() * ( b61 * k1 + b63 * k3 + b64 * k4 + b65 * k5 + b66 * k6 );
}

// Description: Iterates Runge-Kutta with given order depending on the systems settings.
// Type: ODESolver public function
// @param rho: [&Eigen::MatrixXcd] Input (density-) matrix
// @param s: [&System] Class providing set of system functions
// @param t: [double] Time to iterate at
// @param dir: [int] Time direction. Solver order can differ in both t and tau direction, depending on the systems settings. Default value is DIR_T
// @return: [Eigen::MatrixXcd] rho at time t+t_step
MatType ODESolver::iterate( const MatType &rho, System_Parent &s, const double t, std::vector<SaveState> &savedStates, const int dir = DIR_T ) {
    int order = s.getSolverRungeKuttaOrder( dir );
    if ( order == 4 ) {
        return iterateRungeKutta4( rho, s, t, savedStates );
    }
    return iterateRungeKutta5( rho, s, t, savedStates );
}

// Description: Saves a tuple of a complex (density-)matrix and time, ensuring times and matrices don't get mixed up
// Type: ODESolver private function
// @param mat: [&Eigen::MatrixXcd] Matrix to save
// @param t: [double] Corresponding time
// @return: [void]
void ODESolver::saveState( const MatType &mat, const double t, std::vector<SaveState> &savedStates ) {
    savedStates.emplace_back( SaveState( mat, t ) );
}

// Description: Saves a tuple of a complex (Hamilton-)matrix and time, ensuring times and matrices don't get mixed up
// Type: ODESolver private function
// @param mat: [&Eigen::MatrixXcd] Matrix to save
// @param t: [double] Corresponding time
// @return: [void]
void ODESolver::saveHamilton( const MatType &mat, const double t ) {
    savedHamiltons.emplace_back( SaveState( mat, t ) );
}

// Description: Checks wether or not to save the current matrix for later calculations
// Type: ODESolver private function
// @param s: [&System] Class providing set of system functions
// @param curIt: [&int] Iterator variable to use for iteration check
// @return: [bool] True if current matrix should be saved, else false. Saving depends on the iteration skip for tau direction and wether or not they are needed later.
bool ODESolver::queueNow( System_Parent &s, int &curIt ) {
    if ( ( s.calculate_spectrum() || s.calculate_g2() ) && curIt % s.getIterationSkip() == 0 ) {
        curIt = 1;
        return true;
    }
    curIt++;
    return false;
}

// Description: Resets all temporary variables. Currently: out, akf_mat, track variables
// Type: ODESolver private function
// @param s: [&System] Class providing set of system functions
// @return: [int] dimensions of temporary variables
int ODESolver::reset( System_Parent &s ) {
    track_gethamilton_read = 0;
    track_gethamilton_write = 0;
    track_gethamilton_calc = 0;
    track_gethamilton_calcattempt = 0;
    dim = (int)( s.parameters.iterations_t_max / s.getIterationSkip() ) + 10;
    out.clear();
    out.reserve( (int)( std::ceil( s.parameters.iterations_w_resolution ) + 5 ) );
    for ( int w = 0; w < s.parameters.iterations_w_resolution; w++ ) {
        out.push_back( 0 );
    }
    return dim;
}

// Description: Calculates the normal t-direction via solving the von-Neumann equation for rho. May save some of the density matrices for later uses. Logs the calculation and outputs progress.
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @return: [bool] True if calculations are sucessfull, else false
bool ODESolver::calculate_t_direction( System_Parent &s ) {
    Timer &rkTimer = createTimer( "RungeKutta-Main-Loop" );
    ProgressBar progressbar = ProgressBar( s.parameters.iterations_t_max, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    rkTimer.start();

    logs.level2( "Calculating t-direction from {} to {} at stepsize {}... ", s.getTimeborderStart(), s.getTimeborderEnd(), s.getTimeStep() );
    MatType rho = s.getRho0();
    saveState( rho, s.getTimeborderStart(), savedStates );
    s.expectationValues( rho, s.getTimeborderStart() );

    // Main Time Loop
    for ( double t_t = s.getTimeborderStart() + s.getTimeStep(); t_t < s.getTimeborderEnd(); t_t += s.getTimeStep() ) {
        // Runge-Kutta iteration
        rho = iterate( rho, s, t_t, savedStates );
        // Save Rho for tau-direction
        saveState( rho, t_t, savedStates );
        // Expectation Values
        s.expectationValues( rho, t_t );
        // Progress and time output
        rkTimer.iterate();
        outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, "T-Direction: " );
        // Divergent bruder?
        if ( !s.traceValid( rho, t_t ) ) {
            t_t = s.getTimeborderEnd() + s.getTimeStep();
            return false;
        }
    }
    rkTimer.end();
    logs.level2( "Done! Saved {} states, cached {} hamilton matrices.\n", savedStates.size(), savedHamiltons.size() );
    return true;
}
template <typename T>
MatType ODESolver::iterate_definite_integral( const MatType &rho, T rungefunction, const double t, const double step ) {
    MatType rk1 = rungefunction( rho, t );
    MatType rk2 = rungefunction( rho + step * 0.5 * rk1, t + step * 0.5 );
    MatType rk3 = rungefunction( rho + step * 0.5 * rk2, t + step * 0.5 );
    MatType rk4 = rungefunction( rho + step * rk3, t + step );
    // Dichtematrix
    return rho + step / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
}

// Description: Integrates rho from t0 to t1 via Rungefunc.
// Type: ODESolver public static function
// @return: [vector<SaveState>] Vector of save state tuples (matrix, time)
std::vector<SaveState> ODESolver::calculate_definite_integral_vec( MatType rho, std::function<MatType( const MatType &, const double )> const &rungefunction, const double t0, const double t1, const double step ) { //std::function<MatrixXcd( const MatrixXcd &, const double )>
    std::vector<SaveState> ret;
    //ret.reserve( std::ceil( std::abs( t1 - t0 ) / std::abs( step ) ) );
    ret.emplace_back( SaveState( rho, t0 ) );
    if ( step > 0 )
        for ( double t = t0; t < t1; t += step ) {
            rho = iterate_definite_integral( rho, rungefunction, t, step );
            ret.emplace_back( SaveState( rho, t ) );
        }
    else if ( step < 0 )
        for ( double t = t0; t > t1; t += step ) {
            rho = iterate_definite_integral( rho, rungefunction, t, step );
            ret.emplace_back( SaveState( rho, t ) );
        }
    return ret;
}

// Description: Integrates rho from t0 to t1 via Rungefunc.
// Type: ODESolver public static function
// @return: [SaveState] Save state tuple (matrix, time)
SaveState ODESolver::calculate_definite_integral( MatType rho, std::function<MatType( const MatType &, const double )> const &rungefunction, const double t0, const double t1, const double step ) { //std::function<MatrixXcd( const MatrixXcd &, const double )>
    //ret.reserve( std::ceil( std::abs( t1 - t0 ) / std::abs( step ) ) );
    if ( step > 0 )
        for ( double t = t0; t < t1; t += step ) {
            rho = iterate_definite_integral( rho, rungefunction, t, step );
        }
    else if ( step < 0 )
        for ( double t = t0; t > t1; t += step ) {
            rho = iterate_definite_integral( rho, rungefunction, t, step );
        }
    return SaveState( rho, t1 );
}

// Desciption: Function to calculate the number of iterations used for tau direction calculations
// Type: ODESolver private function
// @param s: [&System] Class providing set of system functions
// @return: [int] Number of tau-direction iterations
int ODESolver::getIterationNumberTau( System_Parent &s ) {
    int num = 0;
    // Tau Direction Iteration steps
    for ( int i = 0; i < (int)savedStates.size(); i += s.getIterationSkip() ) {
        double t_t = getTimeAt( i );
        for ( double t_tau = t_t + s.getTimeStep(); t_tau < s.getTimeborderEnd(); t_tau += s.getTimeStep() ) { // t + +s.getTimeStep()
            num++;
        }
    }
    return num;
}

// Description: Function to calculate the number of iterations used for spectru calculations
// Type: ODESolver private function
// @param s: [&System] Class providing set of system functions
// @return: [int] Number of spectrum iterations
int ODESolver::getIterationNumberSpectrum( System_Parent &s ) {
    int num = 0;
    // Spectrum steps
    for ( int spec_w = 0; spec_w < s.parameters.iterations_w_resolution; spec_w++ ) {
        num++;
    }
    return num;
}

// Description: Function to extract the time corresponding to a certain iteration of saved states
// Type: ODESolver private function
// @param i: [int] Iteration number
// @return: [double] Time corresponding to iteration number i
double ODESolver::getTimeAt( int i ) {
    return savedStates.at( i ).t;
}

// Description: Function to extract the matrix corresponding to a certain iteration of saved states
// Type: ODESolver private function
// @param i: [int] Iteration number
// @return: [Eigen::MatrixXcd] (Density-) Matrix corresponding to iteration number i
MatType ODESolver::getRhoAt( int i ) {
    return savedStates.at( i ).mat;
}

// Description: Calculates the G1(tau) function. Uses akf_mat temporary variable to save the tau-direction expectation values. Calculates <b^+(t) * b(t+tau)> via quantum regression theorem. Logs and outputs progress.
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator: [&Eigen::MatrixXcd] Creator operator (adjunct of annihilator)
// @param op_annihilator: [&Eigen::MatrixXcd] Annihilator operator
// @return: [bool] True if calculations were sucessfull, else false
bool ODESolver::calculate_g1( System_Parent &s, const MatType &op_creator, const MatType &op_annihilator, DenseMat &cache, std::string purpose ) {
    if ( !( (int)savedStates.size() > 0 ) ) {
        logs( "Need to calculate t-direction first!\n" );
        return false;
    }
    reset( s );
    Timer &timer = createTimer( "RungeKutta-G1-Loop (" + purpose + ")" );
    int totalIterations = getIterationNumberTau( s );
    ProgressBar progressbar = ProgressBar( totalIterations, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    progressbar.strBarStart = "|";
    progressbar.strBarEnd = "|";
    timer.start();
    logs.level2( "Calculating G1(tau)... purpose: {}, saving to matrix of size {}x{}... ", purpose, cache.rows(), cache.cols() );
    std::string progressstring = "G1("+purpose+"): ";
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < (int)savedStates.size(); i += s.getIterationSkip() ) {
        std::vector<SaveState> past_rhos;
        past_rhos.reserve( (int)( ( savedStates.size() - i ) / s.getIterationSkip() ) );
        int k = i / s.getIterationSkip();
        double t_t = getTimeAt( i );
        MatType rho_tau = s.dgl_calc_rhotau( getRhoAt( i ), op_annihilator, t_t );
        saveState( rho_tau, t_t, past_rhos );

        cache( k, 0 ) = s.dgl_expectationvalue<MatType, dcomplex>( rho_tau, op_creator, t_t );
        int j = 1;
        int curIt_tau = 1;
        for ( double t_tau = t_t + s.getTimeStep(); t_tau < s.getTimeborderEnd(); t_tau += s.getTimeStep() ) { // t + +s.getTimeStep()
            rho_tau = iterate( rho_tau, s, t_tau, past_rhos, DIR_TAU );
            saveState( rho_tau, t_tau, past_rhos );
            timer.iterate();
            if ( queueNow( s, curIt_tau ) ) {
                cache( k, j ) = s.dgl_expectationvalue<MatType, dcomplex>( rho_tau, op_creator, t_tau );
                j++; // equivalent to s.parameters.akf_vecIndex
            }
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring );
        }
    }
    timer.end();
    logs.level2( "G1 ({}): Attempts w/r: {}, Write: {}, Read: {}, Calc: {}. Done!\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );
    return true;
}

// Description: Calculates the G2(tau=0) function. Calculates <b^+(t) * b^+(t) * b(t) * b(t)> / <b^+(t) * b(t)>^2 . Logs and outputs progress. Saves resulting function.
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator: [&Eigen::MatrixXcd] Creator operator (adjunct of annihilator)
// @param op_annihilator: [&Eigen::MatrixXcd] Annihilator operator
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false
/*bool ODESolver::calculate_g2_0( System_Parent &s, const MatType &op_creator, const MatType &op_annihilator, std::string fileOutputName = "g2(0).txt" ) {
    if ( !( (int)savedStates.size() > 0 ) ) {
        logs( "Need to calculate t-direction first!\n" );
        return false;
    }

    Timer &timer = createTimer( "G2-0-Loop" );
    int totalIterations = (int)savedStates.size() / s.getIterationSkip();
    ProgressBar progressbar = ProgressBar( totalIterations, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    progressbar.strBarStart = "|";
    progressbar.strBarEnd = "|";
    timer.start();
    logs.level2( "Calculating G2(0)... " );

    std::vector<dcomplex> g2Values;
    g2Values.reserve( totalIterations );
    for ( int i = 0; i < (int)savedStates.size(); i += s.getIterationSkip() ) {
        double t_t = getTimeAt( i );
        MatType rho = getRhoAt( i );
        MatType M1 = op_creator * op_creator * op_annihilator * op_annihilator;
        MatType M2 = op_creator * op_annihilator;
        g2Values.emplace_back( s.dgl_expectationvalue<MatType,dcomplex>( rho, M1, t_t ) / std::pow( s.dgl_expectationvalue<MatType,dcomplex>( rho, M2, t_t ), 2 ) );
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

// Description: Calculates the Eberly-Wódkiewicz spectrum. Uses precalculated values from akf_mat
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false
bool ODESolver::calculate_spectrum( System_Parent &s, const MatType &op_creator, const MatType &op_annihilator, std::string fileOutputName = "spectrum.txt" ) {
    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
    s.command(ODESolver::CHANGE_TO_SINGLETHREADED_SUBPROGRAM);
    // Cache matrix. Saves complex double entries of G1(t,tau)
    DenseMat akf_mat = DenseMat::Zero( dim, dim );
    // Calculate G1(t,tau) with given operator matrices
    calculate_g1( s, op_creator, op_annihilator, akf_mat, "spectrum" );
    // Create Timer and Progressbar for the spectrum loop
    Timer &timer = createTimer( "Spectrum-Loop" );
    int totalIterations = getIterationNumberSpectrum( s );
    ProgressBar progressbar = ProgressBar( totalIterations, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    progressbar.strBarStart = "|";
    progressbar.strBarEnd = "|";
    timer.start();
    logs.level2( "Calculating spectrum... Calculating frequencies... " );
    //Calculate frequencies:
    std::vector<double> spectrum_frequency_w;
    for ( int w = 0; w < s.parameters.iterations_w_resolution; w++ ) {
        spectrum_frequency_w.push_back( s.parameters.spectrum_frequency_center - ( s.parameters.spectrum_frequency_range ) + w / ( (double)s.parameters.iterations_w_resolution ) * ( 2. * ( s.parameters.spectrum_frequency_range ) ) );
    }
    logs.level2( "Done, calculating fourier transform via direct integral... " );
    // Calculate main fourier transform integral
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads ) //reduction(+:spectrum.out)
    for ( int spec_w = 0; spec_w < s.parameters.iterations_w_resolution; spec_w++ ) {
        std::vector<dcomplex> expfunc;
        expfunc.reserve( savedStates.size() / s.getIterationSkip() + 1 );
        for ( int spec_tau = 0; spec_tau < ( s.getTimeborderEnd() - s.getTimeborderStart() - 0.0 ) / ( s.getTimeStep() * s.getIterationSkip() ); spec_tau++ ) {
            expfunc.emplace_back( std::exp( -1i * spectrum_frequency_w.at( spec_w ) * (double)(spec_tau)*s.getTimeStep() * (double)( s.getIterationSkip() ) ) );
        }
        for ( int i = 0; i < (int)savedStates.size(); i += s.getIterationSkip() ) {
            int k = i / s.getIterationSkip();
            double t_t = getTimeAt( i );
            for ( int spec_tau = 0; spec_tau < ( s.getTimeborderEnd() - s.getTimeborderStart() - t_t ) / ( s.getTimeStep() * s.getIterationSkip() ); spec_tau++ ) {
                out.at( spec_w ) += expfunc.at( spec_tau ) * akf_mat( k, spec_tau );
            }
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "AKF-Spectrum: " );
        }
        timer.iterate();
    }
    // Final output and timer end
    outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "AKF-Spectrum", PROGRESS_FORCE_OUTPUT );
    timer.end();
    // Save output
    logs.level2( "Done, saving to {}... ", fileOutputName );
    std::string filepath = s.parameters.subfolder + fileOutputName;
    FILE *spectrumfile = std::fopen( filepath.c_str(), "w" );
    if ( !spectrumfile ) {
        logs.level2( "Failed to open outputfile for spectrum!\n" );
        return false;
    }
    for ( int spec_w = 0; spec_w < s.parameters.iterations_w_resolution; spec_w++ ) {
        fmt::print( spectrumfile, "{0:.15e}\t{1:.15e}\n", spectrum_frequency_w[spec_w], real( out[spec_w] * s.getTimeStep() * s.getTimeStep() * (double)( s.getIterationSkip() * s.getIterationSkip() ) ) );
    }
    std::fclose( spectrumfile );
    logs.level2( "Done!\n" );
    // Sucessfull
    return true;
}

bool ODESolver::calculate_g2( System_Parent &s, const MatType &op_creator_1, const MatType &op_annihilator_1, const MatType &op_creator_2, const MatType &op_annihilator_2, DenseMat &cache, std::string purpose ) {
    // Ensuring T-Direction was calculated first
    if ( !( (int)savedStates.size() > 0 ) ) {
        logs( "Need to calculate t-direction first!\n" );
        return false;
    }
    // Create Timer and Progresbar
    Timer &timer = createTimer( "RungeKutta-G2-Loop (" + purpose + ")" );
    int totalIterations = getIterationNumberTau( s );
    ProgressBar progressbar = ProgressBar( totalIterations, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    progressbar.strBarStart = "|";
    progressbar.strBarEnd = "|";
    timer.start();
    logs.level2( "Calculating G2(tau)... purpose: {}, saving to matrix of size {}x{}... ", purpose, cache.rows(), cache.cols() );
    MatType evalOperator = op_creator_1*op_creator_2;
    std::string progressstring = "G2("+purpose+"): ";
    // Main G2 Loop
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < (int)savedStates.size(); i += s.getIterationSkip() ) {
        // Create and reserve past rho's vector
        std::vector<SaveState> past_rhos;
        past_rhos.reserve( (int)( ( savedStates.size() - i ) / s.getIterationSkip() ) );
        // Get index incrementing by 1 from global index
        int k = i / s.getIterationSkip();
        // Get Time from saved State
        double t_t = getTimeAt( i );
        // Calculate rho_tau
        MatType rho_tau = s.dgl_calc_rhotau_2( getRhoAt( i ), op_annihilator_2, op_creator_1, t_t );
        saveState( rho_tau, t_t, past_rhos );

        cache( k, 0 ) = s.dgl_expectationvalue<MatType, dcomplex>( rho_tau, evalOperator, t_t );
        int j = 1;
        int curIt_tau = 1;
        for ( double t_tau = t_t + s.getTimeStep(); t_tau < s.getTimeborderEnd(); t_tau += s.getTimeStep() ) { // t + +s.getTimeStep()
            rho_tau = iterate( rho_tau, s, t_tau, past_rhos, DIR_TAU );
            saveState( rho_tau, t_tau, past_rhos );
            timer.iterate();
            if ( queueNow( s, curIt_tau ) ) {
                cache( k, j ) = s.dgl_expectationvalue<MatType, dcomplex>( rho_tau, evalOperator, t_tau );
                j++; 
            }
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, progressstring );
        }
    }
    timer.end();
    logs.level2( "G2 ({}): Attempts w/r: {}, Write: {}, Read: {}, Calc: {}. Done!\n", purpose, track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_read, track_gethamilton_calc );
    return true;
}

bool ODESolver::calculate_advanced_photon_statistics( System_Parent &s, const MatType &op_creator_1, const MatType &op_annihilator_1, const MatType &op_creator_2, const MatType &op_annihilator_2, std::string fileOutputName ) {
    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
    s.command(ODESolver::CHANGE_TO_SINGLETHREADED_SUBPROGRAM);
    // Cache matrix. Saves complex double entries of G1(t,tau)
    DenseMat akf_mat_11 = DenseMat::Zero( dim, dim );
    DenseMat akf_mat_22 = DenseMat::Zero( dim, dim );
    DenseMat akf_mat_12 = DenseMat::Zero( dim, dim );
    // Calculate G1(t,tau) with given operator matrices
    calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_1, op_annihilator_1, akf_mat_11, "Mode 11" );
    calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_2, op_annihilator_2, akf_mat_22, "Mode 22" );
    calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_2, op_annihilator_2, akf_mat_12, "Mode 12" );
    // Create Timer and Progressbar for the spectrum loop
    Timer &timer = createTimer( "Advanced photon Statistics" );
    timer.start();
    logs.level2( "Calculating advanced photon statistics... Calculating Concurrence... " );
    
    logs.level2( "Done!\n" );
    // Sucessfull
    return true;
}