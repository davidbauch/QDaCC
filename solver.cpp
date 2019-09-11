#pragma once
#include "global.h"
#include "solver.h"

ODESolver::ODESolver( const System &s ) {
    logs.level2( "Creating ODESolver Class... " );
    int dim = reset( s );
    savedStates.clear();
    savedStates.reserve( dim );
    logs.level2( "Done!\n" );
}

// TODO: implement hamilton matrix saving (tuple t,H)
MatrixXcd ODESolver::getHamilton( const System &s, const double t ) {
    return s.dgl_getHamilton( t );
}

MatrixXcd ODESolver::iterateRungeKutta4( const MatrixXcd &rho, const System &s, const double t ) {
    /* Verschiedene H's fuer k1-4 ausrechnen */
    //double wallt = omp_get_wtime();
    MatrixXcd H_calc_k1 = getHamilton( s, t );
    MatrixXcd H_calc_k23 = getHamilton( s, t + s.getTimeStep() * 0.5 );
    MatrixXcd H_calc_k4 = getHamilton( s, t + s.getTimeStep() );
    //rk_Hamiltons.add(0, (omp_get_wtime()-wallt)/numThreads);
    /* k1-4 ausrechnen */
    //wallt = omp_get_wtime();
    MatrixXcd rk1 = s.dgl_rungeFunction( rho, H_calc_k1, t );
    MatrixXcd rk2 = s.dgl_rungeFunction( rho + s.getTimeStep() * 0.5 * rk1, H_calc_k23, t + s.getTimeStep() * 0.5 );
    MatrixXcd rk3 = s.dgl_rungeFunction( rho + s.getTimeStep() * 0.5 * rk2, H_calc_k23, t + s.getTimeStep() * 0.5 );
    MatrixXcd rk4 = s.dgl_rungeFunction( rho + s.getTimeStep() * rk3, H_calc_k4, t + s.getTimeStep() );
    //rk_kMatrices.add(0, (omp_get_wtime()-wallt)/numThreads);
    /* Dichtematrix */
    return rho + s.getTimeStep() / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
}

MatrixXcd ODESolver::iterateRungeKutta5( const MatrixXcd &rho, const System &s, const double t ) {
    // Verschiedene H's fuer k1-6 ausrechnen
    //double wallt = omp_get_wtime();
    MatrixXcd H_calc_k1 = getHamilton( s, t );
    MatrixXcd H_calc_k2 = getHamilton( s, t + a2 * s.getTimeStep() );
    MatrixXcd H_calc_k3 = getHamilton( s, t + a3 * s.getTimeStep() );
    MatrixXcd H_calc_k4 = getHamilton( s, t + a4 * s.getTimeStep() );
    MatrixXcd H_calc_k5 = getHamilton( s, t + a5 * s.getTimeStep() );
    MatrixXcd H_calc_k6 = getHamilton( s, t + a6 * s.getTimeStep() );
    //rk_Hamiltons.add(0, (omp_get_wtime()-wallt)/numThreads);
    // k1-4 ausrechnen
    //wallt = omp_get_wtime();
    MatrixXcd k1 = s.dgl_rungeFunction( rho, H_calc_k1, t );
    MatrixXcd k2 = s.dgl_rungeFunction( rho + s.getTimeStep() * b11 * k1, H_calc_k2, t + a2 * s.getTimeStep() );
    MatrixXcd k3 = s.dgl_rungeFunction( rho + s.getTimeStep() * ( b21 * k1 + b22 * k2 ), H_calc_k3, t + a3 * s.getTimeStep() );
    MatrixXcd k4 = s.dgl_rungeFunction( rho + s.getTimeStep() * ( b31 * k1 + b32 * k2 + b33 * k3 ), H_calc_k4, t + a4 * s.getTimeStep() );
    MatrixXcd k5 = s.dgl_rungeFunction( rho + s.getTimeStep() * ( b41 * k1 + b42 * k2 + b43 * k3 + b44 * k4 ), H_calc_k5, t + a5 * s.getTimeStep() );
    MatrixXcd k6 = s.dgl_rungeFunction( rho + s.getTimeStep() * ( b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4 + b55 * k5 ), H_calc_k6, t + a6 * s.getTimeStep() );
    //rk_kMatrices.add(0, (omp_get_wtime()-wallt)/numThreads);
    // Dichtematrix
    return rho + s.getTimeStep() * ( b61 * k1 + b63 * k3 + b64 * k4 + b65 * k5 + b66 * k6 );
}

MatrixXcd ODESolver::iterate( const MatrixXcd &rho, const System &s, const double t, const int dir = DIR_T ) {
    int order = s.getSolverRungeKuttaOrder( dir );
    if ( order == 4 ) {
        return iterateRungeKutta4( rho, s, t );
    } //else if (order == 5) {
    return iterateRungeKutta5( rho, s, t );
    //}
}

void ODESolver::saveState( const MatrixXcd &mat, const double t ) {
    savedStates.emplace_back( SaveState( mat, t ) );
}

bool ODESolver::queueNow( const System &s, int &curIt ) {
    if ( ( s.calculate_spectrum() || s.calculate_g2() ) && curIt % s.getIterationSkip( DIR_TAU ) == 0 ) {
        curIt = 1;
        return true;
    }
    curIt++;
    return false;
}

int ODESolver::reset( const System &s ) {
    int dim = (int)( s.parameters.iterations_t_max / s.getIterationSkip( DIR_TAU ) ) + 10;
    out.clear();
    out.reserve( (int)( std::ceil( s.parameters.spectrum_frequency_iterations ) + 5 ) );
    akf_mat = MatrixXcd::Zero( dim, dim );
    for ( int w = 0; w < s.parameters.spectrum_frequency_iterations; w++ ) {
        out.push_back( 0 );
    }
    return dim;
}

// Numerical functions:
// Calculate T-Direction, output to files via system.expectationValues and save matrices for further calculations
bool ODESolver::calculate_t_direction( System &s ) {
    logs.level2( "Calculating t-direction from {} to {} at stepsize {}... ", s.getTimeborderStart(), s.getTimeborderEnd(), s.getTimeStep() );
    Timer &rkTimer = createTimer( "RungeKutta-Main-Loop" );
    ProgressBar progressbar = ProgressBar( s.parameters.iterations_total_max, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    rkTimer.start();

    int curIt = 1;
    MatrixXcd rho = s.getRho0();
    saveState( rho, s.getTimeborderStart() );
    s.expectationValues( rho, s.getTimeborderStart() ); //maybe: neue funktion für output, das dann in system, eig redundant, da nie t direction ausgerechnet wird ohne sie auch auszugeben.

    // Main Time Loop
    for ( double t_t = s.getTimeborderStart() + s.getTimeStep(); t_t < s.getTimeborderEnd(); t_t += s.getTimeStep() ) {
        // Runge-Kutta iteration
        rho = iterate( rho, s, t_t );
        // Save Rho for tau-direction
        if ( queueNow( s, curIt ) )
            saveState( rho, t_t );
        // Expectation Values
        s.expectationValues( rho, t_t );
        // Progress and time output
        rkTimer.iterate();
        outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_total_max, "T-Direction: " );
        // Divergent bruder?
        if ( !s.traceValid( rho, t_t ) ) {
            t_t = s.getTimeborderEnd() + s.getTimeStep();
            break;
        }
    }
    rkTimer.end();
    logs.level2( "Done! Saved {} states.\n", savedStates.size() );
    return true;
}
// Helperfunctions
int ODESolver::getIterationNumberTau( const System &s ) {
    int num = 0;
    // Tau Direction Iteration steps
    for ( int i = 0; i < (int)savedStates.size(); i++ ) {
        double t_t = getTimeAt( i );
        for ( double t_tau = t_t + s.getTimeStep(); t_tau < s.getTimeborderEnd(); t_tau += s.getTimeStep() ) { // t + +s.getTimeStep()
            num++;
        }
    }
    return num;
}
int ODESolver::getIterationNumberSpectrum( const System &s ) {
    int num = 0;
    // Spectrum steps
    for ( int spec_w = 0; spec_w < s.parameters.spectrum_frequency_iterations; spec_w++ ) {
        num++;
    }
    return num;
}
double ODESolver::getTimeAt( int i ) {
    return savedStates.at( i ).t;
}
MatrixXcd ODESolver::getRhoAt( int i ) {
    return savedStates.at( i ).rho;
}
// Calculate tau-direction for given operators b,b^+ and precalculated rho.
bool ODESolver::calculate_g1( const System &s, const MatrixXcd &op_creator, const MatrixXcd &op_annihilator ) {
    if ( !( (int)savedStates.size() > 0 ) ) {
        logs( "Need to calculate t-direction first!\n" );
        return false;
    }

    Timer &timer = createTimer( "RungeKutta-Tau-Loop" );
    int totalIterations = getIterationNumberTau( s );
    ProgressBar progressbar = ProgressBar( totalIterations, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    progressbar.strBarStart = "|";
    progressbar.strBarEnd = "|";
    timer.start();

#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < (int)savedStates.size(); i++ ) {
        double t_t = getTimeAt( i );
        MatrixXcd rho_tau = s.dgl_calc_rhotau( getRhoAt( i ), op_annihilator, t_t );
        akf_mat( i, 0 ) = s.dgl_expectationvalue( rho_tau, op_creator, t_t );
        int j = 1;
        int curIt_tau = 1;
        for ( double t_tau = t_t + s.getTimeStep(); t_tau < s.getTimeborderEnd(); t_tau += s.getTimeStep() ) { // t + +s.getTimeStep()
            rho_tau = iterate( rho_tau, s, t_tau, DIR_TAU );
            timer.iterate();
            if ( queueNow( s, curIt_tau ) ) {
                akf_mat( i, j ) = s.dgl_expectationvalue( rho_tau, op_creator, t_tau );
                j++; // equivalent to s.parameters.akf_vecIndex
            }
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "AKF-Tau: " );
        }
    }
    timer.end();
    return true;
}

bool ODESolver::calculate_g2_0( const System &s, const MatrixXcd &op_creator, const MatrixXcd &op_annihilator, std::string fileOutputName = "g2(0).txt" ) {
    if ( !( (int)savedStates.size() > 0 ) ) {
        logs( "Need to calculate t-direction first!\n" );
        return false;
    }

    Timer &timer = createTimer( "G2-0-Loop" );
    int totalIterations = (int)savedStates.size();
    ProgressBar progressbar = ProgressBar( totalIterations, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    progressbar.strBarStart = "|";
    progressbar.strBarEnd = "|";
    timer.start();

    std::vector<std::complex<double>> g2Values;
    g2Values.reserve( totalIterations );
    for ( int i = 0; i < totalIterations; i++ ) {
        double t_t = getTimeAt( i );
        MatrixXcd rho = getRhoAt( i );
        MatrixXcd M1 = op_creator * op_creator * op_annihilator * op_annihilator;
        MatrixXcd M2 = op_creator * op_annihilator;
        g2Values.emplace_back( s.dgl_expectationvalue( rho, M1, t_t ) / std::pow( s.dgl_expectationvalue( rho, M2, t_t ), 2 ) );
        timer.iterate();
        outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "G2: " );
    }
    timer.end();
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
    return true;
}

// Calculate Spectrum from previously calculated akf_mat
bool ODESolver::calculate_spectrum( const System &s, std::string fileOutputName = "spectrum.txt" ) {
    Timer &timer = createTimer( "Spectrum-Loop" );
    int totalIterations = getIterationNumberSpectrum( s );
    ProgressBar progressbar = ProgressBar( totalIterations, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    progressbar.strBarStart = "|";
    progressbar.strBarEnd = "|";
    timer.start();

    //Calculate frequencies:
    std::vector<double> spectrum_frequency_w;
    for ( int w = 0; w < s.parameters.spectrum_frequency_iterations; w++ ) {
        spectrum_frequency_w.push_back( s.parameters.spectrum_frequency_center - ( s.parameters.spectrum_frequency_range ) + w / ( (double)s.parameters.iterations_t_max ) * ( ( s.parameters.iterations_skips_tau - 1 ) * ( 1 - s.parameters.iterations_skips_w ) + 1 ) * ( 2. * ( s.parameters.spectrum_frequency_range ) ) );
    }
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads ) //reduction(+:spectrum.out)
    for ( int spec_w = 0; spec_w < s.parameters.spectrum_frequency_iterations; spec_w++ ) {
        std::vector<std::complex<double>> expfunc;
        expfunc.reserve( savedStates.size() );
        for ( int spec_tau = 0; spec_tau < ( s.getTimeborderEnd() - s.getTimeborderStart() - 0.0 ) / ( s.getTimeStep() * s.getIterationSkip( DIR_TAU ) ); spec_tau++ ) {
            expfunc.emplace_back( std::exp( -1i * spectrum_frequency_w.at( spec_w ) * (double)(spec_tau)*s.getTimeStep() * (double)( s.getIterationSkip( DIR_TAU ) ) ) );
        }
        for ( int i = 0; i < (int)savedStates.size(); i++ ) {
            double t_t = getTimeAt( i );
            for ( int spec_tau = 0; spec_tau < ( s.getTimeborderEnd() - s.getTimeborderStart() - t_t ) / ( s.getTimeStep() * s.getIterationSkip( DIR_TAU ) ); spec_tau++ ) {
                out.at( spec_w ) += expfunc.at( spec_tau ) * akf_mat( i, spec_tau );
            }
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "AKF-Spectrum: " );
        }
        timer.iterate();
    }
    outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "AKF-Spectrum", PROGRESS_FORCE_OUTPUT );
    timer.end();
    std::string filepath = s.parameters.subfolder + fileOutputName;
    FILE *spectrumfile = std::fopen( filepath.c_str(), "w" );
    if ( !spectrumfile ) {
        logs.level2( "Failed to open outputfile for spectrum!\n" );
        return false;
    }
    for ( int spec_w = 0; spec_w < s.parameters.spectrum_frequency_iterations; spec_w++ ) {
        fmt::print( spectrumfile, "{0:.15e}\t{1:.15e}\n", spectrum_frequency_w[spec_w], real( out[spec_w] * s.getTimeStep() * s.getTimeStep() * (double)( s.getIterationSkip( DIR_TAU ) * s.getIterationSkip( DIR_TAU ) ) ) );
    }
    std::fclose( spectrumfile );
    return true;
}