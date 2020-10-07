#include "solver/solver_ode.h"

ODESolver::ODESolver( System &s ) {
    logs.level2( "Creating ODESolver Class... " );
    reset( s );
    savedStates.clear();
    savedStates.reserve( dim );
    logs.level2( "Done!\n" );
}

MatType ODESolver::getHamilton( System &s, const double t, bool use_saved_hamiltons ) {
    if ( s.parameters.numerics_use_saved_hamiltons ) {////|| use_saved_hamiltons ) {
        if ( savedHamiltons.size() == 0 || t > savedHamiltons.back().t ) {
            track_gethamilton_calcattempt++;
            // Calculate new Hamilton for t, save it, return it
            saveHamilton( s.dgl_getHamilton( t ), t );
            track_gethamilton_write++;
            return savedHamiltons.back().mat;
        } else {
            long unsigned int approx = std::floor( t / s.parameters.t_step / 2.0 );
            if ( !(approx > savedHamiltons.size() || approx < 0) ) {
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
    }
    track_gethamilton_calc++;
    return s.dgl_getHamilton( t );
}



void ODESolver::saveState( const MatType &mat, const double t, std::vector<SaveState> &savedStates ) {
    savedStates.emplace_back( SaveState( mat, t ) );
}

void ODESolver::saveHamilton( const MatType &mat, const double t ) {
    savedHamiltons.emplace_back( SaveState( mat, t ) );
}

bool ODESolver::queueNow( System &s, int &curIt ) {
    bool queue = s.parameters.numerics_calculate_spectrum_H || s.parameters.numerics_calculate_spectrum_V || s.parameters.numerics_calculate_g2;
    if ( queue && curIt % s.parameters.iterations_t_skip == 0 ) {
        curIt = 1;
        return true;
    }
    curIt++;
    return false;
}

int ODESolver::reset( System &s ) {
    track_gethamilton_read = 0;
    track_gethamilton_write = 0;
    track_gethamilton_calc = 0;
    track_gethamilton_calcattempt = 0;
    dim = (int)std::ceil((s.parameters.t_end-s.parameters.t_start)/s.parameters.t_step)+10;//(int)( s.parameters.iterations_t_max / s.parameters.iterations_t_skip ) + 10;
    return dim;
}

bool ODESolver::calculate_t_direction( System &s ) {
    MatType rho = s.operatorMatrices.rho;
    
    Timer &rkTimer = createTimer( "RungeKutta-Main-Loop" );
    ProgressBar progressbar = ProgressBar( s.parameters.iterations_t_max, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    rkTimer.start();
    logs.level2( "Calculating t-direction from {} to {} at stepsize {}... ", s.parameters.t_start, s.parameters.t_end, s.parameters.t_step );
    
    // Calculate Time evolution on time vector timestamps.
    calculate_runge_kutta( rho, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );

    // Finalize
    outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, "T-Direction: ", PROGRESS_FORCE_OUTPUT );
    rkTimer.end();
    logs.level2( "Done! Saved {} states, cached {} hamilton matrices.\n", savedStates.size(), savedHamiltons.size() );
    
    
    // Interpolate Matrices?
    std::vector<SaveState> interpolate_savedstates;
    if ( s.parameters.numerics_interpolate_outputs ) {
        interpolate_savedstates = Solver::calculate_smooth_curve( savedStates, s.parameters.t_start, s.parameters.t_end, std::max(s.parameters.iterations_t_max*5,2500), s.parameters.output_handlerstrings );
    } else {
        interpolate_savedstates = savedStates;
    }

    // Calculate expectation values
    for ( long unsigned int i = 0; i < interpolate_savedstates.size(); i++) {
        s.traceValid( interpolate_savedstates.at(i).mat, interpolate_savedstates.at(i).t );
        s.expectationValues( interpolate_savedstates.at(i).mat, interpolate_savedstates.at(i).t, interpolate_savedstates );
    }
    return true;
}

int ODESolver::getIterationNumberTau( System &s ) {
    int num = 0;
    // Tau Direction Iteration steps
    for ( int i = 0; i < (int)savedStates.size(); i += s.parameters.iterations_t_skip ) {
        double t_t = getTimeAt( i );
        for ( double t_tau = t_t + s.parameters.t_step; t_tau < s.parameters.t_end; t_tau += s.parameters.t_step ) { // t + +s.parameters.t_step
            num++;
        }
    }
    return num;
}

// Description: Function to calculate the number of iterations used for spectru calculations
// Type: ODESolver private function
// @param s: [&System] Class providing set of system functions
// @return: [int] Number of spectrum iterations

int ODESolver::getIterationNumberSpectrum( System &s ) {
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
// @return: [MatType] (Density-) Matrix corresponding to iteration number i

MatType ODESolver::getRhoAt( int i ) {
    return savedStates.at( i ).mat;
}

// Description: Calculates the Eberly-Wódkiewicz spectrum. Uses precalculated values from akf_mat
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false
bool ODESolver::calculate_spectrum( System &s, const MatType &op_creator, const MatType &op_annihilator, std::string fileOutputName, const int cache_index ) {
    // Send system command to change to single core mainprogram, because this memberfunction is already using multithreading
    s.command( Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
    // Cache matrix. Saves complex double entries of G1(t,tau)
    Dense akf_mat = Dense::Zero( dim, dim );
    // Calculate G1(t,tau) with given operator matrices
    calculate_g1( s, op_creator, op_annihilator, akf_mat, "spectrum" );
    // Create Timer and Progressbar for the spectrum loop
    Timer &timer = createTimer( "Spectrum-Loop" );
    int totalIterations = getIterationNumberSpectrum( s );
    ProgressBar progressbar = ProgressBar( totalIterations );
    timer.start();
    logs.level2( "Calculating spectrum... Calculating frequencies... " );
    //Calculate frequencies:
    std::vector<double> spectrum_frequency_w;
    std::vector<Scalar> out;
    for ( int w = 0; w < s.parameters.iterations_w_resolution; w++ ) {
        spectrum_frequency_w.push_back( s.parameters.spectrum_frequency_center - ( s.parameters.spectrum_frequency_range ) + w / ( (double)s.parameters.iterations_w_resolution ) * ( 2. * ( s.parameters.spectrum_frequency_range ) ) );
        out.push_back( 0 );
    }
    logs.level2( "Done, calculating fourier transform via direct integral... " );
    // Calculate main fourier transform integral
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int spec_w = 0; spec_w < s.parameters.iterations_w_resolution; spec_w++ ) {
        std::vector<Scalar> expfunc;
        expfunc.reserve( dim / s.parameters.iterations_wtau_skip ); ////savedStates.size() / s.parameters.iterations_t_skip + 1 );
        for ( int spec_tau = 0; spec_tau < dim; spec_tau+=s.parameters.iterations_wtau_skip ) { ////( s.parameters.t_end - s.parameters.t_start - 0.0 ) / ( s.parameters.t_step * s.parameters.iterations_t_skip )
            expfunc.emplace_back( std::exp( -1i * spectrum_frequency_w.at( spec_w ) * (double)(spec_tau)*s.parameters.t_step ) );//// (double)( s.parameters.iterations_t_skip ) ) );
        }
        for ( long unsigned int i = 0; i < dim; i+=s.parameters.iterations_t_skip) {////= s.parameters.iterations_t_skip ) { ////savedStates.size()
            ////int k = i / s.parameters.iterations_t_skip;
            double t_t = ((double)i)*s.parameters.t_step; ////getTimeAt( i );
            int dim2 = (int)std::floor(( s.parameters.t_end - s.parameters.t_start - t_t ) / ( s.parameters.t_step ));
            for ( int j = 0; j < dim2; j+=s.parameters.iterations_wtau_skip ) { ////( s.parameters.t_end - s.parameters.t_start - t_t ) / ( s.parameters.t_step * s.parameters.iterations_t_skip )
                out.at( spec_w ) += expfunc.at( j/s.parameters.iterations_wtau_skip ) * akf_mat( i, j ) / s.parameters.t_step/s.parameters.t_step / ( (double)s.parameters.iterations_t_skip*s.parameters.iterations_wtau_skip ); ////akf_mat( k, spec_tau );  
            }
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "Spectrum: " );
        }
        timer.iterate();
    }
    // Final output and timer end
    outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "Spectrum", PROGRESS_FORCE_OUTPUT );
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
        fmt::print( spectrumfile, "{0:.15e}\t{1:.15e}\n", spectrum_frequency_w[spec_w], real( out[spec_w] * s.parameters.t_step * s.parameters.t_step * (double)( s.parameters.iterations_t_skip * s.parameters.iterations_t_skip ) ) );
    }
    std::fclose( spectrumfile );
    logs.level2( "Done!\n" );
    // Check if we need to save calculated g-functions
    if ( s.parameters.numerics_calculate_g2 ) {
        logs.level2( "Saving g1 to temporary matrix {}", cache_index );
        if ( cache_index == 1 )
            cache1 = akf_mat;
        else if ( cache_index == 2 )
            cache2 = akf_mat;
        logs.level2( " succesfull\n" );
    }
    // Sucessfull
    return true;
}


// Description: Calculates G2 photon statistics
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator(annihilator)_X: set of operator matrices and their conjugated forms
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false

bool ODESolver::calculate_advanced_photon_statistics( System &s, const MatType &op_creator_1, const MatType &op_annihilator_1, const MatType &op_creator_2, const MatType &op_annihilator_2, std::string fileOutputName ) {
    //bool calculate_full_g2_of_tau = !s.parameters.numerics_use_simplified_g2;
    bool output_full_g2 = false;
    bool output_full_ind = s.parameters.numerics_calculate_timeresolution_indistinguishability;
    bool calculate_concurrence_with_g2_of_zero = s.parameters.numerics_use_simplified_g2;
    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
    s.command( Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
    int pbsize = (int)savedStates.size() / s.parameters.iterations_t_skip;
    // Progress
    ProgressBar progressbar = ProgressBar( pbsize );
    // Cache matrix. Saves complex double entries of G1(t,tau). Will be created even if they are not needed
    Dense akf_mat_g2_11 = Dense::Zero( dim, dim );
    Dense akf_mat_g2_22 = Dense::Zero( dim, dim );
    Dense akf_mat_g2_12 = Dense::Zero( dim, dim );
    // Calculate G2(t,tau) with given operator matrices
    if ( !calculate_concurrence_with_g2_of_zero ) {
        logs.level2( "Calculating 3 seperate G2(t,tau) matrices...\n" );
        if ( s.parameters.numerics_calculate_g2_H || s.parameters.numerics_calculate_g2_C )
            calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_1, op_annihilator_1, akf_mat_g2_11, "Mode 11" );
        if ( s.parameters.numerics_calculate_g2_V || s.parameters.numerics_calculate_g2_C )
            calculate_g2( s, op_creator_2, op_annihilator_2, op_creator_2, op_annihilator_2, akf_mat_g2_22, "Mode 22" );
        if ( s.parameters.numerics_calculate_g2_C )
            calculate_g2( s, op_creator_1, op_annihilator_2, op_creator_1, op_annihilator_2, akf_mat_g2_12, "Mode 12" );
    }
    // Modified Dimension Step. If Upscaling of G1/G2 is active, step for matrices are i instead of t_skip
    int mat_step = (s.parameters.numerics_stretch_correlation_grid ? 1 : s.parameters.iterations_t_skip);
    logs.level2("Using mat_step = {}\n",mat_step);
    // Modified T-Step
    double deltaT = s.parameters.t_step * ( 1.0 * mat_step );
    // Calculating G2(tau)
    Timer &timer_g2 = createTimer( "G2(tau integral" );
    timer_g2.start();
    logs.level2( "Calculating G2(tau)... " );
    // Prepare output vector
    std::vector<Scalar> g2_11, g2_22, g2_12;
    std::vector<Scalar> g2_11_zero, g2_22_zero, g2_12_zero;
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        g2_11.push_back( 0 );
        g2_22.push_back( 0 );
        g2_12.push_back( 0 );
        g2_11_zero.push_back( 0 );
        g2_22_zero.push_back( 0 );
        g2_12_zero.push_back( 0 );
    }
    MatType M1, M2, M3, rho;
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        int k = i / mat_step;
        if ( !calculate_concurrence_with_g2_of_zero ) {
            for ( int t = 0; t < dim; t++ ) {
                g2_11.at( k ) += akf_mat_g2_11( k, t ) * s.parameters.t_step;
                g2_22.at( k ) += akf_mat_g2_22( k, t ) * s.parameters.t_step;
                g2_12.at( k ) += akf_mat_g2_12( k, t ) * s.parameters.t_step;
            }
        }
        double t_t = getTimeAt( i );
        rho = getRhoAt( i );
        if ( s.parameters.numerics_calculate_g2_H ) {
            M1 = op_creator_1 * op_creator_1 * op_annihilator_1 * op_annihilator_1;
            M2 = op_creator_1 * op_annihilator_1;
            g2_11_zero.at( k ) = s.dgl_expectationvalue<MatType, Scalar>( rho, M1, t_t ); // / std::pow( s.dgl_expectationvalue<MatType, Scalar>( rho, M2, t_t ), 2 );
        }
        if ( s.parameters.numerics_calculate_g2_V ) {
            M1 = op_creator_2 * op_creator_2 * op_annihilator_2 * op_annihilator_2;
            M2 = op_creator_2 * op_annihilator_2;
            g2_22_zero.at( k ) = s.dgl_expectationvalue<MatType, Scalar>( rho, M1, t_t ); // / std::pow( s.dgl_expectationvalue<MatType, Scalar>( rho, M2, t_t ), 2 );
        }
        if ( s.parameters.numerics_calculate_g2_C || s.parameters.numerics_calculate_g2_V || s.parameters.numerics_calculate_g2_H ) {
            M1 = op_creator_1 * op_creator_1 * op_annihilator_2 * op_annihilator_2;
            M2 = op_creator_1 * op_annihilator_1;
            M3 = op_creator_2 * op_creator_2;
            g2_12_zero.at( k ) = s.dgl_expectationvalue<MatType, Scalar>( rho, M1, t_t ); // / ( s.dgl_expectationvalue<MatType, Scalar>( rho, M2, t_t ) * s.dgl_expectationvalue<MatType, Scalar>( rho, M3, t_t ) );
        }
        outputProgress( s.parameters.output_handlerstrings, timer_g2, progressbar, pbsize, "G2(tau): " );
        timer_g2.iterate();
    }
    timer_g2.end();
    // Calculating Concurrence
    // Prepare Concurrence output vector
    std::vector<Scalar> rho_11, rho_22, rho_12;
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        rho_11.push_back( 0 );
        rho_22.push_back( 0 );
        rho_12.push_back( 0 );
    }
    if ( s.parameters.numerics_calculate_g2_C ) {
        Timer &timer_c = createTimer( "Concurrence" );
        timer_c.start();
        logs.level2( "Calculating Concurrence... " );
        // Main integral
        logs.level2( "Concurrence integral timestep: {}\n", deltaT );
        //#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
        for ( long unsigned int T = 0; T < savedStates.size(); T += mat_step ) {
            int t = T / mat_step;
            for ( int i = 0; i < T; i+= mat_step ) { ////+= s.parameters.iterations_t_skip
                ////int k = i / s.parameters.iterations_t_skip;
                // Calculate for t
                if ( !calculate_concurrence_with_g2_of_zero ) {
                    for ( int tau = 0; tau < t - i; tau++ ) {
                        rho_11.at( t ) += akf_mat_g2_11( i, tau ) * s.parameters.t_step * deltaT;
                        rho_22.at( t ) += akf_mat_g2_22( i, tau ) * s.parameters.t_step * deltaT;
                        rho_12.at( t ) += akf_mat_g2_12( i, tau ) * s.parameters.t_step * deltaT;
                    }
                } else {
                    int k = i / mat_step;
                    rho_11.at( t ) += g2_11_zero.at( k ) * deltaT;
                    rho_22.at( t ) += g2_22_zero.at( k ) * deltaT;
                    rho_12.at( t ) += g2_12_zero.at( k ) * deltaT;
                }
                outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, pbsize, "Concurrence: " );
            }
            timer_c.iterate();
        }
        // Final output and timer end
        outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, pbsize, "Concurrence", PROGRESS_FORCE_OUTPUT );
        timer_c.end();
    }
    // Calculating indistinguishability
    Dense akf_mat_g1_1 = Dense::Zero( dim, dim );
    Dense akf_mat_g1_2 = Dense::Zero( dim, dim );
    std::vector<double> p1, p2;
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        p1.push_back( 0 );
        p2.push_back( 0 );
    }
    if ( !calculate_concurrence_with_g2_of_zero ) {
        M1 = op_creator_1 * op_annihilator_1;
        M2 = op_creator_2 * op_annihilator_2;
        if ( !cache1.isZero() )
            akf_mat_g1_1 = cache1;
        else if ( s.parameters.numerics_calculate_g2_H ) //TODO: hier vieleicht weg, nur g1 benutzen wenn auch das spektrum ausgerechnet wurde, sonst g1=0
            calculate_g1( s, op_creator_1, op_annihilator_1, akf_mat_g1_1, "g1_1" );
        if ( !cache2.isZero() )
            akf_mat_g1_2 = cache2;
        else if ( s.parameters.numerics_calculate_g2_V )
            calculate_g1( s, op_creator_2, op_annihilator_2, akf_mat_g1_2, "g1_2" );
        Timer &timer = createTimer( "Indistinguishability" );
        timer.start();
        long unsigned int tstart = output_full_ind ? 0 : savedStates.size() - mat_step;
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
        for ( long unsigned int T = tstart; T < savedStates.size(); T += mat_step ) {
            Scalar top_1 = 0;
            Scalar top_2 = 0;
            Scalar bottom_1 = 0;
            Scalar bottom_2 = 0;
            for ( int i = 0; i < T; i += mat_step ) {
                int k = i / mat_step;
                auto rho = getRhoAt( i );
                double t = getTimeAt( i );
                for ( int j = 0; j < T - i; j += 1 ) { //&& j + i < (int)savedStates.size() ////mat_step
                    int l = j / mat_step;
                    auto rho_tau = getRhoAt( i + j );
                    double t_tau = getTimeAt( i + j );
                    Scalar gpop_1 = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M1, t_tau );
                    Scalar gpop_2 = s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M2, t_tau );
                    Scalar gbot_1 = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator_1, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator_1, t );
                    Scalar gbot_2 = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator_2, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator_2, t );
                    top_1 += ( gpop_1 + akf_mat_g2_11( i, j ) - akf_mat_g1_1( i, j ) * std::conj( akf_mat_g1_1( i, j ) ) ); //k,j nicht k, l
                    bottom_1 += 2.0 * gpop_1 - gbot_1 * std::conj( gbot_1 );
                    //fmt::print( "what contributes: gpop {}, g2_11 {}, g1_1 {}, gbot {}, TOP {}, BOT {}\n", gpop_1, akf_mat_g2_11( k, l ), std::pow( std::abs( akf_mat_g1_1( k, l ) ), 2.0 ), std::pow( std::abs( gbot_1 ), 2.0 ), top_1, bottom_1 );
                    top_2 += ( gpop_2 + akf_mat_g2_22( i, j ) - akf_mat_g1_2( i, j ) * std::conj( akf_mat_g1_2( i, j ) ) ); //T,j nicht k,l
                    bottom_2 += 2.0 * gpop_2 - gbot_2 * std::conj( gbot_2 );
                }
                // If we dont output the full time resoluted indistinguishability, we instead output the time resoluted integral for T = T_max
                if ( !output_full_ind ) {
                    p1.at( i / mat_step ) = ( 1.0 - std::abs( top_1 / bottom_1 ) );
                    p2.at( i / mat_step ) = ( 1.0 - std::abs( top_2 / bottom_2 ) );
                    outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability (Simplified): " );
                    timer.iterate();
                }
            }
            p1.at( T / mat_step ) = ( 1.0 - std::real( top_1 / bottom_1 ) );
            p2.at( T / mat_step ) = ( 1.0 - std::real( top_2 / bottom_2 ) );
            outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability (Full): " );
            timer.iterate();
        }
        // Final output and timer end
        outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability: ", PROGRESS_FORCE_OUTPUT );
        timer.end();
    }
    // Save output
    logs.level2( "Done, saving to {}... ", fileOutputName );
    std::string filepath = s.parameters.subfolder + fileOutputName;
    FILE *concurrencefile = std::fopen( filepath.c_str(), "w" );
    if ( !concurrencefile ) {
        logs.level2( "Failed to open outputfile for advanced photon statistics!\n" );
        return false;
    }
    fmt::print( concurrencefile, "Time\tMatrixelement_11\tMatrixelement_22\tConcurrence\tG2_11(tau)\tG2_22(tau)\tG2_12(tau)\tG2(tau=0)_11\tG2(tau=0)_22\tG2(tau=0)_12\tI_1\tI_2\n" );
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        int k = i / mat_step;
        fmt::print( concurrencefile, "{:.5e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", getTimeAt( i ), std::abs( rho_11.at( k ) / ( rho_11.at( k ) + rho_22.at( k ) ) ), std::abs( rho_22.at( k ) / ( rho_11.at( k ) + rho_22.at( k ) ) ), std::abs( 2.0 * std::abs( rho_12.at( k ) ) / ( rho_11.at( k ) + rho_22.at( k ) ) ), std::real( g2_11.at( k ) ), std::real( g2_22.at( k ) ), std::real( g2_12.at( k ) ), std::real( g2_11_zero.at( k ) ), std::real( g2_22_zero.at( k ) ), std::real( g2_12_zero.at( k ) ), p1.at( k ), p2.at( k ) );
    }
    std::fclose( concurrencefile );
    logs.level2( "Done!\n" );
    logs( "Final Concurrence: {}\n", 2.0 * std::abs( rho_12.back() / ( rho_11.back() + rho_22.back() ) ) );
    if ( !calculate_concurrence_with_g2_of_zero )
        logs( "Final Indistinguishability: {} H, {} V\n", p1.back(), p2.back() );
    //if ()
    //logs( "Final G2(0): {}\n", std::real( 2.0 * std::abs( rho_12.back() ) / ( rho_11.back() + rho_22.back() ) ) );

    if ( output_full_g2 ) {
        logs.level2( "Saving full, non integrated G2(t,tau) matrices... " );
        FILE *ganz = std::fopen( ( s.parameters.subfolder + "g2(ttau)_matrix.txt" ).c_str(), "w" );
        fmt::print( ganz, "t\ttau\tG2_11(t,tau)\tG2_22(t,tau)\tG2_12(t,tau)\n" );
        for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
            int k = i / mat_step;
            for ( long unsigned int j = 0; j < savedStates.size(); j += mat_step ) {
                int l = j / mat_step;
                fmt::print( ganz, "{0:.5e}\t{1:.5e}\t{2:.8e}\t{3:.8e}\t{4:.8e}\n", getTimeAt( i ), getTimeAt( j ), std::abs( akf_mat_g2_11( k, l ) ), std::abs( akf_mat_g2_22( k, l ) ), std::abs( akf_mat_g2_12( k, l ) ) );
            }
            fmt::print( ganz, "\n" );
        }
        std::fclose( ganz );
        logs.level2( "Done!\n" );
    }
    // Sucessfull
    return true;
}