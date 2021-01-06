#include "system/system.h"

System::System( const std::vector<std::string> &input ) : Operatable() {
    // Set Name of this system.
    name = "Biexciton (4NS)";
    Log::L2( "Creating System Class for '{}'\n", name );
    // Initialize all subclasses with the input vector
    parameters = Parameters( input );
    operatorMatrices = OperatorMatrices( parameters );
    operatorMatrices.outputOperators( parameters );
    // Log the operator matrix base
    parameters.log( operatorMatrices.base );
    // Create all possible file outputs
    fileoutput = FileOutput( {"densitymatrix.txt", "atomicinversion.txt", "photonpopulation.txt"}, parameters, operatorMatrices );
    // Initialize / Adjust the remaining system class
    terminate_message = global_message_normaltermination;
    Timer &timer_systeminit = Timers::create( "System Initialization", true, false );
    Log::L2( "System initialization...\n" );
    timer_systeminit.start();
    if ( !init_system() ) {
        Log::L2( "System initialization failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    Log::L2( "Successful! Elapsed time is {}ms\n", timer_systeminit.getWallTime( TIMER_MILLISECONDS ) );
    timer_systeminit.end();
}

bool System::init_system() {
    // Single chirp for single atomic level
    Chirp::Inputs chirpinputs( parameters.t_start, parameters.t_end, parameters.t_step, parameters.chirp_type, parameters.numerics_order_highest );
    chirpinputs.add( parameters.chirp_t, parameters.chirp_y, parameters.chirp_ddt );
    chirp = Chirp( chirpinputs );
    if ( parameters.chirp_total != 0 )
        chirp.fileOutput( parameters.subfolder + "chirp.txt" );

    // Arbitrary number of pulses onto single atomic level.
    Pulse::Inputs pulseinputs_H( parameters.t_start, parameters.t_end, parameters.t_step, parameters.numerics_order_highest );
    pulseinputs_H.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "H" );
    pulseinputs_H.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "+", 1.0 / std::sqrt( 2.0 ) );
    pulseinputs_H.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "-", 1.0 / std::sqrt( 2.0 ) );
    pulse_H = Pulse( pulseinputs_H );
    Pulse::Inputs pulseinputs_V( parameters.t_start, parameters.t_end, parameters.t_step, parameters.numerics_order_highest );
    pulseinputs_V.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "V" );
    pulseinputs_V.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "+", 1i / std::sqrt( 2.0 ) );
    pulseinputs_V.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "-", -1i / std::sqrt( 2.0 ) );
    pulse_V = Pulse( pulseinputs_V );
    if ( parameters.pulse_amp.at( 0 ) != 0 ) {
        Pulse::fileOutput( parameters.subfolder + "pulse.txt", {pulse_H, pulse_V} );
    }

    if ( parameters.numerics_use_saved_coefficients )
        savedCoefficients.reserve( parameters.numerics_saved_coefficients_max_size );

    // Output Phonon functions if phonons are active //TODO: initialize_polaron_frame_functions
    if ( parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
        initialize_path_integral_functions();
    } else {
        initialize_polaron_frame_functions();
    }
    // Time Transformation
    timeTrafoMatrix = ( Dense( 1i * operatorMatrices.H_0 ).exp() ).sparseView(); //.pruned();
    // Check time trafo
    Sparse ttrafo = ( Dense( 1i * operatorMatrices.H_0 * 125E-12 ).exp() * operatorMatrices.H_used * Dense( -1i * operatorMatrices.H_0 * 125E-12 ).exp() ).sparseView();
    Sparse temp = dgl_timetrafo( operatorMatrices.H_used, 125E-12 ) - ttrafo;
    if ( std::abs( Dense( temp ).sum() / parameters.p_omega_atomic_G_H ) >= 0.0001 ) {
        Log::L1( "Unitary timetransformation error is bigger than 0.01% of atomic transition energy!\n\n" );
    }

    return true;
}

Sparse System::dgl_rungeFunction( const Sparse &rho, const Sparse &H, const double t, std::vector<SaveState> &past_rhos ) {
    Sparse ret = -1i * dgl_kommutator( H, rho );
    // Photon losses
    if ( parameters.p_omega_cavity_loss != 0.0 ) {
        ret += parameters.p_omega_cavity_loss / 2.0 * dgl_lindblad( rho, operatorMatrices.photon_annihilate_H, operatorMatrices.photon_create_H ); //BUGFIX: kappa/2.0
        ret += parameters.p_omega_cavity_loss / 2.0 * dgl_lindblad( rho, operatorMatrices.photon_annihilate_V, operatorMatrices.photon_create_V ); //BUGFIX: kappa/2.0
    }
    // Pure Dephasing
    if ( parameters.p_omega_pure_dephasing != 0.0 ) {
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_ground * rho * operatorMatrices.atom_state_H + operatorMatrices.atom_state_H * rho * operatorMatrices.atom_state_ground );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_ground * rho * operatorMatrices.atom_state_V + operatorMatrices.atom_state_V * rho * operatorMatrices.atom_state_ground );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_H * rho * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_biexciton * rho * operatorMatrices.atom_state_H );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_V * rho * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_biexciton * rho * operatorMatrices.atom_state_V );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_H * rho * operatorMatrices.atom_state_V + operatorMatrices.atom_state_V * rho * operatorMatrices.atom_state_H );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_ground * rho * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_biexciton * rho * operatorMatrices.atom_state_ground );
    }
    // Radiative decay
    if ( parameters.p_omega_decay != 0.0 ) {
        ret += parameters.p_omega_decay / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H, operatorMatrices.atom_sigmaplus_G_H );
        ret += parameters.p_omega_decay / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V, operatorMatrices.atom_sigmaplus_G_V );
        ret += parameters.p_omega_decay / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B, operatorMatrices.atom_sigmaplus_H_B );
        ret += parameters.p_omega_decay / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B, operatorMatrices.atom_sigmaplus_V_B );
    }
    if ( parameters.numerics_phonon_approximation_order != PHONON_PATH_INTEGRAL && parameters.p_phonon_T >= 0.0 ) {
        ret += dgl_phonons_pmeq( rho, t, past_rhos );
    }
    return ret; //.pruned();
}

Sparse System::dgl_timetrafo( const Sparse &A, const double t ) {
    Sparse ret = A;
    if ( parameters.numerics_use_interactionpicture == 1 ) {
        // TIMETRANSFORMATION_ANALYTICAL
        if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) {
            std::vector<Eigen::Triplet<Scalar>> ret_v;
            for ( int k = 0; k < A.outerSize(); ++k ) {
                for ( Sparse::InnerIterator it( A, k ); it; ++it ) {
                    // Convert row/col into respective photon numbers / atomic state
                    int i = it.row() % 4;
                    int photon_H_i = std::floor( it.row() / ( 4.0 * ( parameters.p_max_photon_number + 1 ) ) );
                    int photon_V_i = ( (int)std::floor( it.row() / 4.0 ) ) % ( parameters.p_max_photon_number + 1 );
                    int j = it.col() % 4;
                    int photon_H_j = std::floor( it.col() / ( 4.0 * ( parameters.p_max_photon_number + 1 ) ) );
                    int photon_V_j = ( (int)std::floor( it.col() / 4.0 ) ) % ( parameters.p_max_photon_number + 1 );
                    // Add new Value to triplet list
                    Scalar factor = std::exp( 1i * t * ( parameters.p_omega_atomic_G_H * ( delta( i, 1 ) - delta( j, 1 ) ) + parameters.p_omega_atomic_G_V * ( delta( i, 2 ) - delta( j, 2 ) ) + parameters.p_omega_atomic_B * ( delta( i, 3 ) - delta( j, 3 ) ) + parameters.p_omega_cavity_H * ( photon_H_i - photon_H_j ) + parameters.p_omega_cavity_V * ( photon_V_i - photon_V_j ) ) );
                    //Scalar factor = std::exp( 1i * t * ( parameters.p_omega_atomic_G_H / 2.0 * ( delta( i, 1 ) - delta( i, 0 ) - delta( j, 1 ) + delta( j, 0 ) ) + parameters.p_omega_atomic_G_V / 2.0 * ( delta( i, 2 ) - delta( i, 0 ) - delta( j, 2 ) + delta( j, 0 ) ) + parameters.p_omega_atomic_H_B / 2.0 * ( delta( i, 3 ) - delta( i, 1 ) - delta( j, 3 ) + delta( j, 1 ) ) + parameters.p_omega_atomic_V_B / 2.0 * ( delta( i, 3 ) - delta( i, 2 ) - delta( j, 3 ) + delta( j, 2 ) ) + parameters.p_omega_cavity_H * ( photon_H_i - photon_H_j ) + parameters.p_omega_cavity_V * ( photon_V_i - photon_V_j ) ) );
                    ret_v.emplace_back( it.row(), it.col(), it.value() * factor );
                }
            }
            // Generate new Matrix from triplet list
            ret.setFromTriplets( ret_v.begin(), ret_v.end() );
        }
        // TIMETRANSFORMATION_MATRIXEXPONENTIAL
        else if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_MATRIXEXPONENTIAL ) {
            Sparse U = ( Dense( 1i * operatorMatrices.H_0 * t ).exp() ).sparseView();
            ret = U * A * U.adjoint();
        }
    }
    return ret;
}

Sparse System::dgl_chirp( const double t ) {
    if ( parameters.chirp_total == 0 )
        return Sparse( parameters.maxStates, parameters.maxStates );
    return ( 2.0 * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_H + operatorMatrices.atom_state_V ) * chirp.get( t ); //Experimental; corrected chirp
}

Sparse System::dgl_pulse( const double t ) {
    if ( parameters.pulse_amp.at( 0 ) == 0 ) {
        return Sparse( parameters.maxStates, parameters.maxStates );
    }
    return ( ( operatorMatrices.atom_sigmaminus_G_H + operatorMatrices.atom_sigmaminus_H_B ) * std::conj( pulse_H.get( t ) ) + ( operatorMatrices.atom_sigmaminus_G_V + operatorMatrices.atom_sigmaminus_V_B ) * std::conj( pulse_V.get( t ) ) + ( operatorMatrices.atom_sigmaplus_G_H + operatorMatrices.atom_sigmaplus_H_B ) * pulse_H.get( t ) + ( operatorMatrices.atom_sigmaplus_G_V + operatorMatrices.atom_sigmaplus_V_B ) * pulse_V.get( t ) );
    //return 0.5 * ( operatorMatrices.atom_sigmaplus_G_H * pulse_H.get( t ) + operatorMatrices.atom_sigmaminus_G_H * std::conj( pulse_H.get( t ) ) + operatorMatrices.atom_sigmaplus_G_V * pulse_V.get( t ) + operatorMatrices.atom_sigmaminus_G_V * std::conj( pulse_V.get( t ) ) );
}

Scalar System::dgl_raman_population_increment( const std::vector<SaveState> &past_rhos, const char mode, const Scalar before, const double t ) {
    Scalar ret = 0;
    double chirpcorrection = chirp.get( t );
    double w1 = ( mode == 'h' ? parameters.p_omega_atomic_G_H + chirpcorrection : parameters.p_omega_atomic_G_V + chirpcorrection );
    double w2 = ( mode == 'h' ? parameters.p_omega_atomic_H_B + chirpcorrection : parameters.p_omega_atomic_V_B + chirpcorrection );
    double wc = ( mode == 'h' ? parameters.p_omega_cavity_H : parameters.p_omega_cavity_V );
    double sigma1 = parameters.p_omega_pure_dephasing + parameters.p_omega_decay;
    double sigma2 = parameters.p_omega_pure_dephasing + 3. * parameters.p_omega_decay;
    Sparse op = operatorMatrices.atom_sigmaminus_G_B * ( mode == 'h' ? operatorMatrices.photon_create_H : operatorMatrices.photon_create_V );
    double A = std::exp( -parameters.p_omega_cavity_loss * parameters.t_step );
    Scalar B, R;
#pragma omp parallel for ordered schedule( dynamic ) shared( past_rhos ) num_threads( parameters.numerics_phonons_maximum_threads )
    for ( long unsigned int i = 0; i < past_rhos.size(); i++ ) {
        //for ( SaveState savestate : past_rhos ) {
        double tdd = past_rhos.at( i ).t;
        B = std::exp( -1i * ( w2 - wc - 0.5i * ( parameters.p_omega_cavity_loss + sigma2 ) ) * ( t - tdd ) ) - std::exp( -1i * ( w1 - wc - 0.5i * ( parameters.p_omega_cavity_loss + sigma1 ) ) * ( t - tdd ) );
        R = dgl_expectationvalue<Sparse, Scalar>( past_rhos.at( i ).mat, op, tdd ) * ( mode == 'h' ? std::conj( pulse_H.get( tdd ) ) : std::conj( pulse_V.get( tdd ) ) );
        //fmt::print("t = {}, tau = {}, A = {}, B = {}, R = {}\n",t,tdd,A,B,R);
#pragma omp critical
        ret += B * R;
    }
    return A * before + parameters.t_step * parameters.t_step * ret;
}

void System::expectationValues( const Sparse &rho, const double t, const std::vector<SaveState> &past_rhos ) {
    // Calculate reused expectation values
    double exp_photon_H = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.photon_n_H, t ) );
    double exp_photon_V = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.photon_n_V, t ) );
    // Calculating photon emission probability through photon expectation values
    photonemissionprob_integral_H += exp_photon_H * parameters.t_step * parameters.p_omega_cavity_loss;
    photonemissionprob_integral_V += exp_photon_V * parameters.t_step * parameters.p_omega_cavity_loss;
    // Calculating raman photon expectation values and emission probability
    if ( parameters.numerics_output_raman_population ) {
        ramanphotonpopulation_integral_H = dgl_raman_population_increment( past_rhos, 'h', ramanphotonpopulation_integral_H, t );
        ramanphotonpopulation_integral_V = dgl_raman_population_increment( past_rhos, 'v', ramanphotonpopulation_integral_V, t );
        ramanphotonemissionprob_integral_H += std::real( 2.0 * parameters.p_omega_coupling * ramanphotonpopulation_integral_H ) * parameters.t_step * parameters.p_omega_cavity_loss;
        ramanphotonemissionprob_integral_V += std::real( 2.0 * parameters.p_omega_coupling * ramanphotonpopulation_integral_V ) * parameters.t_step * parameters.p_omega_cavity_loss;
    }
    // Output expectation values
    double ev_operators_electronic_H = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.atom_state_H, t ) );
    double ev_operators_electronic_V = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.atom_state_V, t ) );
    fmt::print( fileoutput.fp_atomicinversion, "{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\n", t, std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.atom_state_ground, t ) ), ev_operators_electronic_H, ev_operators_electronic_V, std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.atom_state_biexciton, t ) ) );
    fmt::print( fileoutput.fp_photonpopulation, "{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}", t, exp_photon_H, exp_photon_V, photonemissionprob_integral_H, photonemissionprob_integral_V );
    if ( parameters.numerics_output_raman_population ) {
        fmt::print( fileoutput.fp_photonpopulation, "\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}", std::real( ramanphotonpopulation_integral_H ), std::real( ramanphotonpopulation_integral_V ), ramanphotonemissionprob_integral_H, ramanphotonemissionprob_integral_V );
    }

    if ( parameters.numerics_output_electronic_emission ) {
        electronic_emissionprob_integral_H += ev_operators_electronic_H * parameters.t_step * parameters.p_omega_decay;
        electronic_emissionprob_integral_V += ev_operators_electronic_V * parameters.t_step * parameters.p_omega_decay;
        fmt::print( fileoutput.fp_photonpopulation, "\t{:.5e}\t{:.5e}", electronic_emissionprob_integral_H, electronic_emissionprob_integral_V );
    }
    fmt::print( fileoutput.fp_photonpopulation, "\n" );

    if ( !parameters.output_no_dm ) {
        fmt::print( fileoutput.fp_densitymatrix, "{:.5e}\t", t ); //, rho.nonZeros(), rho.rows() * rho.cols() - rho.nonZeros() );
        if ( parameters.output_full_dm ) {
            for ( int i = 0; i < parameters.maxStates; i++ )
                for ( int j = 0; j < parameters.maxStates; j++ ) {
                    fmt::print( fileoutput.fp_densitymatrix, "{:.10e}\t", std::real( rho.coeff( i, j ) ) );
                }
        } else
            for ( int j = 0; j < parameters.maxStates; j++ )
                fmt::print( fileoutput.fp_densitymatrix, "{:.5e}\t", std::real( rho.coeff( j, j ) ) );
        fmt::print( fileoutput.fp_densitymatrix, "\n" );
    }
}

Sparse System::dgl_getHamilton( const double t ) {
    double local_b = ( parameters.numerics_phonon_approximation_order != PHONON_PATH_INTEGRAL ? 1.0 * parameters.p_phonon_b : 1.0 );
    return dgl_timetrafo( local_b * ( operatorMatrices.H_used + dgl_pulse( t ) ) + dgl_chirp( t ), t );
}

bool System::command( unsigned int index ) {
    // The only reason for classes using this system class to set the main programs maximum threads to 1 is, if the usual T-direction is already done.
    if ( index == Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM ) {
        parameters.numerics_phonons_maximum_threads = 1;
        Log::L2( "Set maximum number of Threads for primary calculations to {}\n", parameters.numerics_phonons_maximum_threads );
        if ( parameters.numerics_use_saved_coefficients ) {
            // Sort after t
            Log::L2( "Sorting saved coefficients by t... " );
            std::sort( savedCoefficients.begin(), savedCoefficients.end(), Save_State_sort_t );
            // Sort after tau
            Log::L2( "Done, sorting saved coefficients chunkwise by tau... " );
            int current_start = 0;
            int i;
            for ( i = 0; i < (int)savedCoefficients.size(); i++ ) {
                if ( savedCoefficients.at( i ).t > savedCoefficients.at( current_start ).t || i >= (int)savedCoefficients.size() ) {
                    std::sort( savedCoefficients.begin() + current_start, savedCoefficients.begin() + i, Save_State_sort_tau );
                    current_start = i;
                }
            }
            // Sort last entries
            std::sort( savedCoefficients.begin() + current_start, savedCoefficients.begin() + i, Save_State_sort_tau );
            // Just in case, test coefficients:
            //for ( i = 0; i < (int)savedCoefficients.size() - 1; i++ ) {
            //    if ( ( savedCoefficients.at( i ).t == savedCoefficients.at( i + 1 ).t && ( savedCoefficients.at( i ).tau >= savedCoefficients.at( i + 1 ).tau ) ) || savedCoefficients.at( i ).t > savedCoefficients.at( i + 1 ).t ) {
            //        Log::L2( "Coefficient mismatch after sorting! i = {}, t = {}, tau = {}, i+1 = {}, t = {}, tau = {}\n", i, savedCoefficients.at( i ).t, savedCoefficients.at( i ).tau, i + 1, savedCoefficients.at( i + 1 ).t, savedCoefficients.at( i + 1 ).tau );
            //    }
            //Log::L2( "T = {}, tau 0 {}\n", savedCoefficients.at( i ).t, savedCoefficients.at( i ).tau );
            //}
            Log::L2( "Done!\n" );
        }
    }
    return true;
}

bool System::exit_system( const int failure ) {
    pulse_H.log();
    pulse_V.log();
    chirp.log();
    Log::L2( "Coefficients: Attempts w/r: {}, Write: {}, Calc: {}, Read: {}, Read-But-Not-Equal: {}. Done!\n", track_getcoefficient_calcattempt, track_getcoefficient_write, track_getcoefficient_calculate, track_getcoefficient_read, track_getcoefficient_read_but_unequal );
    Log::L2( "Number of approx+/- adjustments: {}\n", globaltries );
    Log::L1( "Maximum RAM used: {} MB\n", getPeakRSS() / 1024 / 1024 );
    fileoutput.close();
    return true;
}

bool System::traceValid( Sparse &rho, double t_hit, bool force ) {
    double trace = std::real( getTrace<Scalar>( rho ) );
    parameters.trace.emplace_back( trace );
    if ( trace < 0.99 || trace > 1.01 || force ) {
        if ( force )
            fmt::print( "{} {} -> trace check failed at t = {} with trace(rho) = {}\n", PREFIX_ERROR, global_message_error_divergent, t_hit, trace );
        terminate_message = global_message_error_divergent;
        parameters.numerics_calculate_spectrum = 0;
        FILE *fp_trace = std::fopen( ( parameters.subfolder + "trace.txt" ).c_str(), "w" );
        for ( int i = 0; i < (int)parameters.trace.size() && parameters.t_step * 1.0 * i < t_hit; i++ ) {
            fmt::print( fp_trace, "{:.10e} {:.15e}\n", parameters.t_step * 1.0 * ( i + 1 ), parameters.trace.at( i ) );
        }
        std::fclose( fp_trace );
        return false;
    } else {
        return true;
    }
}