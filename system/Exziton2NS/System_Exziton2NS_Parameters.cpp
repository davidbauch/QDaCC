#include "../parameters.cpp"

class Parameters : public Parameters_Parent {
   public:
    // Mandatory parameters, move to parent class:
    // Time variables
    double t_start, t_end, t_step;
    // System Dimensions
    int maxStates, p_max_photon_number;
    // Starting state:
    int p_initial_state;
    bool startCoherent;

    // Non mandatory parameters, dependant on system chosen:
    // System Parameterss
    double p_omega_atomic;
    double p_omega_cavity;
    double p_omega_coupling;
    double p_omega_cavity_loss;
    double p_omega_pure_dephasing;
    double p_omega_decay;

    // Calculated System properties:
    double init_detuning, max_detuning, init_rabifrequenz, max_rabifrequenz;

    // Chirp and Pulse properties:
    std::vector<double> pulse_center, pulse_amp, pulse_omega, pulse_sigma;
    std::vector<std::string> pulse_type;
    double chirp_total;
    std::vector<double> chirp_t, chirp_y, chirp_ddt;
    std::string chirp_type;

    // Runtime parameters and other stuff
    int iterations_t_max;
    int iterations_total_max;
    std::vector<double> trace;
    bool output_full_dm;

    // AKF & Spectrum
    int numerics_maximum_threads, spectrum_frequency_iterations; //akf_everyXIt //akf_skip_omega
    double akf_deltaWmax, spectrum_frequency_center, spectrum_frequency_range;

    Parameters(){};
    Parameters( const std::vector<std::string> &arguments ) : Parameters_Parent() {
        init( arguments );
    }

    bool parseInput( const std::vector<std::string> &arguments ) {
        int index = 1;
        bool legacy = false;
        // Looking for legacy input, usefull because the python handler doesnt understand named inputs :(
        if ( ( index = vec_find_str( "-legacy", arguments ) ) != -1 ) {
            t_start = 0.0;
            t_end = getNextInput<double>( arguments, "t_end", ++index );
            t_step = getNextInput<double>( arguments, "t_step", index );
            p_omega_atomic = getNextInput<double>( arguments, "p_omega_atomic", index );
            p_omega_cavity = getNextInput<double>( arguments, "p_omega_cavity", index );
            p_omega_coupling = getNextInput<double>( arguments, "p_omega_coupling", index );
            p_omega_cavity_loss = getNextInput<double>( arguments, "p_omega_cavity_loss", index );
            p_omega_pure_dephasing = getNextInput<double>( arguments, "p_omega_pure_dephasing", index );
            p_omega_decay = getNextInput<double>( arguments, "p_omega_decay", index );

            chirp_t = getNextInputVector<double>( arguments, "chirp_time", index );
            chirp_y = getNextInputVector<double>( arguments, "chirp_yval", index );
            chirp_ddt = getNextInputVector<double>( arguments, "chirp_diff", index );
            chirp_type = getNextInputString( arguments, "chirp_file_type", index );

            pulse_center.emplace_back( getNextInput<double>( arguments, "pulse_center", index ) );
            pulse_amp.emplace_back( getNextInput<double>( arguments, "pulse_amp", index ) );
            pulse_omega.emplace_back( getNextInput<double>( arguments, "pulse_omega", index ) );
            pulse_sigma.emplace_back( getNextInput<double>( arguments, "pulse_sigma", index ) );
            pulse_type.emplace_back( getNextInputString( arguments, "pulse_type", index ) );

            p_max_photon_number = getNextInput<int>( arguments, "p_max_photon_number", index );
            p_initial_state = getNextInput<int>( arguments, "p_initial_state", index );

            numerics_calculate_spectrum = getNextInput<double>( arguments, "numerics_calculate_spectrum", index );
            iterations_skips_tau = getNextInput<int>( arguments, "iterations_skips_tau", index );
            spectrum_frequency_center = getNextInput<double>( arguments, "spectrum_frequency_center", index );
            spectrum_frequency_range = getNextInput<double>( arguments, "spectrum_frequency_range", index );
            iterations_skips_w = getNextInput<int>( arguments, "iterations_skips_w", index );

            numerics_use_interactionpicture = getNextInput<int>( arguments, "numerics_use_interactionpicture", index );
            numerics_use_rwa = getNextInput<int>( arguments, "numerics_use_rwa", index );
            numerics_maximum_threads = getNextInput<int>( arguments, "numerics_maximum_threads", index );
            double temp_rk_type = getNextInput<double>( arguments, "rungekuttatype", index );
            output_advanced_log = getNextInput<int>( arguments, "output_advanced_log", index );
            subfolder = getNextInputString( arguments, "subfolder", index );
            numerics_order_t = (int)( temp_rk_type / 10 );   // 4 or 5
            numerics_order_tau = ( (int)temp_rk_type ) % 10; // 4 or 5
            legacy = true;
        }

        // Look for --time, if not found, standard values are used (t0 = 0, t1 = 1ns, deltaT = auto)
        if ( ( index = vec_find_str( "--time", arguments ) ) != -1 ) {
            t_start = 0.0;
            t_end = getNextInput<double>( arguments, "t_end", ++index );
            t_step = getNextInput<double>( arguments, "t_step", index );
        } else if ( !legacy ) {
            t_start = 0.0;
            t_end = convertParam<double>( "1.0ns" );
            t_step = -1;
        }
        // Look for single parameter corrections
        if ( ( index = vec_find_str( "--tstart", arguments ) ) != -1 ) {
            t_start = getNextInput<double>( arguments, "t_start single", ++index );
        }
        if ( ( index = vec_find_str( "--tend", arguments ) ) != -1 ) {
            t_end = getNextInput<double>( arguments, "t_end single", ++index );
        }
        if ( ( index = vec_find_str( "--tstep", arguments ) ) != -1 ) {
            t_step = getNextInput<double>( arguments, "t_step single", ++index );
        }

        // Look for --system, if not found, standard system is used (g=66mueV, k=66mueV, p_omega_pure_dephasing = 3mueV, p_omega_decay = 1mueV)
        if ( ( index = vec_find_str( "--system", arguments ) ) != -1 ) {
            p_omega_atomic = getNextInput<double>( arguments, "p_omega_atomic", ++index );
            p_omega_cavity = getNextInput<double>( arguments, "p_omega_cavity", index );
            p_omega_coupling = getNextInput<double>( arguments, "p_omega_coupling", index );
            p_omega_cavity_loss = getNextInput<double>( arguments, "p_omega_cavity_loss", index );
            p_omega_pure_dephasing = getNextInput<double>( arguments, "p_omega_pure_dephasing", index );
            p_omega_decay = getNextInput<double>( arguments, "p_omega_decay", index );
        } else if ( !legacy ) {
            p_omega_atomic = convertParam<double>( "1.52eV" );
            p_omega_cavity = convertParam<double>( "1.52eV" );
            p_omega_coupling = convertParam<double>( "66mueV" );
            p_omega_cavity_loss = convertParam<double>( "66mueV" );
            p_omega_pure_dephasing = convertParam<double>( "3mueV" );
            p_omega_decay = convertParam<double>( "1mueV" );
        }
        // Look for single parameter corrections
        if ( ( index = vec_find_str( "--we", arguments ) ) != -1 ) {
            p_omega_atomic = getNextInput<double>( arguments, "p_omega_atomic single", ++index );
        }
        if ( ( index = vec_find_str( "--wc", arguments ) ) != -1 ) {
            p_omega_cavity = getNextInput<double>( arguments, "p_omega_cavity single", ++index );
        }
        if ( ( index = vec_find_str( "--coupling", arguments ) ) != -1 ) {
            p_omega_coupling = getNextInput<double>( arguments, "p_omega_coupling single", ++index );
        }
        if ( ( index = vec_find_str( "--kappa", arguments ) ) != -1 ) {
            p_omega_cavity_loss = getNextInput<double>( arguments, "p_omega_cavity_loss single", ++index );
        }
        if ( ( index = vec_find_str( "--gammapure", arguments ) ) != -1 ) {
            p_omega_pure_dephasing = getNextInput<double>( arguments, "p_omega_pure_dephasing single", ++index );
        }
        if ( ( index = vec_find_str( "--gamma", arguments ) ) != -1 ) {
            p_omega_decay = getNextInput<double>( arguments, "p_omega_decay single", ++index );
        }

        // Look for --chirp, if not found, standard system is used (no chirp, everything zero)
        if ( ( index = vec_find_str( "--chirp", arguments ) ) != -1 ) {
            chirp_t = getNextInputVector<double>( arguments, "chirp_time", ++index );
            chirp_y = getNextInputVector<double>( arguments, "chirp_yval", index );
            chirp_ddt = getNextInputVector<double>( arguments, "chirp_diff", index );
            chirp_type = getNextInputString( arguments, "chirp_file_type", index );
        } else if ( !legacy ) {
            chirp_t = {t_start, t_end};
            chirp_y = {0.0, 0.0};
            chirp_ddt = {0.0, 0.0};
            chirp_type = "monotone";
        }
        // Look for single parameter corrections //TODO: redundant
        if ( ( index = vec_find_str( "--chirpT", arguments ) ) != -1 ) {
            chirp_t = getNextInputVector<double>( arguments, "chirp_time single", ++index );
        }
        if ( ( index = vec_find_str( "--chirpY", arguments ) ) != -1 ) {
            chirp_y = getNextInputVector<double>( arguments, "chirp_yval single", ++index );
        }
        if ( ( index = vec_find_str( "--chirpDDT", arguments ) ) != -1 ) {
            chirp_ddt = getNextInputVector<double>( arguments, "chirp_ddt single", ++index );
        }
        if ( ( index = vec_find_str( "--chirpType", arguments ) ) != -1 ) {
            chirp_type = getNextInputString( arguments, "chirp_type single", ++index );
        }

        // Look for --pulse, if not found, standard system is used (no pulse, everything zero)
        if ( ( index = vec_find_str( "--pulse", arguments ) ) != -1 ) {
            // Single pulse:
            if ( arguments.at( index + 1 ).at( 0 ) == '[' ) {
                pulse_center = getNextInputVector<double>( arguments, "pulse_center", ++index );
                pulse_amp = getNextInputVector<double>( arguments, "pulse_amp", index );
                pulse_omega = getNextInputVector<double>( arguments, "pulse_omega", index );
                pulse_sigma = getNextInputVector<double>( arguments, "pulse_sigma", index );
                pulse_type = getNextInputVectorString( arguments, "pulse_type", index );
            } else {
                pulse_center.emplace_back( getNextInput<double>( arguments, "pulse_center", ++index ) );
                pulse_amp.emplace_back( getNextInput<double>( arguments, "pulse_amp", index ) );
                pulse_omega.emplace_back( getNextInput<double>( arguments, "pulse_omega", index ) );
                pulse_sigma.emplace_back( getNextInput<double>( arguments, "pulse_sigma", index ) );
                pulse_type.emplace_back( getNextInputString( arguments, "pulse_type", index ) );
            }
        } else if ( !legacy ) {
            pulse_center.emplace_back( 0.0 );
            pulse_amp.emplace_back( 0.0 );
            pulse_omega.emplace_back( 0.0 );
            pulse_sigma.emplace_back( 1.0 );
            pulse_type.emplace_back( "cw" );
        }
        if ( ( index = vec_find_str( "-pulse", arguments ) ) != -1 ) {
            pulse_amp.emplace_back( 1.0 );
            pulse_omega.emplace_back( p_omega_atomic );
            pulse_sigma.emplace_back( convertParam<double>( "20ps" ) );
            pulse_center.emplace_back( 10.0 * pulse_sigma.back() );
            pulse_type.emplace_back( "gauss_pi" );
            // Look for single parameter corrections
            if ( ( index = vec_find_str( "--pulseCenter", arguments ) ) != -1 ) {
                pulse_center.back() = getNextInput<double>( arguments, "pulse_center single", ++index );
            }
            if ( ( index = vec_find_str( "--pulseAmp", arguments ) ) != -1 ) {
                pulse_amp.back() = getNextInput<double>( arguments, "pulse_amp single", ++index );
            }
            if ( ( index = vec_find_str( "--pulseFreq", arguments ) ) != -1 ) {
                pulse_omega.back() = getNextInput<double>( arguments, "pulse_omega single", ++index );
            }
            if ( ( index = vec_find_str( "--pulseSigma", arguments ) ) != -1 ) {
                pulse_sigma.back() = getNextInput<double>( arguments, "pulse_sigma single", ++index );
            }
            if ( ( index = vec_find_str( "--pulseType", arguments ) ) != -1 ) {
                pulse_type.back() = getNextInputString( arguments, "pulse_type single", ++index );
            }
        }

        // Look for --dimensions, if not found, standard system is used (maxphotons = 0, starting state = |g,0>)
        if ( ( index = vec_find_str( "--dimensions", arguments ) ) != -1 ) {
            p_max_photon_number = getNextInput<int>( arguments, "p_max_photon_number", ++index );
            p_initial_state = getNextInput<int>( arguments, "p_initial_state", index );
        } else {
            p_max_photon_number = 1;
            p_initial_state = 0;
        }
        // Look for single parameter corrections
        if ( ( index = vec_find_str( "--maxPhotons", arguments ) ) != -1 ) {
            p_max_photon_number = getNextInput<int>( arguments, "p_max_photon_number single", ++index );
        }
        if ( ( index = vec_find_str( "--initState", arguments ) ) != -1 ) {
            p_initial_state = getNextInput<int>( arguments, "p_initial_state single", ++index );
        }

        // Look for (-RK4), -RK5, (-RK4T), (-RK4Tau), -RK5T, -RK5Tau
        numerics_order_t = 4;
        numerics_order_tau = 4;
        if ( ( index = vec_find_str( "-RK5", arguments ) ) != -1 ) {
            numerics_order_t = 5;
            numerics_order_tau = 5;
        }
        if ( ( index = vec_find_str( "-RK5T", arguments ) ) != -1 ) {
            numerics_order_t = 5;
        }
        if ( ( index = vec_find_str( "-RK5Tau", arguments ) ) != -1 ) {
            numerics_order_tau = 5;
        }
        numerics_order_highest = numerics_order_t;
        if ( numerics_order_tau > numerics_order_highest )
            numerics_order_highest = numerics_order_tau;

        // Look for --spectrum, if not found, no spectrum is evaluated
        if ( ( index = vec_find_str( "--spectrum", arguments ) ) != -1 ) {
            numerics_calculate_spectrum = 1;
            iterations_skips_tau = getNextInput<int>( arguments, "iterations_skips_tau", ++index );
            spectrum_frequency_center = getNextInput<double>( arguments, "spectrum_frequency_center", index );
            spectrum_frequency_range = getNextInput<double>( arguments, "spectrum_frequency_range", index );
            iterations_skips_w = getNextInput<int>( arguments, "iterations_skips_w", index );
        } else if ( !legacy ) {
            numerics_calculate_spectrum = 0;
            iterations_skips_tau = 1;
            spectrum_frequency_center = p_omega_cavity;
            spectrum_frequency_range = 0.0;
            iterations_skips_w = 1;
        }
        if ( ( index = vec_find_str( "-spectrum", arguments ) ) != -1 ) {
            numerics_calculate_spectrum = 1;
            iterations_skips_tau = 1;
            spectrum_frequency_center = p_omega_cavity;
            spectrum_frequency_range = ( p_omega_coupling + p_omega_cavity_loss + p_omega_decay + p_omega_pure_dephasing ) * 10.0;
            iterations_skips_w = 1;
        }
        // Look for single parameter corrections
        if ( ( index = vec_find_str( "--specTauSkip", arguments ) ) != -1 ) {
            iterations_skips_tau = getNextInput<int>( arguments, "iterations_skips_tau single", ++index );
        }
        if ( ( index = vec_find_str( "--specCenter", arguments ) ) != -1 ) {
            spectrum_frequency_center = getNextInput<double>( arguments, "spectrum_frequency_center single", ++index );
        }
        if ( ( index = vec_find_str( "--specRange", arguments ) ) != -1 ) {
            spectrum_frequency_range = getNextInput<double>( arguments, "spectrum_frequency_range single", ++index );
        }
        if ( ( index = vec_find_str( "--specWSkip", arguments ) ) != -1 ) {
            iterations_skips_w = getNextInput<int>( arguments, "iterations_skips_w single", ++index );
        }

        if ( ( index = vec_find_str( "-g2", arguments ) ) != -1 ) {
            numerics_calculate_g2 = true;
        } else {
            numerics_calculate_g2 = false;
        }

        // Look for other parameters
        if ( ( index = vec_find_str( "-noInteractionpic", arguments ) ) != -1 ) {
            numerics_use_interactionpicture = 0;
        } else if ( !legacy ) {
            numerics_use_interactionpicture = 1;
        }
        if ( ( index = vec_find_str( "-noRWA", arguments ) ) != -1 ) {
            numerics_use_rwa = 0;
        } else if ( !legacy ) {
            numerics_use_rwa = 1;
        }
        if ( ( index = vec_find_str( "--Threads", arguments ) ) != -1 ) {
            numerics_maximum_threads = getNextInput<int>( arguments, "numerics_maximum_threads", ++index );
        } else if ( !legacy ) {
            numerics_maximum_threads = 1;
        }
        if ( ( index = vec_find_str( "-noHandler", arguments ) ) != -1 ) {
            output_handlerstrings = 0;
        } else {
            output_handlerstrings = 1;
        }
        output_operators = 0;
        if ( ( index = vec_find_str( "-output_operators", arguments ) ) != -1 ) {
            output_operators = 2;
        }
        if ( ( index = vec_find_str( "-outputHamiltons", arguments ) ) != -1 ) {
            output_operators = 1;
        }
        if ( ( index = vec_find_str( "-outputOperatorsStop", arguments ) ) != -1 ) {
            output_operators = 3;
        }
        numerics_order_timetrafo = TIMETRANSFORMATION_PRECALCULATED;
        if ( ( index = vec_find_str( "-timeTrafoMatrixExponential", arguments ) ) != -1 ) {
            numerics_order_timetrafo = TIMETRANSFORMATION_MATRIXEXPONENTIAL;
        }
        if ( ( index = vec_find_str( "-startCoherent", arguments ) ) != -1 ) {
            startCoherent = true;
        } else {
            startCoherent = false;
        }
        if ( ( index = vec_find_str( "-fullDM", arguments ) ) != -1 ) {
            output_full_dm = true;
        } else {
            output_full_dm = false;
        }
        
        subfolder = arguments.back();
        return true;
    }

    bool adjustInput() {
        // Calculate/Recalculate some parameters:
        // Adjust pulse area if pulse_type is "gauss_pi"
        for ( int i = 0; i < (int)pulse_amp.size(); i++ )
            if ( pulse_type.at( i ).compare( "gauss_pi" ) == 0 ) {
                pulse_amp.at( i ) = M_PI / ( std::sqrt( 2.0 * M_PI ) * pulse_sigma.at( i ) );
            }

        // Calculate Rabi frequencies:
        init_rabifrequenz = rabiFrequency( p_omega_atomic - p_omega_cavity, p_omega_coupling, (int)( p_initial_state / 2.0 ) );
        max_rabifrequenz = rabiFrequency( p_omega_atomic - p_omega_cavity, p_omega_coupling, (int)( p_initial_state / 2.0 ) + 1 );

        // Calculate minimum step necessary to resolve Rabi-oscillation if step=-1
        if ( t_step == -1 ) {
            if ( numerics_use_rwa )
                t_step = 2. / 8. * M_PI / std::max( std::max( init_rabifrequenz, max_rabifrequenz ), p_omega_coupling + p_omega_cavity_loss + p_omega_decay + p_omega_pure_dephasing );
            if ( !numerics_use_rwa )
                t_step = 2. / 8. * M_PI / std::max( std::max( init_rabifrequenz, max_rabifrequenz ), p_omega_atomic + p_omega_cavity );
        }

        // Calculate the maximum dimensions for operator matrices (max states)
        maxStates = 2 * ( p_max_photon_number + 1 );

        // Calculate stuff for RK
        iterations_t_max = (int)std::ceil( ( t_end - t_start ) / t_step ); // FIXME: rendundant

        // Mandatory: rescale chirp ddt into chirp/ps
        for ( long unsigned i = 0; i < chirp_ddt.size(); i++ )
            chirp_ddt.at( i ) = chirp_ddt.at( i ) * 1E12;
        // Calculate values for chirp:
        chirp_total = vec_max( chirp_y );

        // Initial and maximum system detuning (not taking into account custom chirps)
        init_detuning = ( p_omega_atomic - p_omega_cavity );
        max_detuning = ( init_detuning + chirp_total > init_detuning ) ? init_detuning + chirp_total : init_detuning;

        // Adjust/calculate frequency range for spectrum
        spectrum_frequency_iterations = iterations_t_max / iterations_skips_w;
        if ( spectrum_frequency_center == -1 )
            spectrum_frequency_center = p_omega_cavity;
        if ( spectrum_frequency_range == -1 )
            spectrum_frequency_range = ( std::abs( max_detuning ) + p_omega_coupling + p_omega_cavity_loss / 2. ) * 3.;

        // Calculate total number of iterations necessary //FIXME: redundant
        iterations_total_max = iterations_t_max;
        if ( numerics_calculate_spectrum ) {
            int curIt = 1;
            for ( double ii = t_start + t_step; ii < t_end; ii += t_step ) {
                if ( curIt % iterations_skips_tau == 0 ) {
                    for ( double iii = ii + t_step; iii < t_end; iii += t_step ) {
                        iterations_total_max++;
                    }
                    curIt = 1;
                } else
                    curIt += 1;
            }
        }
        trace.reserve( iterations_t_max + 5 );
        return true;
    }

    void log() {
        logs.wrapInBar( "Parameters" );
        logs( "\n" );
        logs.wrapInBar( "Time borders and t_step", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Timeborder left t_start = {} s\n", t_start );
        logs( "Timeborder right t_end = {} s\n", t_end );
        logs( "Timeborder t_step delta t = {} s\n", t_step );
        logs( "Time iterations (main loop) = {}\n", iterations_t_max );
        logs( "Total time iterations = {}\n\n", iterations_total_max ); //, (int)ceil(iterations_t_max/2.*iterations_t_max/((double)akf_everyXIt)) );
        logs.wrapInBar( "System Parameters", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Energy Level difference |g><g| - |e><e| = {} Hz -> {} eV -> {} nm\n", p_omega_atomic, Hz_to_eV( p_omega_atomic ), Hz_to_wavelength( p_omega_atomic ) );
        logs( "Cavity Frequency w_c = {} Hz -> {} eV -> {} nm\n", p_omega_cavity, Hz_to_eV( p_omega_cavity ), Hz_to_wavelength( p_omega_cavity ) );
        logs( "Coupling strengh g = {} Hz -> {} mueV\n", p_omega_coupling, Hz_to_eV( p_omega_coupling ) * 1E6 );
        logs( "Photon loss rate k = {} Hz -> {} mueV -> Q = {:.2f}\n", p_omega_cavity_loss, Hz_to_eV( p_omega_cavity_loss ) * 1E6, p_omega_cavity / p_omega_cavity_loss );
        logs( "Atomic dephasing rate gamma_pure = {} Hz -> {} mueV\n", p_omega_pure_dephasing, Hz_to_eV( p_omega_pure_dephasing ) * 1E6 );
        logs( "RAD rate gamma = {} Hz -> {} mueV\n", p_omega_decay, Hz_to_eV( p_omega_decay ) * 1E6 );
        logs( "Initial state rho0 = |{},{}> with maximum number of {} photons\n\n", ( p_initial_state % 2 == 0 ? "g" : "e" ), (int)std::floor( p_initial_state / 2 ), p_max_photon_number );
        logs.wrapInBar( "Excitation Pulse", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        if ( pulse_amp.size() > 0 && pulse_center.at( 0 ) != -1 && pulse_amp.at( 0 ) != 0 ) {
            for ( int i = 0; i < (int)pulse_amp.size(); i++ ) {
                logs( "Exiting system at t_0 = {} with amplitude {} ({}meV), frequency {}eV ({}) and FWHM {}\n", pulse_center.at( i ), pulse_amp.at( i ), Hz_to_eV( pulse_amp.at( i ) ) * 1E3, pulse_omega.at( i ), Hz_to_eV( pulse_omega.at( i ) ), pulse_sigma.at( i ) * ( 2 * std::sqrt( 2 * std::log( 2 ) ) ) );
                logs( "Used pulse_type - " + pulse_type.at( i ) + "\n" );
            }
        } else
            logs( "Not using pulse to exite system\n\n" );
        logs.wrapInBar( "Energy Chirp", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        if ( chirp_total != 0 ) {
            for ( int i = 0; i < (int)chirp_t.size() - 1; i++ ) {
                //TODO: das falsch
                logs( "Chirp between t0 = {}ps and t1 = {}ps: {}mueV -> average rate {}mueV/ps\n", chirp_t.at( i ), chirp_t.at( i + 1 ), chirp_y.at( i + 1 ) - chirp_y.at( i ), ( chirp_y.at( i + 1 ) - chirp_y.at( i ) ) * 1e12 / ( chirp_t.at( i + 1 ) - chirp_t.at( i ) ) );
            }
            if ( chirp_type.compare( "none" ) != 0 )
                logs( "\nChirpfile of type '" + chirp_type + "' is used!" );
        } else
            logs( "Not using chirp" );
        logs( "\n\n" );

        logs.wrapInBar( "Caluclated Frequencies" );
        logs( "\nInitial system detuning = {} Hz -> {} mueV\n", init_detuning, Hz_to_eV( init_detuning ) * 1E6 );
        logs( "Maximum system detuning = {} Hz -> {} mueV\n", max_detuning, Hz_to_eV( max_detuning ) * 1E6 );
        logs( "Initial Rabi Frequencies = {} Hz -> {} mueV\n", init_rabifrequenz, Hz_to_eV( init_rabifrequenz ) * 1E6 );
        logs( "Maximum Rabi Frequencies = {} Hz -> {} mueV\n\n", max_rabifrequenz, Hz_to_eV( max_rabifrequenz ) * 1E6 );
        int works = 1;
        if ( ( init_rabifrequenz != 0.0 ) && ( 3. * t_step > 2. * M_PI / init_rabifrequenz ) )
            works = 0;
        else if ( max_rabifrequenz != 0.0 && 3. * t_step > 2. * M_PI / max_rabifrequenz )
            works = 0;
        if ( !works ) {
            fmt::print( "{} WARNING: Step may be too small to resolve predicted oscillation: dT needed vs dT: {:.10e} < {:.10e}\n", PREFIX_WARNING, 2. / 3. * M_PI / std::max( init_rabifrequenz, max_rabifrequenz ), t_step );
            logs( "WARNING: Step may be too small to resolve predicted oscillation: \n-> delta T needed: {:.10e} \n-> delta T used: {:.10e}\n\n", 2. / 3. * M_PI / std::max( init_rabifrequenz, max_rabifrequenz ), t_step );
        }
        logs.wrapInBar( "Spectrum" );
        logs( "\nCenter Frequency: {} Hz -> {} eV\n", spectrum_frequency_center, Hz_to_eV( spectrum_frequency_center ) );
        logs( "Frequency Range: +/- {} Hz -> +/- {} mueV\n", spectrum_frequency_range, Hz_to_eV( spectrum_frequency_range ) * 1E6 );
        logs( "Skipping {} steps in grid\n\n", iterations_skips_tau - 1 );

        logs.wrapInBar( "Numerics" );
        logs( "\nOrder of Runge-Kutta used: Time: RK{}, Spectrum: RK{}\n", numerics_order_t, numerics_order_tau );
        logs( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
        logs( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture == 1 ) ? "YES" : "NO" ) );
        logs( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == TIMETRANSFORMATION_PRECALCULATED ) ? "Analytic" : "Matrix Exponential" ) );
        logs( "Threads used for AFK? - {}\n", numerics_maximum_threads );
        logs( "\n" );
        logs.wrapInBar( "Program Log:", LOG_SIZE_FULL, LOG_LEVEL_2 );
        logs( "\n" );
        logs.level2( "OutputHandlerStrings: {}\n", output_handlerstrings );
    }

    static void help() {
        fmt::print( "--help, -help, -h\tThis screen\n" );
        fmt::print( "--time [start] [end] [step]\n\t--tstart [start]\n\t--tend [end]\n\t--tstep [step]\n" );
        fmt::print( "--system [p_omega_atomic] [p_omega_cavity] [p_omega_coupling] [kappa] [gammapure] [gamma] else standard values (1.52eV,1.52eV,66mueV,66mueV,3mueV,1mueV) are used\n\t--we [p_omega_atomic]\n\t--wc [p_omega_cavity]\n\t--coupling [p_omega_coupling]\n\t--kappa [kappa]\n\t--gamma [p_omega_decay]\n\t--gammapure [p_omega_pure_dephasing]\n" );
        fmt::print( "--chirp ['[Array Time]'] ['[Array Y]'] ['[Array d/dt]'] [type]\n\t--chirpT ['[Array Time]']\n\t--chirpY ['[Array Y]']\n\t--chirpDDT ['[Array d/dt]']\n\t--chirpType [type] where type = monotone, hermite, linear, spline\n" );
        fmt::print( "--pulse [Center] [Amplitude] [Frequency] [Sigma] [Type]\n\t-pulse for standard pulse\n\t--pulseCenter [Center]\n\t--pulseAmp [Amplitude]\n\t--pulseFreq [Frequency]\n\t--pulseSigma [Sigma]\n\t--pulseType [Type] where Type = cw, gauss, gauss_pi\n" );
        fmt::print( "--dimensions [maximum Photons] [Initial state]\n\t--maxPhotons [maximum Photons]\n\t--initState [Initial state], has to be smaller than (2*n+1)\n" );
        fmt::print( "--spectrum [Tau Skips] [Center] [Range] [Omega Skips] enables spectrum\n\t-spectrum enables spectrum centered at cavity\n\t--specTauSkip [Iterations skips (int)]\n\t--specCenter [Center]\n\t--specRange [Range]\n\t--specWSkip [Iteration skips (int)]\n" );
        fmt::print( "-RK5 enables Runge Kutta of order 5 for T and Tau direction\n\t-RK5T enables Runge Kutta of order 5 for T direction\n\t-RK5Tau enables Runge Kutta of order 5 for Tau direction\n" );
        fmt::print( "-noInteractionpic disables Interaction picture - enabled by default\n-noRWA disables rotating wave approximation - enabled by default\n-timeTrafoMatrixExponential enables Time Transformation via Matrix exponential - disabled by default\n-startCoherent enables starting with a coherent state. Starting state is then ground state with alpha = initState\n-fullDM enables full output of densitymatrix including all offdiagonal terms.\n" );
        fmt::print( "--Threads [number] number of threads to use for both AKF and Spectrum integral calculation\n" );
        fmt::print( "Additional commands:\n\t-advLog Enables advanced logging\n\t-noHandler disables handler strings and enables loadbar output (for console)\n\t-output_operators, -outputHamiltons, -outputOperatorsStop Enables output of matrices (requires -advLog)" );
    }
};