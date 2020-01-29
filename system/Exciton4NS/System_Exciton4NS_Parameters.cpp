#include "../../global.h"
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
    double p_omega_atomic_G_V;
    double p_omega_atomic_G_H;
    double p_omega_atomic_V_B;
    double p_omega_atomic_H_B;
    double p_omega_atomic_B;
    double p_omega_cavity_V;
    double p_omega_cavity_H;
    double p_omega_coupling;
    double p_omega_cavity_loss;
    double p_omega_pure_dephasing;
    double p_omega_decay;
    double p_phonon_b, p_phonon_alpha, p_phonon_wcutoff, p_phonon_T, p_phonon_tcutoff;
    bool p_phonon_adjust;
    double p_deltaE;
    double p_biexciton_bindingenergy;

    // Calculated System properties:
    double init_detuning_G_H, init_detuning_G_V, init_detuning_H_B, init_detuning_V_B, max_detuning_G_H, max_detuning_G_V, max_detuning_H_B, max_detuning_V_B;
    double init_rabifrequenz_G_H, init_rabifrequenz_G_V, init_rabifrequenz_H_B, init_rabifrequenz_V_B, max_rabifrequenz_G_H, max_rabifrequenz_G_V, max_rabifrequenz_H_B, max_rabifrequenz_V_B;
    double init_detuning, max_detuning, init_rabifrequenz, max_rabifrequenz;

    // Chirp and Pulse properties:
    std::vector<double> pulse_center, pulse_amp, pulse_omega, pulse_sigma, pulse_omega_chirp;
    std::vector<std::string> pulse_type;
    std::vector<std::string> pulse_pol;
    double chirp_total;
    std::vector<double> chirp_t, chirp_y, chirp_ddt;
    std::string chirp_type;

    // AKF & Spectrum
    int numerics_maximum_threads, iterations_w_resolution, numerics_phonon_approximation_1, numerics_phonon_approximation_2;
    double akf_deltaWmax, spectrum_frequency_center, spectrum_frequency_range;
    bool numerics_calculate_spectrum_H, numerics_calculate_spectrum_V;
    bool numerics_use_saved_coefficients;
    bool numerics_output_raman_population;
    int numerics_phonons_maximum_threads;
    bool numerics_use_saved_hamiltons;
    long unsigned int numerics_saved_coefficients_cutoff; // True: Only save last few coefficients (only viable for T-direction, not for G1/2)
    long unsigned int numerics_saved_coefficients_max_size;
    int logfilecounter;

    Parameters(){};
    Parameters( const std::vector<std::string> &arguments ) : Parameters_Parent() {
        init( arguments );
    }

    int index_to_state( const char mode = 'E', const int state = 0 ) {
        if ( mode == 'E' )
            return state % 4;
        if ( mode == 'H' )
            return (int)std::floor( state / ( 4 * ( p_max_photon_number + 1 ) ) );
        if ( mode == 'V' )
            return ( (int)std::floor( state / 4 ) ) % ( p_max_photon_number + 1 );
        return -1;
    }

    bool parseInput( const std::vector<std::string> &arguments ) {
        Parse_Parameters params;
        // Look for --time, if not found, standard values are used (t0 = 0, t1 = 1ns, deltaT = auto)
        params = Parse_Parameters( arguments, {"--time", "--tstart", "--tend", "--tstep"}, {3, 1, 1, 1}, "Timeparameters" );
        t_start = params.get<double>( {0, 3}, "0.0" );
        t_end = params.get<double>( {1, 4}, "1.0ns" );
        t_step = params.get<double>( {2, 5}, "-1" );

        // Look for --system, if not found, standard system is used (g=66mueV, k=66mueV, p_omega_pure_dephasing = 3mueV, p_omega_decay = 1mueV)
        params = Parse_Parameters( arguments, {"--system", "--we", "--wcH", "--wcV", "--coupling", "--kappa", "--gammapure", "--gamma", "--deltaE", "--excitonBindEnergy"}, {9, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Systemparameters" );
        p_omega_atomic_G_H = params.get<double>( {0, 9}, "1.366eV" );
        p_omega_cavity_V = params.get<double>( {1, 10}, "1.366eV" );
        p_omega_cavity_H = params.get<double>( {2, 11}, "1.366eV" );
        p_omega_coupling = params.get<double>( {3, 12}, "66mueV" );
        p_omega_cavity_loss = params.get<double>( {4, 13}, "66mueV" );
        p_omega_pure_dephasing = params.get<double>( {5, 14}, "3mueV" );
        p_omega_decay = params.get<double>( {6, 15}, "1mueV" );
        p_deltaE = params.get<double>( {7, 16}, "2mueV" );
        p_biexciton_bindingenergy = params.get<double>( {8, 17}, "3meV" );

        // Look for --chirp, if not found, standard system is used (no chirp, everything zero)
        params = Parse_Parameters( arguments, {"--chirp", "--chirpT", "--chirpY", "--chirpDDT", "--chirpType"}, {4, 1, 1, 1, 1}, "Chirpparameters" );
        chirp_t = convertParam<double>( str_to_vec( params.get( {0, 4}, "[0,1E-9]" ) ) );
        chirp_y = convertParam<double>( str_to_vec( params.get( {1, 5}, "[0,0]" ) ) );
        chirp_ddt = convertParam<double>( str_to_vec( params.get( {2, 6}, "[0,0]" ) ) );
        chirp_type = params.get( {3, 7}, "monotone" );

        // Look for --pulse, if not found, standard system is used (no pulse, everything zero)
        params = Parse_Parameters( arguments, {"--pulse", "--pulseCenter", "--pulseAmp", "--pulseFreq", "--pulseSigma", "--pulseType", "--pulsePol", "--pulseChirp", "-pulse"}, {7, 1, 1, 1, 1, 1, 1, 1, 1}, "Pulseparameters" );
        if ( params.get( 0, "Empty" ).at( 0 ) == '[' ) {
            pulse_center = convertParam<double>( str_to_vec( params.get( {0, 7}, "[0]" ) ) );
            pulse_amp = convertParam<double>( str_to_vec( params.get( {1, 8}, "[0]" ) ) );
            pulse_omega = convertParam<double>( str_to_vec( params.get( {2, 9}, "[0]" ) ) );
            pulse_sigma = convertParam<double>( str_to_vec( params.get( {3, 10}, "[1]" ) ) );
            pulse_type = str_to_vec( params.get( {4, 11}, "[cw]" ) );
            pulse_pol = str_to_vec( params.get( {5, 12}, "[H]" ) );
            pulse_omega_chirp = convertParam<double>( str_to_vec( params.get( {6, 13}, "[0]" ) ) );
        } else if ( params.get( 14 ) || params.get( 0 ) ) {
            pulse_center.emplace_back( params.get<double>( {0, 7}, "100ps" ) );
            pulse_amp.emplace_back( params.get<double>( {1, 8}, "4" ) );
            pulse_omega.emplace_back( params.get<double>( {2, 9}, "1.366eV" ) );
            pulse_sigma.emplace_back( params.get<double>( {3, 10}, "20ps" ) );
            pulse_type.emplace_back( params.get( {4, 11}, "gauss_pi" ) );
            pulse_pol.emplace_back( params.get( {4, 12}, "H" ) );
            pulse_omega_chirp.emplace_back( params.get<double>( {4, 13}, "0" ) );
        } else {
            pulse_center.emplace_back( convertParam<double>( "0" ) );
            pulse_amp.emplace_back( convertParam<double>( "0" ) );
            pulse_omega.emplace_back( convertParam<double>( "0" ) );
            pulse_sigma.emplace_back( convertParam<double>( "1" ) );
            pulse_type.emplace_back( "cw" );
            pulse_pol.emplace_back( "H" );
            pulse_omega_chirp.emplace_back( convertParam<double>( "0" ) );
        }

        // Look for --dimensions, if not found, standard system is used (maxphotons = 0, starting state = |g,0>)
        params = Parse_Parameters( arguments, {"--maxPhotons", "--initState"}, {1, 3}, "Initial State parameters" );
        p_max_photon_number = params.get<int>( 0, "2" );
        p_initial_state = instr( "ghvb", params.get( 1, "b" ) ) + 4 * ( p_max_photon_number + 1 ) * params.get<int>( 2, "0" ) + 4 * params.get<int>( 3, "0" );

        // Look for --spectrum, if not found, no spectrum is evaluated
        params = Parse_Parameters( arguments, {"--specTauRes", "--specCenter", "--specRange", "--specWRes", "-spectrum", "-spectrumH", "-spectrumV"}, {1, 1, 1, 1, 1, 1, 1}, "Spectrum Parameters" );
        iterations_tau_resolution = params.get<int>( 0, "250" );
        spectrum_frequency_center = params.get<double>( 1, std::to_string( p_omega_cavity_H ) );
        spectrum_frequency_range = params.get<double>( 2, std::to_string( ( p_omega_coupling + p_omega_cavity_loss + p_omega_decay + p_omega_pure_dephasing ) * 10.0 ) );
        iterations_w_resolution = params.get<int>( 3, "500" );
        numerics_calculate_spectrum_H = params.get( 4 ) || params.get( 5 );
        numerics_calculate_spectrum_V = params.get( 4 ) || params.get( 6 );

        // Look for (-RK4), -RK5, (-RK4T), (-RK4Tau), -RK5T, -RK5Tau
        params = Parse_Parameters( arguments, {"-g2", "-RK5", "-RK5T", "-RK5Tau", "-noInteractionpic", "-noRWA", "--Threads", "-noHandler", "-outputOperators", "-outputHamiltons", "-outputOperatorsStop", "-timeTrafoMatrixExponential", "-startCoherent", "-fullDM", "-scale", "-disableMatrixCaching", "-disableHamiltonCaching", "-disableMainProgramThreading", "-noRaman", "-g2s", "--lfc"}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Other parameters" );
        numerics_calculate_g2 = params.get( 0 ) || params.get( 19 );
        numerics_order_t = ( params.get( 1 ) || params.get( 2 ) ) ? 5 : 4;
        numerics_order_tau = ( params.get( 1 ) || params.get( 3 ) ) ? 5 : 4;
        numerics_order_highest = numerics_order_t;
        if ( numerics_order_tau > numerics_order_highest )
            numerics_order_highest = numerics_order_tau;
        numerics_use_interactionpicture = params.get( 4 ) ? 0 : 1;
        numerics_use_rwa = params.get( 5 ) ? 0 : 1;
        numerics_maximum_threads = params.get<int>( 6, "1" );
        if ( numerics_maximum_threads == -1 )
            numerics_maximum_threads = omp_get_max_threads();
        output_handlerstrings = params.get( 7 ) ? 0 : 1;
        output_operators = params.get( 8 ) ? 2 : ( params.get( 9 ) ? 1 : ( params.get( 10 ) ? 3 : 0 ) );
        numerics_order_timetrafo = params.get( 11 ) ? TIMETRANSFORMATION_MATRIXEXPONENTIAL : TIMETRANSFORMATION_ANALYTICAL;
        startCoherent = params.get( 12 );
        output_full_dm = params.get( 13 );
        scale_parameters = params.get( 14 );
        numerics_use_saved_coefficients = !params.get( 15 );
        numerics_use_saved_hamiltons = !params.get( 16 );
        numerics_phonons_maximum_threads = ( !numerics_use_saved_coefficients || !params.get( 17 ) ) ? numerics_maximum_threads : 1;
        numerics_output_raman_population = !params.get( 18 );
        numerics_use_simplified_g2 = params.get( 19 );
        logfilecounter = params.get<int>( 20, "-1" );

        // Phonon Parameters
        params = Parse_Parameters( arguments, {"--phonons", "--temperature", "-phonons", "--phononorder", "-noMarkov", "-phononcoeffs", "-noPhononAdjust"}, {5, 1, 1, 1, 1, 1, 1} );
        p_phonon_alpha = params.get<double>( 0, "0.03E-24" );
        p_phonon_wcutoff = params.get<double>( 1, "1meV" );
        p_phonon_tcutoff = params.get<double>( 2, "4ps" );
        p_phonon_T = params.get( 6 ) ? params.get<double>( 5, "3" ) : params.get<double>( {3, 5}, "0" );
        numerics_phonon_approximation_2 = params.get<int>( {4, 7}, std::to_string( PHONON_APPROXIMATION_BACKWARDS_INTEGRAL ) ); // Second Markov / Transformation method
        numerics_phonon_approximation_1 = params.get( 8 ) ? 0 : 1;                                                              // First Markov
        output_coefficients = params.get( 9 ) ? 1 : 0;
        p_phonon_adjust = !params.get( 10 );

        subfolder = arguments.back();
        return true;
    }

    bool scaleInputs( const double scaling ) {
        // Adjust normal parameters: time is multiplid by scaling, frequency divided
        t_start = scaleVariable( t_start, scaling );
        t_end = scaleVariable( t_end, scaling );
        t_step = scaleVariable( t_step, scaling );
        p_omega_atomic_G_H = scaleVariable( p_omega_atomic_G_H, 1.0 / scaling );
        p_omega_atomic_G_V = scaleVariable( p_omega_atomic_G_V, 1.0 / scaling );
        p_omega_atomic_H_B = scaleVariable( p_omega_atomic_H_B, 1.0 / scaling );
        p_omega_atomic_V_B = scaleVariable( p_omega_atomic_V_B, 1.0 / scaling );
        p_omega_cavity_H = scaleVariable( p_omega_cavity_H, 1.0 / scaling );
        p_omega_cavity_V = scaleVariable( p_omega_cavity_V, 1.0 / scaling );
        p_omega_coupling = scaleVariable( p_omega_coupling, 1.0 / scaling );
        p_omega_cavity_loss = scaleVariable( p_omega_cavity_loss, 1.0 / scaling );
        p_omega_pure_dephasing = scaleVariable( p_omega_pure_dephasing, 1.0 / scaling );
        p_omega_decay = scaleVariable( p_omega_decay, 1.0 / scaling );
        // Adjust chirp and pulse
        for ( int i = 0; i < (int)chirp_t.size(); i++ ) {
            chirp_t.at( i ) = scaleVariable( chirp_t.at( i ), scaling );
            chirp_y.at( i ) = scaleVariable( chirp_y.at( i ), 1.0 / scaling );
            chirp_ddt.at( i ) = scaleVariable( chirp_ddt.at( i ), 1.0 / scaling );
        }
        for ( int i = 0; i < (int)pulse_center.size(); i++ ) {
            pulse_center.at( i ) = scaleVariable( pulse_center.at( i ), scaling );
            pulse_amp.at( i ) = scaleVariable( pulse_amp.at( i ), 1.0 / scaling );
            pulse_omega.at( i ) = scaleVariable( pulse_omega.at( i ), 1.0 / scaling );
            pulse_sigma.at( i ) = scaleVariable( pulse_sigma.at( i ), scaling );
        }
        // Adjusting spectrum
        spectrum_frequency_center = scaleVariable( spectrum_frequency_center, 1.0 / scaling );
        spectrum_frequency_range = scaleVariable( spectrum_frequency_range, 1.0 / scaling );
        return true;
    }

    double scaleVariable( const double variable, const double scaling ) {
        if ( scale_parameters ) {
            return variable * scaling;
        }
        return variable;
    }

    double getIdealTimestep() {
        if ( numerics_use_rwa )
            return 2. / 8. * M_PI / std::max( std::max( init_rabifrequenz, max_rabifrequenz ), p_omega_coupling + p_omega_cavity_loss + p_omega_decay + p_omega_pure_dephasing + ( vec_max( pulse_amp ) > 0 ? std::abs( p_omega_atomic_G_H - vec_max( pulse_omega ) ) : 0 ) );
        if ( !numerics_use_rwa )
            return 2. / 8. * M_PI / std::max( std::max( init_rabifrequenz, max_rabifrequenz ), p_omega_atomic_G_H + p_omega_cavity_H );
        return 1E-13;
    }

    bool adjustInput() {
        // Calculate/Recalculate some parameters:
        // Adjust pulse area if pulse_type is "gauss_pi"
        for ( int i = 0; i < (int)pulse_amp.size(); i++ )
            if ( pulse_type.at( i ).compare( "gauss_pi" ) == 0 ) {
                if ( pulse_amp.at( i ) < 500 )
                    pulse_amp.at( i ) = pulse_amp.at( i ) * M_PI / ( std::sqrt( 2.0 * M_PI * pulse_sigma.at( i ) * std::sqrt( std::pow( pulse_omega_chirp.at( i ) / pulse_sigma.at( i ), 2.0 ) + std::pow( pulse_sigma.at( i ), 2.0 ) ) ) ) / 2.0; //https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.241306
                pulse_type.at( i ) = "gauss";
            }

        // Calculating remaining atomic frequencies depending on delta E and biexciton binding energy.
        p_omega_atomic_G_V = p_omega_atomic_G_H + p_deltaE / 2.0;
        p_omega_atomic_B = 2.0 * p_omega_atomic_G_H - p_biexciton_bindingenergy;
        p_omega_atomic_G_H = p_omega_atomic_G_H - p_deltaE / 2.0;
        p_omega_atomic_H_B = p_omega_atomic_B - p_omega_atomic_G_H;
        p_omega_atomic_V_B = p_omega_atomic_B - p_omega_atomic_G_V;

        // Calculate Rabi frequencies:
        init_rabifrequenz_G_H = rabiFrequency( p_omega_atomic_G_H - p_omega_cavity_H, p_omega_coupling, index_to_state( 'H', p_initial_state ) );
        init_rabifrequenz_G_V = rabiFrequency( p_omega_atomic_G_V - p_omega_cavity_V, p_omega_coupling, index_to_state( 'V', p_initial_state ) );
        init_rabifrequenz_H_B = rabiFrequency( p_omega_atomic_H_B - p_omega_cavity_H, p_omega_coupling, index_to_state( 'H', p_initial_state ) );
        init_rabifrequenz_V_B = rabiFrequency( p_omega_atomic_V_B - p_omega_cavity_V, p_omega_coupling, index_to_state( 'V', p_initial_state ) );
        max_rabifrequenz_G_H = rabiFrequency( p_omega_atomic_G_H - p_omega_cavity_H, p_omega_coupling, index_to_state( 'H', p_initial_state + 1 ) );
        max_rabifrequenz_G_V = rabiFrequency( p_omega_atomic_G_V - p_omega_cavity_V, p_omega_coupling, index_to_state( 'V', p_initial_state + 1 ) );
        max_rabifrequenz_H_B = rabiFrequency( p_omega_atomic_H_B - p_omega_cavity_H, p_omega_coupling, index_to_state( 'H', p_initial_state + 1 ) );
        max_rabifrequenz_V_B = rabiFrequency( p_omega_atomic_V_B - p_omega_cavity_V, p_omega_coupling, index_to_state( 'V', p_initial_state + 1 ) );
        init_rabifrequenz = vec_max<double>( {init_rabifrequenz_G_H, init_rabifrequenz_G_V, init_rabifrequenz_H_B, init_rabifrequenz_V_B} );
        max_rabifrequenz = vec_max<double>( {max_rabifrequenz_G_H, max_rabifrequenz_G_V, max_rabifrequenz_H_B, max_rabifrequenz_V_B} );

        // Calculate minimum step necessary to resolve Rabi-oscillation if step=-1
        if ( t_step == -1 ) {
            t_step = std::min( 1E-13, getIdealTimestep() );
            t_step = std::max( std::numeric_limits<double>::epsilon(), t_step );
        }

        // Calculate the maximum dimensions for operator matrices (max states)
        maxStates = 4 * ( p_max_photon_number + 1 ) * ( p_max_photon_number + 1 ); // 4 Electronic states x N+1 Photonic states * 2 for H,V modes

        // Calculate stuff for RK
        iterations_t_max = (int)std::ceil( ( t_end - t_start ) / t_step );
        iterations_t_skip = std::max( 1.0, std::ceil( iterations_t_max / iterations_tau_resolution ) );

        // Mandatory: rescale chirp ddt into chirp/ps
        for ( long unsigned i = 0; i < chirp_ddt.size(); i++ )
            chirp_ddt.at( i ) = chirp_ddt.at( i ) * 1E12;
        // Calculate values for chirp:
        chirp_total = vec_max( chirp_y );

        // Initial and maximum system detuning (not taking into account custom chirps)
        init_detuning_G_H = ( p_omega_atomic_G_H - p_omega_cavity_H );
        init_detuning_G_V = ( p_omega_atomic_G_V - p_omega_cavity_V );
        init_detuning_H_B = ( p_omega_atomic_H_B - p_omega_cavity_H );
        init_detuning_V_B = ( p_omega_atomic_V_B - p_omega_cavity_V );
        max_detuning_G_H = ( init_detuning_G_H + chirp_total > init_detuning_G_H ) ? init_detuning_G_H + chirp_total : init_detuning_G_H;
        max_detuning_G_V = ( init_detuning_G_V + chirp_total > init_detuning_G_V ) ? init_detuning_G_V + chirp_total : init_detuning_G_V;
        max_detuning_H_B = ( init_detuning_H_B + chirp_total > init_detuning_H_B ) ? init_detuning_H_B + chirp_total : init_detuning_H_B;
        max_detuning_V_B = ( init_detuning_V_B + chirp_total > init_detuning_V_B ) ? init_detuning_V_B + chirp_total : init_detuning_V_B;
        init_detuning = vec_max<double>( {init_detuning_G_H, init_detuning_G_V, init_detuning_H_B, init_detuning_V_B} );
        max_detuning = vec_max<double>( {max_detuning_G_H, max_detuning_G_V, max_detuning_H_B, max_detuning_V_B} );

        // Adjust/calculate frequency range for spectrum
        if ( spectrum_frequency_center == -1 )
            spectrum_frequency_center = p_omega_cavity_H;
        if ( spectrum_frequency_range == -1 )
            spectrum_frequency_range = ( std::abs( max_detuning ) + p_omega_coupling + p_omega_cavity_loss / 2. ) * 3.;

        // Calculate phonon stuff
        p_phonon_b = 1.0;
        if ( p_phonon_T > 0 ) {
            double integral = 0;
            double stepsize = 0.01 * p_phonon_wcutoff;
            for ( double w = stepsize; w < 10 * p_phonon_wcutoff; w += stepsize ) {
                integral += stepsize * ( p_phonon_alpha * w * std::exp( -w * w / 2.0 / p_phonon_wcutoff / p_phonon_wcutoff ) / std::tanh( 1.0545718E-34 * w / 2.0 / 1.3806488E-23 / p_phonon_T ) );
            }
            p_phonon_b = std::exp( -0.5 * integral );
            if ( p_phonon_adjust ) {
                p_omega_pure_dephasing = convertParam<double>( "1mueV" ) * p_phonon_T;
                p_omega_decay *= p_phonon_b * p_phonon_b;
            }
        }
        numerics_calculate_spectrum = numerics_calculate_spectrum_H || numerics_calculate_spectrum_V;
        numerics_saved_coefficients_max_size = (int)( ( t_end - t_start ) / t_step * 2.0 * ( p_phonon_tcutoff / t_step ) ) + 10;
        trace.reserve( iterations_t_max + 5 );

        numerics_saved_coefficients_cutoff = ( numerics_calculate_spectrum || numerics_calculate_g2 ) ? 0 : ( p_phonon_tcutoff / t_step ) * 5;

        return true;
    }

    /*
        New Parameter Log:
        System Parameters
        - Base Energies
        - Couplings
        - Init State
        - Detunings
        - Pulse
        -- Different Pulses
        - Chirp
        -- Different Chirps
        Phonon Parameters
        Photon Statistics
        - G1
        -- Spectrum
        - G2
        -- all other photon statistics and if they are calculated
        Numerical Parameters
        - Time
        - Used stuff
        unterpunkte mit
        = Unterpunkt
        == UnterUnterpunk (z.b. chirp parameter)
        markieren für auswertunsgskritp
    */
    void log( const std::vector<std::string> &info ) {
        logs.wrapInBar( "System Parameters" );
        logs( "Version: 1.2.3\n\n" );

        logs.wrapInBar( "Base Energies", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Biexciton binding energy: {:.8e} Hz - {:.8} meV\n", p_biexciton_bindingenergy, Hz_to_eV( p_biexciton_bindingenergy ) * 1E3 );
        logs( "Exciton fine structure splitting: {:.8e} Hz - {:.8} meV\n", p_biexciton_bindingenergy, Hz_to_eV( p_biexciton_bindingenergy ) * 1E3 );
        logs( "Transition Energy H - G: {:.8e} Hz - {:.8} eV - {:.8} nm\n", p_omega_atomic_G_H, Hz_to_eV( p_omega_atomic_G_H ), Hz_to_wavelength( p_omega_atomic_G_H ) );
        logs( "Transition Energy V - G: {:.8e} Hz - {:.8} eV - {:.8} nm\n", p_omega_atomic_G_V, Hz_to_eV( p_omega_atomic_G_V ), Hz_to_wavelength( p_omega_atomic_G_V ) );
        logs( "Transition Energy B - H: {:.8e} Hz - {:.8} eV - {:.8} nm\n", p_omega_atomic_H_B, Hz_to_eV( p_omega_atomic_H_B ), Hz_to_wavelength( p_omega_atomic_H_B ) );
        logs( "Transition Energy B - V: {:.8e} Hz - {:.8} eV - {:.8} nm\n", p_omega_atomic_V_B, Hz_to_eV( p_omega_atomic_V_B ), Hz_to_wavelength( p_omega_atomic_V_B ) );
        logs( "Biexciton Resonance: {:.8e} Hz - {:.8} eV - {:.8} nm\n", p_omega_atomic_B, Hz_to_eV( p_omega_atomic_B ), Hz_to_wavelength( p_omega_atomic_B ) );
        logs( "Two Photon Resonance: {:.8e} Hz - {:.8} eV - {:.8} nm\n", p_omega_atomic_B / 2.0, Hz_to_eV( p_omega_atomic_B / 2.0 ), Hz_to_wavelength( p_omega_atomic_B / 2.0 ) );
        logs( "Cavity Energy H: {:.8e} Hz - {:.8} eV - {:.8} nm\n", p_omega_cavity_H, Hz_to_eV( p_omega_cavity_H ), Hz_to_wavelength( p_omega_cavity_H ) );
        logs( "Cavity Energy V: {:.8e} Hz - {:.8} eV - {:.8} nm\n\n", p_omega_cavity_V, Hz_to_eV( p_omega_cavity_V ), Hz_to_wavelength( p_omega_cavity_V ) );

        logs.wrapInBar( "Coupling Parameters", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Coupling strengh g: {:.8e} Hz - {:.8} mueV\n", p_omega_coupling, Hz_to_eV( p_omega_coupling ) * 1E6 );
        logs( "Photon loss rate k: {:.8e} Hz - {:.8} mueV - Q = {:.2f}\n", p_omega_cavity_loss, Hz_to_eV( p_omega_cavity_loss ) * 1E6, p_omega_cavity_H / p_omega_cavity_loss );
        logs( "Atomic dephasing rate gamma_pure: {:.8e} Hz - {:.8} mueV\n", p_omega_pure_dephasing, Hz_to_eV( p_omega_pure_dephasing ) * 1E6 );
        logs( "RAD rate gamma: {:.8e} Hz - {:.8} mueV\n\n", p_omega_decay, Hz_to_eV( p_omega_decay ) * 1E6 );

        logs.wrapInBar( "Initial System Parameters", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Initial state rho0 = |{0}><{0}| with maximum number of {1} photons\n", info.at( p_initial_state ), p_max_photon_number );
        logs( "Initial detuning X_H - G - H = {:.8e} Hz - {:.8} mueV, Maximum detuning X_H - G - H = {:.8e} Hz - {:.8} mueV\n", init_detuning_G_H, Hz_to_eV( init_detuning_G_H ) * 1E6, max_detuning_G_H, Hz_to_eV( max_detuning_G_H ) * 1E6 );
        logs( "Initial detuning X_V - G - V = {:.8e} Hz - {:.8} mueV, Maximum detuning X_V - G - V = {:.8e} Hz - {:.8} mueV\n", init_detuning_G_V, Hz_to_eV( init_detuning_G_V ) * 1E6, max_detuning_G_V, Hz_to_eV( max_detuning_G_V ) * 1E6 );
        logs( "Initial detuning B - X_H - H = {:.8e} Hz - {:.8} mueV, Maximum detuning B - X_H - H = {:.8e} Hz - {:.8} mueV\n", init_detuning_H_B, Hz_to_eV( init_detuning_H_B ) * 1E6, max_detuning_H_B, Hz_to_eV( max_detuning_H_B ) * 1E6 );
        logs( "Initial detuning B - X_V - V = {:.8e} Hz - {:.8} mueV, Maximum detuning B - X_V - V = {:.8e} Hz - {:.8} mueV\n", init_detuning_V_B, Hz_to_eV( init_detuning_V_B ) * 1E6, max_detuning_V_B, Hz_to_eV( max_detuning_V_B ) * 1E6 );

        logs( "Initial Rabi Frequency X_H - G = {:.8e} Hz - {:.8} mueV, Maximum Rabi Frequency X_H - B = {:.8e} Hz - {:.8} mueV\n", init_rabifrequenz_G_H, Hz_to_eV( init_rabifrequenz_G_H ) * 1E6, max_rabifrequenz_G_H, Hz_to_eV( max_rabifrequenz_G_H ) * 1E6 );
        logs( "Initial Rabi Frequency X_V - G = {:.8e} Hz - {:.8} mueV, Maximum Rabi Frequency X_V - G = {:.8e} Hz - {:.8} mueV\n", init_rabifrequenz_G_V, Hz_to_eV( init_rabifrequenz_G_V ) * 1E6, max_rabifrequenz_G_V, Hz_to_eV( max_rabifrequenz_G_V ) * 1E6 );
        logs( "Initial Rabi Frequency B - X_H = {:.8e} Hz - {:.8} mueV, Maximum Rabi Frequency B - X_H = {:.8e} Hz - {:.8} mueV\n", init_rabifrequenz_H_B, Hz_to_eV( init_rabifrequenz_H_B ) * 1E6, max_rabifrequenz_H_B, Hz_to_eV( max_rabifrequenz_H_B ) * 1E6 );
        logs( "Initial Rabi Frequency B - X_V = {:.8e} Hz - {:.8} mueV, Maximum Rabi Frequency B - X_V = {:.8e} Hz - {:.8} mueV\n\n", init_rabifrequenz_V_B, Hz_to_eV( init_rabifrequenz_V_B ) * 1E6, max_rabifrequenz_V_B, Hz_to_eV( max_rabifrequenz_V_B ) * 1E6 );

        logs.wrapInBar( "Pulse", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        if ( pulse_amp.size() > 0 && pulse_center.at( 0 ) != -1 && pulse_amp.at( 0 ) != 0 ) {
            for ( int i = 0; i < (int)pulse_amp.size(); i++ ) {
                logs( "Exiting system at t_0 = {:.8e}\nAmplitude {:.10e} ({:.8} meV, {:.2f} pi)\nFrequency {:.10e} ({:.8}eV)\nFWHM {:.8e}\nPulseChirp: {:.8e}\n", pulse_center.at( i ), pulse_amp.at( i ), Hz_to_eV( pulse_amp.at( i ) ) * 1E3, pulse_amp.at( i ) / ( M_PI / ( std::sqrt( 2.0 * M_PI ) * pulse_sigma.at( i ) ) / 2.0 ), pulse_omega.at( i ), Hz_to_eV( pulse_omega.at( i ) ), pulse_sigma.at( i ) * ( 2 * std::sqrt( 2 * std::log( 2 ) ) ), pulse_omega_chirp.at( i ) );
                logs( "Used pulse_type - {}\nUsed pulse_pol on mode {}", pulse_type.at( i ), pulse_pol.at( i ) );
            }
            logs( "\n\n" );
        } else
            logs( "Not using pulse to exite system\n\n" );
        logs.wrapInBar( "Chirp", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        if ( chirp_total != 0 ) {
            double total = 0;
            int current = 0;
            for ( int i = 0; i < (int)chirp_t.size() - 1; i++ ) {
                if ( chirp_y.at( i + 1 ) - chirp_y.at( i ) != 0.0 ) {
                    //logs.inBar( "Chirp " + toStr( current++ ) );
                    logs( "Chirp {0:} between t0 = {1:.8e} ps\nChirp {0:} t1 = {2:.8e} ps\nTotal Chirp {0:}: {3:.8} mueV\n-> average rate Chirp {0:}: {4:.8} mueV/ps\n",current, chirp_t.at( i ), chirp_t.at( i + 1 ), Hz_to_eV( chirp_y.at( i + 1 ) - chirp_y.at( i ) ) * 1E6, Hz_to_eV( ( chirp_y.at( i + 1 ) - chirp_y.at( i ) ) ) * 1E6 / 1E12 / ( chirp_t.at( i + 1 ) - chirp_t.at( i ) ) );
                    total += chirp_y.at( i + 1 ) - chirp_y.at( i );
                    current++;
                }
            }
            if ( chirp_type.compare( "none" ) != 0 )
                logs( "\nChirpfile of type '" + chirp_type + "' is used!\n" );
            logs( "Total Chirp = {:.8} mueV\n\n", Hz_to_eV( total ) * 1E6 );
        } else
            logs( "Not using chirp\n\n" );

        logs.wrapInBar( "Photon Statistics" );
        logs( "\n" );
        logs.wrapInBar( "G1 Settings", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Anticipated tau-grid resolution is {}x{} resulting in {} skips per timestep\n", iterations_tau_resolution, iterations_tau_resolution, iterations_t_skip );
        if ( numerics_calculate_spectrum_H || numerics_calculate_spectrum_V ) {
            logs( "Calcluating spectrum for modes {}\n", ( numerics_calculate_spectrum_H && numerics_calculate_spectrum_V ? "cavities H and V" : fmt::format( "cavity {}", numerics_calculate_spectrum_H ? "H" : "V" ) ) );
            logs( "Center Frequency: {:.8e} Hz - {:.8} eV\n", spectrum_frequency_center, Hz_to_eV( spectrum_frequency_center ) );
            logs( "Frequency Range: +/- {:.8e} Hz - +/- {:.8} mueV\n", spectrum_frequency_range, Hz_to_eV( spectrum_frequency_range ) * 1E6 );
            logs( "Maximum w-vector resolution is {}\n", iterations_w_resolution );
        } else {
            logs( "Not calculating spectrum\n" );
        }
        logs( "\n" );
        logs.wrapInBar( "G2 Settings", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Calculate advanced statistics? - {}\n", numerics_calculate_g2 ? fmt::format( "YES{}", ( numerics_use_simplified_g2 ? " (simplified)" : "" ) ) : "NO" );

        logs.wrapInBar( "Phonons", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        if ( p_phonon_T ) {
            std::vector<std::string> approximations = {"Transformation integral via d/dt chi = -i/hbar*[H,chi] + d*chi/dt onto interaction picture chi(t-tau)", "Transformation Matrix U(t,tau)=exp(-i/hbar*H_DQ_L(t)*tau) onto interaction picture chi(t-tau)", "No Transformation, only interaction picture chi(t-tau)", "Analytical Lindblad formalism"};
            logs( "Temperature = {}k\nCutoff energy = {}meV\nCutoff Time = {}ps\nAlpha = {}\n<B> = {}\nFirst Markov approximation used? (rho(t) = rho(t-tau)) - {}\nTransformation approximation used: {} - {}\n\n", p_phonon_T, Hz_to_eV( p_phonon_wcutoff ) * 1E3, p_phonon_tcutoff * 1E12, p_phonon_alpha, p_phonon_b, ( numerics_phonon_approximation_1 == 1 ? "Yes" : "No" ), numerics_phonon_approximation_2, approximations.at( numerics_phonon_approximation_2 ) );
        } else {
            logs( "Not using phonons\n\n" );
        }

        logs.wrapInBar( "Numerical Parameters" );
        logs( "\n" );
        logs.wrapInBar( "Time", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Timeborder start: {:.8e} s - {:.2f} ps\n", t_start, t_start * 1E12 );
        logs( "Timeborder end: {:.8e} s - {:.2f} ps\n", t_end, t_end * 1E12 );
        logs( "Timeborder delta: {:.8e} s - {:.2f} fs \n", t_step, t_step * 1E15 );
        logs( "Ideal time delta for this calculation: {:.8e} s - {:.2f} fs, minimum possible: {:.8e} s - {:.2f} fs\n", getIdealTimestep(), getIdealTimestep() * 1E15, std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::epsilon() * 1E15 );
        int works = 1;
        if ( ( init_rabifrequenz != 0.0 ) && ( 3. * t_step > 2. * M_PI / init_rabifrequenz ) )
            works = 0;
        else if ( max_rabifrequenz != 0.0 && 3. * t_step > 2. * M_PI / max_rabifrequenz )
            works = 0;
        if ( !works ) {
            fmt::print( "{} WARNING: Step may be too small to resolve predicted oscillation: dT needed vs dT: {:.10e} < {:.10e}\n", PREFIX_WARNING, 2. / 3. * M_PI / std::max( init_rabifrequenz, max_rabifrequenz ), t_step );
            logs( "WARNING: Step may be too small to resolve predicted oscillation: \n- delta T needed: {:.10e} \n- delta T used: {:.10e}\n", 2. / 3. * M_PI / std::max( init_rabifrequenz, max_rabifrequenz ), t_step );
        }
        logs( "Time iterations (main loop) = {}\n\n", iterations_t_max );

        logs.wrapInBar( "Settings", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        logs( "Order of Runge-Kutta used: Time: RK{}, Spectrum: RK{}\n", numerics_order_t, numerics_order_tau );
        logs( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
        logs( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture == 1 ) ? "YES" : "NO" ) );
        logs( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) ? "Analytic" : "Matrix Exponential" ) );
        logs( "Threads used for primary calculations - {}\nThreads used for Secondary calculations - {}\n", numerics_phonons_maximum_threads, numerics_maximum_threads );
        logs( "Used scaling for parameters? - {}\n", ( scale_parameters ? std::to_string( scale_value ) : "no" ) );
        if ( p_phonon_T )
            logs( "Cache Phonon Coefficient Matrices? - {}\n", ( numerics_use_saved_coefficients ? fmt::format( "Yes (maximum {} matrices saved)", ( numerics_saved_coefficients_cutoff > 0 ) ? numerics_saved_coefficients_cutoff : numerics_saved_coefficients_max_size ) : "No" ) );
        logs( "\n\n" );
        if ( logfilecounter >= 0 )
            logs( "Logfile ident number: {}\n\n", logfilecounter );
        logs.wrapInBar( "Program Log:", LOG_SIZE_FULL, LOG_LEVEL_2 );
        logs( "\n" );
        logs.level2( "OutputHandlerStrings: {}\n", output_handlerstrings );
    }

    static void help() {
        fmt::print( "--help, -help, -h\tThis screen\n" );
        fmt::print( "--time [start] [end] [step]\n\t--tstart [start]\n\t--tend [end]\n\t--tstep [step]\n" );
        fmt::print( "--system [p_omega_atomic] [p_omega_cavity_H] [p_omega_cavity_V] [p_omega_coupling] [kappa] [gammapure] [gamma] [deltaE] [bexcitonbinding] else standard values (1.366eV,1.366eV,1.366eV,66mueV,66mueV,3mueV,1mueV,0,3meV) are used\n\t--we [1.366eV] First (H-polarized) exciton energy\n\t--wcH [1.366eV] Cavity Energy H\n\t--wcV [1.366eV] Cavity Energy V\n\t--coupling [66mueV] QD-Lightfield coupling. Same for H and V.\n\t--kappa [66mueV] Cavity Decay\n\t--gamma [1mueV] Radiative Decay\n\t--gammapure [3mueV] Pure Dephasing\n\t--deltaE [0eV] H and V energy difference\n\t--excitonBindEnergy [3meV] Biexciton binding energy\n" );
        fmt::print( "--chirp ['[Array Time]'] ['[Array Y]'] ['[Array d/dt]'] [type]\n\t--chirpT ['[Array Time]']\n\t--chirpY ['[Array Y]']\n\t--chirpDDT ['[Array d/dt]']\n\t--chirpType [type] where type = monotone, hermite, linear, spline\n" );
        fmt::print( "--pulse [Center] [Amplitude] [Frequency] [Sigma] [Type]\n\t-pulse for standard pulse\n\t--pulseCenter [Center]\n\t--pulseAmp [Amplitude]\n\t--pulseFreq [Frequency]\n\t--pulseSigma [Sigma]\n\t--pulseType [Type] where Type = cw, gauss, gauss_pi\n" );
        fmt::print( "--dimensions [maximum Photons] [Initial state]\n\t--maxPhotons [maximum Photons]\n\t--initState [Initial state], has to be smaller than (2*n+1)\n" );
        fmt::print( "--spectrum [Tau Resolution] [Center] [Range] [Omega Resolution] enables spectrum\n\t-spectrum enables spectrum centered at cavity\n\t--specTauRes [Grid resolution (int)] standard is 1000\n\t--specCenter [Center]\n\t--specRange [Range]\n\t--specWRes [w vector resolution (int)] standard is 1000\n" );
        fmt::print( "--phonons [Alpha] [W Cutoff] [t Cutoff ] [Temperature] Enables phonons with custom settings\n\t--temperature [temperature] Enables phonons with standard values and set temperature\n\t-phonons Enables phonons with standard values at T=3k\n\t--phononorder Sets the order of approximation to use\n\t-noMarkov disables first markov approximation\n\t-phononcoeffs Enables output of additional phonon coefficients\n" );
        fmt::print( "-g2 enables calculation of full G2 statistics\n\t-g2s enables calculation of simple G2 statistics\n-RK5 enables Runge Kutta of order 5 for T and Tau direction\n\t-RK5T enables Runge Kutta of order 5 for T direction\n\t-RK5Tau enables Runge Kutta of order 5 for Tau direction\n" );
        fmt::print( "-noInteractionpic disables Interaction picture - enabled by default\n-noRWA disables rotating wave approximation - enabled by default\n-timeTrafoMatrixExponential enables Time Transformation via Matrix exponential - disabled by default\n-startCoherent enables starting with a coherent state. Starting state is then ground state with alpha = initState\n-fullDM enables full output of densitymatrix including all offdiagonal terms.\n" );
        fmt::print( "--Threads [number] number of threads to use for both AKF and Spectrum integral calculation\n" );
        fmt::print( "Additional commands:\n\t-advLog Enables advanced logging\n\t-noHandler disables handler strings and enables loadbar output (for console)\n\t-output_operators, -outputHamiltons, -outputOperatorsStop Enables output of matrices (requires -advLog)\n\t-noRaman disables slow Raman calculations" );
    }
};