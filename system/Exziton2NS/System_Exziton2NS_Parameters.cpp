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
    double p_phonon_b, p_phonon_alpha, p_phonon_wcutoff, p_phonon_T, p_phonon_tcutoff;

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
    int numerics_maximum_threads, spectrum_frequency_iterations, numerics_phonon_approximation_1, numerics_phonon_approximation_2; //akf_everyXIt //akf_skip_omega
    double akf_deltaWmax, spectrum_frequency_center, spectrum_frequency_range;

    Parameters(){};
    Parameters( const std::vector<std::string> &arguments ) : Parameters_Parent() {
        init( arguments );
    }

    bool parseInput( const std::vector<std::string> &arguments ) {
        Parse_Parameters params;
        // Look for --time, if not found, standard values are used (t0 = 0, t1 = 1ns, deltaT = auto)
        params = Parse_Parameters( arguments, {"--time", "--tstart", "--tend", "--tstep"}, {3, 1, 1, 1}, "Timeparameters" );
        t_start = params.get<double>( {0, 3}, "0.0" );
        t_end = params.get<double>( {1, 4}, "1.0ns" );
        t_step = params.get<double>( {2, 5}, "-1" );

        // Look for --system, if not found, standard system is used (g=66mueV, k=66mueV, p_omega_pure_dephasing = 3mueV, p_omega_decay = 1mueV)
        params = Parse_Parameters( arguments, {"--system", "--we", "--wc", "--coupling", "--kappa", "--gammapure", "--gamma"}, {6, 1, 1, 1, 1, 1, 1}, "Systemparameters" );
        p_omega_atomic = params.get<double>( {0, 6}, "1.52eV" );
        p_omega_cavity = params.get<double>( {1, 7}, "1.52eV" );
        p_omega_coupling = params.get<double>( {2, 8}, "66mueV" );
        p_omega_cavity_loss = params.get<double>( {3, 9}, "66mueV" );
        p_omega_pure_dephasing = params.get<double>( {4, 10}, "3mueV" );
        p_omega_decay = params.get<double>( {5, 11}, "1mueV" );

        // Look for --chirp, if not found, standard system is used (no chirp, everything zero)
        params = Parse_Parameters( arguments, {"--chirp", "--chirpT", "--chirpY", "--chirpDDT", "--chirpType"}, {4, 1, 1, 1, 1}, "Chirpparameters" );
        chirp_t = convertParam<double>( str_to_vec( params.get( {0, 4}, "[0,1E-9]" ) ) );
        chirp_y = convertParam<double>( str_to_vec( params.get( {1, 5}, "[0,0]" ) ) );
        chirp_ddt = convertParam<double>( str_to_vec( params.get( {2, 6}, "[0,0]" ) ) );
        chirp_type = params.get( {3, 7}, "monotone" );

        // Look for --pulse, if not found, standard system is used (no pulse, everything zero)
        params = Parse_Parameters( arguments, {"--pulse", "--pulseCenter", "--pulseAmp", "--pulseFreq", "--pulseSigma", "--pulseType", "-pulse"}, {5, 1, 1, 1, 1, 1, 1}, "Pulseparameters" );
        if ( params.get( 0, "Empty" ).at( 0 ) == '[' ) {
            pulse_center = convertParam<double>( str_to_vec( params.get( {0, 5}, "[0]" ) ) );
            pulse_amp = convertParam<double>( str_to_vec( params.get( {1, 6}, "[0]" ) ) );
            pulse_omega = convertParam<double>( str_to_vec( params.get( {2, 7}, "[0]" ) ) );
            pulse_sigma = convertParam<double>( str_to_vec( params.get( {3, 8}, "[1]" ) ) );
            pulse_type = str_to_vec( params.get( {4, 9}, "[cw]" ) );
        } else if ( params.get( 10 ) || params.get( 0 ) ) {
            pulse_center.emplace_back( params.get<double>( {0, 5}, "100ps" ) );
            pulse_amp.emplace_back( params.get<double>( {1, 6}, "2" ) );
            pulse_omega.emplace_back( params.get<double>( {2, 7}, "1.52eV" ) );
            pulse_sigma.emplace_back( params.get<double>( {3, 8}, "20ps" ) );
            pulse_type.emplace_back( params.get( {4, 9}, "gauss_pi" ) );
        } else {
            pulse_center.emplace_back( convertParam<double>( "0" ) );
            pulse_amp.emplace_back( convertParam<double>( "0" ) );
            pulse_omega.emplace_back( convertParam<double>( "0" ) );
            pulse_sigma.emplace_back( convertParam<double>( "1" ) );
            pulse_type.emplace_back( "cw" );
        }

        // Look for --dimensions, if not found, standard system is used (maxphotons = 0, starting state = |g,0>)
        params = Parse_Parameters( arguments, {"--dimensions", "--maxPhotons", "--initState"}, {2, 1, 1}, "Initial State parameters" );
        p_max_photon_number = params.get<int>( {0, 2}, "1" );
        p_initial_state = params.get<int>( {1, 3}, "0" );

        // Look for --spectrum, if not found, no spectrum is evaluated
        params = Parse_Parameters( arguments, {"--spectrum", "--specTauRes", "--specCenter", "--specRange", "--specWRes", "-spectrum"}, {4, 1, 1, 1, 1, 1}, "Spectrum Parameters" );
        numerics_calculate_spectrum = params.get( 0 ) || params.get( 8 ) ? 1 : 0;
        iterations_tau_resolution = params.get<int>( {0, 4}, "1" );
        spectrum_frequency_center = params.get<double>( {1, 5}, std::to_string( p_omega_cavity ) );
        spectrum_frequency_range = params.get<double>( {2, 6}, std::to_string( ( p_omega_coupling + p_omega_cavity_loss + p_omega_decay + p_omega_pure_dephasing ) * 10.0 ) );
        iterations_w_resolution = params.get<int>( {3, 7}, "1" );

        // Look for (-RK4), -RK5, (-RK4T), (-RK4Tau), -RK5T, -RK5Tau
        params = Parse_Parameters( arguments, {"-g2", "-RK5", "-RK5T", "-RK5Tau", "-noInteractionpic", "-noRWA", "--Threads", "-noHandler", "-outputOperators", "-outputHamiltons", "-outputOperatorsStop", "-timeTrafoMatrixExponential", "-startCoherent", "-fullDM", "-scale"}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Other parameters" );
        numerics_calculate_g2 = params.get( 0 );
        numerics_order_t = ( params.get( 1 ) || params.get( 2 ) ) ? 5 : 4;
        numerics_order_tau = ( params.get( 1 ) || params.get( 3 ) ) ? 5 : 4;
        numerics_order_highest = numerics_order_t;
        if ( numerics_order_tau > numerics_order_highest )
            numerics_order_highest = numerics_order_tau;
        numerics_use_interactionpicture = params.get( 4 ) ? 0 : 1;
        numerics_use_rwa = params.get( 5 ) ? 0 : 1;
        numerics_maximum_threads = params.get<int>( 6, "1" );
        output_handlerstrings = params.get( 7 ) ? 0 : 1;
        output_operators = params.get( 8 ) ? 2 : ( params.get( 9 ) ? 1 : ( params.get( 10 ) ? 3 : 0 ) );
        numerics_order_timetrafo = params.get( 11 ) ? TIMETRANSFORMATION_MATRIXEXPONENTIAL : TIMETRANSFORMATION_PRECALCULATED;
        startCoherent = params.get( 12 );
        output_full_dm = params.get( 13 );
        scale_parameters = params.get( 14 );

        // Phonon Parameters
        params = Parse_Parameters( arguments, {"--phonons", "--temperature", "-phonons", "--phononorder", "-noMarkov", "-phononcoeffs"}, {5, 1, 1, 1, 1, 1} );
        p_phonon_alpha = params.get<double>( 0, "0.03E-24" );
        p_phonon_wcutoff = params.get<double>( 1, "1meV" );
        p_phonon_tcutoff = params.get<double>( 2, "4ps" );
        p_phonon_T = params.get( 6 ) ? 3.0 : params.get<double>( {3, 5}, "0" );
        numerics_phonon_approximation_2 = params.get<int>( {4, 7}, std::to_string( PHONON_APPROXIMATION_BACKWARDS_INTEGRAL ) ); // Second Markov / Transformation method
        numerics_phonon_approximation_1 = params.get( 8 ) ? 0 : 1;                                                              // First Markov
        output_coefficients = params.get( 9 ) ? 1 : 0;

        subfolder = arguments.back();
        return true;
    }

    bool scaleInputs( const double scaling ) {
        // Adjust normal parameters: time is multiplid by scaling, frequency divided
        t_start = scaleVariable( t_start, scaling );
        t_end = scaleVariable( t_end, scaling );
        t_step = scaleVariable( t_step, scaling );
        p_omega_atomic = scaleVariable( p_omega_atomic, 1.0 / scaling );
        p_omega_cavity = scaleVariable( p_omega_cavity, 1.0 / scaling );
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

    bool adjustInput() {
        // Calculate/Recalculate some parameters:
        // Adjust pulse area if pulse_type is "gauss_pi"
        for ( int i = 0; i < (int)pulse_amp.size(); i++ )
            if ( pulse_type.at( i ).compare( "gauss_pi" ) == 0 ) {
                pulse_amp.at( i ) = pulse_amp.at( i ) * M_PI / ( std::sqrt( 2.0 * M_PI ) * pulse_sigma.at( i ) );
                pulse_type.at( i ) = "gauss";
            }

        // Calculate Rabi frequencies:
        init_rabifrequenz = rabiFrequency( p_omega_atomic - p_omega_cavity, p_omega_coupling, (int)( p_initial_state / 2.0 ) );
        max_rabifrequenz = rabiFrequency( p_omega_atomic - p_omega_cavity, p_omega_coupling, (int)( p_initial_state / 2.0 ) + 1 );

        // Calculate minimum step necessary to resolve Rabi-oscillation if step=-1
        if ( t_step == -1 ) {
            if ( numerics_use_rwa )
                t_step = 2. / 8. * M_PI / std::max( std::max( init_rabifrequenz, max_rabifrequenz ), p_omega_coupling + p_omega_cavity_loss + p_omega_decay + p_omega_pure_dephasing + ( vec_max( pulse_amp ) > 0 ? std::abs( p_omega_atomic - vec_max( pulse_omega ) ) : 0 ) );
            if ( !numerics_use_rwa )
                t_step = 2. / 8. * M_PI / std::max( std::max( init_rabifrequenz, max_rabifrequenz ), p_omega_atomic + p_omega_cavity );
        }

        // Calculate the maximum dimensions for operator matrices (max states)
        maxStates = 2 * ( p_max_photon_number + 1 );

        // Calculate stuff for RK
        iterations_t_max = (int)std::ceil( ( t_end - t_start ) / t_step );

        //FIXME Adjust tau resolution grid points to be at 1000 or the maximum t iterations in case that this number is smaller than 1000
        //iterations_tau_resolution = std::max( iterations_t_max, iterations_tau_resolution );

        // Mandatory: rescale chirp ddt into chirp/ps
        for ( long unsigned i = 0; i < chirp_ddt.size(); i++ )
            chirp_ddt.at( i ) = chirp_ddt.at( i ) * 1E12;
        // Calculate values for chirp:
        chirp_total = vec_max( chirp_y );

        // Initial and maximum system detuning (not taking into account custom chirps)
        init_detuning = ( p_omega_atomic - p_omega_cavity );
        max_detuning = ( init_detuning + chirp_total > init_detuning ) ? init_detuning + chirp_total : init_detuning;

        // Adjust/calculate frequency range for spectrum
        spectrum_frequency_iterations = iterations_t_max / iterations_w_resolution;
        if ( spectrum_frequency_center == -1 )
            spectrum_frequency_center = p_omega_cavity;
        if ( spectrum_frequency_range == -1 )
            spectrum_frequency_range = ( std::abs( max_detuning ) + p_omega_coupling + p_omega_cavity_loss / 2. ) * 3.;

        // Calculate total number of iterations necessary
        iterations_total_max = iterations_t_max;
        if ( numerics_calculate_spectrum ) {
            int curIt = 1;
            for ( double ii = t_start + t_step; ii < t_end; ii += t_step ) {
                if ( curIt % iterations_tau_resolution == 0 ) {
                    for ( double iii = ii + t_step; iii < t_end; iii += t_step ) {
                        iterations_total_max++;
                    }
                    curIt = 1;
                } else
                    curIt += 1;
            }
        }

        // Calculate phonon stuff
        p_phonon_b = 1.0;
        if ( p_phonon_T > 0 ) {
            double integral = 0;
            double stepsize = 0.01 * p_phonon_wcutoff;
            for ( double w = stepsize; w < 10 * p_phonon_wcutoff; w += stepsize ) {
                integral += stepsize * ( p_phonon_alpha * w * std::exp( -w * w / 2.0 / p_phonon_wcutoff / p_phonon_wcutoff ) / std::tanh( 1.0545718E-34 * w / 2.0 / 1.3806488E-23 / p_phonon_T ) );
            }
            p_phonon_b = std::exp( -0.5 * integral );
            //p_omega_pure_dephasing = convertParam<double>( "1mueV" ) * p_phonon_T;
            //p_omega_decay *= p_phonon_b*p_phonon_b;
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
                logs( "Exiting system at t_0 = {}\nAmplitude {} ({}meV)\nFrequency {} ({}eV)\nFWHM {}\n", pulse_center.at( i ), pulse_amp.at( i ), Hz_to_eV( pulse_amp.at( i ) ) * 1E3, pulse_omega.at( i ), Hz_to_eV( pulse_omega.at( i ) ), pulse_sigma.at( i ) * ( 2 * std::sqrt( 2 * std::log( 2 ) ) ) );
                logs( "Used pulse_type - " + pulse_type.at( i ) + "\n" );
            }
            logs( "\n" );
        } else
            logs( "Not using pulse to exite system\n\n" );
        logs.wrapInBar( "Energy Chirp", LOG_SIZE_HALF, LOG_LEVEL_1, LOG_BAR_1 );
        if ( chirp_total != 0 ) {
            for ( int i = 0; i < (int)chirp_t.size() - 1; i++ ) {
                logs( "Chirp between t0 = {}ps\nt1 = {}ps\nTotal Chirp: {}mueV\n-> average rate {}mueV/ps\n", chirp_t.at( i ), chirp_t.at( i + 1 ), Hz_to_eV( chirp_y.at( i + 1 ) - chirp_y.at( i ) ) * 1E6, Hz_to_eV( ( chirp_y.at( i + 1 ) - chirp_y.at( i ) ) ) * 1E6 / 1E12 / ( chirp_t.at( i + 1 ) - chirp_t.at( i ) ) );
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
        logs( "Maximum tau-grid resolution is {}x{}\nMaximum w-vector resolution is {}\n\n", iterations_tau_resolution, iterations_tau_resolution, iterations_w_resolution );

        logs.wrapInBar( "Phonons" );
        std::vector<std::string> approximations = {"Transformation integral via d/dt chi = -i/hbar*[H,chi] + d*chi/dt onto interaction picture chi(t-tau)", "Transformation Matrix U(t,tau)=exp(-i/hbar*H_DQ_L(t)*tau) onto interaction picture chi(t-tau)", "No Transformation, only interaction picture chi(t-tau)", "Analytical Lindblad formalism"};
        logs( "\nTemperature = {}k\nCutoff energy = {}meV\nCutoff Time = {}ps\nAlpha = {}\n<B> = {}\nFirst Markov approximation used? (rho(t) = rho(t-tau)) - {}\nTransformation approximation used: {} - {}\n\n", p_phonon_T, Hz_to_eV( p_phonon_wcutoff ) * 1E3, p_phonon_tcutoff * 1E12, p_phonon_alpha, p_phonon_b, ( numerics_phonon_approximation_1 == 1 ? "Yes" : "No" ), numerics_phonon_approximation_2, approximations.at( numerics_phonon_approximation_2 ) );

        logs.wrapInBar( "Numerics" );
        logs( "\nOrder of Runge-Kutta used: Time: RK{}, Spectrum: RK{}\n", numerics_order_t, numerics_order_tau );
        logs( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
        logs( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture == 1 ) ? "YES" : "NO" ) );
        logs( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == TIMETRANSFORMATION_PRECALCULATED ) ? "Analytic" : "Matrix Exponential" ) );
        logs( "Threads used for AFK? - {}\n", numerics_maximum_threads );
        logs( "Used scaling for parameters? - {}\n", ( scale_parameters ? std::to_string( scale_value ) : "no" ) );
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
        fmt::print( "--spectrum [Tau Resolution] [Center] [Range] [Omega Resolution] enables spectrum\n\t-spectrum enables spectrum centered at cavity\n\t--specTauRes [Grid resolution (int)] standard is 1000\n\t--specCenter [Center]\n\t--specRange [Range]\n\t--specWRes [w vector resolution (int)] standard is 1000\n" );
        fmt::print( "--phonons [Alpha] [W Cutoff] [t Cutoff ] [Temperature] Enables phonons with custom settings\n\t--temperature [temperature] Enables phonons with standard values and set temperature\n\t-phonons Enables phonons with standard values at T=3k\n\t-phononorder Sets the order of approximation to use\n\t-noMarkov disables first markov approximation\n\t-phononcoeffs Enables output of additional phonon coefficients\n" );
        fmt::print( "-g2 enables calculation of G2(tau=0)\n-RK5 enables Runge Kutta of order 5 for T and Tau direction\n\t-RK5T enables Runge Kutta of order 5 for T direction\n\t-RK5Tau enables Runge Kutta of order 5 for Tau direction\n" );
        fmt::print( "-noInteractionpic disables Interaction picture - enabled by default\n-noRWA disables rotating wave approximation - enabled by default\n-timeTrafoMatrixExponential enables Time Transformation via Matrix exponential - disabled by default\n-startCoherent enables starting with a coherent state. Starting state is then ground state with alpha = initState\n-fullDM enables full output of densitymatrix including all offdiagonal terms.\n" );
        fmt::print( "--Threads [number] number of threads to use for both AKF and Spectrum integral calculation\n" );
        fmt::print( "Additional commands:\n\t-advLog Enables advanced logging\n\t-noHandler disables handler strings and enables loadbar output (for console)\n\t-output_operators, -outputHamiltons, -outputOperatorsStop Enables output of matrices (requires -advLog)" );
    }
};