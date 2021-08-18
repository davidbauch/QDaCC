#include "system/parameters.h"

namespace std {
double max( const Parameter &p1, const Parameter &p2 ) {
    return max( p2.get(), p1.get() );
}
double min( const Parameter &p1, const Parameter &p2 ) {
    return min( p2.get(), p1.get() );
}
}; // namespace std

Parameters::Parameters( const std::vector<std::string> &arguments ) {
    numerics_use_interactionpicture = 0;
    numerics_use_rwa = 0;
    numerics_order_timetrafo = TIMETRANSFORMATION_MATRIXEXPONENTIAL;
    output_advanced_log = 0;
    output_handlerstrings = 0;
    output_operators = 0;
    iterations_t_skip = 1;
    iterations_tau_resolution = 1;
    iterations_w_resolution = 1;
    scale_parameters = false;
    scale_value = 1E12;
    iterations_t_max = 1;
    output_full_dm = false;
    output_no_dm = false;
    maxStates = 0;
    numerics_calculate_till_converged = false;

    // Parsing input:
    Timer &timer_parseInput = Timers::create( "Parsing parameters", true, false );
    Log::wrapInBar( "Conversion of input variables", Log::BAR_SIZE_FULL, Log::LEVEL_2, Log::BAR_0 );
    Log::L2( "\n" );
    Log::L2( "Parsing input variables... " );
    timer_parseInput.start();
    if ( !parseInput( arguments ) ) {
        Log::L2( "Parsing input variables failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    timer_parseInput.end();
    Log::L2( "successful. Elapsed time is {}ms\n", timer_parseInput.getWallTime( TIMER_MILLISECONDS ) );

    // Scaling inputs:
    if ( scale_parameters ) {
        Log::L2( "Rescaling parameters to {}... ", scale_value );
        scaleInputs( scale_value );
        Log::L2( "Done!\n" );
    }

    // Adjusting inputs:
    Timer &timer_adjustInput = Timers::create( "Adjusting parameters", true, false );
    Log::L2( "Adjusting input variables... " );
    timer_adjustInput.start();
    if ( !adjustInput() ) {
        Log::L2( "Adjusting input variables failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    Log::L2( "successful. Elapsed time is {}ms\n", timer_adjustInput.getWallTime( TIMER_MILLISECONDS ) );
    timer_adjustInput.end();
}

bool Parameters::parseInput( const std::vector<std::string> &arguments ) {
    //Parse_Parameters params;
    // Look for --time, if not found, standard values are used (t0 = 0, t1 = 1ns, deltaT = auto)
    t_start = get_parameter<double>( "--time", "tstart" );
    t_end = get_parameter<double>( "--time", "tend" );
    t_step = get_parameter<double>( "--time", "tstep" );

    // Runge Kutta Parameters
    numerics_rk_order = get_parameter<double>( "--rk", "rkorder" );
    numerics_rk_tol = get_parameter<double>( "--rk", "rktol" );
    numerics_rk_stepdelta = get_parameter<double>( "--rk", "rkstepdelta" );
    numerics_rk_stepmin = get_parameter<double>( "--rk", "rkstepmin" );
    numerics_rk_stepmax = get_parameter<double>( "--rk", "rkstepmax" );
    numerics_rk_interpolate = !get_parameter_passed( "-rknointerpolate" );
    numerics_rk_usediscrete_timesteps = numerics_rk_stepdelta > 0 ? true : false;

    inputstring_electronic = get_parameter( "--S", "SE" );
    inputstring_photonic = get_parameter( "--S", "SO" );
    inputstring_pulse = get_parameter( "--S", "SP" );
    inputstring_chirp = get_parameter( "--S", "SC" );
    inputstring_spectrum = get_parameter( "--G", "GS" );
    inputstring_indist = get_parameter( "--G", "GI" );
    inputstring_conc = get_parameter( "--G", "GC" );
    inputstring_gfunc = get_parameter( "--G", "GF" );
    inputstring_wigner = get_parameter( "--G", "GW" );

    p_omega_coupling = get_parameter<double>( "--system", "coupling" );
    p_omega_cavity_loss = get_parameter<double>( "--system", "kappa" );
    p_omega_pure_dephasing = get_parameter<double>( "--system", "gammapure" );
    p_omega_decay = get_parameter<double>( "--system", "gamma" );
    p_initial_state_s = get_parameter( "--R" );

    // Look for --spectrum, if not found, no spectrum is evaluated
    iterations_tau_resolution = get_parameter<int>( "--spectrum", "gridres" );
    numerics_use_interactionpicture = get_parameter_passed( "-noInteractionpic" ) ? 0 : 1;
    numerics_use_rwa = get_parameter_passed( "-noRWA" ) ? 0 : 1;
    numerics_maximum_threads = get_parameter<int>( "--Threads" );
    if ( numerics_maximum_threads == -1 )
        numerics_maximum_threads = omp_get_max_threads();
    output_handlerstrings = get_parameter_passed( "-noHandler" ) ? 0 : 1;
    output_operators = get_parameter_passed( "-outputOp" ) ? 2 : ( get_parameter_passed( "-outputHamiltons" ) ? 1 : ( get_parameter_passed( "-outputOpStop" ) ? 3 : 0 ) );
    numerics_order_timetrafo = get_parameter_passed( "-timeTrafoMatrixExponential" ) ? TIMETRANSFORMATION_MATRIXEXPONENTIAL : TIMETRANSFORMATION_ANALYTICAL;
    startCoherent = false; //TODO://get_parameter_passed( "-startCoherent" ) || ( p_initial_state_electronic.front() == 'c' );
    output_full_dm = get_parameter_passed( "-fullDM" );
    output_no_dm = get_parameter_passed( "-noDM" );
    scale_parameters = get_parameter_passed( "-scale" );
    numerics_use_saved_coefficients = !get_parameter_passed( "-disableMatrixCaching" );
    numerics_use_saved_hamiltons = !get_parameter_passed( "-disableHamiltonCaching" );
    numerics_phonons_maximum_threads = ( !numerics_use_saved_coefficients || !get_parameter_passed( "-disableMainProgramThreading" ) ) ? numerics_maximum_threads : 1;
    numerics_output_raman_population = get_parameter_passed( "-raman" );
    logfilecounter = convertParam<int>( splitline( get_parameter( "--lfc" ), ',' ) );
    numerics_calculate_timeresolution_indistinguishability = get_parameter_passed( "-timedepInd" );
    numerics_output_electronic_emission = get_parameter_passed( "-oElec" );
    numerics_stretch_correlation_grid = false; //FIXME: Doesnt work right now
    numerics_interpolate_outputs = get_parameter_passed( "-interpolate" );

    // Phonon Parameters
    p_phonon_alpha = get_parameter<double>( "--phonons", "phononalpha" );
    p_phonon_wcutoff = get_parameter<double>( "--phonons", "phononwcutoff" );
    p_phonon_tcutoff = get_parameter<double>( "--phonons", "phonontcutoff" );
    p_phonon_T = get_parameter<double>( "--phonons", "temperature" );
    numerics_phonon_approximation_order = get_parameter<int>( "--phonons", "phononorder" );
    numerics_phonon_approximation_markov1 = get_parameter_passed( "-noMarkov" ) ? 0 : 1; // First Markov
    output_coefficients = get_parameter_passed( "-phononcoeffs" ) ? 1 : 0;
    p_phonon_adjust = !get_parameter_passed( "-noPhononAdjust" );
    p_phonon_pure_dephasing = convertParam<double>( "1mueV" );
    // Path Integral Parameters
    p_phonon_nc = get_parameter<int>( "--pathintegral", "NC" );
    numerics_pathintegral_stepsize_iterator = get_parameter<double>( "--pathintegral", "iteratorStepsize" );
    numerics_pathintegral_squared_threshold = get_parameter<double>( "--pathintegral", "squaredThreshold" );
    numerics_pathintegral_sparse_prune_threshold = get_parameter<double>( "--pathintegral", "sparsePruneThreshold" );
    numerics_pathintegral_dynamiccutoff_iterations_max = get_parameter<double>( "--pathintegral", "iteratorStepsize" );
    numerics_pathintegral_docutoff_propagator = get_parameter_passed( "-cutoffPropagator" );
    numerics_pathint_partially_summed = true;
    t_step_pathint = get_parameter<double>( "--pathintegral", "tstepPath" );

    kb = 1.3806488E-23;   // J/K, scaling needs to be for energy
    hbar = 1.0545718E-34; // J/s, scaling will be 1

    parse_system();

    subfolder = arguments.back();
    return true;
}

//TODO: this is broken. and unncesesary actually. remove.
bool Parameters::scaleInputs( const double scaling ) {
    // Adjust normal parameters: time is multiplid by scaling, frequency divided
    t_start.setScale( scaling, Parameter::SCALE_TIME );
    t_end.setScale( scaling, Parameter::SCALE_TIME );
    t_step.setScale( scaling, Parameter::SCALE_TIME );
    numerics_rk_stepdelta.setScale( scaling, Parameter::SCALE_TIME );

    //TODO: inputs scalen! pulse un chirp auch!
    //p_omega_atomic_G_H.setScale( scaling, Parameter::SCALE_ENERGY );
    //p_omega_atomic_G_V.setScale( scaling, Parameter::SCALE_ENERGY );
    //p_omega_atomic_H_B.setScale( scaling, Parameter::SCALE_ENERGY );
    //p_omega_atomic_V_B.setScale( scaling, Parameter::SCALE_ENERGY );
    //p_deltaE.setScale( scaling, Parameter::SCALE_ENERGY );
    //p_biexciton_bindingenergy.setScale( scaling, Parameter::SCALE_ENERGY );
    //p_omega_cavity_H.setScale( scaling, Parameter::SCALE_ENERGY );
    //p_omega_cavity_V.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_coupling.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_cavity_loss.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_pure_dephasing.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_decay.setScale( scaling, Parameter::SCALE_ENERGY );
    // Adjust chirp and pulse
    //for ( int i = 0; i < (int)chirp_t.size(); i++ ) {
    //    chirp_t.at( i ).setScale( scaling, Parameter::SCALE_TIME );
    //    chirp_y.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
    //    chirp_ddt.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
    //}
    //for ( int i = 0; i < (int)pulse_center.size(); i++ ) {
    //    pulse_center.at( i ).setScale( scaling, Parameter::SCALE_TIME );
    //    if ( pulse_type.at( i ).compare( "gauss_pi" ) != 0 )
    //        pulse_amp.at( i ) = scaleVariable( pulse_amp.at( i ), 1.0 / scaling );
    //    pulse_omega.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
    //    pulse_sigma.at( i ).setScale( scaling, Parameter::SCALE_TIME );
    //    pulse_omega_chirp.at( i ).setScale( scaling * scaling, Parameter::SCALE_ENERGY );
    //}
    // Adjusting spectrum
    //TODO: neuen input saclen!
    //spectrum_frequency_center.setScale( scaling, Parameter::SCALE_ENERGY );
    //spectrum_frequency_range.setScale( scaling, Parameter::SCALE_ENERGY );
    // Phonons
    p_phonon_wcutoff.setScale( scaling, Parameter::SCALE_ENERGY );
    p_phonon_tcutoff.setScale( scaling, Parameter::SCALE_TIME );
    p_phonon_alpha.setScale( scaling * scaling, Parameter::SCALE_ENERGY );
    p_phonon_pure_dephasing.setScale( scaling, Parameter::SCALE_ENERGY );
    kb.setScale( scaling, Parameter::SCALE_ENERGY );
    return true;
}

bool Parameters::adjustInput() {
    Log::L2( "Adjusting Inputs...\n" );

    // For threadsafety
    if ( numerics_rk_order > 5 )
        numerics_use_saved_hamiltons = false;

    // Calculate/Recalculate some parameters:
    // Adjust pulse area if pulse_type is "gauss_pi"
    for ( auto &[name, mat] : input_pulse ) {
        for ( int i = 0; i < mat.string_v["Type"].size(); i++ ) {
            auto pos = mat.string_v["Type"][i].find( "_pi" );
            if ( pos != std::string::npos ) {
                if ( mat.string_v["Type"][i].find( "gauss" ) != std::string::npos )
                    mat.numerical_v["Amplitude"][i] = mat.numerical_v["Amplitude"][i] * M_PI / ( std::sqrt( 2.0 * M_PI * mat.numerical_v["Width"][i] * std::sqrt( std::pow( mat.numerical_v["Chirp"][i] / mat.numerical_v["Width"][i], 2.0 ) + std::pow( mat.numerical_v["Width"][i], 2.0 ) ) ) ) / 2.0; //https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.241306
                else if ( mat.string_v["Type"][i].find( "cutoff" ) != std::string::npos )
                    mat.numerical_v["Amplitude"][i] = mat.numerical_v["Amplitude"][i] * M_PI / ( std::sqrt( 2.0 * M_PI * mat.numerical_v["Width"][i] * mat.numerical_v["Width"][i] ) ) / 2.0; //https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.241306
                mat.string_v["Type"][i].erase( pos, 3 );
            }
        }
    }
    for ( auto &[name, mat] : input_chirp ) {
        if ( mat.string["Type"].compare( "sine" ) != 0 ) {
            for ( long unsigned i = 0; i < mat.numerical_v["ddt"].size(); i++ )
                mat.numerical_v["ddt"][i] = mat.numerical_v["ddt"][i] * 1E12;
        }
    }

    // Calculate minimum step necessary to resolve Rabi-oscillation if step=-1
    if ( t_step == -1 ) {
        t_step = 1E-13; //std::min( scaleVariable( 1E-13, scale_value ), getIdealTimestep() );
        //t_step = std::max( std::numeric_limits<double>::epsilon(), t_step );
    }
    if ( t_end == -1 ) {
        // If this is given, we calculate the t-direction until 99% ground state poulation is reached after any pulses.
        numerics_calculate_till_converged = true;
        for ( auto &[name, mat] : input_pulse ) {
            for ( auto &t : mat.numerical_v["Center"] ) {
                t_end = std::max( t_end.get(), 2.0 * t );
            }
        }
        if ( t_end == -1 )
            t_end = 10E-12;
        Log::L2( "Calculate till at least {} and adjust accordingly to guarantee convergence.\n", t_end );
    }

    // Calculate stuff for RK
    iterations_t_max = (int)std::ceil( ( t_end - t_start ) / ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? t_step_pathint : t_step ) );
    iterations_t_skip = std::max( 1.0, std::ceil( iterations_t_max / iterations_tau_resolution ) );

    // No phonon adjust if pathintegral is chosen
    if ( numerics_phonon_approximation_order == 5 ) {
        p_phonon_adjust = false;
    }
    // Calculate phonon stuff
    p_phonon_b = 1.0;
    if ( p_phonon_T >= 0 ) {
        double integral = 0;
        double stepsize = 0.01 * p_phonon_wcutoff;
        for ( double w = stepsize; w < 10 * p_phonon_wcutoff; w += stepsize ) {
            integral += stepsize * ( p_phonon_alpha * w * std::exp( -w * w / 2.0 / p_phonon_wcutoff / p_phonon_wcutoff ) / std::tanh( hbar * w / 2.0 / kb / p_phonon_T ) );
        }
        p_phonon_b = std::exp( -0.5 * integral );
        if ( p_phonon_adjust ) {
            p_omega_pure_dephasing = p_phonon_pure_dephasing * p_phonon_T;
            p_omega_decay = p_omega_decay * p_phonon_b * p_phonon_b;
        }
    }
    numerics_saved_coefficients_max_size = (int)( ( t_end - t_start ) / t_step * 2.0 * ( p_phonon_tcutoff / t_step ) ) + 10;
    trace.reserve( iterations_t_max + 5 );

    numerics_saved_coefficients_cutoff = 0; //( numerics_calculate_spectrum || numerics_calculate_g2 ) ? 0 : ( p_phonon_tcutoff / t_step ) * 5;
    Log::L2( "Adjusting Inputs Done!\n" );
    return true;
}

void Parameters::parse_system() {
    // Generate the input variables for the electronic system:
    auto levels = splitline( inputstring_electronic, ';' );
    for ( std::string &level : levels ) {
        auto conf = splitline( level, ':' );
        input_s conf_s;
        conf_s.numerical["Energy"] = convertParam<Parameter>( conf[1] );           // Energy
        conf_s.string_v["CoupledTo"] = splitline( conf[2], ',' );                  // Coupled to Levels
        conf_s.numerical["DecayScaling"] = convertParam<Parameter>( conf[3] );     // Decay Scaling, Per Mode
        conf_s.numerical["DephasingScaling"] = convertParam<Parameter>( conf[4] ); // Dephasing Scaling
        conf_s.numerical["PhononCoupling"] = convertParam<Parameter>( conf[5] );   // Phonon Coupling
        input_electronic[conf[0]] = conf_s;
    }
    auto cavities = splitline( inputstring_photonic, ';' );
    for ( std::string &cavity : cavities ) {
        auto conf = splitline( cavity, ':' );
        input_s conf_s;
        conf_s.numerical["Energy"] = convertParam<Parameter>( conf[1] );                              // Energy
        conf_s.numerical["MaxPhotons"] = convertParam<Parameter>( conf[2] );                          // Maximum Photons
        conf_s.string_v["CoupledTo"] = splitline( conf[3], ',' );                                     // Coupled to Transitions
        conf_s.numerical_v["CouplingScaling"] = convertParam<Parameter>( splitline( conf[4], ',' ) ); // Coupling Scaling, per transition INTO cavity
        //conf_s.numerical_v["BackCouplingScaling"] = convertParam<Parameter>( splitline( conf[5], ',' ) ); // BackCoupling Scaling, per transition from cavity back into the electronic system
        conf_s.numerical["DecayScaling"] = convertParam<Parameter>( conf[5] ); // Decay Scaling, for all transitions
        input_photonic[conf[0]] = conf_s;
    }

    auto pulses = splitline( inputstring_pulse, ';' );
    for ( std::string &pulse : pulses ) {
        auto conf = splitline( pulse, ':' );
        input_s conf_s;
        conf_s.string_v["CoupledTo"] = splitline( conf[1], ',' );                               // Coupled to Transitions
        conf_s.numerical_v["Amplitude"] = convertParam<Parameter>( splitline( conf[2], ',' ) ); // Pulse Amp
        conf_s.numerical_v["Frequency"] = convertParam<Parameter>( splitline( conf[3], ',' ) ); // Frequency
        conf_s.numerical_v["Width"] = convertParam<Parameter>( splitline( conf[4], ',' ) );     // Width
        conf_s.numerical_v["Center"] = convertParam<Parameter>( splitline( conf[5], ',' ) );    // Center
        conf_s.numerical_v["Chirp"] = convertParam<Parameter>( splitline( conf[6], ',' ) );     // Chirp
        conf_s.string_v["Type"] = splitline( conf[7], ',' );                                    // Type
        input_pulse[conf[0]] = conf_s;
    }
    auto chirps = splitline( inputstring_chirp, ';' );
    for ( std::string &chirp : chirps ) {
        auto conf = splitline( chirp, ':' );
        input_s conf_s;
        conf_s.string_v["CoupledTo"] = splitline( conf[1], ',' );                               // Coupled to Transitions
        conf_s.numerical_v["AmpFactor"] = convertParam<Parameter>( splitline( conf[2], ',' ) ); // Amplitude Scaling for coupled_to
        conf_s.numerical_v["Amplitude"] = convertParam<Parameter>( splitline( conf[3], ',' ) ); // Amplitudes
        conf_s.numerical_v["Times"] = convertParam<Parameter>( splitline( conf[4], ',' ) );     // "Times"
        conf_s.numerical_v["ddt"] = convertParam<Parameter>( splitline( conf[5], ',' ) );       // "d/dt"
        conf_s.string["Type"] = conf[6];                                                        // Type
        input_chirp[conf[0]] = conf_s;
    }
    for ( std::string &spectrum : splitline( inputstring_spectrum, ';' ) ) {
        auto conf = splitline( spectrum, ':' );
        input_s conf_s;
        conf_s.string_v["Modes"] = splitline( conf[0], ',' );                                // Modes to calculate Spectrum for. Single modes can again be split with "+", meaning a+b;a to calculate for a+b and a seperately
        conf_s.numerical_v["Center"] = convertParam<Parameter>( splitline( conf[1], ',' ) ); // Center
        conf_s.numerical_v["Range"] = convertParam<Parameter>( splitline( conf[2], ',' ) );  // Range
        conf_s.numerical_v["resW"] = convertParam<Parameter>( splitline( conf[3], ',' ) );   // Resolution for w
        input_correlation["Spectrum"] = conf_s;
    }
    for ( std::string &indist : splitline( inputstring_indist, ';' ) ) {
        auto conf = splitline( indist, ':' );
        input_s conf_s;
        conf_s.string_v["Modes"] = splitline( conf[0], ',' ); // Modes to calculate Indistinguishgability for. Single modes can again be split with "+", meaning a+b;a to calculate for a+b and a seperately
        input_correlation["Indist"] = conf_s;
    }
    for ( std::string &conc : splitline( inputstring_conc, ';' ) ) {
        auto conf = splitline( conc, ':' );
        input_s conf_s;
        conf_s.string_v["Modes"] = splitline( conf[0], ',' ); // Modes to calculate Concurrence for
        input_correlation["Conc"] = conf_s;
    }
    for ( std::string &g_func : splitline( inputstring_gfunc, ';' ) ) {
        auto conf = splitline( g_func, ':' );
        auto n = conf.size();
        input_s conf_s;
        conf_s.string_v["Modes"] = splitline( conf[0], ',' );                                                                                                               // Modes to calculate G1/G2 functions for
        conf_s.numerical_v["Order"] = convertParam<Parameter>( splitline( conf[1], ',' ) );                                                                                 // 1 or 2
        conf_s.numerical_v["Integrated"] = convertParam<Parameter>( n > 2 ? splitline( conf[2], ',' ) : std::vector<std::string>( conf_s.string_v["Modes"].size(), "2" ) ); // 0 or 1 or 2 for false/true/both
        input_correlation["GFunc"] = conf_s;
    }
    for ( std::string &wigner : splitline( inputstring_wigner, ';' ) ) {
        auto conf = splitline( wigner, ':' );
        auto n = conf.size();
        input_s conf_s;
        conf_s.string_v["Modes"] = splitline( conf[0], ',' );                                                                                                         // Modes to calculate Wigner function for
        conf_s.numerical_v["X"] = convertParam<Parameter>( splitline( conf[1], ',' ) );                                                                               // -X to X
        conf_s.numerical_v["Y"] = convertParam<Parameter>( n > 2 ? splitline( conf[2], ',' ) : splitline( conf[1], ',' ) );                                           // -Y to Y
        conf_s.numerical_v["Res"] = convertParam<Parameter>( n > 3 ? splitline( conf[3], ',' ) : std::vector<std::string>( conf_s.numerical_v["X"].size(), "100" ) ); // Resolution
        conf_s.numerical_v["Skip"] = convertParam<Parameter>( n > 4 ? splitline( conf[4], ',' ) : std::vector<std::string>( conf_s.numerical_v["X"].size(), "1" ) );  // Skips in t-direction
        input_correlation["Wigner"] = conf_s;
    }
}

void Parameters::log( const std::vector<std::string> &info ) {
    Log::L1( "                                                                                                                                                      \n                                                                                                                  #################                   \n                                                                                                              ..  :::::::::::::::::  ..               \n                                                                                                            ...       -      :+:      ...             \n                                                                                                           ...       .+:     :+:       ...            \n                                                                                                          ...        =+=     :+:        ...           \n                                                                                                         ....       :+++-    :+:         ...          \n         .:----.       ::::::.        :-:          .:---:            .:---:              :---:.          ...        :-+-:    :+:         ...          \n       +@@@@@@@@@*.   =@@@@@@@@@*:   .@@@.       =%@@@@@@@=         %@@%@@@@+          =@@@%@@@*        ....         .+:     :+:          ...         \n     .%@@*:   :*@@@.  =@@#   :+@@@+  .@@@.      %@@%-   :+-         =:   =@@@.        -@@#.  *@@#       ....         .+:     :+:          ...         \n     *@@#       #@@*  =@@#     :@@@: .@@@.     *@@%.                   .:#@@#         %@@=   .@@@.      ...          .+:     :+:          ...         \n     %@@+       +@@#  =@@#      @@@= .@@@.     %@@*         .....    @@@@@@*:         @@@-    @@@:      ...          .+:     :+:          ...         \n     #@@*       #@@*  =@@#     .@@@: .@@@.     #@@#        :@@@@@%   ...:+@@@+        %@@-   .@@@.      ....         .+:     :+:          ...         \n     -@@@-     +@@@:  =@@#    :%@@*  .@@@.     -@@@+     :: ......        %@@#   :-.  *@@*   +@@#       ....         .+:     :+:          ...         \n      -%@@@##%@@@@-   =@@@%%%@@@%=   .@@@%%%%%- :%@@@%#%@@+        =@%#*#@@@%:  *@@%   #@@%*%@@#.        ...         .+:    .-+:.        ....         \n        :=+***==%@@#+ .=++++==-.      =+++++++:   :=+**+=:          -=+**+=:    .++-    :+***=:          ...         .+:    =+++-        ...          \n                 :+#%.                                                                                    ...        .+:    .++=        ...           \n                                                                                                           ...       .+:     -+:        ...           \n      -##########################################################################################           ...      .+:      =        ..             \n      .::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::             ..  ::::-::::::::::::  ...              \n                                                                                                                  #################                   \n                                                                                                                                                      \n                                                                                                                                                      \n" );
    Log::wrapInBar( "System Parameters" );
    Log::L1( "Version: {} ({})\n\n", GLOBAL_PROGRAM_VERSION, GLOBAL_PROGRAM_LASTCHANGE );

    Log::wrapInBar( "Electronic Configuration", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    for ( auto &[name, mat] : input_electronic ) {
        Log::L1( "Electronic State: {} with energy {:.8e} Hz - {:.8e} eV - {:.8f} nm\n", name, mat.numerical["Energy"], mat.numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ), mat.numerical["Energy"].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
        for ( auto i = 0; i < mat.string_v["CoupledTo"].size(); i++ ) {
            if ( mat.string_v["CoupledTo"][i].front() == '-' ) break;
            Log::L1( " - Coupled to {} with transition energy {:.8e} Hz - {:.8e} eV\n", mat.string_v["CoupledTo"][i], std::abs( mat.numerical["Energy"] - input_electronic[mat.string_v["CoupledTo"][i]].numerical["Energy"] ), std::abs( mat.numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ) - input_electronic[mat.string_v["CoupledTo"][i]].numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ) ) );
        }
        Log::L1( " - Decay Scaling: {}\n", mat.numerical["DecayScaling"] );
        Log::L1( " - Dephasing Scaling: {}\n", mat.numerical["DephasingScaling"] );
    }
    Log::L1( "\n" );
    if ( input_electronic.size() > 0 ) {
        Log::wrapInBar( "Photonic Configuration", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
        for ( auto &[name, mat] : input_photonic ) {
            Log::L1( "Cavity: {} with energy {:.8e} Hz - {:.8e} eV - {:.8f} nm\n", name, mat.numerical["Energy"], mat.numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ), mat.numerical["Energy"].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
            for ( auto i = 0; i < mat.string_v["CoupledTo"].size(); i++ ) {
                Log::L1( " - Coupled to electronic transition {} \n", mat.string_v["CoupledTo"][i] );
                Log::L1( " - - Coupling Scaling: {}\n", mat.numerical_v["CouplingScaling"][i] );
            }
            Log::L1( " - Cavity Q-Factor: {:.0f}\n", mat.numerical["Energy"] / p_omega_cavity_loss );
            Log::L1( " - Decay Scaling: {}\n", mat.numerical["DecayScaling"] );
        }
        Log::L1( "\n" );
    }

    Log::wrapInBar( "Coupling Parameters", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Coupling strengh g: {:.8e} Hz - {:.8} mueV\n", p_omega_coupling, p_omega_coupling.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Photon loss rate k: {:.8e} Hz - {:.8} mueV\n", p_omega_cavity_loss, p_omega_cavity_loss.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Atomic dephasing rate gamma_pure: {:.8e} Hz - {:.8} mueV\n", p_omega_pure_dephasing, p_omega_pure_dephasing.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "RAD rate gamma: {:.8e} Hz - {:.8} mueV\n\n", p_omega_decay, p_omega_decay.getSI( Parameter::UNIT_ENERGY_MUEV ) );

    Log::wrapInBar( "Initial System Parameters", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Initial state rho0 = {}\n\n", info.at( p_initial_state ) );

    Log::wrapInBar( "Pulse", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_pulse.size() > 0 ) {
        for ( auto &[name, mat] : input_pulse ) {
            Log::L1( "Pulse {}:\n", name );
            Log::L1( " - Coupled to Transitions: " );
            for ( auto &mode : mat.string_v["CoupledTo"] )
                Log::L1( "{} ", mode );
            Log::L1( "\n" );
            for ( int i = 0; i < mat.numerical_v["Amplitude"].size(); i++ ) {
                Log::L1( " - Single Pulse:\n" );
                Log::L1( " - - Amplitude: {} Hz - {:.8} mueV\n", mat.numerical_v["Amplitude"][i], mat.numerical_v["Amplitude"][i].getSI( Parameter::UNIT_ENERGY_MUEV ) );
                Log::L1( " - - Frequency: {} Hz - {:.8} mueV\n", mat.numerical_v["Frequency"][i], mat.numerical_v["Frequency"][i].getSI( Parameter::UNIT_ENERGY_EV ) );
                Log::L1( " - - Width: {} s - {:.8} ps\n", mat.numerical_v["Width"][i], mat.numerical_v["Width"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Center: {} s - {:.8} ps\n", mat.numerical_v["Center"][i], mat.numerical_v["Center"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Chirp: {}\n", mat.numerical_v["Chirp"][i] );
                Log::L1( " - - Type: {}\n", mat.string_v["Type"][i] );
            }
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using any pulses to exite or drive the system.\n\n" );
    }
    Log::wrapInBar( "Chirp", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_chirp.size() > 0 ) {
        for ( auto &[name, mat] : input_chirp ) {
            Log::L1( "Chirp: {}\n", name );
            Log::L1( " - Coupled to States:\n" );
            for ( int i = 0; i < mat.string_v["CoupledTo"].size(); i++ )
                Log::L1( " - - {} with scaling {}\n", mat.string_v["CoupledTo"][i], mat.numerical_v["AmpFactor"][i] );
            for ( int i = 0; i < mat.numerical_v["Amplitude"].size(); i++ ) {
                Log::L1( " - Chirp Point {}:\n", i );
                Log::L1( " - - Amplitude: {} Hz - {} meV:\n", mat.numerical_v["Amplitude"][i], mat.numerical_v["Amplitude"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                Log::L1( " - - Time: {} s - {} ps:\n", mat.numerical_v["Times"][i], mat.numerical_v["Times"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Derivative DDT: {}\n", mat.numerical_v["ddt"][i] );
            }
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using any electronic chirps to shift the system.\n\n" );
    }

    Log::wrapInBar( "Phonons", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( p_phonon_T >= 0 ) {
        std::vector<std::string> approximations = {"Transformation integral via d/dt chi = -i/hbar*[H,chi] + d*chi/dt onto interaction picture chi(t-tau)", "Transformation Matrix U(t,tau)=exp(-i/hbar*H_DQ_L(t)*tau) onto interaction picture chi(t-tau)", "No Transformation, only interaction picture chi(t-tau)", "Analytical Lindblad formalism", "Mixed", "Path Integral"};
        Log::L1( "Temperature = {} k\nCutoff energy = {} meV\nCutoff Time = {} ps\nAlpha = {}\n<B> = {}\nFirst Markov approximation used? (rho(t) = rho(t-tau)) - {}\nTransformation approximation used: {} - {}\n", p_phonon_T, p_phonon_wcutoff.getSI( Parameter::UNIT_ENERGY_MEV ), p_phonon_tcutoff * 1E12, p_phonon_alpha, p_phonon_b, ( numerics_phonon_approximation_markov1 == 1 ? "Yes" : "No" ), numerics_phonon_approximation_order, approximations.at( numerics_phonon_approximation_order ) );
        // Pathintegral
        if ( numerics_phonon_approximation_order == 5 ) {
            Log::L1( " - Path Integral Settings:\n" );
            Log::L1( " - Backsteps NC: {}\n", p_phonon_nc );
            Log::L1( " - Iterator Stepsize: {}\n", numerics_pathintegral_stepsize_iterator );
            Log::L1( " - Thresholds: Squared({}), SparsePrune({}), CutoffIterations({}), PropagatorMapping({})\n", numerics_pathintegral_squared_threshold, numerics_pathintegral_sparse_prune_threshold, numerics_pathintegral_dynamiccutoff_iterations_max, numerics_pathintegral_docutoff_propagator );
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using phonons\n\n" );
    }

    Log::wrapInBar( "Numerical Parameters" );
    Log::L1( "\n" );
    Log::wrapInBar( "Time", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Timeborder start: {:.8e} s - {:.2f} ps\n", t_start, t_start * 1E12 );
    Log::L1( "Timeborder end: {:.8e} s - {:.2f} ps{}\n", t_end, t_end * 1E12, numerics_calculate_till_converged ? " (variable time end at 99.9\% convergence)" : "" );
    Log::L1( "Timeborder delta: {:.8e} s - {:.2f} fs \n", t_step, t_step * 1E15 );
    if ( numerics_phonon_approximation_order == 5 ) {
        Log::L1( "Timeborder delta path integral: {:.8e} s - {:.2f} ps\n", t_step_pathint, t_step_pathint * 1E12 );
    }
    //Log::L1( "Ideal time delta for this calculation: {:.8e} s - {:.2f} fs, minimum possible: {:.8e} s - {:.2f} fs\n", getIdealTimestep(), getIdealTimestep() * 1E15, std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::epsilon() * 1E15 );
    //int works = 1;
    //if ( ( init_rabifrequenz != 0.0 ) && ( 3. * t_step > 2. * M_PI / init_rabifrequenz ) )
    //    works = 0;
    //else if ( max_rabifrequenz != 0.0 && 3. * t_step > 2. * M_PI / max_rabifrequenz )
    //    works = 0;
    //if ( !works ) {
    //    if ( output_handlerstrings )
    //        fmt::print( "{} WARNING: Step may be too small to resolve predicted oscillation: dT needed vs dT: {:.10e} < {:.10e}\n", PREFIX_WARNING, 2. / 3. * M_PI / std::max( init_rabifrequenz, max_rabifrequenz ), t_step );
    //    Log::L1( "WARNING: Step may be too small to resolve predicted oscillation: \n- delta T needed: {:.10e} \n- delta T used: {:.10e}\n", 2. / 3. * M_PI / std::max( init_rabifrequenz, max_rabifrequenz ), t_step );
    //}
    Log::L1( "Time iterations (main loop) = {}\n\n", iterations_t_max );

    Log::wrapInBar( "G-Function Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_correlation.size() > 0 ) {
        Log::L1( "Anticipated tau-grid resolution is {}x{} resulting in {} skips per timestep\n", iterations_tau_resolution, iterations_tau_resolution, iterations_t_skip );
        Log::L1( "Calculating:\n" );
        for ( auto &[name, mat] : input_correlation ) {
            Log::L1( " - {} on mode(s) ", name );
            for ( auto &mode : mat.string_v["Modes"] )
                Log::L1( "{} ", mode );
            Log::L1( "\n" );
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using any G1 or G2 correlation functions.\n\n" );
    }
    Log::L1( "\n" );

    Log::wrapInBar( "Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Solver used: RK{}{}\n", numerics_rk_order, numerics_rk_order != 45 ? "" : fmt::format( " (Tolerance: {}, Stepdelta: {}, Steplimits: [{},{}])", numerics_rk_tol, numerics_rk_stepdelta, numerics_rk_stepmin, numerics_rk_stepmax ) );
    Log::L1( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
    Log::L1( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture == 1 ) ? "YES" : "NO" ) );
    Log::L1( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) ? "Analytic" : "Matrix Exponential" ) );
    Log::L1( "Threads used for primary calculations - {}\nThreads used for Secondary calculations - {}\nThreads used by Eigen: {}\n", numerics_phonons_maximum_threads, numerics_maximum_threads, Eigen::nbThreads() );
    Log::L1( "Used scaling for parameters? - {}\n", ( scale_parameters ? std::to_string( scale_value ) : "no" ) );
    if ( p_phonon_T )
        Log::L1( "Cache Phonon Coefficient Matrices? - {}\n", ( numerics_use_saved_coefficients ? fmt::format( "Yes (maximum {} matrices saved)", ( numerics_saved_coefficients_cutoff > 0 ) ? numerics_saved_coefficients_cutoff : numerics_saved_coefficients_max_size ) : "No" ) );
    if ( numerics_interpolate_outputs )
        Log::L1( "WARNING: Temporal outputs are interpolated!\n" );
    Log::L1( "\n" );
    for ( int i = 0; i < logfilecounter.size(); i++ ) {
        if ( logfilecounter[i] >= 0 )
            Log::L1( "Logfile ident number {}: {}\n", i, logfilecounter[i] );
    }
    Log::L1( "\n" );

    Log::wrapInBar( "Program Log:", Log::BAR_SIZE_FULL, Log::LEVEL_2 );
    Log::L1( "\n" );
    Log::L2( "OutputHandlerStrings: {}\n", output_handlerstrings );
}