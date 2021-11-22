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
    Log::L2( "[System] Parsing input variables...\n" );
    timer_parseInput.start();
    if ( !parseInput( arguments ) ) {
        Log::L2( "[System] Parsing input variables failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    timer_parseInput.end();
    Log::L2( "[System] Successful. Elapsed time is {}ms\n", timer_parseInput.getWallTime( Timers::MILLISECONDS ) );

    // Scaling inputs:
    if ( scale_parameters ) {
        Log::L2( "[System] Rescaling parameters to {}...\n", scale_value );
        scaleInputs( scale_value );
        Log::L2( "[System] Done!\n" );
    }

    // Adjusting inputs:
    Timer &timer_adjustInput = Timers::create( "Adjusting parameters", true, false );
    Log::L2( "[System] Adjusting input variables...\n" );
    timer_adjustInput.start();
    if ( !adjustInput() ) {
        Log::L2( "[System] Adjusting input variables failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    Log::L2( "[System] Successful. Elapsed time is {}ms\n", timer_adjustInput.getWallTime( Timers::MILLISECONDS ) );
    timer_adjustInput.end();
}

bool Parameters::parseInput( const std::vector<std::string> &arguments ) {
    //Parse_Parameters params;
    // Look for --time, if not found, standard values are used (t0 = 0, t1 = 1ns, deltaT = auto)
    t_start = QDLC::CommandlineArguments::get_parameter<double>( "--time", "tstart" );
    t_end = QDLC::CommandlineArguments::get_parameter<double>( "--time", "tend" );
    t_step = QDLC::CommandlineArguments::get_parameter<double>( "--time", "tstep" );

    // Runge Kutta Parameters
    numerics_rk_order = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rkorder" );
    numerics_rk_tol = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rktol" );
    numerics_rk_stepdelta = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepdelta" );
    numerics_rk_stepmin = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepmin" );
    numerics_rk_stepmax = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepmax" );
    numerics_rk_interpolate = !QDLC::CommandlineArguments::get_parameter_passed( "-rknointerpolate" );
    numerics_rk_usediscrete_timesteps = numerics_rk_stepdelta > 0 ? true : false;

    inputstring_electronic = QDLC::CommandlineArguments::get_parameter( "--S", "SE" );
    inputstring_photonic = QDLC::CommandlineArguments::get_parameter( "--S", "SO" );
    inputstring_pulse = QDLC::CommandlineArguments::get_parameter( "--S", "SP" );
    inputstring_chirp = QDLC::CommandlineArguments::get_parameter( "--S", "SC" );
    inputstring_spectrum = QDLC::CommandlineArguments::get_parameter( "--G", "GS" );
    inputstring_indist = QDLC::CommandlineArguments::get_parameter( "--G", "GI" );
    inputstring_conc = QDLC::CommandlineArguments::get_parameter( "--G", "GC" );
    inputstring_gfunc = QDLC::CommandlineArguments::get_parameter( "--G", "GF" );
    inputstring_wigner = QDLC::CommandlineArguments::get_parameter( "--G", "GW" );
    inputstring_correlation_resolution = QDLC::CommandlineArguments::get_parameter( "--G", "grid" );

    p_omega_coupling = QDLC::CommandlineArguments::get_parameter<double>( "--system", "coupling" );
    p_omega_cavity_loss = QDLC::CommandlineArguments::get_parameter<double>( "--system", "kappa" );
    p_omega_pure_dephasing = QDLC::CommandlineArguments::get_parameter<double>( "--system", "gammapure" );
    p_omega_decay = QDLC::CommandlineArguments::get_parameter<double>( "--system", "gamma" );
    p_initial_state_s = QDLC::CommandlineArguments::get_parameter( "--R" );

    // Look for --spectrum, if not found, no spectrum is evaluated
    iterations_tau_resolution = QDLC::CommandlineArguments::get_parameter<int>( "--spectrum", "gridres" );
    numerics_use_interactionpicture = QDLC::CommandlineArguments::get_parameter_passed( "-noInteractionpic" ) ? false : true;
    numerics_use_rwa = QDLC::CommandlineArguments::get_parameter_passed( "-noRWA" ) ? 0 : 1;
    numerics_maximum_threads = QDLC::CommandlineArguments::get_parameter<int>( "--Threads" );
    if ( numerics_maximum_threads == -1 )
        numerics_maximum_threads = omp_get_max_threads();
    output_handlerstrings = QDLC::CommandlineArguments::get_parameter_passed( "-noHandler" ) ? 0 : 1;
    output_operators = QDLC::CommandlineArguments::get_parameter_passed( "-outputOp" ) ? 2 : ( QDLC::CommandlineArguments::get_parameter_passed( "-outputHamiltons" ) ? 1 : ( QDLC::CommandlineArguments::get_parameter_passed( "-outputOpStop" ) ? 3 : 0 ) );
    numerics_order_timetrafo = QDLC::CommandlineArguments::get_parameter_passed( "-timeTrafoMatrixExponential" ) ? TIMETRANSFORMATION_MATRIXEXPONENTIAL : TIMETRANSFORMATION_ANALYTICAL;
    output_full_dm = QDLC::CommandlineArguments::get_parameter_passed( "-fullDM" );
    output_no_dm = QDLC::CommandlineArguments::get_parameter_passed( "-noDM" );
    scale_parameters = QDLC::CommandlineArguments::get_parameter_passed( "-scale" ); // MHMHMH
    numerics_use_saved_coefficients = !QDLC::CommandlineArguments::get_parameter_passed( "-disableMatrixCaching" );
    numerics_use_saved_hamiltons = !QDLC::CommandlineArguments::get_parameter_passed( "-disableHamiltonCaching" );
    numerics_use_function_caching = !QDLC::CommandlineArguments::get_parameter_passed( "-disableFunctionCaching" );
    numerics_phonons_maximum_threads = ( !numerics_use_saved_coefficients || !QDLC::CommandlineArguments::get_parameter_passed( "-disableMainProgramThreading" ) ) ? numerics_maximum_threads : 1;
    numerics_output_raman_population = QDLC::CommandlineArguments::get_parameter_passed( "-raman" ); // DEPRECATED
    logfilecounter = QDLC::Misc::convertParam<int>( QDLC::String::splitline( QDLC::CommandlineArguments::get_parameter( "--lfc" ), ',' ) );
    numerics_calculate_timeresolution_indistinguishability = QDLC::CommandlineArguments::get_parameter_passed( "-timedepInd" ); //DEPRECATED
    numerics_stretch_correlation_grid = false;                                                                                  //FIXME: Doesnt work right now //DEPRECATED
    numerics_interpolate_outputs = QDLC::CommandlineArguments::get_parameter_passed( "-interpolate" );

    // Phonon Parameters
    p_phonon_alpha = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phononalpha" );
    p_phonon_wcutoff = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phononwcutoff" );
    p_phonon_tcutoff = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phonontcutoff" );
    p_phonon_T = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "temperature" );
    numerics_phonon_approximation_order = QDLC::CommandlineArguments::get_parameter<int>( "--phonons", "phononorder" );
    numerics_phonon_approximation_markov1 = QDLC::CommandlineArguments::get_parameter_passed( "-noMarkov" ) ? 0 : 1; // First Markov
    numerics_phonon_nork45 = QDLC::CommandlineArguments::get_parameter_passed( "-noPhononRK45" );                    // Disables RK45 for phonon backwards integral; use if detunings are high. //TODO: andersrum, RK backwards integral soll rkorder nehmen, dann mit nork45 für backwards integral deaktivieren!
    output_coefficients = QDLC::CommandlineArguments::get_parameter_passed( "-phononcoeffs" ) ? 1 : 0;
    p_phonon_adjust = !QDLC::CommandlineArguments::get_parameter_passed( "-noPhononAdjust" );
    p_phonon_pure_dephasing = QDLC::Misc::convertParam<double>( "1mueV" );
    // Path Integral Parameters
    p_phonon_nc = QDLC::CommandlineArguments::get_parameter<int>( "--pathintegral", "NC" );
    numerics_pathintegral_stepsize_iterator = QDLC::CommandlineArguments::get_parameter<double>( "--pathintegral", "iteratorStepsize" );
    numerics_pathintegral_squared_threshold = QDLC::CommandlineArguments::get_parameter<double>( "--pathintegral", "squaredThreshold" );
    numerics_pathintegral_sparse_prune_threshold = QDLC::CommandlineArguments::get_parameter<double>( "--pathintegral", "sparsePruneThreshold" );
    numerics_pathintegral_dynamiccutoff_iterations_max = QDLC::CommandlineArguments::get_parameter<double>( "--pathintegral", "iteratorStepsize" );
    numerics_pathintegral_docutoff_propagator = QDLC::CommandlineArguments::get_parameter_passed( "-cutoffPropagator" );
    numerics_pathint_partially_summed = true;
    t_step_pathint = QDLC::CommandlineArguments::get_parameter<double>( "--pathintegral", "tstepPath" );

    kb = 1.3806488E-23;   // J/K, scaling needs to be for energy
    hbar = 1.0545718E-34; // J/s, scaling will be 1

    parse_system();

    subfolder = arguments.back();
    return true;
}

//TODO: this is broken. and unncesesary actually. remove. oder vernünftig alles scalen.
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
    Log::L2( "[System] Adjusting Inputs...\n" );

    if (output_handlerstrings)
        Timers::toggleHandler();

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
                    mat.numerical_v["Amplitude"][i] = mat.numerical_v["Amplitude"][i] * QDLC::Math::PI / ( std::sqrt( 2.0 * QDLC::Math::PI * mat.numerical_v["Width"][i] * std::sqrt( std::pow( mat.numerical_v["Chirp"][i] / mat.numerical_v["Width"][i], 2.0 ) + std::pow( mat.numerical_v["Width"][i], 2.0 ) ) ) ) / 2.0; //https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.241306
                else if ( mat.string_v["Type"][i].find( "cutoff" ) != std::string::npos )
                    mat.numerical_v["Amplitude"][i] = mat.numerical_v["Amplitude"][i] * QDLC::Math::PI / ( std::sqrt( 2.0 * QDLC::Math::PI * mat.numerical_v["Width"][i] * mat.numerical_v["Width"][i] ) ) / 2.0; //https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.241306
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
    if ( t_step < 0 ) {
        t_step = 1E-13; //std::min( scaleVariable( 1E-13, scale_value ), getIdealTimestep() );
        //t_step = std::max( std::numeric_limits<double>::epsilon(), t_step );
    }
    if ( t_end < 0 ) {
        // If this is given, we calculate the t-direction until 99% ground state poulation is reached after any pulses.
        numerics_calculate_till_converged = true;
        for ( auto &[name, mat] : input_pulse ) {
            for ( auto &t : mat.numerical_v["Center"] ) {
                t_end = std::max( t_end.get(), 2.0 * t );
            }
        }
        if ( t_end < 0 )
            t_end = 10E-12;
        Log::L2( "[System] Calculate till at least {} and adjust accordingly to guarantee convergence.\n", t_end );
    }

    if ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
        numerics_use_saved_hamiltons = false;
    }

    // Calculate stuff for RK
    iterations_t_max = (int)std::ceil( ( t_end - t_start ) / ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? t_step_pathint : t_step ) );
    iterations_t_skip = std::max( 1.0, std::ceil( iterations_t_max / iterations_tau_resolution ) );

    // Build dt vector. Use standard if not specified otherwise for all calculations. Path integral cannot use other timestep than the original.
    {
        auto &settings = input_correlation_resolution.count( "Modified" ) ? input_correlation_resolution["Modified"] : input_correlation_resolution["Standard"];
        double skip = input_correlation_resolution.count( "Modified" ) == 0 ? 1.0*iterations_t_skip : 1.0;
        double t_t = 0;
        int current = 0;
        grid_values.emplace_back( t_start );
        grid_value_indices[t_start] = 0;
        while ( t_t < t_end ) {
            if ( t_t > settings.numerical_v["Time"][current] and current < settings.numerical_v["Time"].size() )
                current++;
            grid_steps.emplace_back( settings.numerical_v["Delta"][current]*skip );
            t_t += grid_steps.back();
            grid_values.emplace_back( t_t );
            grid_value_indices[t_t] = grid_values.size()-1;
        }
        //std::cout << "Values for "<<mode<<": " << t_values[mode] << std::endl;
    }

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
    Log::L2( "[System] Adjusting Inputs Done!\n" );
    return true;
}

void Parameters::parse_system() {
    // Generate the input variables for the electronic system:
    auto levels = QDLC::String::splitline( inputstring_electronic, ';' );
    for ( std::string &level : levels ) {
        auto conf = QDLC::String::splitline( level, ':' );
        input_s conf_s;
        conf_s.numerical["Energy"] = QDLC::Misc::convertParam<Parameter>( conf[1] );           // Energy
        conf_s.string_v["CoupledTo"] = QDLC::String::splitline( conf[2], ',' );                // Coupled to Levels
        conf_s.numerical["DecayScaling"] = QDLC::Misc::convertParam<Parameter>( conf[3] );     // Decay Scaling, Per Mode
        conf_s.numerical["DephasingScaling"] = QDLC::Misc::convertParam<Parameter>( conf[4] ); // Dephasing Scaling
        conf_s.numerical["PhononCoupling"] = QDLC::Misc::convertParam<Parameter>( conf[5] );   // Phonon Coupling
        input_electronic[conf[0]] = conf_s;
    }
    auto cavities = QDLC::String::splitline( inputstring_photonic, ';' );
    for ( std::string &cavity : cavities ) {
        auto conf = QDLC::String::splitline( cavity, ':' );
        input_s conf_s;
        conf_s.numerical["Energy"] = QDLC::Misc::convertParam<Parameter>( conf[1] );                                            // Energy
        conf_s.numerical["MaxPhotons"] = QDLC::Misc::convertParam<Parameter>( conf[2] );                                        // Maximum Photons
        conf_s.string_v["CoupledTo"] = QDLC::String::splitline( conf[3], ',' );                                                 // Coupled to Transitions
        conf_s.numerical_v["CouplingScaling"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[4], ',' ) ); // Coupling Scaling, per transition INTO cavity
        //conf_s.numerical_v["BackCouplingScaling"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[5], ',' ) ); // BackCoupling Scaling, per transition from cavity back into the electronic system
        conf_s.numerical["DecayScaling"] = QDLC::Misc::convertParam<Parameter>( conf[5] ); // Decay Scaling, for all transitions
        input_photonic[conf[0]] = conf_s;
    }

    //TODO: erst pulsetype übergeben, dann parameter parsen. verschiedene types können dann auch verschieden viele paramter haben.
    auto pulses = QDLC::String::splitline( inputstring_pulse, ';' );
    for ( std::string &pulse : pulses ) {
        auto conf = QDLC::String::splitline( pulse, ':' );
        input_s conf_s;
        conf_s.string_v["CoupledTo"] = QDLC::String::splitline( conf[1], ',' );                                                                                                                                 // Coupled to Transitions
        conf_s.numerical_v["Amplitude"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) );                                                                                       // Pulse Amp
        conf_s.numerical_v["Frequency"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) );                                                                                       // Frequency
        conf_s.numerical_v["Width"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[4], ',' ) );                                                                                           // Width
        conf_s.numerical_v["Center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[5], ',' ) );                                                                                          // Center
        conf_s.numerical_v["Chirp"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[6], ',' ) );                                                                                           //TODO: move one down so it becomes optional                                                                                        // Chirp
        conf_s.string_v["Type"] = QDLC::String::splitline( conf[7], ',' );                                                                                                                                      // Type
        conf_s.numerical_v["SuperAmp"] = conf.size() > 8 ? QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[8], ',' ) ) : std::vector<Parameter>( conf_s.numerical_v["Center"].size(), 2.0 ); // Optional: SuperGaussian Amplitude
        input_pulse[conf[0]] = conf_s;
    }
    auto chirps = QDLC::String::splitline( inputstring_chirp, ';' );
    for ( std::string &chirp : chirps ) {
        auto conf = QDLC::String::splitline( chirp, ':' );
        input_s conf_s;
        conf_s.string_v["CoupledTo"] = QDLC::String::splitline( conf[1], ',' );                                           // Coupled to Transitions
        conf_s.numerical_v["AmpFactor"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) ); // Amplitude Scaling for coupled_to
        conf_s.numerical_v["Amplitude"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) ); // Amplitudes
        conf_s.numerical_v["Times"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[4], ',' ) );     // "Times"
        conf_s.numerical_v["ddt"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[5], ',' ) );       // "d/dt"
        conf_s.string["Type"] = conf[6];                                                                                  // Type
        input_chirp[conf[0]] = conf_s;
    }
    for ( std::string &spectrum : QDLC::String::splitline( inputstring_spectrum, ';' ) ) {
        auto conf = QDLC::String::splitline( spectrum, ':' );
        input_s conf_s;
        conf_s.string_v["Modes"] = QDLC::String::splitline( conf[0], ',' );                                            // Modes to calculate Spectrum for. Single modes can again be split with "+", meaning a+b;a to calculate for a+b and a seperately
        conf_s.numerical_v["Center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) ); // Center
        conf_s.numerical_v["Range"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) );  // Range
        conf_s.numerical_v["resW"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) );   // Resolution for w
        input_correlation["Spectrum"] = conf_s;
    }
    for ( std::string &indist : QDLC::String::splitline( inputstring_indist, ';' ) ) {
        auto conf = QDLC::String::splitline( indist, ':' );
        input_s conf_s;
        conf_s.string_v["Modes"] = QDLC::String::splitline( conf[0], ',' ); // Modes to calculate Indistinguishgability for. Single modes can again be split with "+", meaning a+b;a to calculate for a+b and a seperately
        input_correlation["Indist"] = conf_s;
    }
    for ( std::string &conc : QDLC::String::splitline( inputstring_conc, ';' ) ) {
        auto conf = QDLC::String::splitline( conc, ':' );
        input_s conf_s;
        conf_s.string_v["Modes"] = QDLC::String::splitline( conf[0], ',' ); // Modes to calculate Concurrence for
        input_correlation["Conc"] = conf_s;
    }
    for ( std::string &g_func : QDLC::String::splitline( inputstring_gfunc, ';' ) ) {
        auto conf = QDLC::String::splitline( g_func, ':' );
        auto n = conf.size();
        input_s conf_s;
        conf_s.string_v["Modes"] = QDLC::String::splitline( conf[0], ',' );                                                                                                                           // Modes to calculate G1/G2 functions for
        conf_s.numerical_v["Order"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) );                                                                                 // 1 or 2
        conf_s.numerical_v["Integrated"] = QDLC::Misc::convertParam<Parameter>( n > 2 ? QDLC::String::splitline( conf[2], ',' ) : std::vector<std::string>( conf_s.string_v["Modes"].size(), "2" ) ); // 0 or 1 or 2 for false/true/both
        input_correlation["GFunc"] = conf_s;
    }
    for ( std::string &wigner : QDLC::String::splitline( inputstring_wigner, ';' ) ) {
        auto conf = QDLC::String::splitline( wigner, ':' );
        auto n = conf.size();
        input_s conf_s;
        conf_s.string_v["Modes"] = QDLC::String::splitline( conf[0], ',' );                                                                                                                     // Modes to calculate Wigner function for
        conf_s.numerical_v["X"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) );                                                                               // -X to X
        conf_s.numerical_v["Y"] = QDLC::Misc::convertParam<Parameter>( n > 2 ? QDLC::String::splitline( conf[2], ',' ) : QDLC::String::splitline( conf[1], ',' ) );                             // -Y to Y
        conf_s.numerical_v["Res"] = QDLC::Misc::convertParam<Parameter>( n > 3 ? QDLC::String::splitline( conf[3], ',' ) : std::vector<std::string>( conf_s.numerical_v["X"].size(), "100" ) ); // Resolution
        conf_s.numerical_v["Skip"] = QDLC::Misc::convertParam<Parameter>( n > 4 ? QDLC::String::splitline( conf[4], ',' ) : std::vector<std::string>( conf_s.numerical_v["X"].size(), "1" ) );  // Skips in t-direction
        input_correlation["Wigner"] = conf_s;
    }
    for ( std::string &gconf : QDLC::String::splitline( inputstring_correlation_resolution, ';' ) ) { // Redundant, only the latest will be used.
        auto single = QDLC::String::splitline( gconf, ':' );
        input_s conf_s;
        std::vector<Parameter> times, dts;
        for ( int i = 0; i < single.size(); i++ ) {
            auto cur = QDLC::String::splitline( single[i], '-' );
            times.emplace_back( QDLC::Misc::convertParam<Parameter>( cur[0] ) );
            dts.emplace_back( QDLC::Misc::convertParam<Parameter>( cur[1] ) );
        }
        conf_s.numerical_v["Time"] = times;
        conf_s.numerical_v["Delta"] = dts;
        input_correlation_resolution["Modified"] = conf_s;
    }
    {
        input_s conf_s;
        conf_s.numerical_v["Time"] = { t_end };
        conf_s.numerical_v["Delta"] = { t_step };
        input_correlation_resolution["Standard"] = conf_s;
    }
}

void Parameters::log( const Dense &initial_state_vector_ket ) {
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
    Log::L1( "Initial state rho0 = [{}]\n\n", initial_state_vector_ket.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ) );

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
                Log::L1( " - - Frequency: {} Hz - {:.8} eV\n", mat.numerical_v["Frequency"][i], mat.numerical_v["Frequency"][i].getSI( Parameter::UNIT_ENERGY_EV ) );
                Log::L1( " - - Width: {} s - {:.8} ps\n", mat.numerical_v["Width"][i], mat.numerical_v["Width"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Center: {} s - {:.8} ps\n", mat.numerical_v["Center"][i], mat.numerical_v["Center"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Chirp: {}\n", mat.numerical_v["Chirp"][i] );
                Log::L1( " - - Type: {}{}\n", mat.string_v["Type"][i], mat.string_v["Type"][i] == "gauss" ? fmt::format( " (Gaussian Amplitude: {})", mat.numerical_v["SuperAmp"][i] ) : "" );
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
        std::vector<std::string> approximations = { "Transformation integral via d/dt chi = -i/hbar*[H,chi] + d*chi/dt onto interaction picture chi(t-tau)", "Transformation Matrix U(t,tau)=exp(-i/hbar*H_DQ_L(t)*tau) onto interaction picture chi(t-tau)", "No Transformation, only interaction picture chi(t-tau)", "Analytical Lindblad formalism", "Mixed", "Path Integral" };
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
    Log::L1( "Time iterations (main loop) = {}\n\n", iterations_t_max );

    Log::wrapInBar( "G-Function Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_correlation.size() > 0 ) {
        Log::L1( "Tau-grid resolution is {}x{}\n", grid_values.size(), grid_values.size() );
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
    if ( numerics_rk_order == 45 and numerics_phonon_nork45 )
        Log::L1( "Will NOT use RK45 for the phonon backwards integral!\n" );
    Log::L1( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
    Log::L1( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture ) ? "YES" : "NO" ) );
    Log::L1( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) ? "Analytic" : "Matrix Exponential" ) );
    Log::L1( "Threads used for primary calculations - {}\nThreads used for Secondary calculations - {}\nThreads used by Eigen: {}\n", numerics_phonons_maximum_threads, numerics_maximum_threads, Eigen::nbThreads() );
    Log::L1( "Used scaling for parameters? - {}\n", ( scale_parameters ? std::to_string( scale_value ) : "no" ) );
    if ( p_phonon_T )
        Log::L1( "Cache Phonon Coefficient Matrices? - {}\n", ( numerics_use_saved_coefficients ? fmt::format( "Yes (maximum {} matrices saved)", ( numerics_saved_coefficients_cutoff > 0 ) ? numerics_saved_coefficients_cutoff : numerics_saved_coefficients_max_size ) : "No" ) );
    if ( numerics_interpolate_outputs )
        Log::L1( "WARNING: Temporal outputs are interpolated!\n" );
    if ( !numerics_use_function_caching )
        Log::L1( "NOT using function caching.\n" );
    if ( !numerics_use_saved_hamiltons )
        Log::L1( "NOT using Hamilton caching.\n" );
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