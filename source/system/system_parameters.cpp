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
    scale_parameters = false;
    scale_value = 1E12;
    iterations_t_max = 1;
    maxStates = 0;
    numerics_calculate_till_converged = false;

    // Parsing input:
    Timer &timer_parseInput = Timers::create( "Parsing parameters", true, false );
    Log::wrapInBar( "Conversion of input variables", Log::BAR_SIZE_FULL, Log::LEVEL_2, Log::BAR_0 );
    LOG2( "\n" );
    LOG2( "[System] Parsing input variables...\n" );

    timer_parseInput.start();
    parse_input( arguments );
    timer_parseInput.end();
    LOG2( "[System] Successful. Elapsed time is {}ms\n", timer_parseInput.getWallTime( Timers::MILLISECONDS ) );

    // Scaling inputs:
    if ( scale_parameters ) {
        LOG2( "[System] Rescaling parameters to {}...\n", scale_value );
        scale_inputs( scale_value );
    }

    // Adjusting inputs:
    Timer &timer_adjustInput = Timers::create( "Adjusting parameters", true, false );
    LOG2( "[System] Adjusting input variables...\n" );

    timer_adjustInput.start();
    adjust_input();
    timer_adjustInput.end();
    LOG2( "[System] Successful. Elapsed time is {}ms\n", timer_adjustInput.getWallTime( Timers::MILLISECONDS ) );
}

void Parameters::parse_input( const std::vector<std::string> &arguments ) {
    // Parse_Parameters params;
    //  Look for --time, if not found, standard values are used (t0 = 0, t1 = 1ns, deltaT = auto)
    t_start = QDLC::CommandlineArguments::get_parameter<double>( "--time", "tstart" );
    t_end = QDLC::CommandlineArguments::get_parameter<double>( "--time", "tend" );
    t_step = QDLC::CommandlineArguments::get_parameter<double>( "--time", "tstep" );
    numerics_groundstate_string = QDLC::CommandlineArguments::get_parameter( "--groundstate" );

    // Runge Kutta Parameters
    numerics_rk_order = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rkorder" );
    inputstring_rk45_config = QDLC::CommandlineArguments::get_parameter( "--rk", "rktol" );
    numerics_rk_stepdelta = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepdelta" );
    numerics_rk_stepmin = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepmin" );
    numerics_rk_stepmax = QDLC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepmax" );
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
    inputstring_raman = QDLC::CommandlineArguments::get_parameter( "--G", "GR" );
    inputstring_correlation_resolution = QDLC::CommandlineArguments::get_parameter( "--G", "grid" );
    inputstring_detector = QDLC::CommandlineArguments::get_parameter( "--detector" );
    inputstring_densitymatrix_config = QDLC::CommandlineArguments::get_parameter( "--DMconfig" );
    inputstring_SPconf = QDLC::CommandlineArguments::get_parameter( "--SPconfig" );

    p_omega_coupling = QDLC::CommandlineArguments::get_parameter<double>( "--system", "coupling" );
    p_omega_cavity_loss = QDLC::CommandlineArguments::get_parameter<double>( "--system", "kappa" );
    p_omega_pure_dephasing = QDLC::CommandlineArguments::get_parameter<double>( "--system", "gammapure" );
    p_omega_decay = QDLC::CommandlineArguments::get_parameter<double>( "--system", "gamma" );
    p_initial_state_s = QDLC::CommandlineArguments::get_parameter( "--R" );

    grid_resolution = QDLC::CommandlineArguments::get_parameter<int>( "--G", "gridres" );
    numerics_use_interactionpicture = not QDLC::CommandlineArguments::get_parameter_passed( "-noInteractionpic" );
    numerics_use_rwa = not QDLC::CommandlineArguments::get_parameter_passed( "-noRWA" );
    numerics_maximum_threads = QDLC::CommandlineArguments::get_parameter<int>( "--Threads" );
    if ( numerics_maximum_threads == -1 )
        numerics_maximum_threads = omp_get_max_threads();
    output_handlerstrings = QDLC::CommandlineArguments::get_parameter_passed( "-handler" );
    output_operators = QDLC::CommandlineArguments::get_parameter_passed( "-outputOp" ) ? 2 : ( QDLC::CommandlineArguments::get_parameter_passed( "-outputHamiltons" ) ? 1 : ( QDLC::CommandlineArguments::get_parameter_passed( "-outputOpStop" ) ? 3 : 0 ) );
    output_eigenvalues = QDLC::CommandlineArguments::get_parameter_passed( "-outputEigs" );
    numerics_order_timetrafo = QDLC::CommandlineArguments::get_parameter_passed( "-timeTrafoMatrixExponential" ) ? TIMETRANSFORMATION_MATRIXEXPONENTIAL : TIMETRANSFORMATION_ANALYTICAL;
    scale_parameters = QDLC::CommandlineArguments::get_parameter_passed( "-scale" ); // MHMHMH
    numerics_use_saved_coefficients = !QDLC::CommandlineArguments::get_parameter_passed( "-disableMatrixCaching" );
    numerics_use_saved_hamiltons = !QDLC::CommandlineArguments::get_parameter_passed( "-disableHamiltonCaching" );
    numerics_use_function_caching = !QDLC::CommandlineArguments::get_parameter_passed( "-disableFunctionCaching" );
    numerics_force_caching = false; // If true, even if any saving was disabled internally (not by the user), the matrices will still be cached.
    numerics_phonons_maximum_threads = ( !numerics_use_saved_coefficients || !QDLC::CommandlineArguments::get_parameter_passed( "-disableMainProgramThreading" ) ) ? numerics_maximum_threads : 1;
    logfilecounter = QDLC::Misc::convertParam<double>( QDLC::String::splitline( QDLC::CommandlineArguments::get_parameter( "--lfc" ), ',' ) );
    numerics_interpolate_outputs = QDLC::CommandlineArguments::get_parameter_passed( "-interpolate" );
    s_numerics_interpolate = QDLC::CommandlineArguments::get_parameter( "--interpolateOrder" );
    numerics_output_rkerror = QDLC::CommandlineArguments::get_parameter_passed( "-error" );
    output_detector_transformations = QDLC::CommandlineArguments::get_parameter_passed( "-outputDetector" );

    // Phonon Parameters
    p_phonon_alpha = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phononalpha" );
    p_phonon_wcutoff = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phononwcutoff" );
    p_phonon_wcutoffdelta = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phononwcutoffdelta" );
    p_phonon_tcutoff = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phonontcutoff" );
    p_phonon_T = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "temperature" );
    numerics_phonon_approximation_order = QDLC::CommandlineArguments::get_parameter<int>( "--phonons", "phononorder" );
    numerics_phonon_approximation_markov1 = QDLC::CommandlineArguments::get_parameter_passed( "-noMarkov" ) ? 0 : 1; // First Markov
    numerics_phonon_nork45 = !QDLC::CommandlineArguments::get_parameter_passed( "-usePhononRK45" );                  // Enables. RK45 for phonon backwards integral; use if detunings are low, otherwise expensive.
    output_coefficients = QDLC::CommandlineArguments::get_parameter_passed( "-phononcoeffs" ) ? 1 : 0;
    output_path = QDLC::CommandlineArguments::get_parameter_passed( "-oPath" ) ? 1 : 0;
    p_phonon_adjust = !QDLC::CommandlineArguments::get_parameter_passed( "-noPhononAdjust" );
    p_phonon_pure_dephasing = QDLC::Misc::convertParam<double>( "1mueV" );

    // Phonon Quantum Dot Paramters. These are used to overwrite the wcutoff and alpha values if needed.
    p_phonon_qd_de = QDLC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDDe" );
    p_phonon_qd_dh = QDLC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDDh" );
    p_phonon_qd_rho = QDLC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDrho" );
    p_phonon_qd_cs = QDLC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDcs" );
    p_phonon_qd_ratio = QDLC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDratio" );
    p_phonon_qd_ae = QDLC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDae" );

    // Path Integral Parameters
    p_phonon_nc = QDLC::CommandlineArguments::get_parameter<int>( "--pathintegral", "NC" );
    numerics_subiterator_stepsize = QDLC::CommandlineArguments::get_parameter<double>( "--pathintegral", "iteratorStepsize" );
    numerics_pathintegral_squared_threshold = QDLC::CommandlineArguments::get_parameter<double>( "--numericalpathintegral", "squaredThreshold" );
    numerics_pathintegral_sparse_prune_threshold = QDLC::CommandlineArguments::get_parameter<double>( "--numericalpathintegral", "sparsePruneThreshold" );
    numerics_pathintegral_dynamiccutoff_iterations_max = QDLC::CommandlineArguments::get_parameter<double>( "--numericalpathintegral", "cutoffADM" ); // FIXME ????
    numerics_pathintegral_sparse_to_dense_threshold = QDLC::CommandlineArguments::get_parameter<double>( "--numericalpathintegral", "denseTensorThreshold" );
    numerics_pathintegral_force_dense = QDLC::CommandlineArguments::get_parameter_passed( "-pathIntForceDense" );
    numerics_pathintegral_docutoff_propagator = QDLC::CommandlineArguments::get_parameter_passed( "-cutoffPropagator" );
    numerics_pathint_partially_summed = !QDLC::CommandlineArguments::get_parameter_passed( "-disablePSPath" );
    t_step_pathint = QDLC::CommandlineArguments::get_parameter<double>( "--pathintegral", "tstepPath" );

    kb = 1.3806488E-23;   // J/K, scaling needs to be for energy
    hbar = 1.0545718E-34; // J/s, scaling will be 1

    parse_system();

    working_directory = arguments.back();
}

// TODO: this is broken. and unncesesary actually. remove. oder vernünftig alles scalen.
void Parameters::scale_inputs( const double scaling ) {
    // Adjust normal parameters: time is multiplid by scaling, frequency divided
    t_start.setScale( scaling, Parameter::SCALE_TIME );
    t_end.setScale( scaling, Parameter::SCALE_TIME );
    t_step.setScale( scaling, Parameter::SCALE_TIME );
    numerics_rk_stepdelta.setScale( scaling, Parameter::SCALE_TIME );

    // TODO: inputs scalen! pulse un chirp auch!
    // p_omega_atomic_G_H.setScale( scaling, Parameter::SCALE_ENERGY );
    // p_omega_atomic_G_V.setScale( scaling, Parameter::SCALE_ENERGY );
    // p_omega_atomic_H_B.setScale( scaling, Parameter::SCALE_ENERGY );
    // p_omega_atomic_V_B.setScale( scaling, Parameter::SCALE_ENERGY );
    // p_deltaE.setScale( scaling, Parameter::SCALE_ENERGY );
    // p_biexciton_bindingenergy.setScale( scaling, Parameter::SCALE_ENERGY );
    // p_omega_cavity_H.setScale( scaling, Parameter::SCALE_ENERGY );
    // p_omega_cavity_V.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_coupling.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_cavity_loss.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_pure_dephasing.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_decay.setScale( scaling, Parameter::SCALE_ENERGY );
    // Adjust chirp and pulse
    // for ( int i = 0; i < (int)chirp_t.size(); i++ ) {
    //    chirp_t.at( i ).setScale( scaling, Parameter::SCALE_TIME );
    //    chirp_y.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
    //    chirp_ddt.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
    //}
    // for ( int i = 0; i < (int)pulse_center.size(); i++ ) {
    //    pulse_center.at( i ).setScale( scaling, Parameter::SCALE_TIME );
    //    if ( pulse_type.at( i ).compare( "gauss_pi" ) != 0 )
    //        pulse_amp.at( i ) = scaleVariable( pulse_amp.at( i ), 1.0 / scaling );
    //    pulse_omega.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
    //    pulse_sigma.at( i ).setScale( scaling, Parameter::SCALE_TIME );
    //    pulse_omega_chirp.at( i ).setScale( scaling * scaling, Parameter::SCALE_ENERGY );
    //}
    // Adjusting spectrum
    // TODO: neuen input saclen!
    // spectrum_frequency_center.setScale( scaling, Parameter::SCALE_ENERGY );
    // spectrum_frequency_range.setScale( scaling, Parameter::SCALE_ENERGY );
    // Phonons
    p_phonon_wcutoff.setScale( scaling, Parameter::SCALE_ENERGY );
    p_phonon_tcutoff.setScale( scaling, Parameter::SCALE_TIME );
    p_phonon_alpha.setScale( scaling * scaling, Parameter::SCALE_ENERGY );
    p_phonon_pure_dephasing.setScale( scaling, Parameter::SCALE_ENERGY );
    kb.setScale( scaling, Parameter::SCALE_ENERGY );
}

void Parameters::adjust_input() {
    LOG2( "[System] Adjusting Inputs...\n" );

    if ( output_handlerstrings )
        Timers::toggleHandler();

    // For threadsafety
    if ( numerics_rk_order > 5 )
        numerics_use_saved_hamiltons = false;

    // Calculate/Recalculate some parameters:
    // Adjust pulse data
    // TODO: für complexe ampltiude is_imag in Parameters = true -> dann mit 1i multiplizeiren in get() funktion
    for ( auto &[name, mat] : input_pulse ) {
        mat.numerical_v["Amplitude"] = QDLC::Misc::convertParam<Parameter>( mat.string_v["Amplitude"] );
        // Set all optional parameters to default
        mat.numerical_v["ChirpRate"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );
        mat.numerical_v["GaussAmp"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 2.0 );
        mat.numerical_v["SUPERDelta"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );
        mat.numerical_v["SUPERFreq"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );
        mat.numerical_v["CutoffDelta"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );

        for ( int i = 0; i < mat.string_v["Amplitude"].size(); i++ ) {
            auto type_params = QDLC::String::splitline( mat.string_v["Type"][i], '+' );
            // Extract optional Parameters
            for ( const auto &type : type_params ) {
                if ( type.find( "cw" ) != std::string::npos )
                    mat.string_v["Type"][i] = "cw";
                else if ( type.find( "gauss" ) != std::string::npos )
                    mat.string_v["Type"][i] = "gauss";
                else {
                    std::string param = QDLC::String::splitline( type, '(' ).back();
                    param.pop_back();
                    if ( type.find( "chirped" ) != std::string::npos ) {
                        mat.numerical_v["ChirpRate"][i] = QDLC::Misc::convertParam<Parameter>( param );
                    } else if ( type.find( "cutoff" ) != std::string::npos ) {
                        mat.numerical_v["CutoffDelta"][i] = QDLC::Misc::convertParam<Parameter>( param );
                    } else if ( type.find( "super" ) != std::string::npos ) {
                        auto splitparam = QDLC::String::splitline( param, '_' );
                        mat.numerical_v["SUPERDelta"][i] = QDLC::Misc::convertParam<Parameter>( splitparam[0] );
                        mat.numerical_v["SUPERFreq"][i] = QDLC::Misc::convertParam<Parameter>( splitparam[1] );
                    } else if ( type.find( "exponent" ) != std::string::npos ) {
                        mat.numerical_v["GaussAmp"][i] = QDLC::Misc::convertParam<Parameter>( param );
                    }
                }
            }
            if ( mat.string_v["Amplitude"][i].find( "pi" ) != std::string::npos ) {
                mat.numerical_v["Amplitude"][i] = mat.numerical_v["Amplitude"][i] * QDLC::Math::PI / ( std::sqrt( 2.0 * QDLC::Math::PI * mat.numerical_v["Width"][i] * std::sqrt( std::pow( mat.numerical_v["ChirpRate"][i] / mat.numerical_v["Width"][i], 2.0 ) + std::pow( mat.numerical_v["Width"][i], 2.0 ) ) ) ) / 2.0; // https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.241306
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
        t_step = 1E-13; // std::min( scaleVariable( 1E-13, scale_value ), getIdealTimestep() );
        // t_step = std::max( std::numeric_limits<double>::epsilon(), t_step );
        LOG2( "[System] Delta t was set to {}\n", t_step );
    }

    if ( t_end >= 0 and ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? ( t_step_pathint > 0 ) : ( t_step > 0 ) ) ) {
        iterations_t_max = (int)std::ceil( ( t_end - t_start ) / ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? t_step_pathint : t_step ) );
        if ( grid_resolution < 1 and t_end >= 0 )
            grid_resolution = iterations_t_max + 1;
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
            t_end = std::max<double>( 1E-12, 10.0 * t_step );
        LOG2( "[System] Calculate till at least {} and adjust accordingly to guarantee convergence. The matrix index used is {}\n", t_end, numerics_groundstate );
    }

    if ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
        numerics_use_saved_hamiltons = false;
    }

    // Calculate stuff for RK

    LOG2( "[System] Maximum t-value for temporal calculations is {}\n", t_end );
    iterations_t_skip = std::max( 1.0, std::ceil( iterations_t_max / grid_resolution ) );

    // Build dt vector. Use standard if not specified otherwise for all calculations. Path integral cannot use other timestep than the original.
    if ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? ( t_step_pathint > 0 ) : ( t_step > 0 ) ) {
        input_correlation_resolution["Standard"].numerical_v["Time"] = { t_end };
        input_correlation_resolution["Standard"].numerical_v["Delta"] = { numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? t_step_pathint : t_step };
        auto &settings = input_correlation_resolution.contains( "Modified" ) ? input_correlation_resolution["Modified"] : input_correlation_resolution["Standard"];
        double skip = input_correlation_resolution.contains( "Modified" ) ? 1.0 : 1.0 * iterations_t_skip;
        LOG2( "[System] Iteration Skip for Grid is {}.\n", skip );
        double t_t = 0;
        int current = 0;
        grid_values.clear();
        grid_steps.clear();
        grid_value_indices.clear();
        grid_values.emplace_back( t_start );
        grid_value_indices[t_start] = 0;
        LOG2( "[System] Initial Timestep Limit is {} at a timestep of {}.\n", settings.numerical_v["Time"][current], settings.numerical_v["Delta"][current] * skip );
        while ( t_t < t_end ) {
            if ( t_t > settings.numerical_v["Time"][current] and current < settings.numerical_v["Time"].size() ) {
                current++;
                LOG2( "[System] New Timestep Limit is {} at a timestep of {}.\n", settings.numerical_v["Time"][current], settings.numerical_v["Delta"][current] * skip );
            }
            grid_steps.emplace_back( settings.numerical_v["Delta"][current] * skip );
            t_t += grid_steps.back();
            grid_values.emplace_back( t_t );
            grid_value_indices[t_t] = grid_values.size() - 1;
        }
        // std::cout << "Values for "<<mode<<": " << t_values[mode] << std::endl;
        LOG2( "[System] Setting correlation grid resolution to {0}x{0} for a t_end = {1}\n", grid_values.size(), t_end );
    } else {
        LOG2( "[System] Not setting time vector because timestep is negative!\n" );
    }

    // Set interpolation order:
    {
        auto orders = QDLC::String::splitline( s_numerics_interpolate, ',' );
        std::map<std::string, int> methods = { { "monotone", 3 }, { "linear", 0 } };
        std::string method_time = orders.front();
        std::string method_tau = orders.size() > 1 ? orders.back() : "linear";
        numerics_interpolate_method_time = methods[method_time];
        numerics_interpolate_method_tau = methods[method_tau];
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
            double J;
            if ( p_phonon_qd_ae == 0.0 )
                J = p_phonon_alpha * w * std::exp( -w * w / 2.0 / p_phonon_wcutoff / p_phonon_wcutoff );
            else
                J = w * hbar * std::pow( p_phonon_qd_de * std::exp( -w * w * p_phonon_qd_ae * p_phonon_qd_ae / ( 4. * p_phonon_qd_cs * p_phonon_qd_cs ) ) - p_phonon_qd_dh * std::exp( -w * w * p_phonon_qd_ae / p_phonon_qd_ratio * p_phonon_qd_ae / p_phonon_qd_ratio / ( 4. * p_phonon_qd_cs * p_phonon_qd_cs ) ), 2. ) / ( 4. * 3.1415 * 3.1415 * p_phonon_qd_rho * std::pow( p_phonon_qd_cs, 5. ) );

            integral += stepsize * ( J / std::tanh( hbar * w / 2.0 / kb / p_phonon_T ) );
        }
        p_phonon_b = std::exp( -0.5 * integral );
        if ( p_phonon_adjust ) {
            p_omega_pure_dephasing = p_phonon_pure_dephasing * p_phonon_T;
            p_omega_decay = p_omega_decay * p_phonon_b * p_phonon_b;
        }
        if ( numerics_rk_order >= 45 and numerics_use_saved_coefficients )
            numerics_force_caching = true;
    }
    // numerics_saved_coefficients_max_size = (int)( ( t_end - t_start ) / t_step * 2.0 * ( p_phonon_tcutoff / t_step ) ) + 10;
    trace.reserve( iterations_t_max + 5 );
}

void Parameters::parse_system() {
    // Generate the input variables for the electronic system:
    for ( auto levels = QDLC::String::splitline( inputstring_electronic, ';' ); const std::string &level : levels ) {
        auto conf = QDLC::String::splitline( level, ':' );
        input_s conf_s;
        conf_s.numerical["Energy"] = QDLC::Misc::convertParam<Parameter>( conf[1] );           // Energy
        conf_s.string_v["CoupledTo"] = QDLC::String::splitline( conf[2], ',' );                // Coupled to Levels
        conf_s.numerical["DecayScaling"] = QDLC::Misc::convertParam<Parameter>( conf[3] );     // Decay Scaling, Per Mode
        conf_s.numerical["DephasingScaling"] = QDLC::Misc::convertParam<Parameter>( conf[4] ); // Dephasing Scaling
        conf_s.numerical["PhononCoupling"] = QDLC::Misc::convertParam<Parameter>( conf[5] );   // Phonon Coupling
        input_electronic[conf[0]] = conf_s;
    }
    for ( auto cavities = QDLC::String::splitline( inputstring_photonic, ';' ); const std::string &cavity : cavities ) {
        auto conf = QDLC::String::splitline( cavity, ':' );
        input_s conf_s;
        conf_s.numerical["Energy"] = QDLC::Misc::convertParam<Parameter>( conf[1] );                                            // Energy
        conf_s.numerical["MaxPhotons"] = QDLC::Misc::convertParam<Parameter>( conf[2] );                                        // Maximum Photons
        conf_s.string_v["CoupledTo"] = QDLC::String::splitline( conf[3], ',' );                                                 // Coupled to Transitions
        conf_s.numerical_v["CouplingScaling"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[4], ',' ) ); // Coupling Scaling, per transition INTO cavity
        // conf_s.numerical_v["BackCouplingScaling"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[5], ',' ) ); // BackCoupling Scaling, per transition from cavity back into the electronic system
        conf_s.numerical["DecayScaling"] = QDLC::Misc::convertParam<Parameter>( conf[5] ); // Decay Scaling, for all transitions
        input_photonic[conf[0]] = conf_s;
    }

    // --SP 'p:GX:1pi,5pi:1.5eV,1eV:4ps,2ps:20ps,35ps:gauss:'
    // Type is cw or superposition (chained with '+', e.g. 'gauss+chirped(1E-24)' ) of gauss,cutoff,chirped(rate),super(delta),exponent(exponent)
    // p:TYPE:...parameters...
    int pindex = 0;
    for ( auto pulses = QDLC::String::splitline( inputstring_pulse, ';' ); const std::string &pulse : pulses ) {
        auto conf = QDLC::String::splitline( pulse, ':' );
        input_s conf_s;
        conf_s.string_v["CoupledTo"] = QDLC::String::splitline( conf[1], ',' );                                           // Coupled to Transitions
        conf_s.string_v["Amplitude"] = QDLC::String::splitline( conf[2], ',' );                                           // Pulse Amp
        conf_s.numerical_v["Frequency"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) ); // Frequency
        conf_s.numerical_v["Width"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[4], ',' ) );     // Width
        conf_s.numerical_v["Center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[5], ',' ) );    // Center
        conf_s.string_v["Type"] = QDLC::String::splitline( conf[6], ',' );                                                // Type
        // Non-Mandatory values
        // conf_s.numerical_v["Chirp"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[7], ',' ) );                                                                                           // TODO: move one down so it becomes optional                                                                                        // Chirp
        // conf_s.numerical_v["SuperAmp"] = conf.size() > 8 ? QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[8], ',' ) ) : std::vector<Parameter>( conf_s.numerical_v["Center"].size(), 2.0 ); // Optional: SuperGaussian Amplitude
        // For counting purposes:
        conf_s.numerical["PulseIndex"] = pindex;
        pindex += 2;
        input_pulse[conf[0]] = conf_s;
    }
    // TODO: hier auch erst Type übergeben, dann kann man verschiedene chirps (z.b. die modulation) besser implementieren.
    for ( auto chirps = QDLC::String::splitline( inputstring_chirp, ';' ); const std::string &chirp : chirps ) {
        auto conf = QDLC::String::splitline( chirp, ':' );
        input_s conf_s;
        conf_s.string_v["CoupledTo"] = QDLC::String::splitline( conf[1], ',' );                                           // Coupled to Transitions
        conf_s.numerical_v["AmpFactor"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) ); // Amplitude Scaling for coupled_to
        conf_s.numerical_v["Amplitude"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) ); // Amplitudes
        conf_s.numerical_v["Times"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[4], ',' ) );     // "Times"
        if ( conf[5].find( "," ) != std::string::npos ) {
            conf_s.numerical_v["ddt"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[5], ',' ) ); // "d/dt"
            conf_s.string["Type"] = conf[6];                                                                            // Type
        } else {
            conf_s.string["Type"] = conf[5];                                                               // Type
            conf_s.numerical_v["ddt"] = std::vector<Parameter>( conf_s.numerical_v["Times"].size(), 0.0 ); // "d/dt"
        }
        input_chirp[conf[0]] = conf_s;
    }
    for ( const std::string &spectrum : QDLC::String::splitline( inputstring_spectrum, ';' ) ) {
        auto conf = QDLC::String::splitline( spectrum, ':' );
        input_s conf_s;
        conf_s.string_v["Modes"] = QDLC::String::splitline( conf[0], ',' );                                                                                                 // Modes to calculate Spectrum for. Single modes can again be split with "+", meaning a+b;a to calculate for a+b and a seperately
        conf_s.numerical_v["Center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) );                                                      // Center
        conf_s.numerical_v["Range"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) );                                                       // Range
        conf_s.numerical_v["resW"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) );                                                        // Resolution for w
        conf_s.string_v["Normalize"] = conf.size() > 4 ? QDLC::String::splitline( conf[4], ',' ) : std::vector<std::string>( conf_s.numerical_v["Range"].size(), "False" ); // Normalize?
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
    for ( std::string &raman : QDLC::String::splitline( inputstring_raman, ';' ) ) {
        auto conf = QDLC::String::splitline( raman, ':' );
        auto n = conf.size();
        input_s conf_s;
        conf_s.string_v["SourceModes"] = QDLC::String::splitline( conf[0], ',' );
        conf_s.string_v["RamanMode"] = QDLC::String::splitline( conf[1], ',' );
        conf_s.string_v["OpMode"] = QDLC::String::splitline( conf[2], ',' );
        conf_s.string_v["PMode"] = QDLC::String::splitline( conf[3], ',' );
        input_correlation["Raman"] = conf_s;
    }
    // Correlation Grid
    if ( std::find( inputstring_correlation_resolution.begin(), inputstring_correlation_resolution.end(), ':' ) != inputstring_correlation_resolution.end() ) {
        auto single = QDLC::String::splitline( inputstring_correlation_resolution, ';' );
        input_s conf_s;
        std::vector<Parameter> times, dts;
        for ( int i = 0; i < single.size(); i++ ) {
            auto cur = QDLC::String::splitline( single[i], ':' );
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
    {
        input_s conf_s;
        if ( inputstring_SPconf == "inherit" and conf_s.numerical_v["Center"].size() > 0 ) {
            LOG2( "[System-Parameters] Inheriting Pulse Fourier Configuration from Spectrum.\n" );
            conf_s.numerical["Center"] = conf_s.numerical_v["Center"][0];
            conf_s.numerical["Range"] = conf_s.numerical_v["Range"][0];
            conf_s.numerical["Res"] = conf_s.numerical_v["resW"][0];
            conf_s.numerical["dt"] = numerics_subiterator_stepsize;
        } else if ( inputstring_SPconf.size() > 0 ) {
            LOG2( "[System-Parameters] Generating Pulse Fourier Configuration from Parameters using {}.\n", inputstring_SPconf );
            auto conf = QDLC::String::splitline( inputstring_SPconf, ':' );
            conf_s.numerical["Center"] = conf.size() > 0 ? QDLC::Misc::convertParam<Parameter>( conf[0] ) : Parameter( 0.0 );
            conf_s.numerical["Range"] = conf.size() > 1 ? QDLC::Misc::convertParam<Parameter>( conf[1] ) : Parameter( 0.0 );
            conf_s.numerical["Res"] = conf.size() > 2 ? QDLC::Misc::convertParam<Parameter>( conf[2] ) : Parameter( 0.0 );
            conf_s.numerical["dt"] = conf.size() > 3 ? QDLC::Misc::convertParam<Parameter>( conf[3] ) : Parameter( numerics_subiterator_stepsize );
            input_conf["PulseConf"] = conf_s;
        }
    }
    // Detector Input. Split temporal and spectral filtering with ";"
    {
        input_s conf_s;
        if ( !( inputstring_detector == "none" ) ) {
            LOG2( "[System-Parameters] Setting up detector using {}...\n", inputstring_detector );
            auto conf_sep = QDLC::String::splitline( inputstring_detector, ';' );
            if ( conf_sep[0] != "none" ) {
                auto conf = QDLC::String::splitline( conf_sep[0], ':' );
                conf_s.numerical_v["time_center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) );
                conf_s.numerical_v["time_range"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[0], ',' ) );
                conf_s.numerical_v["time_power_amplitude"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) );
                LOG2( "[System-Parameters] Adding Temporal Detector mask using center = {}, range = {} and power_amp = {}.\n", conf_s.numerical_v["time_center"], conf_s.numerical_v["time_range"], conf_s.numerical_v["time_power_amplitude"] );
            } else {
                conf_s.numerical_v["time_center"] = {};
                conf_s.numerical_v["time_range"] = {};
                conf_s.numerical_v["time_power_amplitude"] = {};
            }
            if ( conf_sep.size() > 1 and conf_sep[1] != "none" ) {
                auto conf = QDLC::String::splitline( conf_sep[1], ':' );
                conf_s.numerical_v["spectral_range"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[0], ',' ) );
                conf_s.numerical_v["spectral_center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) );
                conf_s.numerical_v["spectral_power_amplitude"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) );
                conf_s.numerical_v["spectral_number_points"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) );
                LOG2( "[System-Parameters] Adding Spectral Detector mask using center = {}, range = {}, power_amp = {} and {} points for the Fourier Transformation.\n", conf_s.numerical_v["spectral_center"], conf_s.numerical_v["spectral_range"], conf_s.numerical_v["spectral_power_amplitude"], conf_s.numerical_v["spectral_number_points"] );
            } else {
                conf_s.numerical_v["spectral_range"] = {};
                conf_s.numerical_v["spectral_center"] = {};
                conf_s.numerical_v["spectral_power_amplitude"] = {};
                conf_s.numerical_v["spectral_number_points"] = {};
            }
        }
        input_conf["Detector"] = conf_s;
    }
    {
        input_s conf_s;
        LOG2( "[System-Parameters] Setting up density matrix output using {}...\n", inputstring_densitymatrix_config );
        auto conf = QDLC::String::splitline( inputstring_densitymatrix_config, ':' );
        conf_s.string["output_mode"] = conf[0];
        conf_s.string["interaction_picture"] = conf.size() > 1 ? conf[1] : "schroedinger";
        input_conf["DMconfig"] = conf_s;
    }
    // RK 45 Tolerance Vector
    {
        if ( std::find( inputstring_rk45_config.begin(), inputstring_rk45_config.end(), ':' ) == inputstring_rk45_config.end() ) {
            numerics_rk_tol.emplace_back( 1.0, QDLC::Misc::convertParam<double>( inputstring_rk45_config ) );
            LOG2( "[System-Parameters] Set fixed RK45 Tolerance to {}.\n", std::get<1>( numerics_rk_tol.back() ) );
        } else {
            LOG2( "[System-Parameters] Setting up multiple tolerances...\n" );
            for ( auto partial : QDLC::String::splitline( inputstring_rk45_config, ';' ) ) {
                auto tuple = QDLC::String::splitline( partial, ':' );
                numerics_rk_tol.emplace_back( std::make_tuple( QDLC::Misc::convertParam<double>( tuple.front() ), QDLC::Misc::convertParam<double>( tuple.back() ) ) );
                LOG2( "[System-Parameters] Set RK45 Tolerance to {} (from {}).\n", std::get<1>( numerics_rk_tol.back() ), partial );
            }
        }
    }
}

void Parameters::log( const Dense &initial_state_vector_ket ) {
    // LOG( "                                                                                                                                                      \n                                                                                                                  #################                   \n                                                                                                              ..  :::::::::::::::::  ..               \n                                                                                                            ...       -      :+:      ...             \n                                                                                                           ...       .+:     :+:       ...            \n                                                                                                          ...        =+=     :+:        ...           \n                                                                                                         ....       :+++-    :+:         ...          \n         .:----.       ::::::.        :-:          .:---:            .:---:              :---:.          ...        :-+-:    :+:         ...          \n       +@@@@@@@@@*.   =@@@@@@@@@*:   .@@@.       =%@@@@@@@=         %@@%@@@@+          =@@@%@@@*        ....         .+:     :+:          ...         \n     .%@@*:   :*@@@.  =@@#   :+@@@+  .@@@.      %@@%-   :+-         =:   =@@@.        -@@#.  *@@#       ....         .+:     :+:          ...         \n     *@@#       #@@*  =@@#     :@@@: .@@@.     *@@%.                   .:#@@#         %@@=   .@@@.      ...          .+:     :+:          ...         \n     %@@+       +@@#  =@@#      @@@= .@@@.     %@@*         .....    @@@@@@*:         @@@-    @@@:      ...          .+:     :+:          ...         \n     #@@*       #@@*  =@@#     .@@@: .@@@.     #@@#        :@@@@@%   ...:+@@@+        %@@-   .@@@.      ....         .+:     :+:          ...         \n     -@@@-     +@@@:  =@@#    :%@@*  .@@@.     -@@@+     :: ......        %@@#   :-.  *@@*   +@@#       ....         .+:     :+:          ...         \n      -%@@@##%@@@@-   =@@@%%%@@@%=   .@@@%%%%%- :%@@@%#%@@+        =@%#*#@@@%:  *@@%   #@@%*%@@#.        ...         .+:    .-+:.        ....         \n        :=+***==%@@#+ .=++++==-.      =+++++++:   :=+**+=:          -=+**+=:    .++-    :+***=:          ...         .+:    =+++-        ...          \n                 :+#%.                                                                                    ...        .+:    .++=        ...           \n                                                                                                           ...       .+:     -+:        ...           \n      -##########################################################################################           ...      .+:      =        ..             \n      .::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::             ..  ::::-::::::::::::  ...              \n                                                                                                                  #################                   \n                                                                                                                                                      \n                                                                                                                                                      \n" );
    // LOG( "╭━━━┳━━━┳╮╱╱╭━━━╮╱╱╭━┳━━╮╭━━┳━╮╱╱╭╮╱╱╱╱╱╱╭━━━╮╭━━╮\n┃╭━╮┣╮╭╮┃┃╱╱┃╭━╮┃╱╭╯╭┻┫┣╯╰┫┣┻╮╰╮╱┃┃╱╱╱╱╱╱╰╮╭╮┃┃╭╮┃\n┃┃╱┃┃┃┃┃┃┃╱╱┃┃╱╰╯╭╯╭╯╱┃┃╱╱┃┃╱╰╮╰╮┃╰━┳╮╱╭╮╱┃┃┃┃┃╰╯╰╮\n┃┃╱┃┃┃┃┃┃┃╱╭┫┃╱╭╮┃┃┃╱╱┃┃╱╱┃┃╱╱┃┃┃┃╭╮┃┃╱┃┃╱┃┃┃┃┃╭━╮┃\n┃╰━╯┣╯╰╯┃╰━╯┃╰━╯┃┃┃┃╱╭┫┣╮╭┫┣╮╱┃┃┃┃╰╯┃╰━╯┃╭╯╰╯┣┫╰━╯┣╮\n╰━━╮┣━━━┻━━━┻━━━╯╰╮╰╮╰━━╯╰━━╯╭╯╭╯╰━━┻━╮╭╯╰━━━┻┻━━━┻╯\n╱╱╱╰╯╱╱╱╱╱╱╱╱╱╱╱╱╱╰╮╰╮╱╱╱╱╱╱╭╯╭╯╱╱╱╱╭━╯┃\n╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╰━╯╱╱╱╱╱╱╰━╯╱╱╱╱╱╰━━╯\n" );
    Log::P1( "\033[38;2;128;128;128m          \033[38;2;87;107;143m▄\033[38;2;70;99;149m▄\033[38;2;63;95;152m▓                                                                                                           \033[38;2;63;95;152m▓\033[38;2;71;99;148m▄\033[38;2;89;108;142m▄\033[0m\n\033[38;2;128;128;128m       \033[38;2;92;110;141m▄\033[38;2;45;86;158m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒                                                                                                           ▒▒▒▒\033[38;2;49;88;157m▓\033[38;2;98;112;139m▄         \033[38;2;96;133;70m▌\033[0m\n\033[38;2;128;128;128m      \033[38;2;50;89;156m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;38;83;161m▓                                                                                                           \033[38;2;37;82;161m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;56;91;154m▓        \033[38;2;81;136;43m▓\033[38;2;109;131;94m▄                        \033[38;2;233;181;21m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;235;181;19m▒\033[0m\n\033[38;2;128;128;128m     \033[38;2;37;82;161m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;69;39;94m▓\033[38;2;172;7;15m▓\033[38;2;190;5;5m█\033[38;2;180;23;23m█\033[38;2;171;42;41m▓\033[38;2;162;59;59m▄\033[38;2;153;77;77m▄\033[38;2;144;95;95m▄    \033[38;2;118;102;129m▄\033[38;2;113;90;130m▄\033[38;2;110;82;131m▄\033[38;2;108;75;131m▄\033[38;2;107;75;131m▄\033[38;2;107;75;131m▄\033[38;2;108;77;131m▄\033[38;2;111;85;131m▄\033[38;2;115;95;130m▄                                                                      \033[38;2;118;103;129m▄\033[38;2;114;91;130m▄\033[38;2;111;83;131m▄\033[38;2;108;76;131m▄\033[38;2;107;75;131m▄\033[38;2;107;75;131m▄\033[38;2;108;76;131m▄\033[38;2;111;85;131m▄\033[38;2;115;94;130m▄     \033[38;2;146;92;91m▄\033[38;2;155;73;73m▄\033[38;2;164;55;55m▄\033[38;2;173;37;37m▓\033[38;2;183;19;19m█\033[38;2;191;2;2m▓\033[38;2;165;9;21m▓\033[38;2;55;43;104m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;40;84;160m▓       \033[38;2;81;136;43m▓\033[38;2;78;136;37m▓                        \033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒    \033[38;2;186;157;69m▄\033[38;2;187;157;68m▄\033[0m\n\033[38;2;128;128;128m    \033[38;2;37;82;161m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;96;31;73m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓\033[38;2;182;3;13m▓\033[38;2;147;11;59m▓\033[38;2;116;18;99m▒\033[38;2;90;25;132m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓\033[38;2;92;35;134m▓\033[38;2;103;62;132m▓\033[38;2;115;95;130m▄          \033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;89;26;135m▓\033[38;2;93;38;134m▓\033[38;2;100;55;133m▓\033[38;2;107;75;131m▄\033[38;2;116;99;130m▄          \033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;96;46;134m▌           \033[38;2;110;81;131m▄\033[38;2;98;51;133m▓\033[38;2;89;27;135m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓\033[38;2;95;24;126m▓\033[38;2;123;16;89m▓\033[38;2;156;8;47m▓\033[38;2;189;1;3m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓\033[38;2;75;37;89m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;41;84;160m▓      \033[38;2;81;136;43m▓\033[38;2;78;137;36m▒\033[38;2;86;135;51m▓                      \033[38;2;193;159;62m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒   \033[38;2;230;179;24m▒▒▒▒\033[0m\n\033[38;2;128;128;128m   \033[38;2;49;88;156m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;49;45;109m▒\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓\033[38;2;159;8;43m▓\033[38;2;110;19;107m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;89;27;135m▓\033[38;2;106;70;132m▄        \033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;88;24;135m▓\033[38;2;98;51;133m▓\033[38;2;115;96;130m▄       \033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;96;46;134m▌        \033[38;2;116;98;130m▄\033[38;2;97;47;133m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;90;24;133m▒\033[38;2;129;15;82m▓\033[38;2;182;2;13m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓\033[38;2;192;0;0m▓\033[38;2;41;59;131m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;59;93;153m▌     \033[38;2;81;136;43m▓\033[38;2;78;137;36m▒\033[38;2;78;137;36m▒\033[38;2;101;132;79m▌                     \033[38;2;229;178;26m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒\033[38;2;207;166;48m▌  ▒▒▒▒\033[38;2;204;165;51m▌\033[0m\n\033[38;2;128;128;128m  \033[38;2;98;112;139m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;167;8;19m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;175;4;22m▓\033[38;2;105;20;113m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;90;25;133m▓\033[38;2;144;21;66m▓\033[38;2;183;19;18m█\033[38;2;179;27;26m█\033[38;2;175;34;33m▓\033[38;2;171;42;41m▓\033[38;2;167;50;49m▄\033[38;2;157;55;63m▄\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;93;38;134m▓      ▓▓▓▓▓▓▓▓\033[38;2;102;34;122m▓\033[38;2;159;66;65m▄\033[38;2;164;57;56m▄\033[38;2;167;49;49m▓\033[38;2;171;41;41m▓\033[38;2;176;33;32m▓\033[38;2;180;25;24m▓\033[38;2;171;19;31m█\033[38;2;105;26;114m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;89;24;134m▓\033[38;2;146;11;60m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓\033[38;2;147;14;35m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒     \033[38;2;81;136;43m▓\033[38;2;78;137;36m▒\033[38;2;78;137;36m▒\033[38;2;77;137;36m▒                     \033[38;2;253;191;1m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒\033[38;2;232;179;22m▒ \033[38;2;169;147;86m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒\033[38;2;246;187;8m▒\033[0m\n\033[38;2;128;128;128m  \033[38;2;47;87;157m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;42;57;128m▒\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;164;7;36m▓\033[38;2;89;24;134m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;105;20;113m▓\033[38;2;132;14;78m▓\033[38;2;147;10;58m▓\033[38;2;151;9;53m▓\033[38;2;150;10;55m▓\033[38;2;137;13;71m▓\033[38;2;115;18;100m▓\033[38;2;89;24;134m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;129;15;81m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;183;3;12m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;108;19;109m▓\033[38;2;128;15;83m▓\033[38;2;128;15;83m▓\033[38;2;128;15;83m▓\033[38;2;125;16;87m▓\033[38;2;119;17;95m▓\033[38;2;108;19;109m▓\033[38;2;92;23;130m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;88;24;135m▓\033[38;2;149;18;58m▓\033[38;2;189;6;5m▓\033[38;2;191;4;3m▓\033[38;2;191;3;2m▓\033[38;2;192;1;0m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;110;19;108m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓\033[38;2;172;5;26m▓\033[38;2;90;24;132m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;109;19;108m▓\033[38;2;134;14;75m▓\033[38;2;148;10;57m▓\033[38;2;151;9;53m▓\033[38;2;150;10;55m▓\033[38;2;139;12;69m▓\033[38;2;118;17;96m▓\033[38;2;90;24;133m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;131;14;80m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓\033[38;2;36;74;152m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;59;93;153m▌    \033[38;2;81;136;43m▓\033[38;2;78;137;36m▒\033[38;2;78;137;36m▒▒\033[38;2;81;136;42m▓              \033[38;2;165;146;90m▄\033[38;2;182;155;72m▄\033[38;2;163;145;92m▄   \033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒\033[38;2;251;190;3m▒ \033[38;2;202;164;53m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒\033[0m\n\033[38;2;128;128;128m  \033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;97;30;73m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;186;1;7m▓\033[38;2;90;24;133m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;118;17;96m▓\033[38;2;178;3;18m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓\033[38;2;189;1;4m▓\033[38;2;140;12;67m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;158;8;44m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;183;3;12m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;140;12;67m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓\033[38;2;178;3;18m▓\033[38;2;132;14;77m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;88;24;135m▓\033[38;2;181;3;15m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;110;19;108m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;191;0;2m▓\033[38;2;93;23;128m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;123;16;90m▓\033[38;2;184;2;11m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓\033[38;2;189;1;4m▓\033[38;2;133;14;77m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;156;9;47m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;73;37;91m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;34;81;162m▒ \033[38;2;78;137;36m▒\033[38;2;78;137;36m▒▒▒▒▒▒▒\033[38;2;93;134;64m▌            \033[38;2;230;179;25m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒\033[38;2;211;169;44m▌ \033[38;2;165;145;90m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒▒ \033[38;2;232;180;22m▐▒▒▒▒▒  \033[38;2;201;164;53m▄\033[38;2;242;186;12m▒\033[38;2;242;185;12m▒\033[38;2;195;161;60m▄ \033[38;2;196;161;59m▄\033[38;2;171;149;84m▄\033[0m\n\033[38;2;128;128;128m  \033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;141;16;39m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;145;11;61m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;120;17;93m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓▓▓▓\033[38;2;160;8;41m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;106;20;112m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;183;3;12m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;140;12;67m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓▓\033[38;2;141;12;66m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;138;13;71m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;110;19;108m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;153;9;50m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;118;17;96m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;156;0;0m▓\033[38;2;73;0;0m█\033[38;2;43;0;0m█\033[38;2;50;0;0m█\033[38;2;95;0;0m█\033[38;2;184;0;0m▒\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓\033[38;2;151;10;54m▓\033[38;2;139;12;69m▓\033[38;2;141;12;66m▓\033[38;2;144;11;62m▓\033[38;2;150;10;55m▓\033[38;2;152;9;52m▓\033[38;2;155;9;48m▓\033[38;2;160;7;41m▓\033[38;2;162;7;39m▓\033[38;2;169;6;30m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;109;26;63m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒ \033[38;2;78;137;36m▒\033[38;2;78;137;36m▒▒▒▒▒▒▒▒           \033[38;2;203;165;52m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒ \033[38;2;195;160;60m▐▒▒\033[38;2;227;177;27m▌\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒ \033[38;2;254;192;0m▒▒▒\033[38;2;250;189;4m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;197;162;58m▌\033[38;2;185;156;70m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒\033[38;2;195;160;60m▄\033[38;2;245;186;9m░\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;247;188;7m▒\033[38;2;216;172;39m▒\033[38;2;184;156;70m▄\033[0m\n\033[38;2;128;128;128m \033[38;2;92;110;141m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;151;13;32m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;133;14;78m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;158;8;44m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓▓▓▓\033[38;2;192;0;1m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;92;23;130m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;183;3;12m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;140;12;67m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓▓\033[38;2;166;6;33m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;127;16;83m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;110;19;108m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;140;12;68m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;153;9;49m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;0;0;0m█\033[38;2;0;0;0m█\033[38;2;0;0;0m███\033[38;2;64;1;0m█\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓\033[38;2;182;2;13m▓\033[38;2;145;11;61m▓\033[38;2;146;11;59m▓\033[38;2;148;10;58m▓\033[38;2;149;10;56m▓\033[38;2;150;10;54m▓\033[38;2;152;9;52m▓\033[38;2;153;9;50m▓\033[38;2;154;9;48m▓\033[38;2;158;8;43m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;140;16;40m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒ \033[38;2;78;137;36m▒\033[38;2;78;137;36m▒▒▒▒▒▒▒▒\033[38;2;79;136;38m▓ \033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒▒▒\033[38;2;251;190;3m▒\033[38;2;195;161;59m▄\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;238;183;16m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;193;160;61m▌\033[38;2;225;176;30m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒ ▒▒▒ ▒▒▒\033[38;2;198;161;56m▐▒▒\033[38;2;232;180;23m▒\033[38;2;248;189;6m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒▒▒▒▒▒▒▒\033[38;2;253;192;1m▒\033[38;2;218;172;36m▒\033[0m\n\033[38;2;128;128;128m \033[38;2;92;110;141m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;151;13;32m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;156;9;47m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;111;19;105m▓\033[38;2;192;0;0m▓\033[38;2;186;2;8m▓\033[38;2;131;14;80m▓\033[38;2;131;14;80m▓▓▓▓▓▓\033[38;2;161;7;40m▓\033[38;2;192;0;0m▓\033[38;2;145;11;61m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;123;16;90m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;183;3;12m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;140;12;67m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓\033[38;2;191;0;1m▓\033[38;2;113;18;103m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;150;10;54m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;109;19;108m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;157;8;45m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;110;19;107m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓\033[38;2;138;1;0m▓\033[38;2;53;0;0m█\033[38;2;23;0;0m█\033[38;2;30;0;0m█\033[38;2;75;0;0m█\033[38;2;175;0;0m▒\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓\033[38;2;123;16;90m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;130;15;80m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;140;16;40m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒ \033[38;2;78;137;36m▒\033[38;2;78;137;36m▒▒▒▒▒▒▒▒        \033[38;2;236;182;19m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒ \033[38;2;180;153;75m▐▒▒\033[38;2;243;186;11m▒\033[38;2;251;190;3m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒ \033[38;2;240;184;14m▒▒▒\033[38;2;212;168;42m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒ ▒▒▒▒▒\033[38;2;251;190;3m▒ \033[38;2;202;165;52m▀\033[38;2;210;169;44m▀\033[38;2;248;188;6m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒\033[38;2;232;180;23m▀\033[38;2;200;164;55m▀\033[0m\n\033[38;2;128;128;128m  \033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;141;16;39m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;191;0;1m▓\033[38;2;98;22;121m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;107;20;111m▒\033[38;2;157;9;45m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓\033[38;2;130;14;81m▓\033[38;2;125;16;86m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;89;24;134m▓\033[38;2;181;3;14m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;183;3;12m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;140;12;67m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓\033[38;2;191;0;1m▓\033[38;2;174;4;23m▓\033[38;2;144;11;62m▓\033[38;2;99;22;121m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;96;22;125m▓\033[38;2;191;0;2m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;88;24;134m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓\033[38;2;92;23;130m▓\033[38;2;182;3;13m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓\033[38;2;192;0;0m▓\033[38;2;97;22;123m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;107;20;111m▒\033[38;2;166;6;34m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓\033[38;2;177;4;20m▓\033[38;2;119;17;94m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;92;23;130m▓\033[38;2;186;1;7m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;109;26;64m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒ \033[38;2;78;137;36m▒\033[38;2;78;137;36m▒▒▒▒▒▒▒\033[38;2;89;134;57m▌          \033[38;2;202;165;52m▀\033[38;2;201;164;54m▀   \033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒▒ \033[38;2;221;174;33m▐▒▒\033[38;2;251;190;3m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;236;182;18m▒ ▒▒▒▒▒    \033[38;2;209;167;45m▀\033[0m\n\033[38;2;128;128;128m  \033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;96;30;73m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;181;3;15m▓\033[38;2;97;22;123m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;90;24;133m▓\033[38;2;165;6;35m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;183;3;12m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;140;12;67m▓\033[38;2;187;2;7m▓\033[38;2;99;22;121m▒\033[38;2;99;22;121m▒\033[38;2;97;22;123m▒\033[38;2;89;24;134m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓\033[38;2;96;22;125m▓\033[38;2;178;3;18m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓\033[38;2;115;18;100m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;92;23;129m▒\033[38;2;99;22;121m▒\033[38;2;162;7;38m▓\033[38;2;192;0;0m▓\033[38;2;179;3;17m▓\033[38;2;95;22;126m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;91;23;131m▒\033[38;2;113;18;103m▒\033[38;2;125;15;87m▒\033[38;2;130;14;81m▒\033[38;2;128;15;84m▒\033[38;2;116;18;99m▒\033[38;2;96;22;124m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;90;24;132m▓\033[38;2;170;5;28m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓\033[38;2;73;37;91m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;34;81;162m▒    \033[38;2;80;136;42m▓\033[38;2;78;137;36m▒\033[38;2;78;137;36m▒▒\033[38;2;80;136;39m▓                \033[38;2;234;180;21m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒\033[38;2;252;191;2m▒ \033[38;2;196;162;59m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒\033[38;2;210;168;45m▌ \033[38;2;227;176;28m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒\033[38;2;253;191;1m▒\033[0m\n\033[38;2;128;128;128m  \033[38;2;47;87;157m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;42;57;128m▒\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓\033[38;2;188;1;5m▓\033[38;2;126;15;86m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;109;19;108m▓\033[38;2;179;3;17m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓\033[38;2;183;3;12m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;140;12;67m▓\033[38;2;186;2;7m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;88;24;135m▓\033[38;2;119;32;100m▓\033[38;2;177;21;25m▀\033[38;2;181;23;22m▀\033[38;2;183;19;18m▀\033[38;2;184;16;15m█\033[38;2;186;13;13m█\033[38;2;188;8;7m█\033[38;2;180;6;15m▓\033[38;2;93;23;128m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓\033[38;2;158;8;43m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓\033[38;2;186;1;7m▓\033[38;2;122;16;91m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;104;20;114m▓\033[38;2;178;3;19m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓\033[38;2;36;74;152m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;59;93;153m▌    \033[38;2;81;136;43m▓\033[38;2;78;137;36m▒\033[38;2;78;137;36m▒▒                 \033[38;2;176;151;79m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒\033[38;2;213;170;41m▌ \033[38;2;174;150;80m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒   \033[38;2;246;187;8m▒▒\033[38;2;250;189;4m▒\033[0m\n\033[38;2;128;128;128m  \033[38;2;98;112;139m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;167;8;19m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓\033[38;2;178;3;18m▓\033[38;2;128;15;83m▓\033[38;2;91;23;131m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;88;24;135m▓\033[38;2;114;18;101m▓\033[38;2;160;13;43m▓\033[38;2;185;15;14m█\033[38;2;181;23;22m▀\033[38;2;177;31;30m▀\033[38;2;172;40;39m▀\033[38;2;168;48;47m▀\033[38;2;164;56;55m▀\033[38;2;159;65;65m▀ \033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;110;71;126m▌ ▓▓▓▓▓▓▓▓▓▓▓\033[38;2;91;33;134m▓\033[38;2;106;72;132m▀          \033[38;2;102;60;132m▀\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;133;59;94m▒ \033[38;2;163;58;58m▀\033[38;2;167;49;49m▀\033[38;2;172;41;40m▀\033[38;2;159;35;52m▀\033[38;2;114;37;108m▀\033[38;2;89;26;135m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\033[38;2;106;20;112m▓\033[38;2;156;8;47m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓\033[38;2;147;14;35m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒     \033[38;2;81;136;43m▓\033[38;2;78;137;36m▒\033[38;2;78;137;36m▒\033[38;2;96;133;71m▌                  \033[38;2;244;186;10m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒   ▒▒▒▒▒\033[0m\n\033[38;2;128;128;128m   \033[38;2;49;88;156m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;49;45;109m▒\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓\033[38;2;191;0;1m▓\033[38;2;166;6;34m▓\033[38;2;137;13;72m▓\033[38;2;112;19;103m▓\033[38;2;93;23;129m▒\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓\033[38;2;89;26;135m▓\033[38;2;95;42;134m▀\033[38;2;104;65;132m▀            \033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓\033[38;2;108;76;131m▌ ▓▓▓▓▓\033[38;2;88;24;135m▓\033[38;2;92;34;134m▓\033[38;2;97;47;133m▀\033[38;2;103;65;132m▀                 \033[38;2;103;63;132m▀\033[38;2;96;44;134m▀\033[38;2;91;32;134m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓         \033[38;2;105;68;132m▀\033[38;2;95;44;134m▀\033[38;2;89;26;135m▓\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓▓▓▓\033[38;2;96;22;124m▓\033[38;2;118;17;96m▓\033[38;2;145;11;60m▓\033[38;2;179;3;17m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓\033[38;2;192;0;0m▓\033[38;2;41;59;131m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;59;93;153m▌     \033[38;2;81;136;43m▓\033[38;2;78;137;36m▒\033[38;2;83;135;46m▌                     \033[38;2;199;164;55m▀    \033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒▒▒\033[0m\n\033[38;2;128;128;128m    \033[38;2;37;82;161m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;96;31;73m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓\033[38;2;187;11;11m█\033[38;2;173;39;38m▀\033[38;2;160;64;63m▀ \033[38;2;118;58;112m▐\033[38;2;88;24;135m▓\033[38;2;88;24;135m▓▓▓▓▓▓\033[38;2;95;43;134m▌                                                                                \033[38;2;160;64;64m▀\033[38;2;171;41;41m▀\033[38;2;183;18;18m█\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓\033[38;2;192;0;0m▓▓▓▓▓▓▓\033[38;2;75;37;89m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;41;84;160m▓      \033[38;2;81;136;43m▓\033[38;2;78;137;36m▒                           \033[38;2;250;190;4m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒\033[38;2;241;184;13m▒\033[0m\n\033[38;2;128;128;128m     \033[38;2;37;82;161m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;69;39;94m▒\033[38;2;172;7;15m▓\033[38;2;192;0;0m▓\033[38;2;184;17;16m█\033[38;2;170;43;43m▀                                                                                                      \033[38;2;158;68;67m▀\033[38;2;170;44;44m▀\033[38;2;181;21;21m▀\033[38;2;191;1;1m▓\033[38;2;165;9;21m▓\033[38;2;55;43;104m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;40;84;160m▓       \033[38;2;81;136;43m▓                            \033[38;2;212;170;42m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒▒\033[0m\n\033[38;2;128;128;128m      \033[38;2;50;89;156m▀\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;38;83;161m▓   \033[38;2;102;60;132m▓\033[38;2;99;54;133m▀ \033[38;2;108;77;131m▀\033[38;2;92;35;134m▓   \033[38;2;116;98;130m▄ \033[38;2;114;92;130m▄   \033[38;2;112;88;130m▄\033[38;2;112;86;131m▄   \033[38;2;113;90;130m▄\033[38;2;118;102;129m▄\033[38;2;112;86;131m▄\033[38;2;118;104;129m▄  \033[38;2;94;40;134m▓\033[38;2;107;73;132m▄  \033[38;2;115;96;130m▄  \033[38;2;112;87;130m▄  \033[38;2;116;99;130m▄\033[38;2;118;102;129m▄\033[38;2;112;86;131m▄\033[38;2;117;99;130m▄\033[38;2;115;96;130m▄\033[38;2;112;86;131m▄      \033[38;2;89;27;135m▓    \033[38;2;108;76;131m▐   \033[38;2;114;93;130m▄\033[38;2;114;91;130m▄\033[38;2;112;88;130m▄ \033[38;2;112;86;131m▄\033[38;2;116;97;130m▄  \033[38;2;113;89;130m▄  \033[38;2;112;87;130m▄  \033[38;2;94;39;134m▓   \033[38;2;113;89;130m▄\033[38;2;112;86;131m▄\033[38;2;118;102;129m▄  \033[38;2;96;44;134m▓\033[38;2;105;68;132m▄  \033[38;2;109;80;131m▐\033[38;2;110;81;131m▄  \033[38;2;116;98;130m▄\033[38;2;112;86;131m▄\033[38;2;113;91;130m▄   \033[38;2;114;93;130m▄\033[38;2;114;91;130m▄\033[38;2;112;87;130m▄   \033[38;2;113;89;130m▄\033[38;2;112;86;130m▄   \033[38;2;37;82;161m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒▒\033[38;2;56;91;154m▀        \033[38;2;92;134;64m▌                             \033[38;2;232;180;22m▀\033[38;2;251;190;3m▒\033[38;2;223;176;31m▀\033[0m\n\033[38;2;128;128;128m        \033[38;2;45;86;158m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒   \033[38;2;88;24;135m▓   \033[38;2;88;24;135m▓  \033[38;2;109;78;131m▐\033[38;2;98;49;133m▌ \033[38;2;92;34;134m▓  \033[38;2;115;95;130m▄\033[38;2;108;77;131m▄\033[38;2;107;74;131m▐\033[38;2;88;24;135m▓  \033[38;2;90;30;135m▓  \033[38;2;88;24;135m▓  \033[38;2;91;33;134m▓   \033[38;2;96;44;134m▓  \033[38;2;88;24;135m▓  \033[38;2;99;52;133m▐\033[38;2;104;65;132m▌ \033[38;2;88;24;135m▓ \033[38;2;105;69;132m▐\033[38;2;99;53;133m▌      \033[38;2;106;71;132m▀\033[38;2;91;33;134m▓  \033[38;2;101;58;133m▐\033[38;2;106;71;132m▌  \033[38;2;88;24;135m▓ \033[38;2;98;49;133m▐\033[38;2;103;63;132m▌ \033[38;2;88;24;135m▓  \033[38;2;89;27;135m▓  ▓  \033[38;2;94;39;134m▓   \033[38;2;107;74;131m▄\033[38;2;111;83;131m▀\033[38;2;88;24;135m▓  \033[38;2;95;42;134m▓   \033[38;2;103;63;132m▐\033[38;2;104;66;132m▌ \033[38;2;109;79;131m▐\033[38;2;96;44;134m▌ \033[38;2;105;67;132m▐\033[38;2;98;50;133m▌  \033[38;2;88;24;135m▓ \033[38;2;100;56;133m▐\033[38;2;105;69;132m▌  \033[38;2;100;55;133m▀\033[38;2;105;68;132m▄   \033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒▒\033[38;2;49;88;157m▀\033[0m\n\033[38;2;128;128;128m           \033[38;2;70;99;149m▀\033[38;2;63;95;152m▀    \033[38;2;105;69;132m▀\033[38;2;107;73;132m▀\033[38;2;105;69;132m▀\033[38;2;95;42;134m▀\033[38;2;109;78;131m▄  \033[38;2;103;64;132m▀\033[38;2;107;75;131m▀\033[38;2;106;72;132m▀   \033[38;2;106;71;132m▀\033[38;2;110;82;131m▀\033[38;2;104;64;132m▀  \033[38;2;105;68;132m▀  \033[38;2;103;64;132m▀   \033[38;2;105;68;132m▀   \033[38;2;103;64;132m▀\033[38;2;108;77;131m▀\033[38;2;104;65;132m▀  \033[38;2;110;82;131m▀  \033[38;2;103;64;132m▀  \033[38;2;110;82;131m▀     \033[38;2;106;72;132m▀\033[38;2;107;74;131m▀\033[38;2;107;74;131m▀      \033[38;2;104;65;132m▀    \033[38;2;103;64;132m▀   \033[38;2;104;65;132m▀\033[38;2;109;78;131m▀\033[38;2;103;65;132m▀  \033[38;2;107;74;131m▀   \033[38;2;105;70;132m▀\033[38;2;110;82;131m▀\033[38;2;103;64;132m▀   \033[38;2;105;67;132m▀      \033[38;2;106;71;132m▀\033[38;2;108;77;131m▀\033[38;2;106;70;132m▀   \033[38;2;104;65;132m▀     \033[38;2;110;81;131m▀\033[38;2;106;70;132m▀   \033[38;2;63;95;151m▀\033[38;2;71;99;148m▀\033[0m\n" );
    Log::wrapInBar( "System Parameters" );
    LOG( "Version: {} ({})\n\n", GLOBAL_PROGRAM_VERSION, GLOBAL_PROGRAM_LASTCHANGE );

    Log::wrapInBar( "Electronic Configuration", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    for ( auto &[name, mat] : input_electronic ) {
        LOG( "Electronic State: {} with energy {:.8e} Hz - {:.8e} eV - {:.8f} nm\n", name, mat.numerical["Energy"], mat.numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ), mat.numerical["Energy"].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
        for ( auto i = 0; i < mat.string_v["CoupledTo"].size(); i++ ) {
            if ( mat.string_v["CoupledTo"][i].front() == '-' ) break;
            LOG( " - Coupled to {} with transition energy {:.8e} Hz - {:.8e} eV\n", mat.string_v["CoupledTo"][i], std::abs( mat.numerical["Energy"] - input_electronic[mat.string_v["CoupledTo"][i]].numerical["Energy"] ), std::abs( mat.numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ) - input_electronic[mat.string_v["CoupledTo"][i]].numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ) ) );
        }
        LOG( " - Decay Scaling: {}\n", mat.numerical["DecayScaling"] );
        LOG( " - Dephasing Scaling: {}\n", mat.numerical["DephasingScaling"] );
    }
    LOG( "\n" );
    if ( input_electronic.size() > 0 ) {
        Log::wrapInBar( "Photonic Configuration", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
        for ( auto &[name, mat] : input_photonic ) {
            LOG( "Cavity: {} with energy {:.8e} Hz - {:.8e} eV - {:.8f} nm\n", name, mat.numerical["Energy"], mat.numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ), mat.numerical["Energy"].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
            for ( auto i = 0; i < mat.string_v["CoupledTo"].size(); i++ ) {
                auto [from, to] = QDLC::String::split_pair( mat.string_v["CoupledTo"][i], transition_delimiter );
                double purcell = ( 2.0 * p_omega_coupling * p_omega_coupling * mat.numerical_v["CouplingScaling"][i] * mat.numerical_v["CouplingScaling"][i] / p_omega_cavity_loss / p_omega_decay ) * ( p_omega_cavity_loss * p_omega_cavity_loss / ( std::pow( mat.numerical["Energy"] - ( input_electronic[to].numerical["Energy"] - input_electronic[from].numerical["Energy"] ), 2 ) + p_omega_cavity_loss * p_omega_cavity_loss ) );
                LOG( " - Coupled to electronic transition {} (Purcell Enhancement: {}) \n", mat.string_v["CoupledTo"][i], purcell );
                LOG( " - - Coupling Scaling: {}\n", mat.numerical_v["CouplingScaling"][i] );
            }
            LOG( " - Cavity Q-Factor: {:.0f}\n", mat.numerical["Energy"] / p_omega_cavity_loss );
            LOG( " - Decay Scaling: {}\n", mat.numerical["DecayScaling"] );
        }
        LOG( "\n" );
    }

    Log::wrapInBar( "Coupling Parameters", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    LOG( "Coupling strengh g: {:.8e} Hz - {:.8} mueV\n", p_omega_coupling, p_omega_coupling.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    LOG( "Photon loss rate k: {:.8e} Hz - {:.8} mueV\n", p_omega_cavity_loss, p_omega_cavity_loss.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    LOG( "Atomic dephasing rate gamma_pure: {:.8e} Hz - {:.8} mueV\n", p_omega_pure_dephasing, p_omega_pure_dephasing.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    LOG( "RAD rate gamma: {:.8e} Hz - {:.8} mueV\n\n", p_omega_decay, p_omega_decay.getSI( Parameter::UNIT_ENERGY_MUEV ) );

    Log::wrapInBar( "Initial System Parameters", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    LOG( "Initial state rho0 = [{}]\n\n", initial_state_vector_ket.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ) );

    Log::wrapInBar( "Pulse", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_pulse.size() > 0 ) {
        for ( auto &[name, mat] : input_pulse ) {
            LOG( "Pulse {}:\n", name );
            LOG( " - Coupled to Transitions: " );
            for ( auto &mode : mat.string_v["CoupledTo"] )
                LOG( "{} ", mode );
            LOG( "\n" );
            for ( int i = 0; i < mat.numerical_v["Amplitude"].size(); i++ ) {
                LOG( " - Single Pulse:\n" );
                LOG( " - - Amplitude: {} Hz - {:.8} mueV\n", mat.numerical_v["Amplitude"][i], mat.numerical_v["Amplitude"][i].getSI( Parameter::UNIT_ENERGY_MUEV ) );
                LOG( " - - Frequency: {} Hz - {:.8} eV - {:.8} nm\n", mat.numerical_v["Frequency"][i], mat.numerical_v["Frequency"][i].getSI( Parameter::UNIT_ENERGY_EV ), mat.numerical_v["Frequency"][i].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
                LOG( " - - Width: {} s - {:.8} ps\n", mat.numerical_v["Width"][i], mat.numerical_v["Width"][i].getSI( Parameter::UNIT_TIME_PS ) );
                LOG( " - - Center: {} s - {:.8} ps\n", mat.numerical_v["Center"][i], mat.numerical_v["Center"][i].getSI( Parameter::UNIT_TIME_PS ) );
                if ( QDLC::Math::abs2( mat.numerical_v["ChirpRate"][i] != 0.0 ) )
                    LOG( " - - Chirp: {}\n", mat.numerical_v["ChirpRate"][i] );
                if ( QDLC::Math::abs2( mat.numerical_v["SUPERDelta"][i] != 0.0 ) ) {
                    LOG( " - - SUPER Amplitude: {} - {} meV\n", mat.numerical_v["SUPERDelta"][i], mat.numerical_v["SUPERDelta"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                    LOG( " - - SUPER Frequency: {} - {} meV\n", mat.numerical_v["SUPERFreq"][i], mat.numerical_v["SUPERFreq"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                }
                LOG( " - - Type: {}{}\n", mat.string_v["Type"][i], mat.string_v["Type"][i] == "gauss" ? fmt::format( " (Gaussian Amplitude: {})", mat.numerical_v["GaussAmp"][i] ) : "" );
            }
        }
        LOG( "\n" );
    } else {
        LOG( "Not using any pulses to exite or drive the system.\n\n" );
    }
    Log::wrapInBar( "Chirp", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_chirp.size() > 0 ) {
        for ( auto &[name, mat] : input_chirp ) {
            LOG( "Chirp: {}\n", name );
            LOG( " - Coupled to States:\n" );
            for ( int i = 0; i < mat.string_v["CoupledTo"].size(); i++ )
                LOG( " - - {} with scaling {}\n", mat.string_v["CoupledTo"][i], mat.numerical_v["AmpFactor"][i] );
            for ( int i = 0; i < mat.numerical_v["Amplitude"].size(); i++ ) {
                LOG( " - Chirp Point {}:\n", i );
                LOG( " - - Amplitude: {} Hz - {} meV\n", mat.numerical_v["Amplitude"][i], mat.numerical_v["Amplitude"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                LOG( " - - Time: {} s - {} ps\n", mat.numerical_v["Times"][i], mat.numerical_v["Times"][i].getSI( Parameter::UNIT_TIME_PS ) );
                LOG( " - - Derivative DDT: {}\n", mat.numerical_v["ddt"][i] );
            }
        }
        LOG( "\n" );
    } else {
        LOG( "Not using any electronic chirps to shift the system.\n\n" );
    }

    Log::wrapInBar( "Phonons", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( p_phonon_T >= 0 ) {
        std::vector<std::string> approximations = { "Transformation integral via d/dt chi = -i/hbar*[H,chi] + d*chi/dt onto interaction picture chi(t-tau)", "Transformation Matrix U(t,tau)=exp(-i/hbar*H_DQ_L(t)*tau) onto interaction picture chi(t-tau)", "No Transformation, only interaction picture chi(t-tau)", "Analytical Lindblad formalism", "Mixed", "Path Integral" };
        LOG( "Temperature = {} k\n", p_phonon_T );
        LOG( "Cutoff energy = {} Hz - {} meV\n", p_phonon_wcutoff, p_phonon_wcutoff.getSI( Parameter::UNIT_ENERGY_MEV ) );
        LOG( "Cutoff Time = {} ps\n", p_phonon_tcutoff * 1E12 );
        if ( p_phonon_qd_ae == 0.0 ) {
            LOG( "Alpha = {}\n", p_phonon_alpha );
        } else {
            LOG( "Quantum Dot Parameters:\n" );
            LOG( " - Electron Energy D_e = {} Hz - {} eV\n", p_phonon_qd_de, p_phonon_qd_de.getSI( Parameter::UNIT_ENERGY_EV ) );
            LOG( " - Hole Energy D_h = {} Hz - {} eV\n", p_phonon_qd_dh, p_phonon_qd_dh.getSI( Parameter::UNIT_ENERGY_EV ) );
            LOG( " - Material Density rho = {} kg/m^3\n", p_phonon_qd_rho );
            LOG( " - Material Speed of Sound c_s = {} m/s\n", p_phonon_qd_cs );
            LOG( " - Electron Radius = {} nm\n", 1E9 * p_phonon_qd_ae );
            LOG( " - Hole Radius = {} nm\n", 1E9 * p_phonon_qd_ae / p_phonon_qd_ratio );
        }
        LOG( "<B> = {}\n", p_phonon_b );
        LOG( "First Markov approximation used? (rho(t) = rho(t-tau)) - {}\n", ( numerics_phonon_approximation_markov1 ? "Yes" : "No" ) );
        LOG( "Transformation approximation used: {} - {}\n", numerics_phonon_approximation_order, approximations.at( numerics_phonon_approximation_order ) );
        // Pathintegral
        if ( numerics_phonon_approximation_order == 5 ) {
            LOG( " - Path Integral Settings:\n" );
            LOG( " - Backsteps NC: {}\n", p_phonon_nc );
            LOG( " - Iterator Stepsize: {}\n", numerics_subiterator_stepsize );
            LOG( " - Thresholds: Squared({}), SparsePrune({}), CutoffIterations({}), PropagatorMapping({})\n", numerics_pathintegral_squared_threshold, numerics_pathintegral_sparse_prune_threshold, numerics_pathintegral_dynamiccutoff_iterations_max, numerics_pathintegral_docutoff_propagator );
            LOG( " - Used partially summed algorithm?: {}\n", numerics_pathint_partially_summed ? "Yes" : "No" );
        }
        LOG( "\n" );
    } else {
        LOG( "Not using phonons\n\n" );
    }

    Log::wrapInBar( "Numerical Parameters" );
    LOG( "\n" );
    Log::wrapInBar( "Time", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    LOG( "Timeborder start: {:.8e} s - {:.2f} ps\n", t_start, t_start * 1E12 );
    LOG( "Timeborder end: {:.8e} s - {:.2f} ps{}\n", t_end, t_end * 1E12, numerics_calculate_till_converged ? " (variable time end at 99.9\% convergence)" : "" );
    LOG( "Timeborder delta: {:.8e} s - {:.2f} fs \n", t_step, t_step * 1E15 );
    LOG( "Subiterator delta: {:.8e} s - {:.2f} fs \n", numerics_subiterator_stepsize, numerics_subiterator_stepsize * 1E15 );
    if ( numerics_phonon_approximation_order == 5 ) {
        LOG( "Timeborder delta path integral: {:.8e} s - {:.2f} ps\n", t_step_pathint, t_step_pathint * 1E12 );
    }
    LOG( "Time iterations (main loop) = {}\n\n", iterations_t_max );

    Log::wrapInBar( "G-Function Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_correlation.size() > 0 ) {
        LOG( "Tau-grid resolution is {}\n", numerics_calculate_till_converged ? "to be determined." : fmt::format( "{}x{}", grid_values.size(), grid_values.size() ) );
        LOG( "Calculating:\n" );
        for ( auto &[name, mat] : input_correlation ) {
            LOG( " - {} on mode(s) ", name );
            for ( auto &mode : mat.string_v["Modes"] )
                LOG( "{} ", mode );
            LOG( "\n" );
        }
        LOG( "\n" );
        // Detector Stuff
        if ( input_conf["Detector"].numerical_v["time_center"].size() > 0 ) {
            LOG( "Using {} temporal detection windows with parameters:\n", input_conf["Detector"].numerical_v["time_center"].size() );
            for ( int i = 0; i < input_conf["Detector"].numerical_v["time_center"].size(); i++ ) {
                LOG( " - Temporal Detection Window {}:\n", i );
                LOG( " - - Center: {} - {} ps\n", input_conf["Detector"].numerical_v["time_center"][i], input_conf["Detector"].numerical_v["time_center"][i].getSI( Parameter::UNIT_TIME_PS ) );
                LOG( " - - Sigma: {} - {} ps\n", input_conf["Detector"].numerical_v["time_range"][i], input_conf["Detector"].numerical_v["time_range"][i].getSI( Parameter::UNIT_TIME_PS ) );
                LOG( " - - Power: {}\n", input_conf["Detector"].numerical_v["time_power_amplitude"][i] );
            }
        }
        if ( input_conf["Detector"].numerical_v["spectral_center"].size() > 0 ) {
            LOG( "Using {} spectral detection windows with parameters:\n", input_conf["Detector"].numerical_v["spectral_center"].size() );
            for ( int i = 0; i < input_conf["Detector"].numerical_v["spectral_center"].size(); i++ ) {
                LOG( " - Spectral Detection Window {}:\n", i );
                LOG( " - - Center: {} - {} eV\n", input_conf["Detector"].numerical_v["spectral_center"][i], input_conf["Detector"].numerical_v["spectral_center"][i].getSI( Parameter::UNIT_ENERGY_EV ) );
                LOG( " - - Sigma: {} - {} meV\n", input_conf["Detector"].numerical_v["spectral_range"][i], input_conf["Detector"].numerical_v["spectral_range"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                LOG( " - - Power: {}\n", input_conf["Detector"].numerical_v["spectral_power_amplitude"][i] );
                LOG( " - - FT Points: {}\n", input_conf["Detector"].numerical_v["spectral_number_points"][i] );
            }
        }

    } else {
        LOG( "Not using any G1 or G2 correlation functions.\n\n" );
    }
    LOG( "\n" );

    Log::wrapInBar( "Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    LOG( "Solver used: RK{}{}\n", numerics_rk_order, numerics_rk_order != 45 ? "" : fmt::format( " (Tolerance: {}, Stepdelta: {}, Steplimits: [{},{}])", inputstring_rk45_config, numerics_rk_stepdelta, numerics_rk_stepmin, numerics_rk_stepmax ) );
    if ( numerics_rk_order == 45 and numerics_phonon_nork45 )
        LOG( "Will NOT use RK45 for the phonon backwards integral!\n" );
    LOG( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
    LOG( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture ) ? "YES" : "NO" ) );
    LOG( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) ? "Analytic" : "Matrix Exponential" ) );
    LOG( "Threads used for primary calculations - {}\nThreads used for Secondary calculations - {}\nThreads used by Eigen: {}\n", numerics_phonons_maximum_threads, numerics_maximum_threads, Eigen::nbThreads() );
    LOG( "Used scaling for parameters? - {}\n", ( scale_parameters ? std::to_string( scale_value ) : "no" ) );
    if ( p_phonon_T )
        LOG( "Cache Phonon Coefficient Matrices? - {}\n", ( numerics_use_saved_coefficients ? "Yes" : "No" ) );
    if ( numerics_interpolate_outputs )
        LOG( "WARNING: Temporal outputs are interpolated!\n" );
    if ( !numerics_use_function_caching )
        LOG( "NOT using function caching.\n" );
    if ( !numerics_use_saved_hamiltons )
        LOG( "NOT using Hamilton caching.\n" );
    LOG( "\n" );
    for ( int i = 0; i < logfilecounter.size(); i++ ) {
        LOG( "Logfile ident number {}: {}\n", i, logfilecounter[i] );
    }
    LOG( "\n" );

    Log::wrapInBar( "Program Log:", Log::BAR_SIZE_FULL, Log::LEVEL_2 );
    LOG( "\n" );
}