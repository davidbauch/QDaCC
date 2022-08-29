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
    iterations_t_skip = 1;
    scale_parameters = false;
    scale_value = 1E12;
    iterations_t_max = -1;
    maxStates = 0;
    numerics_calculate_till_converged = false;

    // Parsing input:
    Timer &timer_parseInput = Timers::create( "Parsing parameters", true, false );
    Log::Logger::wrapInBar( "Conversion of input variables", Log::Logger::BAR_SIZE_FULL, Log::Logger::LEVEL_2, Log::Logger::BAR_0 );
    Log::L2( "\n" );
    Log::L2( "[System] Parsing input variables...\n" );

    timer_parseInput.start();
    parse_input( arguments );
    timer_parseInput.end();
    Log::L2( "[System] Successful. Elapsed time is {}ms\n", timer_parseInput.getWallTime( Timers::MILLISECONDS ) );

    // Scaling inputs:
    if ( scale_parameters ) {
        Log::L2( "[System] Rescaling parameters to {}...\n", scale_value );
        scale_inputs( scale_value );
    }

    // Adjusting inputs:
    Timer &timer_adjustInput = Timers::create( "Adjusting parameters", true, false );
    Log::L2( "[System] Adjusting input variables...\n" );

    timer_adjustInput.start();
    pre_adjust_input();
    timer_adjustInput.end();
    Log::L2( "[System] Successful. Elapsed time is {}ms\n", timer_adjustInput.getWallTime( Timers::MILLISECONDS ) );
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
    numerics_maximum_primary_threads = QDLC::CommandlineArguments::get_parameter<int>( "--Threads" );
    if ( numerics_maximum_primary_threads == -1 )
        numerics_maximum_primary_threads = omp_get_max_threads();
    output_handlerstrings = QDLC::CommandlineArguments::get_parameter_passed( "-handler" );
    numerics_order_timetrafo = QDLC::CommandlineArguments::get_parameter_passed( "-timeTrafoMatrixExponential" ) ? TIMETRANSFORMATION_MATRIXEXPONENTIAL : TIMETRANSFORMATION_ANALYTICAL;
    scale_parameters = QDLC::CommandlineArguments::get_parameter_passed( "-scale" ); // MHMHMH
    numerics_use_saved_coefficients = not QDLC::CommandlineArguments::get_parameter_passed( "-disableMatrixCaching" );
    numerics_use_saved_hamiltons = not QDLC::CommandlineArguments::get_parameter_passed( "-disableHamiltonCaching" );
    numerics_use_function_caching = not QDLC::CommandlineArguments::get_parameter_passed( "-disableFunctionCaching" );
    numerics_enable_saving_coefficients = false; // If true, even if any saving was disabled internally (not by the user), the matrices will still be cached.
    numerics_maximum_secondary_threads = ( !numerics_use_saved_coefficients || !QDLC::CommandlineArguments::get_parameter_passed( "-disableMainProgramThreading" ) ) ? numerics_maximum_primary_threads : 1;
    logfilecounter = QDLC::Misc::convertParam<double>( QDLC::String::splitline( QDLC::CommandlineArguments::get_parameter( "--lfc" ), ',' ) );
    numerics_interpolate_outputs = QDLC::CommandlineArguments::get_parameter_passed( "-interpolate" );
    s_numerics_interpolate = QDLC::CommandlineArguments::get_parameter( "--interpolateOrder" );

    // Phonon Parameters
    p_phonon_alpha = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phononalpha" );
    p_phonon_wcutoff = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phononwcutoff" );
    p_phonon_wcutoffdelta = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phononwcutoffdelta" );
    p_phonon_tcutoff = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "phonontcutoff" );
    p_phonon_T = QDLC::CommandlineArguments::get_parameter<double>( "--phonons", "temperature" );
    numerics_phonon_approximation_order = QDLC::CommandlineArguments::get_parameter<int>( "--phonons", "phononorder" );
    numerics_phonon_approximation_markov1 = QDLC::CommandlineArguments::get_parameter_passed( "-noMarkov" ) ? 0 : 1; // First Markov
    numerics_phonon_nork45 = not QDLC::CommandlineArguments::get_parameter_passed( "-usePhononRK45" );               // Enables. RK45 for phonon backwards integral; use if detunings are low, otherwise expensive.
    p_phonon_adjust = not QDLC::CommandlineArguments::get_parameter_passed( "-noPhononAdjust" );
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

    auto output_dict_vec = QDLC::String::splitline( QDLC::CommandlineArguments::get_parameter( "--output" ) );
    // Move all elements into the set
    output_dict = std::set<std::string>( std::make_move_iterator( output_dict_vec.begin() ), std::make_move_iterator( output_dict_vec.end() ) );

    kb = 1.3806488E-23;   // J/K, scaling needs to be for energy
    hbar = 1.0545718E-34; // J/s, scaling will be 1

    numerics_main_direction_done = false;

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

void Parameters::pre_adjust_input() {
    Log::L2( "[System] Adjusting Pre-Mainloop Inputs...\n" );

    if ( output_handlerstrings )
        Timers::toggleHandler();

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

    // Adjust Chirp
    for ( auto &[name, mat] : input_chirp ) {
        if ( mat.string["Type"].compare( "sine" ) != 0 ) {
            for ( long unsigned i = 0; i < mat.numerical_v["ddt"].size(); i++ )
                mat.numerical_v["ddt"][i] = mat.numerical_v["ddt"][i] * 1E12;
        }
    }

    // Set interpolation order:
    auto orders = QDLC::String::splitline( s_numerics_interpolate, ',' );
    std::map<std::string, int> methods = { { "monotone", 3 }, { "linear", 0 } };
    std::string method_time = orders.front();
    std::string method_tau = orders.size() > 1 ? orders.back() : "linear";
    numerics_interpolate_method_time = methods[method_time];
    numerics_interpolate_method_tau = methods[method_tau];

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
            // p_omega_pure_dephasing = p_phonon_pure_dephasing * p_phonon_T;
            p_omega_decay = p_omega_decay * p_phonon_b * p_phonon_b; // TODO: unterschiedliche gammas für unterschiedliche B. macht aber auch kaum was.
        }
        if ( numerics_rk_order >= 45 and numerics_use_saved_coefficients )
            numerics_enable_saving_coefficients = true;
    }

    // Calculate minimum step necessary to resolve Rabi-oscillation if step=-1
    if ( t_step < 0 ) {
        t_step = 1E-13; // std::min( scaleVariable( 1E-13, scale_value ), getIdealTimestep() );
        // t_step = std::max( std::numeric_limits<double>::epsilon(), t_step );
        Log::L2( "[System] Delta t was set to {}\n", t_step );
    }

    // Calculate till converged
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
        Log::L2( "[System] Calculate till at least {} and adjust accordingly to guarantee convergence. The matrix index used is {}\n", t_end, numerics_groundstate );
    }

    // Disable Hamilton Caching and RK45 for PI
    if ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
        // numerics_use_saved_hamiltons = false;
        // Log::L2( "[System] Disabled Caching of Hamilton matrices" );
        if ( numerics_rk_order > 5 ) {
            numerics_rk_order = 4;
            Log::L2( "[System] Adjusted RK order to {}\n", numerics_rk_order );
        }
    }

    // Automatically determin subiterator stepsize
    if ( numerics_subiterator_stepsize < 0 ) {
        if ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
            numerics_subiterator_stepsize = t_step_pathint / 5.0;
            Log::L2( "[System] Setting the subiterator stepsize to {}, according to a Path Integral stepsize of {}\n", numerics_subiterator_stepsize, t_step_pathint );
        } else {
            if ( p_phonon_T < 0 and t_step > 0 ) {
                numerics_subiterator_stepsize = t_step / 5.0;
                Log::L2( "[System] Setting the subiterator stepsize to {}, according to a stepsize of {}\n", numerics_subiterator_stepsize, t_step );
            } else if ( p_phonon_T >= 0 and p_phonon_tcutoff > 0 ) {
                numerics_subiterator_stepsize = p_phonon_tcutoff / 200.0;
                Log::L2( "[System] Setting the subiterator stepsize to {}, according to a phonon cutoff time of {}\n", numerics_subiterator_stepsize, p_phonon_tcutoff );
            }
        }
        if ( numerics_subiterator_stepsize == -1 ) {
            numerics_subiterator_stepsize = 1E-13;
            Log::L2( "[System] Setting the subiterator stepsize to a fixed value of {}\n", numerics_subiterator_stepsize );
        }
    }

    // Set Threads to 1 if L3 logging is enabled
    if ( Log::Logger::max_log_level() == Log::Logger::LEVEL_3 ) {
        numerics_maximum_secondary_threads = 1;
        numerics_maximum_primary_threads = 1;
        Log::L2( "[System] Set maximum threads to 1 for all calculations because deeplogging is enabled.\n" );
    }

    // Reserve Trace Vector
    trace.reserve( iterations_t_max + 5 );
}

void Parameters::post_adjust_input() {
    Log::L2( "[System] Adjusting Post-Mainloop Inputs...\n" );

    // Grid Resolution
    if ( t_end >= 0 and ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? ( t_step_pathint > 0 ) : ( t_step > 0 ) ) ) {
        if ( iterations_t_max < 1 ) {
            iterations_t_max = (int)std::ceil( ( t_end - t_start ) / ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? t_step_pathint : t_step ) );
            Log::L2( "[System] Set iterations_t_max to {}\n", iterations_t_max );
        }
        if ( grid_resolution < 1 and iterations_t_max > 0 ) {
            grid_resolution = iterations_t_max; // was +1
        }
        Log::L2( "[System] Set grid_resolution to {}\n", grid_resolution );
    }

    // adjust grids
    post_adjust_grids();
}

void Parameters::post_adjust_grids() {
    Log::L2( "[System] Adjusting Grids...\n" );

    iterations_t_skip = grid_resolution > 0 ? std::max( 1.0, std::ceil( iterations_t_max / grid_resolution ) ) : 1;
    Log::L2( "[System] Maximum t-value for temporal calculations is {}, skip will be {}\n", t_end, iterations_t_skip );

    // Build dt vector. Use standard if not specified otherwise for all calculations. Path integral cannot use other timestep than the original.
    if ( numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? ( t_step_pathint > 0 ) : ( t_step > 0 ) ) {
        input_correlation_resolution["Standard"].numerical_v["Time"] = { t_end };
        // if ( grid_resolution < 0 ) {
        //     grid_resolution = 500;
        //     Log::L2( "[System] No grid_resolution was set! Set gird_resolution to {}\n", grid_resolution );
        // }
        double time_delta = numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? t_step_pathint : Parameter( ( t_end - t_start ) / ( 1. * grid_resolution ) );
        input_correlation_resolution["Standard"].numerical_v["Delta"] = { time_delta };
        Log::L2( "[System] Initial Grid Timestep is {}.\n", time_delta );
        auto &settings = input_correlation_resolution.contains( "Modified" ) ? input_correlation_resolution["Modified"] : input_correlation_resolution["Standard"];
        double skip = 1.0; // input_correlation_resolution.contains( "Modified" ) ? 1.0 : 1.0 * iterations_t_skip; //FIXME: skip braucht man doch eh nicht mehr, da immer das grid erstellt wird und das grid immer entscheidet
        Log::L2( "[System] Iteration Skip for Grid is {}.\n", skip );
        double t_t = 0;
        int current = 0;
        grid_values.clear();
        grid_steps.clear();
        grid_value_indices.clear();
        grid_values.emplace_back( t_start );
        grid_value_indices[t_start] = 0;
        Log::L2( "[System] Initial Timestep Limit is {} at a timestep of {}.\n", settings.numerical_v["Time"][current], settings.numerical_v["Delta"][current] * skip );
        while ( t_t < t_end ) {
            if ( t_t > settings.numerical_v["Time"][current] and current < settings.numerical_v["Time"].size() ) {
                current++;
                Log::L2( "[System] New Timestep Limit is {} at a timestep of {}.\n", settings.numerical_v["Time"][current], settings.numerical_v["Delta"][current] * skip );
            }
            grid_steps.emplace_back( settings.numerical_v["Delta"][current] * skip );
            t_t += grid_steps.back();
            grid_values.emplace_back( t_t );
            grid_value_indices[t_t] = grid_values.size() - 1;
        }
        // std::cout << "Values for "<<mode<<": " << t_values[mode] << std::endl;
        Log::L2( "[System] Setting correlation grid resolution to {0}x{0} for a t_end = {1}\n", grid_values.size(), t_end );
    } else {
        Log::L2( "[System] Not setting time vector because timestep is negative!\n" );
    }
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
        // For counting purposes:
        conf_s.numerical["PulseIndex"] = pindex;
        pindex += 2;
        input_pulse[conf[0]] = conf_s;
    }
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
        conf_s.string_v["Modes"] = QDLC::String::splitline( conf[0], ',' );                                                                                                                               // Modes to calculate Spectrum for. Single modes can again be split with "+", meaning a+b;a to calculate for a+b and a seperately
        conf_s.numerical_v["Center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) );                                                                                    // Center
        conf_s.numerical_v["Range"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) );                                                                                     // Range
        conf_s.numerical_v["resW"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) );                                                                                      // Resolution for w
        conf_s.numerical_v["Order"] = conf.size() > 4 ? QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[4], ',' ) ) : std::vector<Parameter>( conf_s.numerical_v["Range"].size(), 1 ); // Order (1 or 2)?
        conf_s.string_v["Normalize"] = conf.size() > 5 ? QDLC::String::splitline( conf[5], ',' ) : std::vector<std::string>( conf_s.numerical_v["Range"].size(), "False" );                               // Normalize?
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
        if ( conf.size() > 1 ) {
            // Experimental: Calculate spectrum for all matrix entries. Also outputs the non-normalized matrices
            conf_s.numerical_v["Center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) ); // Center
            conf_s.numerical_v["Range"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) );  // Range
            conf_s.numerical_v["resW"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[3], ',' ) );   // Resolution for w
        }
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
    if ( std::ranges::find( inputstring_correlation_resolution.begin(), inputstring_correlation_resolution.end(), ':' ) != inputstring_correlation_resolution.end() ) {
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
            Log::L2( "[System-Parameters] Inheriting Pulse Fourier Configuration from Spectrum.\n" );
            conf_s.numerical["Center"] = conf_s.numerical_v["Center"][0];
            conf_s.numerical["Range"] = conf_s.numerical_v["Range"][0];
            conf_s.numerical["Res"] = conf_s.numerical_v["resW"][0];
            conf_s.numerical["dt"] = t_end / 200;
        } else if ( inputstring_SPconf.size() > 0 ) {
            Log::L2( "[System-Parameters] Generating Pulse Fourier Configuration from Parameters using {}.\n", inputstring_SPconf );
            auto conf = QDLC::String::splitline( inputstring_SPconf, ':' );
            conf_s.numerical["Center"] = conf.size() > 0 ? QDLC::Misc::convertParam<Parameter>( conf[0] ) : Parameter( 0.0 );
            conf_s.numerical["Range"] = conf.size() > 1 ? QDLC::Misc::convertParam<Parameter>( conf[1] ) : Parameter( 0.0 );
            conf_s.numerical["Res"] = conf.size() > 2 ? QDLC::Misc::convertParam<Parameter>( conf[2] ) : Parameter( 0.0 );
            conf_s.numerical["dt"] = conf.size() > 3 ? QDLC::Misc::convertParam<Parameter>( conf[3] ) : Parameter( t_end / 200 );
            input_conf["PulseConf"] = conf_s;
        }
    }
    // Detector Input. Split temporal and spectral filtering with ";"
    {
        input_s conf_s;
        if ( !( inputstring_detector == "none" ) ) {
            Log::L2( "[System-Parameters] Setting up detector using {}...\n", inputstring_detector );
            auto conf_sep = QDLC::String::splitline( inputstring_detector, ';' );
            if ( conf_sep[0] != "none" ) {
                auto conf = QDLC::String::splitline( conf_sep[0], ':' );
                conf_s.numerical_v["time_range"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[0], ',' ) );
                conf_s.numerical_v["time_center"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[1], ',' ) );
                conf_s.numerical_v["time_power_amplitude"] = QDLC::Misc::convertParam<Parameter>( QDLC::String::splitline( conf[2], ',' ) );
                Log::L2( "[System-Parameters] Adding Temporal Detector mask using center = {}, range = {} and power_amp = {}.\n", conf_s.numerical_v["time_center"], conf_s.numerical_v["time_range"], conf_s.numerical_v["time_power_amplitude"] );
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
                Log::L2( "[System-Parameters] Adding Spectral Detector mask using center = {}, range = {}, power_amp = {} and {} points for the Fourier Transformation.\n", conf_s.numerical_v["spectral_center"], conf_s.numerical_v["spectral_range"], conf_s.numerical_v["spectral_power_amplitude"], conf_s.numerical_v["spectral_number_points"] );
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
        Log::L2( "[System-Parameters] Setting up density matrix output using {}...\n", inputstring_densitymatrix_config );
        auto conf = QDLC::String::splitline( inputstring_densitymatrix_config, ':' );
        conf_s.string["output_mode"] = conf[0];
        conf_s.string["interaction_picture"] = conf.size() > 1 ? conf[1] : "schroedinger";
        input_conf["DMconfig"] = conf_s;
    }
    // RK 45 Tolerance Vector
    {
        if ( std::ranges::find( inputstring_rk45_config.begin(), inputstring_rk45_config.end(), ':' ) == inputstring_rk45_config.end() ) {
            numerics_rk_tol.emplace_back( 1.0, QDLC::Misc::convertParam<double>( inputstring_rk45_config ) );
            Log::L2( "[System-Parameters] Set fixed RK45 Tolerance to {}.\n", std::get<1>( numerics_rk_tol.back() ) );
        } else {
            Log::L2( "[System-Parameters] Setting up multiple tolerances...\n" );
            for ( const auto &partial : QDLC::String::splitline( inputstring_rk45_config, ';' ) ) {
                auto tuple = QDLC::String::splitline( partial, ':' );
                numerics_rk_tol.emplace_back( std::make_tuple( QDLC::Misc::convertParam<double>( tuple.front() ), QDLC::Misc::convertParam<double>( tuple.back() ) ) );
                Log::L2( "[System-Parameters] Set RK45 Tolerance to {} (from {}).\n", std::get<1>( numerics_rk_tol.back() ), partial );
            }
        }
    }
}

void Parameters::log( const Dense &initial_state_vector_ket ) {
    // Log::L1( "                                                                                                                                                      \n                                                                                                                  #################                   \n                                                                                                              ..  :::::::::::::::::  ..               \n                                                                                                            ...       -      :+:      ...             \n                                                                                                           ...       .+:     :+:       ...            \n                                                                                                          ...        =+=     :+:        ...           \n                                                                                                         ....       :+++-    :+:         ...          \n         .:----.       ::::::.        :-:          .:---:            .:---:              :---:.          ...        :-+-:    :+:         ...          \n       +@@@@@@@@@*.   =@@@@@@@@@*:   .@@@.       =%@@@@@@@=         %@@%@@@@+          =@@@%@@@*        ....         .+:     :+:          ...         \n     .%@@*:   :*@@@.  =@@#   :+@@@+  .@@@.      %@@%-   :+-         =:   =@@@.        -@@#.  *@@#       ....         .+:     :+:          ...         \n     *@@#       #@@*  =@@#     :@@@: .@@@.     *@@%.                   .:#@@#         %@@=   .@@@.      ...          .+:     :+:          ...         \n     %@@+       +@@#  =@@#      @@@= .@@@.     %@@*         .....    @@@@@@*:         @@@-    @@@:      ...          .+:     :+:          ...         \n     #@@*       #@@*  =@@#     .@@@: .@@@.     #@@#        :@@@@@%   ...:+@@@+        %@@-   .@@@.      ....         .+:     :+:          ...         \n     -@@@-     +@@@:  =@@#    :%@@*  .@@@.     -@@@+     :: ......        %@@#   :-.  *@@*   +@@#       ....         .+:     :+:          ...         \n      -%@@@##%@@@@-   =@@@%%%@@@%=   .@@@%%%%%- :%@@@%#%@@+        =@%#*#@@@%:  *@@%   #@@%*%@@#.        ...         .+:    .-+:.        ....         \n        :=+***==%@@#+ .=++++==-.      =+++++++:   :=+**+=:          -=+**+=:    .++-    :+***=:          ...         .+:    =+++-        ...          \n                 :+#%.                                                                                    ...        .+:    .++=        ...           \n                                                                                                           ...       .+:     -+:        ...           \n      -##########################################################################################           ...      .+:      =        ..             \n      .::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::             ..  ::::-::::::::::::  ...              \n                                                                                                                  #################                   \n                                                                                                                                                      \n                                                                                                                                                      \n" );
    // Log::L1( "╭━━━┳━━━┳╮╱╱╭━━━╮╱╱╭━┳━━╮╭━━┳━╮╱╱╭╮╱╱╱╱╱╱╭━━━╮╭━━╮\n┃╭━╮┣╮╭╮┃┃╱╱┃╭━╮┃╱╭╯╭┻┫┣╯╰┫┣┻╮╰╮╱┃┃╱╱╱╱╱╱╰╮╭╮┃┃╭╮┃\n┃┃╱┃┃┃┃┃┃┃╱╱┃┃╱╰╯╭╯╭╯╱┃┃╱╱┃┃╱╰╮╰╮┃╰━┳╮╱╭╮╱┃┃┃┃┃╰╯╰╮\n┃┃╱┃┃┃┃┃┃┃╱╭┫┃╱╭╮┃┃┃╱╱┃┃╱╱┃┃╱╱┃┃┃┃╭╮┃┃╱┃┃╱┃┃┃┃┃╭━╮┃\n┃╰━╯┣╯╰╯┃╰━╯┃╰━╯┃┃┃┃╱╭┫┣╮╭┫┣╮╱┃┃┃┃╰╯┃╰━╯┃╭╯╰╯┣┫╰━╯┣╮\n╰━━╮┣━━━┻━━━┻━━━╯╰╮╰╮╰━━╯╰━━╯╭╯╭╯╰━━┻━╮╭╯╰━━━┻┻━━━┻╯\n╱╱╱╰╯╱╱╱╱╱╱╱╱╱╱╱╱╱╰╮╰╮╱╱╱╱╱╱╭╯╭╯╱╱╱╱╭━╯┃\n╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╰━╯╱╱╱╱╱╱╰━╯╱╱╱╱╱╰━━╯\n" );
    Log::Logger::P1( "\033[38;2;128;128;128m                                                                           \033[38;2;165;146;89m▄\033[38;2;165;146;90m▄\033[0m\n\033[38;2;128;128;128m                                                                   \033[38;2;179;151;76m▐\033[38;2;250;189;4m▒\033[38;2;228;178;27m▒\033[38;2;169;148;85m▄    \033[38;2;244;185;10m▒\033[38;2;254;192;0m▒\033[38;2;251;190;3m▒\033[38;2;202;164;53m▄            \033[38;2;186;155;69m▐\033[38;2;212;170;43m▄\033[38;2;223;175;31m▒\033[38;2;235;181;20m▒\033[38;2;246;187;8m▒\033[38;2;250;189;4m▒\033[0m\n\033[38;2;128;128;128m                                                                    \033[38;2;227;176;27m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;239;183;16m▒\033[38;2;170;148;84m▄   \033[38;2;233;179;22m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;253;191;1m▒\033[38;2;208;167;47m▄    \033[38;2;181;153;74m▄\033[38;2;213;170;41m▒\033[38;2;195;161;60m▄   \033[38;2;200;163;55m▄\033[38;2;253;191;1m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;246;187;8m▒\033[0m\n\033[38;2;128;128;128m       \033[38;2;78;102;146m▄\033[38;2;60;94;152m▓                                            \033[38;2;117;129;110m▐\033[38;2;92;134;64m▓\033[38;2;92;134;64m▓▓▓▓▓▓▓▓▓▓▓▓▓ \033[38;2;205;165;49m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;235;181;20m▒\033[38;2;162;144;92m▄  \033[38;2;207;166;48m▐\033[38;2;254;192;0m▒\033[38;2;252;190;2m▒\033[38;2;238;182;17m░\033[38;2;253;191;1m▒\033[38;2;208;167;47m▄  \033[38;2;198;160;57m▐\033[38;2;254;192;0m▒\033[38;2;254;192;0m▒\033[38;2;250;189;4m▒\033[38;2;249;189;5m▒\033[38;2;254;192;0m▒\033[38;2;224;175;31m▒  \033[38;2;212;169;43m▀               \033[38;2;106;116;135m▐\033[38;2;62;95;152m▓\033[38;2;83;105;144m▄\033[0m\n\033[38;2;128;128;128m     \033[38;2;62;95;152m▓\033[38;2;33;80;163m▒\033[38;2;34;81;163m▒\033[38;2;34;81;162m▒                                            \033[38;2;113;130;102m▐\033[38;2;77;138;36m▒\033[38;2;77;138;36m▒▒▒▒▒▒▒▒▒▒▒▒▒   \033[38;2;252;190;2m▒\033[38;2;254;192;0m▒\033[38;2;235;180;20m▒\033[38;2;254;192;0m▒\033[38;2;226;176;28m▒   \033[38;2;251;190;3m▒\033[38;2;244;186;10m▒ \033[38;2;233;180;21m▐\033[38;2;253;191;1m▒\033[38;2;206;166;49m▄ \033[38;2;204;164;50m▐\033[38;2;254;192;0m▒\033[38;2;218;172;36m▒                     \033[38;2;99;113;138m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;34;81;162m▒\033[38;2;75;101;147m▄\033[0m\n\033[38;2;128;128;128m    \033[38;2;43;85;159m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;45;88;159m▒\033[38;2;78;111;150m░\033[38;2;80;105;133m░\033[38;2;92;108;125m░\033[38;2;105;111;119m▄                                            \033[38;2;37;57;68m▓\033[38;2;71;85;88m▌       \033[38;2;63;78;83m▓      \033[38;2;238;183;16m▒\033[38;2;253;191;1m▒\033[38;2;204;164;50m░\033[38;2;223;174;32m▀\033[38;2;254;192;0m▒\033[38;2;215;171;40m▒  \033[38;2;235;181;20m▐\033[38;2;253;191;1m▒\033[38;2;190;157;65m▄ \033[38;2;234;181;20m▐\033[38;2;253;191;1m▒\033[38;2;205;166;50m▄\033[38;2;209;166;45m▐\033[38;2;254;192;0m▒\033[38;2;223;174;32m▒                 \033[38;2;104;112;120m▄ \033[38;2;87;107;129m▒ \033[38;2;73;108;151m░\033[38;2;39;84;161m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;55;91;154m▓\033[0m\n\033[38;2;128;128;128m   \033[38;2;47;87;157m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;41;85;160m▒\033[38;2;86;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;87;117;148m░ \033[38;2;82;106;132m░\033[38;2;89;108;129m▒\033[38;2;94;110;127m░\033[38;2;99;111;124m▄\033[38;2;104;112;122m▄\033[38;2;108;113;119m▄    \033[38;2;109;112;115m▄\033[38;2;101;105;110m▄\033[38;2;97;102;108m▄\033[38;2;96;101;107m▄\033[38;2;98;102;108m▄\033[38;2;103;106;111m▄        \033[38;2;105;108;112m▄\033[38;2;105;108;113m▄\033[38;2;105;108;113m▄▄▄▄▄\033[38;2;105;108;113m▄\033[38;2;108;111;115m▄          \033[38;2;34;50;72m▓\033[38;2;71;81;94m▌      \033[38;2;47;61;80m▓\033[38;2;34;50;72m▓\033[38;2;47;61;80m█      \033[38;2;134;120;76m▄\033[38;2;134;120;76m▄\033[38;2;110;108;100m▄\033[38;2;125;116;84m▄\033[38;2;134;120;76m▄ \033[38;2;186;154;68m▄ \033[38;2;203;164;51m▐\033[38;2;254;192;0m▒\033[38;2;225;175;30m▒  \033[38;2;233;180;21m▀\033[38;2;250;189;4m▒\033[38;2;214;170;36m▀  \033[38;2;103;107;112m▄\033[38;2;98;103;109m▄\033[38;2;96;101;107m▄\033[38;2;97;102;108m▄\033[38;2;103;107;111m▄    \033[38;2;108;113;119m▄\033[38;2;102;112;122m▄ \033[38;2;91;109;128m░\033[38;2;85;107;131m▒\033[38;2;86;87;96m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░\033[38;2;83;114;149m░\033[38;2;35;81;162m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;64;95;151m▌\033[0m\n\033[38;2;128;128;128m  \033[38;2;81;104;144m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;33;80;162m▒\033[38;2;76;110;151m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░\033[38;2;87;117;148m░ \033[38;2;87;93;102m▄\033[38;2;47;60;79m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓▓▓\033[38;2;35;51;72m▓\033[38;2;56;68;84m█\033[38;2;103;107;111m▄    \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓▓▓▓\033[38;2;34;50;72m▓\033[38;2;37;52;74m█\033[38;2;61;72;87m█\033[38;2;103;107;111m▄      \033[38;2;34;50;72m▓\033[38;2;71;81;94m▌     \033[38;2;75;84;96m▀\033[38;2;55;67;84m▀\033[38;2;34;50;72m▓\033[38;2;55;67;84m▀\033[38;2;75;84;96m▀ \033[38;2;227;177;27m▒\033[38;2;231;179;23m▒ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;223;172;30m▐\033[38;2;250;189;4m▒\033[38;2;189;157;66m▄ \033[38;2;249;188;5m▒\033[38;2;247;188;7m▒  \033[38;2;104;105;105m▄\033[38;2;56;68;84m█\033[38;2;35;51;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓▓\033[38;2;34;50;72m▓\033[38;2;38;54;74m█\033[38;2;69;79;92m█ \033[38;2;96;118;141m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░\033[38;2;131;93;73m▐\033[38;2;147;95;64m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░\033[38;2;64;101;154m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;33;80;163m▒\033[0m\n\033[38;2;128;128;128m  \033[38;2;35;81;162m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;39;84;161m▓\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░\033[38;2;98;119;139m░\033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓\033[38;2;34;50;72m▓\033[38;2;81;89;98m░ \033[38;2;65;76;90m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;35;51;73m▓    \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;55;68;84m▌ \033[38;2;67;77;90m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓\033[38;2;98;102;108m▌     \033[38;2;34;50;72m▓\033[38;2;71;81;94m▌      \033[38;2;93;99;107m▐\033[38;2;34;50;72m▓\033[38;2;93;99;107m▌  \033[38;2;249;188;5m▒\033[38;2;254;192;0m▒ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;110;113;116m▄\033[38;2;156;134;56m▒\033[38;2;247;188;5m▒\033[38;2;244;186;10m▒\033[38;2;174;145;60m▄\033[38;2;216;170;23m░\033[38;2;245;186;8m▒ \033[38;2;36;52;73m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;49;62;80m▌ \033[38;2;80;88;98m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;40;55;75m█ \033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;131;93;73m▐\033[38;2;212;101;24m▒\033[38;2;139;94;69m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░\033[38;2;86;116;148m░\033[38;2;34;81;162m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;47;87;157m▓\033[0m\n\033[38;2;128;128;128m \033[38;2;96;112;139m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[38;2;58;97;155m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░░ \033[38;2;86;93;101m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓  \033[38;2;87;93;102m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;95;118;142m░\033[38;2;93;118;144m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;65;75;89m▌ \033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;73;83;94m▌ \033[38;2;91;109;129m▒\033[38;2;91;109;128m░\033[38;2;92;109;128m░\033[38;2;84;102;120m░\033[38;2;34;50;72m▓\033[38;2;58;74;94m▌\033[38;2;95;110;127m░\033[38;2;95;110;126m░\033[38;2;95;110;126m░\033[38;2;95;110;127m░\033[38;2;94;110;127m░\033[38;2;94;110;127m░\033[38;2;71;88;107m▐\033[38;2;34;50;72m▓\033[38;2;70;88;108m▒\033[38;2;91;109;129m▒\033[38;2;90;109;129m▒\033[38;2;119;115;79m░\033[38;2;254;192;0m▒ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;172;145;58m▐\033[38;2;89;115;138m░\033[38;2;89;114;135m░\033[38;2;169;145;52m░\033[38;2;252;191;2m▒\033[38;2;229;178;16m▒\033[38;2;205;163;37m▒ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;59;71;86m▌ \033[38;2;101;105;110m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;35;51;72m▓ \033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;131;93;73m▐\033[38;2;212;101;25m▒\033[38;2;212;101;24m▒\033[38;2;131;93;73m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░\033[38;2;45;88;159m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;33;80;163m▒\033[0m\n\033[38;2;128;128;128m \033[38;2;65;96;151m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[38;2;72;107;152m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░░ \033[38;2;86;93;101m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓  \033[38;2;87;93;102m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;95;118;142m░\033[38;2;93;118;144m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;65;75;89m▌ \033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;73;83;94m▌ \033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;80;109;138m▐\033[38;2;34;50;72m▓\033[38;2;55;77;102m▌\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░\033[38;2;68;93;120m▐\033[38;2;34;50;72m▓\033[38;2;90;92;66m▒\033[38;2;177;150;46m▒\033[38;2;187;156;40m▒\033[38;2;161;141;56m░\033[38;2;218;171;23m▒ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;227;175;25m▐\033[38;2;232;179;14m▒\033[38;2;93;110;112m░\033[38;2;87;117;148m░\033[38;2;91;111;121m░\033[38;2;186;155;42m░\033[38;2;235;180;18m▒ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;59;71;86m▌ \033[38;2;101;105;110m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓\033[38;2;35;51;72m▓ \033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;131;93;73m▐\033[38;2;212;101;25m▒\033[38;2;212;101;25m▒\033[38;2;212;101;24m▒\033[38;2;123;93;78m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;59;98;155m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[0m\n\033[38;2;128;128;128m \033[38;2;45;86;158m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;34;81;162m▒\033[38;2;81;113;150m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░░ \033[38;2;86;93;101m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓  \033[38;2;87;93;102m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;95;118;142m░\033[38;2;93;118;144m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;65;75;89m▌ \033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;73;83;94m▌ \033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;80;109;138m▐\033[38;2;34;50;72m▓\033[38;2;55;77;102m▌\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;77;104;131m▄\033[38;2;74;99;125m▄ \033[38;2;102;110;91m░\033[38;2;107;100;56m▐\033[38;2;34;50;72m▓\033[38;2;154;124;30m▒\033[38;2;165;144;54m░\033[38;2;121;120;80m░\033[38;2;157;139;59m░\033[38;2;189;156;40m░ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;227;175;25m▐\033[38;2;254;192;0m▒\033[38;2;247;188;5m▒\033[38;2;96;107;95m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;96;111;115m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;59;71;86m▌ \033[38;2;105;119;133m░\033[38;2;104;118;133m░\033[38;2;104;118;133m░░░\033[38;2;104;118;133m░\033[38;2;105;119;134m░\033[38;2;87;117;148m░\033[38;2;103;91;90m▐\033[38;2;187;98;39m▒\033[38;2;187;98;39m▒\033[38;2;199;99;32m▒\033[38;2;212;101;25m▒\033[38;2;212;101;25m▒▒\033[38;2;212;101;24m▒\033[38;2;115;92;83m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;68;104;153m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[38;2;68;97;150m▌\033[0m\n\033[38;2;128;128;128m \033[38;2;37;82;161m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;34;80;162m▒\033[38;2;85;115;149m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░░ \033[38;2;86;93;101m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓  \033[38;2;87;93;102m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;95;118;142m░\033[38;2;93;118;144m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;65;75;89m▌ \033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;73;83;94m▌ \033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;80;109;138m▐\033[38;2;34;50;72m▓\033[38;2;55;77;102m▌\033[38;2;74;99;125m▐\033[38;2;8;10;13m█\033[38;2;0;0;0m█\033[38;2;0;0;0m█\033[38;2;1;1;0m█\033[38;2;107;81;3m█\033[38;2;123;108;47m▐\033[38;2;34;50;72m▓\033[38;2;68;93;120m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;117;118;89m░\033[38;2;184;154;43m░\033[38;2;215;170;24m░\033[38;2;141;130;69m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;91;118;145m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;59;71;86m▌ \033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░\033[38;2;111;91;85m▐\033[38;2;212;101;25m▒\033[38;2;212;101;25m▒▒▒▒▒▒\033[38;2;212;100;24m▒\033[38;2;108;91;87m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;72;107;152m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[38;2;58;93;153m▌\033[0m\n\033[38;2;128;128;128m \033[38;2;37;82;161m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;34;80;162m▒\033[38;2;85;115;149m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░░ \033[38;2;86;93;101m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓  \033[38;2;87;93;102m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;95;118;142m░\033[38;2;93;118;144m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;65;75;89m▌ \033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;73;83;94m▌ \033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;80;109;138m▐\033[38;2;34;50;72m▓\033[38;2;55;77;102m▌\033[38;2;74;99;125m▐\033[38;2;8;10;13m█\033[38;2;0;0;0m█\033[38;2;0;0;0m█\033[38;2;1;1;1m█\033[38;2;38;51;64m▀\033[38;2;68;93;120m▐\033[38;2;34;50;72m▓\033[38;2;68;93;120m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;93;118;143m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░\033[38;2;91;118;145m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;59;71;86m▌ \033[38;2;91;118;144m░\033[38;2;91;118;144m░\033[38;2;91;118;144m░░░░\033[38;2;91;118;145m░\033[38;2;87;117;148m░\033[38;2;111;91;85m▐\033[38;2;212;101;25m▒\033[38;2;212;101;25m▒▒▒▒▒▒\033[38;2;212;100;24m▒\033[38;2;109;91;87m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;72;107;152m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[38;2;59;93;153m▌\033[0m\n\033[38;2;128;128;128m \033[38;2;45;86;158m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;34;81;162m▒\033[38;2;81;113;150m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░░ \033[38;2;86;93;101m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓  \033[38;2;87;93;102m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;95;118;142m░\033[38;2;93;118;144m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;65;75;89m▌ \033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;73;83;94m▌ \033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;80;109;138m▐\033[38;2;34;50;72m▓\033[38;2;55;77;102m▌\033[38;2;87;117;148m░\033[38;2;87;117;148m░   ░\033[38;2;68;93;120m▐\033[38;2;34;50;72m▓\033[38;2;68;93;120m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;93;118;143m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░\033[38;2;91;118;145m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;59;71;86m▌ \033[38;2;107;110;113m▐\033[38;2;60;71;87m█\033[38;2;60;71;87m█\033[38;2;60;71;87m██\033[38;2;60;71;87m█ \033[38;2;87;117;148m░\033[38;2;103;91;90m▐\033[38;2;188;98;39m▒\033[38;2;188;98;39m▒\033[38;2;199;99;32m▒\033[38;2;212;101;25m▒\033[38;2;212;101;25m▒▒\033[38;2;212;100;24m▒\033[38;2;116;92;83m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;67;103;153m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[38;2;68;97;150m▌\033[0m\n\033[38;2;128;128;128m \033[38;2;65;96;150m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[38;2;72;107;152m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░░ \033[38;2;86;93;101m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓  \033[38;2;87;93;102m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;95;118;142m░\033[38;2;93;118;144m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;65;75;89m▌ \033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;73;83;94m▌ \033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;80;109;138m▐\033[38;2;34;50;72m▓\033[38;2;55;77;102m▌\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░\033[38;2;68;93;120m▐\033[38;2;34;50;72m▓\033[38;2;68;93;120m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;93;118;143m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░\033[38;2;91;118;145m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;59;71;86m▌ \033[38;2;101;105;110m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;35;51;72m▓ \033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;131;93;73m▐\033[38;2;212;101;25m▒\033[38;2;212;101;25m▒\033[38;2;212;101;24m▒\033[38;2;123;93;78m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;58;98;155m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[0m\n\033[38;2;128;128;128m \033[38;2;98;112;138m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒▒\033[38;2;58;97;156m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░░ \033[38;2;86;93;101m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓  \033[38;2;87;93;102m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;95;118;142m░\033[38;2;93;118;144m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;65;75;89m▌ \033[38;2;100;104;109m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;74;83;94m▌ \033[38;2;87;108;131m░\033[38;2;88;109;130m▀\033[38;2;89;109;130m▀\033[38;2;82;101;121m▀\033[38;2;34;50;72m▓\033[38;2;56;74;95m▌\033[38;2;90;109;129m▀\033[38;2;90;109;129m▀▀▀▀▀\033[38;2;69;87;108m▐\033[38;2;34;50;72m▓\033[38;2;69;88;109m▒\033[38;2;88;108;130m▀\033[38;2;87;108;130m▀\033[38;2;87;108;131m░\033[38;2;86;108;131m░ \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓ \033[38;2;89;114;138m░\033[38;2;86;115;146m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;91;118;145m░ \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;59;71;86m▌ \033[38;2;101;105;110m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;35;51;72m▓ \033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;131;93;73m▐\033[38;2;212;101;25m▒\033[38;2;212;100;24m▒\033[38;2;131;93;73m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░\033[38;2;44;87;159m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;33;80;163m▒\033[0m\n\033[38;2;128;128;128m  \033[38;2;36;81;162m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;38;84;161m▓\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░░░\033[38;2;100;119;138m░\033[38;2;96;101;107m▐\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓\033[38;2;34;50;72m▓\033[38;2;80;88;98m▄\033[38;2;98;104;112m▄\033[38;2;65;75;89m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;35;51;72m▓    \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;55;67;84m▌\033[38;2;100;104;109m▄\033[38;2;66;76;90m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓      \033[38;2;34;50;72m▓\033[38;2;71;81;94m▌      \033[38;2;93;99;107m▐\033[38;2;34;50;72m▓\033[38;2;93;99;107m▌     \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;34;50;72m▓\033[38;2;91;97;104m▄\033[38;2;101;105;110m▄\033[38;2;97;102;108m▄\033[38;2;91;97;105m▄\033[38;2;87;94;103m▄    \033[38;2;35;51;73m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;48;61;80m█\033[38;2;100;104;109m▄\033[38;2;78;86;97m▄\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓\033[38;2;37;53;74m█ \033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░\033[38;2;131;93;73m▐\033[38;2;212;101;24m▒\033[38;2;139;94;69m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░\033[38;2;86;116;148m░\033[38;2;34;80;162m▒\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;48;88;157m▌\033[0m\n\033[38;2;128;128;128m  \033[38;2;83;105;144m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;33;80;162m▒\033[38;2;75;109;151m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░░░░░\033[38;2;87;117;148m░\033[38;2;83;112;142m░\033[38;2;80;107;134m▒  \033[38;2;49;63;81m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓▓▓\033[38;2;36;51;73m▓\033[38;2;60;72;87m▀     \033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓▓▓▓\033[38;2;34;50;72m▓\033[38;2;38;53;74m█\033[38;2;63;74;88m▀     \033[38;2;87;94;103m▐\033[38;2;67;77;91m▓\033[38;2;34;50;72m▓\033[38;2;49;63;81m█\033[38;2;73;83;95m▓     \033[38;2;93;99;107m▐\033[38;2;34;50;72m▓\033[38;2;93;99;107m▌     \033[38;2;60;72;87m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓▓▓▓▓▓\033[38;2;71;80;92m▌    \033[38;2;66;77;90m▀\033[38;2;40;55;75m█\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓▓▓▓\033[38;2;34;50;72m▓\033[38;2;36;52;73m▓\033[38;2;62;73;87m▀ \033[38;2;91;110;130m▄\033[38;2;86;117;147m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;131;93;73m▐\033[38;2;147;95;64m▒\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░░\033[38;2;63;100;154m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;33;80;163m▒\033[0m\n\033[38;2;128;128;128m   \033[38;2;49;88;157m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;40;85;161m▓\033[38;2;86;116;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░\033[38;2;87;117;147m░\033[38;2;80;107;136m▄\033[38;2;87;107;128m▀\033[38;2;97;110;123m▀           \033[38;2;52;65;82m▓\033[38;2;34;50;72m▓\033[38;2;34;50;72m▓\033[38;2;69;79;92m▓\033[38;2;88;94;103m▄\033[38;2;90;95;104m▄                      \033[38;2;55;68;85m▀\033[38;2;34;50;72m▓\033[38;2;37;53;74m█      \033[38;2;93;99;107m▐\033[38;2;34;50;72m▓\033[38;2;93;99;107m▌                                     \033[38;2;91;109;127m▀\033[38;2;91;87;92m▒\033[38;2;81;110;139m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░\033[38;2;87;117;148m░░░\033[38;2;82;114;149m░\033[38;2;35;81;162m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;67;97;150m▌\033[0m\n\033[38;2;128;128;128m    \033[38;2;44;86;158m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;43;87;160m▓\033[38;2;76;109;151m░\033[38;2;80;107;135m▄\033[38;2;88;107;127m▀                 \033[38;2;61;72;87m▀\033[38;2;41;56;76m█\033[38;2;36;52;73m▓\033[38;2;36;52;73m▓\033[38;2;36;51;73m▓                       \033[38;2;66;78;89m▀       \033[38;2;92;100;104m▐\033[38;2;36;53;72m▓\033[38;2;93;100;105m▌                                          \033[38;2;89;108;127m▀\033[38;2;80;107;136m▒\033[38;2;72;107;152m▄\033[38;2;38;83;161m▓\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;57;92;153m▀\033[0m\n\033[38;2;128;128;128m     \033[38;2;65;96;150m▀\033[38;2;33;80;162m▒\033[38;2;34;81;163m▒\033[38;2;34;81;162m▒                                            \033[38;2;104;132;86m▐\033[38;2;76;138;35m▒\033[38;2;76;138;35m▒▒▒▒▒▒▒▒▒\033[38;2;76;138;35m▒\033[38;2;77;138;35m▒\033[38;2;76;138;35m▒\033[38;2;76;138;35m▒                                          \033[38;2;99;113;138m▐\033[38;2;34;81;163m▒\033[38;2;34;81;163m▒\033[38;2;34;81;162m▓\033[38;2;78;103;146m▀\033[0m\n\033[38;2;128;128;128m       \033[38;2;81;104;145m▀\033[38;2;63;95;151m▀                                            \033[38;2;110;131;96m▐\033[38;2;89;135;57m▀\033[38;2;89;135;57m▀▀▀▀▀▀▀▀▀▀▀▀▀                                           \033[38;2;65;96;150m▀\033[0m\n\033[38;2;128;128;128m \033[0m\n                                                   Quantum  Simulation\n" );
    Log::Logger::wrapInBar( "System Parameters" );
    Log::L1( "Version: {} ({})\n\n", GLOBAL_PROGRAM_VERSION, GLOBAL_PROGRAM_LASTCHANGE );

    Log::Logger::wrapInBar( "Electronic Configuration", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
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
        Log::Logger::wrapInBar( "Photonic Configuration", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
        for ( auto &[name, mat] : input_photonic ) {
            Log::L1( "Cavity: {} with energy {:.8e} Hz - {:.8e} eV - {:.8f} nm\n", name, mat.numerical["Energy"], mat.numerical["Energy"].getSI( Parameter::UNIT_ENERGY_EV ), mat.numerical["Energy"].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
            for ( auto i = 0; i < mat.string_v["CoupledTo"].size(); i++ ) {
                auto [from, to] = QDLC::String::split_pair( mat.string_v["CoupledTo"][i], transition_delimiter );
                double purcell = ( 2.0 * p_omega_coupling * p_omega_coupling * mat.numerical_v["CouplingScaling"][i] * mat.numerical_v["CouplingScaling"][i] / p_omega_cavity_loss / p_omega_decay ) * ( p_omega_cavity_loss * p_omega_cavity_loss / ( std::pow( mat.numerical["Energy"] - ( input_electronic[to].numerical["Energy"] - input_electronic[from].numerical["Energy"] ), 2 ) + p_omega_cavity_loss * p_omega_cavity_loss ) );
                Log::L1( " - Coupled to electronic transition {} (Purcell Enhancement: {}) \n", mat.string_v["CoupledTo"][i], purcell );
                Log::L1( " - - Coupling Scaling: {}\n", mat.numerical_v["CouplingScaling"][i] );
            }
            Log::L1( " - Cavity Q-Factor: {:.0f}\n", mat.numerical["Energy"] / p_omega_cavity_loss );
            Log::L1( " - Decay Scaling: {}\n", mat.numerical["DecayScaling"] );
        }
        Log::L1( "\n" );
    }

    Log::Logger::wrapInBar( "Coupling Parameters", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
    Log::L1( "Coupling strengh g: {:.8e} Hz - {:.8} mueV\n", p_omega_coupling, p_omega_coupling.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Photon loss rate k: {:.8e} Hz - {:.8} mueV\n", p_omega_cavity_loss, p_omega_cavity_loss.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Atomic dephasing rate gamma_pure: {:.8e} Hz - {:.8} mueV\n", p_omega_pure_dephasing, p_omega_pure_dephasing.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "RAD rate gamma: {:.8e} Hz - {:.8} mueV\n\n", p_omega_decay, p_omega_decay.getSI( Parameter::UNIT_ENERGY_MUEV ) );

    Log::Logger::wrapInBar( "Initial System Parameters", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
    Log::L1( "Initial state rho0 = [{}]\n\n", initial_state_vector_ket.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ) );

    Log::Logger::wrapInBar( "Pulse", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
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
                Log::L1( " - - Frequency: {} Hz - {:.8} eV - {:.8} nm\n", mat.numerical_v["Frequency"][i], mat.numerical_v["Frequency"][i].getSI( Parameter::UNIT_ENERGY_EV ), mat.numerical_v["Frequency"][i].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
                Log::L1( " - - Width: {} s - {:.8} ps\n", mat.numerical_v["Width"][i], mat.numerical_v["Width"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Center: {} s - {:.8} ps\n", mat.numerical_v["Center"][i], mat.numerical_v["Center"][i].getSI( Parameter::UNIT_TIME_PS ) );
                if ( QDLC::Math::abs2( mat.numerical_v["ChirpRate"][i] != 0.0 ) )
                    Log::L1( " - - Chirp: {}\n", mat.numerical_v["ChirpRate"][i] );
                if ( QDLC::Math::abs2( mat.numerical_v["SUPERDelta"][i] != 0.0 ) ) {
                    Log::L1( " - - SUPER Amplitude: {} - {} meV\n", mat.numerical_v["SUPERDelta"][i], mat.numerical_v["SUPERDelta"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                    Log::L1( " - - SUPER Frequency: {} - {} meV\n", mat.numerical_v["SUPERFreq"][i], mat.numerical_v["SUPERFreq"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                }
                Log::L1( " - - Type: {}{}\n", mat.string_v["Type"][i], mat.string_v["Type"][i] == "gauss" ? fmt::format( " (Gaussian Amplitude: {})", mat.numerical_v["GaussAmp"][i] ) : "" );
            }
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using any pulses to exite or drive the system.\n\n" );
    }
    Log::Logger::wrapInBar( "Chirp", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
    if ( input_chirp.size() > 0 ) {
        for ( auto &[name, mat] : input_chirp ) {
            Log::L1( "Chirp: {}\n", name );
            Log::L1( " - Coupled to States:\n" );
            for ( int i = 0; i < mat.string_v["CoupledTo"].size(); i++ )
                Log::L1( " - - {} with scaling {}\n", mat.string_v["CoupledTo"][i], mat.numerical_v["AmpFactor"][i] );
            for ( int i = 0; i < mat.numerical_v["Amplitude"].size(); i++ ) {
                Log::L1( " - Chirp Point {}:\n", i );
                Log::L1( " - - Amplitude: {} Hz - {} meV\n", mat.numerical_v["Amplitude"][i], mat.numerical_v["Amplitude"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                Log::L1( " - - Time: {} s - {} ps\n", mat.numerical_v["Times"][i], mat.numerical_v["Times"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Derivative DDT: {}\n", mat.numerical_v["ddt"][i] );
            }
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using any electronic chirps to shift the system.\n\n" );
    }

    Log::Logger::wrapInBar( "Phonons", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
    if ( p_phonon_T >= 0 ) {
        std::vector<std::string> approximations = { "Transformation integral via d/dt chi = -i/hbar*[H,chi] + d*chi/dt onto interaction picture chi(t-tau)", "Transformation Matrix U(t,tau)=exp(-i/hbar*H_DQ_L(t)*tau) onto interaction picture chi(t-tau)", "No Transformation, only interaction picture chi(t-tau)", "Analytical Lindblad formalism", "Mixed", "Path Integral" };
        Log::L1( "Temperature = {} k\n", p_phonon_T );
        Log::L1( "Cutoff energy = {} Hz - {} meV\n", p_phonon_wcutoff, p_phonon_wcutoff.getSI( Parameter::UNIT_ENERGY_MEV ) );
        Log::L1( "Cutoff Time = {} ps\n", p_phonon_tcutoff * 1E12 );
        Log::L1( "Phonon Integral Iterator Stepsize = {} ps\n", numerics_subiterator_stepsize * 1E12 );
        if ( p_phonon_qd_ae == 0.0 ) {
            Log::L1( "Alpha = {}\n", p_phonon_alpha );
        } else {
            Log::L1( "Quantum Dot Parameters:\n" );
            Log::L1( " - Electron Energy D_e = {} Hz - {} eV\n", p_phonon_qd_de, p_phonon_qd_de.getSI( Parameter::UNIT_ENERGY_EV ) );
            Log::L1( " - Hole Energy D_h = {} Hz - {} eV\n", p_phonon_qd_dh, p_phonon_qd_dh.getSI( Parameter::UNIT_ENERGY_EV ) );
            Log::L1( " - Material Density rho = {} kg/m^3\n", p_phonon_qd_rho );
            Log::L1( " - Material Speed of Sound c_s = {} m/s\n", p_phonon_qd_cs );
            Log::L1( " - Electron Radius = {} nm\n", 1E9 * p_phonon_qd_ae );
            Log::L1( " - Hole Radius = {} nm\n", 1E9 * p_phonon_qd_ae / p_phonon_qd_ratio );
        }
        Log::L1( "<B> = {}\n", p_phonon_b );
        Log::L1( "First Markov approximation used? (rho(t) = rho(t-tau)) - {}\n", ( numerics_phonon_approximation_markov1 ? "Yes" : "No" ) );
        Log::L1( "Transformation approximation used: {} - {}\n", numerics_phonon_approximation_order, approximations.at( numerics_phonon_approximation_order ) );
        // Pathintegral
        if ( numerics_phonon_approximation_order == 5 ) {
            Log::L1( " - Path Integral Settings:\n" );
            Log::L1( " - Backsteps NC: {}\n", p_phonon_nc );
            Log::L1( " - Iterator Stepsize: {}\n", numerics_subiterator_stepsize );
            Log::L1( " - Thresholds: Squared({}), SparsePrune({}), CutoffIterations({}), PropagatorMapping({})\n", numerics_pathintegral_squared_threshold, numerics_pathintegral_sparse_prune_threshold, numerics_pathintegral_dynamiccutoff_iterations_max, numerics_pathintegral_docutoff_propagator );
            Log::L1( " - Used partially summed algorithm?: {}\n", numerics_pathint_partially_summed ? "Yes" : "No" );
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using phonons\n\n" );
    }

    Log::Logger::wrapInBar( "Numerical Parameters" );
    Log::L1( "\n" );
    Log::Logger::wrapInBar( "Time", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
    Log::L1( "Timeborder start: {:.8e} s - {:.2f} ps\n", t_start, t_start * 1E12 );
    Log::L1( "Timeborder end: {:.8e} s - {:.2f} ps{}\n", t_end, t_end * 1E12, numerics_calculate_till_converged ? " (variable time end at 99.9\% convergence)" : "" );
    Log::L1( "Timeborder delta: {:.8e} s - {:.2f} fs \n", t_step, t_step * 1E15 );
    Log::L1( "Subiterator delta: {:.8e} s - {:.2f} fs \n", numerics_subiterator_stepsize, numerics_subiterator_stepsize * 1E15 );
    if ( numerics_phonon_approximation_order == 5 ) {
        Log::L1( "Timeborder delta path integral: {:.8e} s - {:.2f} ps\n", t_step_pathint, t_step_pathint * 1E12 );
    }
    Log::L1( "Time iterations (main loop) = {}\n\n", iterations_t_max );

    Log::Logger::wrapInBar( "G-Function Settings", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
    if ( input_correlation.size() > 0 ) {
        Log::L1( "Tau-grid resolution is {}\n", numerics_calculate_till_converged ? "to be determined." : fmt::format( "{}x{}", grid_values.size(), grid_values.size() ) );
        Log::L1( "Calculating:\n" );
        for ( auto &[name, mat] : input_correlation ) {
            Log::L1( " - {} on mode(s) ", name );
            for ( auto &mode : mat.string_v["Modes"] )
                Log::L1( "{} ", mode );
            Log::L1( "\n" );
        }
        Log::L1( "\n" );
        // Detector Stuff
        if ( input_conf["Detector"].numerical_v["time_center"].size() > 0 ) {
            Log::L1( "Using {} temporal detection windows with parameters:\n", input_conf["Detector"].numerical_v["time_center"].size() );
            for ( int i = 0; i < input_conf["Detector"].numerical_v["time_center"].size(); i++ ) {
                Log::L1( " - Temporal Detection Window {}:\n", i );
                Log::L1( " - - Center: {} - {} ps\n", input_conf["Detector"].numerical_v["time_center"][i], input_conf["Detector"].numerical_v["time_center"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Sigma: {} - {} ps\n", input_conf["Detector"].numerical_v["time_range"][i], input_conf["Detector"].numerical_v["time_range"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Power: {}\n", input_conf["Detector"].numerical_v["time_power_amplitude"][i] );
            }
        }
        if ( input_conf["Detector"].numerical_v["spectral_center"].size() > 0 ) {
            Log::L1( "Using {} spectral detection windows with parameters:\n", input_conf["Detector"].numerical_v["spectral_center"].size() );
            for ( int i = 0; i < input_conf["Detector"].numerical_v["spectral_center"].size(); i++ ) {
                Log::L1( " - Spectral Detection Window {}:\n", i );
                Log::L1( " - - Center: {} - {} eV\n", input_conf["Detector"].numerical_v["spectral_center"][i], input_conf["Detector"].numerical_v["spectral_center"][i].getSI( Parameter::UNIT_ENERGY_EV ) );
                Log::L1( " - - Sigma: {} - {} meV\n", input_conf["Detector"].numerical_v["spectral_range"][i], input_conf["Detector"].numerical_v["spectral_range"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                Log::L1( " - - Power: {}\n", input_conf["Detector"].numerical_v["spectral_power_amplitude"][i] );
                Log::L1( " - - FT Points: {}\n", input_conf["Detector"].numerical_v["spectral_number_points"][i] );
            }
        }

    } else {
        Log::L1( "Not using any G1 or G2 correlation functions.\n\n" );
    }
    Log::L1( "\n" );

    Log::Logger::wrapInBar( "Settings", Log::Logger::BAR_SIZE_HALF, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
    Log::L1( "Solver used: RK{}{}\n", numerics_rk_order, numerics_rk_order != 45 ? "" : fmt::format( " (Tolerance: {}, Stepdelta: {}, Steplimits: [{},{}])", inputstring_rk45_config, numerics_rk_stepdelta, numerics_rk_stepmin, numerics_rk_stepmax ) );
    if ( numerics_rk_order == 45 and numerics_phonon_nork45 )
        Log::L1( "Will NOT use RK45 for the phonon backwards integral!\n" );
    Log::L1( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
    Log::L1( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture ) ? "YES" : "NO" ) );
    Log::L1( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) ? "Analytic" : "Matrix Exponential" ) );
    Log::L1( "Threads used for primary calculations - {}\nThreads used for Secondary calculations - {}\nThreads used by Eigen: {}\n", numerics_maximum_secondary_threads, numerics_maximum_primary_threads, Eigen::nbThreads() );
    Log::L1( "Used scaling for parameters? - {}\n", ( scale_parameters ? std::to_string( scale_value ) : "no" ) );
    if ( p_phonon_T )
        Log::L1( "Cache Phonon Coefficient Matrices? - {}\n", ( numerics_use_saved_coefficients ? "Yes" : "No" ) );
    if ( numerics_interpolate_outputs )
        Log::L1( "WARNING: Temporal outputs are interpolated!\n" );
    if ( not numerics_use_function_caching )
        Log::L1( "NOT using function caching.\n" );
    if ( not numerics_use_saved_hamiltons )
        Log::L1( "NOT using Hamilton caching.\n" );
    Log::L1( "\n" );
    if ( not output_dict.empty() ) {
        Log::L1( "Additional Outputs:\n" );
        std::ranges::for_each( output_dict, []( const auto &el ) { Log::L1( " - {}\n", el ); } );
    }

    for ( int i = 0; i < logfilecounter.size(); i++ ) {
        Log::L1( "Logfile ident number {}: {}\n", i, logfilecounter[i] );
    }
    Log::L1( "\n" );

    Log::Logger::wrapInBar( "Program Log:", Log::Logger::BAR_SIZE_FULL, Log::Logger::LEVEL_2 );
    Log::L1( "\n" );
}