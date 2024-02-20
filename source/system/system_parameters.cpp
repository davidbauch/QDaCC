#include "system/parameters.h"

using namespace QDACC;

Parameters::Parameters( const std::vector<std::string> &arguments ) {
    // Parsing input:
    Timer &timer_parseInput = Timers::create( "Parsing parameters", true, false ).start();
    Log::Logger::wrapInBar( "Conversion of input variables", Log::BAR_SIZE_FULL, Log::LEVEL_2, Log::BAR_0 );
    Log::L2( "\n" );
    Log::L2( "[System] Parsing input variables...\n" );

    parse_input( arguments );
    timer_parseInput.end();
    Log::L2( "[System] Successful. Elapsed time is {}ms\n", timer_parseInput.getWallTime( Timers::MILLISECONDS ) );

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
    t_start = QDACC::CommandlineArguments::get_parameter<double>( "--time", "tstart" );
    t_end = QDACC::CommandlineArguments::get_parameter<double>( "--time", "tend" );
    t_step = QDACC::CommandlineArguments::get_parameter<double>( "--time", "tstep" );
    numerics_hard_t_max = QDACC::CommandlineArguments::get_parameter<double>( "--hardmax" );
    numerics_groundstate_string = QDACC::CommandlineArguments::get_parameter( "--groundstate" );

    // Runge Kutta Parameters
    numerics_rk_order = QDACC::CommandlineArguments::get_parameter<double>( "--rk", "rkorder" );
    inputstring_rk45_config = QDACC::CommandlineArguments::get_parameter( "--rk", "rktol" );
    numerics_rk_stepdelta = QDACC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepdelta" );
    numerics_rk_stepmin = QDACC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepmin" );
    numerics_rk_stepmax = QDACC::CommandlineArguments::get_parameter<double>( "--rk", "rkstepmax" );
    numerics_rk_usediscrete_timesteps = numerics_rk_stepdelta > 0 ? true : false;

    inputstring_electronic = QDACC::CommandlineArguments::get_parameter( "--S", "SE" );
    inputstring_photonic = QDACC::CommandlineArguments::get_parameter( "--S", "SO" );
    inputstring_pulse = QDACC::CommandlineArguments::get_parameter( "--S", "SP" );
    inputstring_chirp = QDACC::CommandlineArguments::get_parameter( "--S", "SC" );
    inputstring_spectrum = QDACC::CommandlineArguments::get_parameter( "--G", "GS" );
    inputstring_indist = QDACC::CommandlineArguments::get_parameter( "--G", "GI" );
    inputstring_conc = QDACC::CommandlineArguments::get_parameter( "--G", "GC" );
    inputstring_gfunc = QDACC::CommandlineArguments::get_parameter( "--G", "GF" );
    inputstring_timebins = QDACC::CommandlineArguments::get_parameter( "--G", "GT" );
    inputstring_wigner = QDACC::CommandlineArguments::get_parameter( "--G", "GW" );
    inputstring_raman = QDACC::CommandlineArguments::get_parameter( "--G", "GR" );
    inputstring_correlation_resolution = QDACC::CommandlineArguments::get_parameter( "--G", "grid" );
    inputstring_detector_time = QDACC::CommandlineArguments::get_parameter( "--detector", "temporalDetector" );
    inputstring_detector_spectral = QDACC::CommandlineArguments::get_parameter( "--detector", "spectralDetector" );
    inputstring_densitymatrix_config = QDACC::CommandlineArguments::get_parameter( "--DMconfig" );
    inputstring_SPconf = QDACC::CommandlineArguments::get_parameter( "--SPconfig" );

    p_omega_coupling = QDACC::CommandlineArguments::get_parameter<double>( "--system", "coupling" );
    p_omega_cavity_loss = QDACC::CommandlineArguments::get_parameter<double>( "--system", "kappa" );
    p_omega_pure_dephasing = QDACC::CommandlineArguments::get_parameter<double>( "--system", "gammapure" );
    p_omega_decay = QDACC::CommandlineArguments::get_parameter<double>( "--system", "gamma" );
    p_initial_state_s = QDACC::CommandlineArguments::get_parameter( "--R" );

    grid_resolution = QDACC::CommandlineArguments::get_parameter<int>( "--G", "gridres" );
    numerics_use_interactionpicture = not QDACC::CommandlineArguments::get_parameter_passed( "-noInteractionpic" );
    numerics_use_rwa = not QDACC::CommandlineArguments::get_parameter_passed( "-noRWA" );
    numerics_maximum_primary_threads = QDACC::CommandlineArguments::get_parameter<int>( "--Threads" );
    if ( numerics_maximum_primary_threads == -1 ) numerics_maximum_primary_threads = omp_get_max_threads();
    output_handlerstrings = QDACC::CommandlineArguments::get_parameter_passed( "-handler" );
    numerics_order_timetrafo = QDACC::CommandlineArguments::get_parameter_passed( "-timeTrafoMatrixExponential" ) ? QDACC::TransformationOrder::MatrixExponential : QDACC::TransformationOrder::Analytical;
    numerics_use_saved_coefficients = not QDACC::CommandlineArguments::get_parameter_passed( "-disableMatrixCaching" );
    numerics_use_saved_hamiltons = not QDACC::CommandlineArguments::get_parameter_passed( "-disableHamiltonCaching" );
    numerics_use_function_caching = not QDACC::CommandlineArguments::get_parameter_passed( "-disableFunctionCaching" );
    numerics_enable_saving_coefficients = false; // If true, even if any saving was disabled internally (not by the user), the matrices will still be cached.
    numerics_maximum_secondary_threads = ( !numerics_use_saved_coefficients || !QDACC::CommandlineArguments::get_parameter_passed( "-disableMainProgramThreading" ) ) ? numerics_maximum_primary_threads : 1;
    logfilecounter = QDACC::Misc::convertParam<double>( QDACC::String::splitline( QDACC::CommandlineArguments::get_parameter( "--lfc" ), ',' ) );
    numerics_interpolate_outputs = QDACC::CommandlineArguments::get_parameter_passed( "-interpolate" );
    s_numerics_interpolate = QDACC::CommandlineArguments::get_parameter( "--interpolateOrder" );

    // Phonon Parameters
    p_phonon_alpha = QDACC::CommandlineArguments::get_parameter<double>( "--phonons", "phononalpha" );
    p_phonon_ohm = QDACC::CommandlineArguments::get_parameter<double>( "--phonons", "phononohm" );
    p_phonon_wcutoff = QDACC::CommandlineArguments::get_parameter<double>( "--phonons", "phononwcutoff" );
    p_phonon_wcutoffdelta = QDACC::CommandlineArguments::get_parameter<double>( "--phonons", "phononwcutoffdelta" );
    p_phonon_tcutoff = QDACC::CommandlineArguments::get_parameter<double>( "--phonons", "phonontcutoff" );
    p_phonon_T = QDACC::CommandlineArguments::get_parameter<double>( "--phonons", "temperature" );
    numerics_phonon_approximation_order = (QDACC::PhononApproximation)QDACC::CommandlineArguments::get_parameter<int>( "--phonons", "phononorder" );
    numerics_phonon_approximation_markov1 = QDACC::CommandlineArguments::get_parameter_passed( "-noMarkov" ) ? 0 : 1; // First Markov
    numerics_phonon_nork45 = not QDACC::CommandlineArguments::get_parameter_passed( "-usePhononRK45" );               // Enables. RK45 for phonon backwards integral; use if detunings are low, otherwise expensive.
    p_phonon_adjust_rad = QDACC::CommandlineArguments::get_parameter<double>( "--phononAdjust", "pARad" );
    p_phonon_adjust_dep = QDACC::CommandlineArguments::get_parameter<double>( "--phononAdjust", "pAPure" );
    p_phonon_adjust_b = QDACC::CommandlineArguments::get_parameter<double>( "--phononAdjust", "pARescaling" );
    p_phonon_pure_dephasing = QDACC::Misc::convertParam<double>( "1mueV" );

    // Phonon Quantum Dot Paramters. These are used to overwrite the wcutoff and alpha values if needed.
    p_phonon_qd_de = QDACC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDDe" );
    p_phonon_qd_dh = QDACC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDDh" );
    p_phonon_qd_rho = QDACC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDrho" );
    p_phonon_qd_cs = QDACC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDcs" );
    p_phonon_qd_ratio = QDACC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDratio" );
    p_phonon_qd_ae = QDACC::CommandlineArguments::get_parameter<double>( "--quantumdot", "QDae" );

    // Path Integral Parameters
    p_phonon_nc = QDACC::CommandlineArguments::get_parameter<int>( "--pathintegral", "NC" );
    numerics_subiterator_stepsize = QDACC::CommandlineArguments::get_parameter<double>( "--pathintegral", "iteratorStepsize" );
    numerics_pathintegral_squared_threshold = QDACC::CommandlineArguments::get_parameter<double>( "--numericalpathintegral", "squaredThreshold" );
    numerics_pathintegral_sparse_prune_threshold = QDACC::CommandlineArguments::get_parameter<double>( "--numericalpathintegral", "sparsePruneThreshold" );
    numerics_pathintegral_dynamiccutoff_iterations_max = QDACC::CommandlineArguments::get_parameter<double>( "--numericalpathintegral", "cutoffADM" ); // FIXME ????
    numerics_pathintegral_sparse_to_dense_threshold = QDACC::CommandlineArguments::get_parameter<double>( "--numericalpathintegral", "denseTensorThreshold" );
    numerics_pathintegral_force_dense = QDACC::CommandlineArguments::get_parameter_passed( "-pathIntForceDense" );
    numerics_pathintegral_docutoff_propagator = QDACC::CommandlineArguments::get_parameter_passed( "-cutoffPropagator" );
    numerics_pathint_partially_summed = !QDACC::CommandlineArguments::get_parameter_passed( "-disablePSPath" );
    t_step_pathint = QDACC::CommandlineArguments::get_parameter<double>( "--pathintegral", "tstepPath" );
    numerics_pathintegral_use_qrt = QDACC::CommandlineArguments::get_parameter_passed( "-useQRT" );
    numerics_pathintegral_set_couplings_zero = QDACC::CommandlineArguments::get_parameter_passed( "-setCouplingsZero" );

    auto output_dict_vec = QDACC::String::splitline( QDACC::CommandlineArguments::get_parameter( "--output" ), ';' );
    // Move all elements into the set
    output_dict = std::set<std::string>( std::make_move_iterator( output_dict_vec.begin() ), std::make_move_iterator( output_dict_vec.end() ) );

    numerics_custom_expectation_values = QDACC::String::splitline( QDACC::CommandlineArguments::get_parameter( "--expv" ), ';' );

    kb = 1.3806488E-23;   // J/K, scaling needs to be for energy
    hbar = 1.0545718E-34; // J/s, scaling will be 1

    maxStates = 0;
    numerics_calculate_till_converged = false;
    numerics_main_direction_done = false;

    parse_system();

    working_directory = arguments.back();
}

void Parameters::pre_adjust_input() {
    Log::L2( "[System] Adjusting Pre-Mainloop Inputs...\n" );

    if ( output_handlerstrings ) Timers::toggleHandler();

    // Calculate/Recalculate some parameters:
    // Adjust pulse data
    // TODO: für complexe ampltiude is_imag in Parameters = true -> dann mit 1i multiplizeiren in get() funktion
    for ( auto &[name, mat] : input_pulse ) {
        mat.property_set["Amplitude"] = QDACC::Misc::convertParam<Parameter>( mat.string_v["Amplitude"] );
        // Set all optional parameters to default
        mat.property_set["ChirpRate"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );
        mat.property_set["GaussAmp"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 1.0 );
        mat.property_set["SUPERDelta"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );
        mat.property_set["SUPERFreq"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );
        mat.property_set["CutoffDelta"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );
        mat.property_set["Phase"] = std::vector<Parameter>( mat.string_v["Amplitude"].size(), 0.0 );

        for ( int i = 0; i < mat.string_v["Amplitude"].size(); i++ ) {
            const auto type_params = QDACC::String::splitline( mat.string_v["Type"][i], '+' );
            // Extract optional Parameters
            for ( const auto &type : type_params ) {
                if ( type.find( "cw" ) != std::string::npos )
                    mat.string_v["Type"][i] = "cw";
                else if ( type.find( "gauss" ) != std::string::npos )
                    mat.string_v["Type"][i] = "gauss";
                else if ( type.find( "sech" ) != std::string::npos )
                    mat.string_v["Type"][i] = "sech";
                else if ( type.find( "lorentz" ) != std::string::npos )
                    mat.string_v["Type"][i] = "lorentz";
                else {
                    std::string param = QDACC::String::splitline( type, '(' ).back();
                    param.pop_back();
                    if ( type.find( "chirped" ) != std::string::npos ) {
                        mat.property_set["ChirpRate"][i] = QDACC::Misc::convertParam<Parameter>( param );
                    } else if ( type.find( "cutoff" ) != std::string::npos ) {
                        mat.property_set["CutoffDelta"][i] = QDACC::Misc::convertParam<Parameter>( param );
                    } else if ( type.find( "super" ) != std::string::npos ) {
                        auto splitparam = QDACC::String::splitline( param, '_' );
                        mat.property_set["SUPERDelta"][i] = QDACC::Misc::convertParam<Parameter>( splitparam[0] );
                        mat.property_set["SUPERFreq"][i] = QDACC::Misc::convertParam<Parameter>( splitparam[1] );
                    } else if ( type.find( "exponent" ) != std::string::npos ) {
                        mat.property_set["GaussAmp"][i] = QDACC::Misc::convertParam<Parameter>( param );
                    } else if ( type.find( "phase" ) != std::string::npos ) {
                        mat.property_set["Phase"][i] = QDACC::Misc::convertParam<Parameter>( param ) * QDACC::Math::PI;
                    }
                }
            }
            if ( mat.string_v["Amplitude"][i].find( "pi" ) != std::string::npos ) {
                if ( mat.string_v["Type"][i] == "sech" ) {
                    mat.property_set["Amplitude"][i] = mat.property_set["Amplitude"][i] / mat.property_set["Width"][i];
                } else {
                    mat.property_set["Amplitude"][i] = mat.property_set["Amplitude"][i] * QDACC::Math::PI /
                                                       ( std::sqrt( 2.0 * QDACC::Math::PI * mat.property_set["Width"][i] *
                                                                    std::sqrt( std::pow( mat.property_set["ChirpRate"][i] / mat.property_set["Width"][i], 2.0 ) + std::pow( mat.property_set["Width"][i], 2.0 ) ) ) ) /
                                                       2.0; // https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.241306
                }
            }
        }
    }

    // Adjust Chirp
    for ( auto &[name, mat] : input_chirp ) {
        if ( mat.string["Type"].compare( "sine" ) != 0 ) {
            for ( long unsigned i = 0; i < mat.property_set["ddt"].size(); i++ ) mat.get( "ddt", i ) = mat.get( "ddt", i ) * 1E12;
        }
    }

    // Set interpolation order:
    auto orders = QDACC::String::splitline( s_numerics_interpolate, ',' );
    std::map<std::string, int> methods = { { "cubic", 1 }, { "monotone", 2 }, { "linear", 0 } };
    std::string method_time = orders.front();
    std::string method_tau = orders.size() > 1 ? orders.back() : "linear";
    numerics_interpolate_method_time = methods[method_time];
    numerics_interpolate_method_tau = methods[method_tau];

    // No phonon adjust if pathintegral is chosen
    if ( numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
        p_phonon_adjust_b = 0;
        p_phonon_adjust_rad = 0;
        p_phonon_adjust_dep = 0;
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
                J = w * hbar *
                    std::pow( p_phonon_qd_de * std::exp( -w * w * p_phonon_qd_ae * p_phonon_qd_ae / ( 4. * p_phonon_qd_cs * p_phonon_qd_cs ) ) -
                                  p_phonon_qd_dh * std::exp( -w * w * p_phonon_qd_ae / p_phonon_qd_ratio * p_phonon_qd_ae / p_phonon_qd_ratio / ( 4. * p_phonon_qd_cs * p_phonon_qd_cs ) ),
                              2. ) /
                    ( 4. * 3.1415 * 3.1415 * p_phonon_qd_rho * std::pow( p_phonon_qd_cs, 5. ) );
            integral += stepsize * ( J / std::tanh( hbar * w / 2.0 / kb / p_phonon_T ) );
        }
        if ( p_phonon_adjust_b > 0 ) {
            p_phonon_b = std::exp( -0.5 * integral );
            const auto new_value = p_phonon_b * p_phonon_adjust_b;
            Log::L2( "[System] Adjusting phononscaling from <B> = {} to <B> = {}\n", p_phonon_b, new_value );
            p_phonon_b = new_value;
        }
        if ( p_phonon_adjust_rad > 0 ) {
            const auto new_value = p_omega_decay * p_phonon_b * p_phonon_b * p_phonon_adjust_rad;
            Log::L2( "[System] Adjusting radiative decay from gamma = {} to <B>^2*gamma = \n", p_omega_decay, new_value );
            p_omega_decay = new_value; // TODO: unterschiedliche gammas für unterschiedliche B. macht aber auch kaum was.
        }
        if ( p_phonon_adjust_dep > 0 ) {
            const auto new_value = p_phonon_pure_dephasing * p_phonon_T * p_phonon_adjust_dep;
            Log::L2( "[System] Adjusting radiative decay from gamma = {} to <B>^2*gamma = \n", p_phonon_pure_dephasing, new_value );
            p_omega_pure_dephasing = new_value;
        }
        if ( numerics_rk_order >= 45 and numerics_use_saved_coefficients ) numerics_enable_saving_coefficients = true;
    }

    // Calculate minimum step necessary to resolve Rabi-oscillation if step=-1
    if ( t_step < 0 ) {
        Log::L2( "[System] Delta t is automatically determined.\n" );
        const auto N = 1000;
        if ( t_end > 0 )
            t_step = t_end / N;
        else
            t_step = 100E-15;
        // t_step = std::max( std::numeric_limits<double>::epsilon(), t_step );
        Log::L2( "[System] Delta t was set to {}\n", t_step );
    }

    // Calculate till converged
    if ( t_end < 0 ) {
        Log::L2( "[System] t_end is automatically determined.\n" );
        // If this is given, we calculate the t-direction until 99% ground state poulation is reached after any pulses.
        numerics_calculate_till_converged = true;
        for ( auto &[name, mat] : input_pulse ) {
            for ( auto &t : mat.property_set["Center"] ) {
                t_end = std::max( t_end.get(), 2.0 * t );
            }
        }
        if ( t_end < 0 ) t_end = std::max<double>( 1E-12, 10.0 * std::max<double>( t_step, t_step_pathint ) );
        Log::L2( "[System] Calculate till at least {} and adjust accordingly to guarantee convergence. The matrix index used is {}\n", t_end, numerics_groundstate );
    }

    // Disable Hamilton Caching and RK45 for PI
    if ( numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
        // numerics_use_saved_hamiltons = false;
        // Log::L2( "[System] Disabled Caching of Hamilton matrices" );
        if ( numerics_rk_order > 5 ) {
            numerics_rk_order = 4;
            Log::L2( "[System] Adjusted RK order to {}\n", numerics_rk_order );
        }
    }

    // Automatically determin subiterator stepsize
    if ( numerics_subiterator_stepsize < 0 ) {
        if ( numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
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
    if ( Log::Logger::max_log_level() == Log::LEVEL_3 ) {
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
    if ( t_end >= 0 and ( numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ? ( t_step_pathint > 0 ) : ( t_step > 0 ) ) ) {
        if ( iterations_t_max < 1 ) {
            iterations_t_max = (int)std::ceil( ( t_end - t_start ) / ( numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ? t_step_pathint : t_step ) );
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
    if ( input_correlation_resolution.contains( "Modified" ) and grid_values.size() > 1 ) {
        Log::L2( "[System-Prameters] Grid was not again modified because the gridinput is fixed by the user.\n" );
        return;
    }
    if ( numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ? ( t_step_pathint > 0 ) : ( t_step > 0 ) ) {
        input_correlation_resolution["Standard"].property_set["Time"] = { t_end };
        double time_delta = numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ? t_step_pathint : Parameter( ( t_end - t_start ) / ( 1. * grid_resolution ) );
        input_correlation_resolution["Standard"].property_set["Delta"] = { time_delta };
        auto &settings = input_correlation_resolution.contains( "Modified" ) ? input_correlation_resolution["Modified"] : input_correlation_resolution["Standard"];
        Log::L2( "[System] Initial Grid Timestep is {}.\n", settings.property_set["Delta"].front() );
        double t_t = 0;
        int current = 0;
        grid_values.clear();
        grid_steps.clear();
        grid_value_indices.clear();
        grid_values.emplace_back( t_start );
        grid_value_indices[t_start] = 0;
        Log::L2( "[System] Initial Timestep Limit is {} at a timestep of {}.\n", settings.property_set["Time"][current], settings.property_set["Delta"][current] );
        while ( t_t < t_end ) {
            if ( t_t > settings.property_set["Time"][current] and current < settings.property_set["Time"].size() ) {
                current++;
                Log::L2( "[System] New Timestep Limit is {} at a timestep of {}.\n", settings.property_set["Time"][current], settings.property_set["Delta"][current] );
            }
            grid_steps.emplace_back( settings.property_set["Delta"][current] );
            t_t += grid_steps.back();
            grid_values.emplace_back( t_t );
            grid_value_indices[t_t] = grid_values.size() - 1;
        }
        Log::L2( "[System] Setting correlation grid resolution to {0}x{0} for a t_end = {1}{2}\n", grid_values.size(), t_end, input_correlation_resolution.contains( "Modified" ) ? " using a modified grid." : "" );
    } else {
        Log::L2( "[System] Not setting time vector because timestep is negative!\n" );
    }
}

void Parameters::parse_system() {
    // Generate the input variables for the electronic system:
    for ( auto levels = QDACC::String::splitline( inputstring_electronic, ';' ); const std::string &level : levels ) {
        auto conf = QDACC::String::splitline( level, ':' );
        universal_config conf_s;
        conf_s.property["Energy"] = QDACC::Misc::convertParam<Parameter>( conf[1] );           // Energy
        conf_s.string_v["CoupledTo"] = QDACC::String::splitline( conf[2], ',' );               // Coupled to Levels
        conf_s.property["DecayScaling"] = QDACC::Misc::convertParam<Parameter>( conf[3] );     // Decay Scaling, Per Mode
        conf_s.property["DephasingScaling"] = QDACC::Misc::convertParam<Parameter>( conf[4] ); // Dephasing Scaling
        conf_s.property["PhononCoupling"] = QDACC::Misc::convertParam<Parameter>( conf[5] );   // Phonon Coupling
        input_electronic[conf[0]] = conf_s;
    }
    for ( auto cavities = QDACC::String::splitline( inputstring_photonic, ';' ); const std::string &cavity : cavities ) {
        auto conf = QDACC::String::splitline( cavity, ':' );
        universal_config conf_s;
        conf_s.property["Energy"] = QDACC::Misc::convertParam<Parameter>( conf[1] );                                               // Energy
        conf_s.property["MaxPhotons"] = QDACC::Misc::convertParam<Parameter>( conf[2] );                                           // Maximum Photons
        conf_s.string_v["CoupledTo"] = QDACC::String::splitline( conf[3], ',' );                                                   // Coupled to Transitions
        conf_s.property_set["CouplingScaling"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[4], ',' ) ); // Coupling Scaling, per transition INTO cavity
        // conf_s.property_set["BackCouplingScaling"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[5], ',' ) ); // BackCoupling Scaling, per transition from cavity back into the electronic system
        conf_s.property["DecayScaling"] = QDACC::Misc::convertParam<Parameter>( conf[5] ); // Decay Scaling, for all transitions
        input_photonic[conf[0]] = conf_s;
    }

    // --SP 'p:GX:1pi,5pi:1.5eV,1eV:4ps,2ps:20ps,35ps:gauss:'
    // Type is cw or superposition (chained with '+', e.g. 'gauss+chirped(1E-24)' ) of gauss,cutoff,chirped(rate),super(delta),exponent(exponent),phase(angle)
    // p:TYPE:...parameters...
    int pindex = 0;
    for ( auto pulses = QDACC::String::splitline( inputstring_pulse, ';' ); const std::string &pulse : pulses ) {
        auto conf = QDACC::String::splitline( pulse, ':' );
        universal_config conf_s;
        conf_s.string_v["CoupledTo"] = QDACC::String::splitline( conf[1], ',' );                                             // Coupled to Transitions
        conf_s.string_v["Amplitude"] = QDACC::String::splitline( conf[2], ',' );                                             // Pulse Amp
        conf_s.property_set["Frequency"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[3], ',' ) ); // Frequency
        conf_s.property_set["Width"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[4], ',' ) );     // Width
        conf_s.property_set["Center"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[5], ',' ) );    // Center
        conf_s.string_v["Type"] = QDACC::String::splitline( conf[6], ',' );                                                  // Type
        // For counting purposes:
        conf_s.property["PulseIndex"] = pindex;
        pindex += 2;
        input_pulse[conf[0]] = conf_s;
    }
    for ( auto chirps = QDACC::String::splitline( inputstring_chirp, ';' ); const std::string &chirp : chirps ) {
        auto conf = QDACC::String::splitline( chirp, ':' );
        universal_config conf_s;
        conf_s.string_v["CoupledTo"] = QDACC::String::splitline( conf[1], ',' );                                             // Coupled to Transitions
        conf_s.property_set["AmpFactor"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[2], ',' ) ); // Amplitude Scaling for coupled_to
        conf_s.property_set["Amplitude"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[3], ',' ) ); // Amplitudes
        conf_s.property_set["Times"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[4], ',' ) );     // "Times"
        if ( conf[5].find( "," ) != std::string::npos or conf.size() > 6 ) {
            conf_s.property_set["ddt"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[5], ',' ) ); // "d/dt"
            conf_s.string["Type"] = conf[6];                                                                               // Type
        } else {
            conf_s.string["Type"] = conf[5];                                                                 // Type
            conf_s.property_set["ddt"] = std::vector<Parameter>( conf_s.property_set["Times"].size(), 0.0 ); // "d/dt"
        }
        input_chirp[conf[0]] = conf_s;
    }
    for ( const std::string &spectrum : QDACC::String::splitline( inputstring_spectrum, ';' ) ) {
        auto conf = QDACC::String::splitline( spectrum, ':' );
        universal_config conf_s;
        conf_s.string_v["Modes"] = QDACC::String::splitline( conf[0], ',' ); // Modes to calculate Spectrum for. Single modes can again be split with "+", meaning a+b;a to calculate for a+b and a seperately
        conf_s.property_set["Center"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[1], ',' ) );                                                                                     // Center
        conf_s.property_set["Range"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[2], ',' ) );                                                                                      // Range
        conf_s.property_set["resW"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[3], ',' ) );                                                                                       // Resolution for w
        conf_s.property_set["Order"] = conf.size() > 4 ? QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[4], ',' ) ) : std::vector<Parameter>( conf_s.property_set["Range"].size(), 1 ); // Order (1 or 2)?
        conf_s.string_v["Normalize"] = conf.size() > 5 ? QDACC::String::splitline( conf[5], ',' ) : std::vector<std::string>( conf_s.property_set["Range"].size(), "False" );                                 // Normalize?
        input_correlation["Spectrum"].emplace_back( conf_s );
    }
    for ( std::string &indist : QDACC::String::splitline( inputstring_indist, ';' ) ) {
        auto conf = QDACC::String::splitline( indist, ':' );
        universal_config conf_s;
        conf_s.string_v["Modes"] = QDACC::String::splitline( conf[0], ',' ); // Modes to calculate Indistinguishgability for. Single modes can again be split with "+", meaning a+b;a to calculate for a+b and a seperately
        input_correlation["Indist"].emplace_back( conf_s );
    }
    for ( std::string &conc : QDACC::String::splitline( inputstring_conc, ';' ) ) {
        auto conf = QDACC::String::splitline( conc, ':' );
        universal_config conf_s;
        conf_s.string_v["Modes"] = QDACC::String::splitline( conf[0], ',' ); // Modes to calculate Concurrence for
        conf_s.string_v["Order"] = conf.size() > 1 ? QDACC::String::splitline( conf[1], ',' ) : std::vector<std::string>( conf_s.string_v["Modes"].size(), "full" );
        if ( conf.size() > 2 ) {
            conf_s.string_v["Method"] = conf.size() > 2 ? QDACC::String::splitline( conf[2], ',' ) : std::vector<std::string>( conf_s.string_v["Modes"].size(), "wootters" ); // Can be "wooters" or "seidelmann"
        }
        if ( conf.size() > 3 ) {
            // Experimental: Calculate spectrum for all matrix entries. Also outputs the non-normalized matrices
            conf_s.property_set["Center"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[1], ',' ) ); // Center
            conf_s.property_set["Range"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[2], ',' ) );  // Range
            conf_s.property_set["resW"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[3], ',' ) );   // Resolution for w
        }
        input_correlation["Conc"].emplace_back( conf_s );
    }
    for ( std::string &g_func : QDACC::String::splitline( inputstring_gfunc, ';' ) ) {
        auto conf = QDACC::String::splitline( g_func, ':' );
        auto n = conf.size();
        universal_config conf_s;
        conf_s.string_v["Modes"] = QDACC::String::splitline( conf[0], ',' );                                                                                    // Modes to calculate G1/G2 functions for
        conf_s.property_set["Order"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[1], ',' ) );                                        // 1 or 2 or 3
        conf_s.string_v["Integrated"] = n > 2 ? QDACC::String::splitline( conf[2], ',' ) : std::vector<std::string>( conf_s.string_v["Modes"].size(), "time" ); // time,matrix,both for false/true/both
        input_correlation["GFunc"].emplace_back( conf_s );
    }
    for ( std::string &g_func : QDACC::String::splitline( inputstring_timebins, ';' ) ) {
        auto conf = QDACC::String::splitline( g_func, ':' );
        auto n = conf.size();
        universal_config conf_s;
        conf_s.string_v["Modes"] = QDACC::String::splitline( conf[0], ',' );                                             // Modes to calculate G1/G2 functions for
        conf_s.property_set["Order"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[1], ',' ) ); // 2 or 3
        conf_s.property_set["Start"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[2], ',' ) );
        conf_s.property_set["BinLength"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[3], ',' ) );
        conf_s.property_set["Deph"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[4], ',' ) );
        conf_s.string_v["Integrated"] = n > 5 ? QDACC::String::splitline( conf[5], ',' ) : std::vector<std::string>( conf_s.string_v["Modes"].size(), "time" ); // time,matrix,both for false/true/both
        input_correlation["GFuncTimeBins"].emplace_back( conf_s );
    }
    for ( std::string &wigner : QDACC::String::splitline( inputstring_wigner, ';' ) ) {
        auto conf = QDACC::String::splitline( wigner, ':' );
        auto n = conf.size();
        universal_config conf_s;
        conf_s.string_v["Modes"] = QDACC::String::splitline( conf[0], ',' );                                                                                            // Modes to calculate Wigner function for
        conf_s.property_set["X"] = QDACC::Misc::convertParam<Parameter>( QDACC::String::splitline( conf[1], ',' ) );                                                    // -X to X
        conf_s.property_set["Y"] = QDACC::Misc::convertParam<Parameter>( n > 2 ? QDACC::String::splitline( conf[2], ',' ) : QDACC::String::splitline( conf[1], ',' ) ); // -Y to Y
        conf_s.property_set["Res"] = QDACC::Misc::convertParam<Parameter>( n > 3 ? QDACC::String::splitline( conf[3], ',' ) : std::vector<std::string>( conf_s.property_set["X"].size(), "100" ) ); // Resolution
        conf_s.property_set["Skip"] = QDACC::Misc::convertParam<Parameter>( n > 4 ? QDACC::String::splitline( conf[4], ',' ) : std::vector<std::string>( conf_s.property_set["X"].size(), "1" ) );  // Skips in t-direction
        input_correlation["Wigner"].emplace_back( conf_s );
    }
    for ( std::string &raman : QDACC::String::splitline( inputstring_raman, ';' ) ) {
        auto conf = QDACC::String::splitline( raman, ':' );
        auto n = conf.size();
        universal_config conf_s;
        conf_s.string_v["SourceModes"] = QDACC::String::splitline( conf[0], ',' );
        conf_s.string_v["RamanMode"] = QDACC::String::splitline( conf[1], ',' );
        conf_s.string_v["OpMode"] = QDACC::String::splitline( conf[2], ',' );
        conf_s.string_v["PMode"] = QDACC::String::splitline( conf[3], ',' );
        input_correlation["Raman"].emplace_back( conf_s );
    }
    // Correlation Grid
    if ( std::ranges::find( inputstring_correlation_resolution, ':' ) != inputstring_correlation_resolution.end() ) {
        auto single = QDACC::String::splitline( inputstring_correlation_resolution, ';' );
        universal_config conf_s;
        std::vector<Parameter> times, dts;
        for ( const auto &el : single ) {
            Log::L2( "[System-Parameters] Adding Grid Subspace {}\n", el );
            auto cur = QDACC::String::splitline( el, ':' );
            times.emplace_back( QDACC::Misc::convertParam<Parameter>( cur[0] ) );
            dts.emplace_back( QDACC::Misc::convertParam<Parameter>( cur[1] ) );
        }
        conf_s.property_set["Time"] = times;
        conf_s.property_set["Delta"] = dts;
        input_correlation_resolution["Modified"] = conf_s;
    }
    {
        universal_config conf_s;
        conf_s.property_set["Time"] = { t_end };
        conf_s.property_set["Delta"] = { t_step };
        input_correlation_resolution["Standard"] = conf_s;
    }
    {
        universal_config conf_s;
        if ( inputstring_SPconf == "inherit" and conf_s.property_set["Center"].size() > 0 ) {
            Log::L2( "[System-Parameters] Inheriting Pulse Fourier Configuration from Spectrum.\n" );
            conf_s.property["Center"] = conf_s.property_set["Center"][0];
            conf_s.property["Range"] = conf_s.property_set["Range"][0];
            conf_s.property["Res"] = conf_s.property_set["resW"][0];
            conf_s.property["dt"] = t_end / 200;
        } else if ( inputstring_SPconf.size() > 0 ) {
            Log::L2( "[System-Parameters] Generating Pulse Fourier Configuration from Parameters using {}.\n", inputstring_SPconf );
            auto conf = QDACC::String::splitline( inputstring_SPconf, ':' );
            conf_s.property["Center"] = conf.size() > 0 ? QDACC::Misc::convertParam<Parameter>( conf[0] ) : Parameter( 0.0 );
            conf_s.property["Range"] = conf.size() > 1 ? QDACC::Misc::convertParam<Parameter>( conf[1] ) : Parameter( 0.0 );
            conf_s.property["Res"] = conf.size() > 2 ? QDACC::Misc::convertParam<Parameter>( conf[2] ) : Parameter( 0.0 );
            conf_s.property["dt"] = conf.size() > 3 ? QDACC::Misc::convertParam<Parameter>( conf[3] ) : Parameter( t_end / 200 );
            input_conf["PulseConf"] = conf_s;
        }
    }
    // Detector Input.
    {
        if ( inputstring_detector_time != "none" ) {
            universal_config conf_s;
            Log::L2( "[System-Parameters] Setting up detector using {}...\n", inputstring_detector_time );
            auto configs = QDACC::String::splitline( inputstring_detector_time, ';' );
            conf_s.property_set["time_center"] = {};
            conf_s.property_set["time_range"] = {};
            conf_s.property_set["time_amplitude"] = {};
            conf_s.property_set["time_power_amplitude"] = {};
            conf_s.string_v["time_mode"] = {}; // Mode can be either G1,G2 or a specific name like G1-hbd-hb
            for ( const auto &conf_sep : configs ) {
                auto conf = QDACC::String::splitline( conf_sep, ':' );
                conf_s.property_set["time_range"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[0] ) );
                conf_s.property_set["time_center"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[1] ) );
                conf_s.property_set["time_amplitude"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[2] ) );
                conf_s.property_set["time_power_amplitude"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[3] ) );
                conf_s.string_v["time_mode"].emplace_back( conf[4] );
                Log::L2( "[System-Parameters] Adding Temporal Detector mask using center = {}, range = {}, ampltitude {} and power_amp = {} on mode(s) {}.\n", conf[1], conf[0], conf[2], conf[3], conf[4] );
            }
            input_conf["DetectorTime"] = conf_s;
        }
        if ( inputstring_detector_spectral != "none" ) {
            universal_config conf_s;
            Log::L2( "[System-Parameters] Setting up detector using {}...\n", inputstring_detector_spectral );
            auto configs = QDACC::String::splitline( inputstring_detector_spectral, ';' );
            conf_s.property_set["spectral_range"] = {};
            conf_s.property_set["spectral_center"] = {};
            conf_s.property_set["spectral_amplitude"] = {};
            conf_s.property_set["spectral_power_amplitude"] = {};
            conf_s.property_set["spectral_number_points"] = {};
            conf_s.string_v["spectral_mode"] = {};
            for ( const auto &conf_sep : configs ) {
                auto conf = QDACC::String::splitline( conf_sep, ':' );
                conf_s.property_set["spectral_range"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[0] ) );
                conf_s.property_set["spectral_center"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[1] ) );
                conf_s.property_set["spectral_amplitude"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[2] ) );
                conf_s.property_set["spectral_power_amplitude"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[3] ) );
                conf_s.property_set["spectral_number_points"].emplace_back( QDACC::Misc::convertParam<Parameter>( conf[4] ) );
                conf_s.string_v["spectral_mode"].emplace_back( conf[5] );
                Log::L2(
                    "[System-Parameters] Adding Spectral Detector mask using center = {}, range = {}, amplitude {} power_amp = {} and {} points for the Fourier Transformation on "
                    "mode(s) {}.\n",
                    conf[1], conf[0], conf[2], conf[3], conf[4], conf[5] );
            }
            input_conf["DetectorSpectral"] = conf_s;
        }
    }
    {
        universal_config conf_s;
        Log::L2( "[System-Parameters] Setting up density matrix output using {}...\n", inputstring_densitymatrix_config );
        auto conf = QDACC::String::splitline( inputstring_densitymatrix_config, ':' );
        conf_s.string["output_mode"] = conf[0];
        conf_s.string["interaction_picture"] = conf.size() > 1 ? conf[1] : "schroedinger";
        input_conf["DMconfig"] = conf_s;
    }
    // RK 45 Tolerance Vector
    {
        if ( std::ranges::find( inputstring_rk45_config.begin(), inputstring_rk45_config.end(), ':' ) == inputstring_rk45_config.end() ) {
            numerics_rk_tol.emplace_back( 1.0, QDACC::Misc::convertParam<double>( inputstring_rk45_config ) );
            Log::L2( "[System-Parameters] Set fixed RK45 Tolerance to {}.\n", std::get<1>( numerics_rk_tol.back() ) );
        } else {
            Log::L2( "[System-Parameters] Setting up multiple tolerances...\n" );
            for ( const auto &partial : QDACC::String::splitline( inputstring_rk45_config, ';' ) ) {
                auto tuple = QDACC::String::splitline( partial, ':' );
                numerics_rk_tol.emplace_back( std::make_tuple( QDACC::Misc::convertParam<double>( tuple.front() ), QDACC::Misc::convertParam<double>( tuple.back() ) ) );
                Log::L2( "[System-Parameters] Set RK45 Tolerance to {} (from {}).\n", std::get<1>( numerics_rk_tol.back() ), partial );
            }
        }
    }
}

std::string _get_interpolator_name( int index ) {
    const std::vector names = { "Linear", "Quintic Hermite", "Cubic Hermite" };
    return names[index];
}

void Parameters::log( const Dense &initial_state_vector_ket ) {
    Log::L1(
        "         _______                   _____                    _____     _____          \n        /::\\    \\                 /\\    \\                  /\\    \\   /\\    "
        "\\         \n       /::::\\    \\               /::\\    \\                /::\\____\\ /::\\    \\        \n      /::::::\\    \\             /::::\\    \\              "
        "/:::/    //::::\\    \\       \n     /::::::::\\    \\           /::::::\\    \\            /:::/    //::::::\\    \\      \n    /:::/~~\\:::\\    \\         "
        "/:::/\\:::\\    \\          /:::/    //:::/\\:::\\    \\     \n   /:::/    \\:::\\    \\       /:::/  \\:::\\    \\        /:::/    //:::/  \\:::\\    \\    \n  /:::/    "
        "/ \\:::\\    \\     /:::/    \\:::\\    \\      /:::/    //:::/    \\:::\\    \\   \n /:::/____/   \\:::\\____\\   /:::/    / \\:::\\    \\    /:::/    //:::/    / "
        "\\:::\\    \\  \n|:::|    |     |:::|    | /:::/    /   \\:::\\ ___\\  /:::/    //:::/    /   \\:::\\    \\ \n|:::|____|     |:::|____|/:::/____/     \\:::|    "
        "|/:::/____//:::/____/     \\:::\\____\\\n \\:::\\   _\\___/:::/    / \\:::\\    \\     /:::|____|\\:::\\    \\\\:::\\    \\      \\::/    /\n  \\:::\\ |::| /:::/    /   "
        "\\:::\\    \\   /:::/    /  \\:::\\    \\\\:::\\    \\      \\/____/ \n   \\:::\\|::|/:::/    /     \\:::\\    \\ /:::/    /    \\:::\\    \\\\:::\\    \\             \n "
        "   \\::::::::::/    /       \\:::\\    /:::/    /      \\:::\\    \\\\:::\\    \\            \n     \\::::::::/    /         \\:::\\  /:::/    /        \\:::\\    "
        "\\\\:::\\    \\           \n      \\::::::/    /           \\:::\\/:::/    /          \\:::\\    \\\\:::\\    \\          \n       \\::::/____/             \\::::::/    "
        "/            \\:::\\    \\\\:::\\    \\         \n        |::|    |               \\::::/    /              \\:::\\____\\\\:::\\____\\        \n        |::|____|         "
        "       \\::/____/                \\::/    / \\::/    /        \n         ~~                       ~~                       \\/____/   \\/____/         \n                 "
        "                                                                    \n" );
    Log::Logger::wrapInBar( "System Parameters" );
    Log::L1( "Version: {} ({})\n\n", GLOBAL_PROGRAM_VERSION, GLOBAL_PROGRAM_LASTCHANGE );
    #ifdef USE_SPARSE_MATRIX
    Log::L1("Matrix Type: Sparse\n");
    #else
    Log::L1("Matrix Type: Dense\n");
    #endif

    Log::Logger::wrapInBar( "Electronic Configuration", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    for ( auto &[name, mat] : input_electronic ) {
        Log::L1( "Electronic State: {} with energy {:.8e} Hz - {:.8e} eV - {:.8f} nm\n", name, mat.property["Energy"], mat.property["Energy"].getSI( Parameter::UNIT_ENERGY_EV ),
                 mat.property["Energy"].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
        for ( auto i = 0; i < mat.string_v["CoupledTo"].size(); i++ ) {
            if ( mat.string_v["CoupledTo"][i].front() == '-' ) break;
            Log::L1( " - Coupled to {} with transition energy {:.8e} Hz - {:.8e} eV\n", mat.string_v["CoupledTo"][i], std::abs( mat.property["Energy"] - input_electronic[mat.string_v["CoupledTo"][i]].property["Energy"] ),
                     std::abs( mat.property["Energy"].getSI( Parameter::UNIT_ENERGY_EV ) - input_electronic[mat.string_v["CoupledTo"][i]].property["Energy"].getSI( Parameter::UNIT_ENERGY_EV ) ) );
        }
        Log::L1( " - Decay Scaling: {}\n", mat.property["DecayScaling"] );
        Log::L1( " - Dephasing Scaling: {}\n", mat.property["DephasingScaling"] );
    }
    Log::L1( "\n" );
    if ( input_electronic.size() > 0 ) {
        Log::Logger::wrapInBar( "Photonic Configuration", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
        for ( auto &[name, mat] : input_photonic ) {
            Log::L1( "Cavity: {} with energy {:.8e} Hz - {:.8e} eV - {:.8f} nm\n", name, mat.property["Energy"], mat.property["Energy"].getSI( Parameter::UNIT_ENERGY_EV ),
                     mat.property["Energy"].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
            for ( auto i = 0; i < mat.string_v["CoupledTo"].size(); i++ ) {
                auto [from, to] = QDACC::String::split_pair( mat.string_v["CoupledTo"][i], transition_delimiter );
                double purcell = ( 2.0 * p_omega_coupling * p_omega_coupling * mat.property_set["CouplingScaling"][i] * mat.property_set["CouplingScaling"][i] / p_omega_cavity_loss / p_omega_decay ) *
                                 ( p_omega_cavity_loss * p_omega_cavity_loss /
                                   ( std::pow( mat.property["Energy"] - ( input_electronic[to].property["Energy"] - input_electronic[from].property["Energy"] ), 2 ) + p_omega_cavity_loss * p_omega_cavity_loss ) );
                Log::L1( " - Coupled to electronic transition {} (Purcell Enhancement: {}) \n", mat.string_v["CoupledTo"][i], purcell );
                Log::L1( " - - Coupling Scaling: {}\n", mat.property_set["CouplingScaling"][i] );
            }
            Log::L1( " - Cavity Q-Factor: {:.0f}\n", mat.property["Energy"] / p_omega_cavity_loss );
            Log::L1( " - Decay Scaling: {}\n", mat.property["DecayScaling"] );
        }
        Log::L1( "\n" );
    }

    Log::Logger::wrapInBar( "Coupling Parameters", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Coupling strengh g: {:.8e} Hz - {:.8} mueV\n", p_omega_coupling, p_omega_coupling.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Photon loss rate k: {:.8e} Hz - {:.8} mueV\n", p_omega_cavity_loss, p_omega_cavity_loss.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Atomic dephasing rate gamma_pure: {:.8e} Hz - {:.8} mueV\n", p_omega_pure_dephasing, p_omega_pure_dephasing.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "RAD rate gamma: {:.8e} Hz - {:.8} mueV\n\n", p_omega_decay, p_omega_decay.getSI( Parameter::UNIT_ENERGY_MUEV ) );

    Log::Logger::wrapInBar( "Initial System Parameters", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    // Log::L1( "Initial state rho0 = [{}]\n\n", initial_state_vector_ket.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ) );
    Log::L1( "Matrix Initial State = {}\n", p_initial_state_s );
    Log::L1( "Matrix Groundstate Index = {}\n\n", numerics_groundstate );

    Log::Logger::wrapInBar( "Pulse", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_pulse.size() > 0 ) {
        for ( auto &[name, mat] : input_pulse ) {
            Log::L1( "Pulse {}:\n", name );
            Log::L1( " - Coupled to Transitions: " );
            for ( auto &mode : mat.string_v["CoupledTo"] ) Log::L1( "{} ", mode );
            Log::L1( "\n" );
            for ( int i = 0; i < mat.property_set["Amplitude"].size(); i++ ) {
                Log::L1( " - Single Pulse:\n" );
                Log::L1( " - - Amplitude: {} Hz - {:.8} mueV\n", mat.property_set["Amplitude"][i], mat.property_set["Amplitude"][i].getSI( Parameter::UNIT_ENERGY_MUEV ) );
                Log::L1( " - - Frequency: {} Hz - {:.8} eV - {:.8} nm\n", mat.property_set["Frequency"][i], mat.property_set["Frequency"][i].getSI( Parameter::UNIT_ENERGY_EV ),
                         mat.property_set["Frequency"][i].getSI( Parameter::UNIT_WAVELENGTH_NM ) );
                Log::L1( " - - Width: {} s - {:.8} ps - {:.8} mueV bandwith\n", mat.property_set["Width"][i], mat.property_set["Width"][i].getSI( Parameter::UNIT_TIME_PS ),
                         Parameter( 1.0 / mat.property_set["Width"][i] ).getSI( Parameter::UNIT_ENERGY_MUEV ) );
                Log::L1( " - - Center: {} s - {:.8} ps\n", mat.property_set["Center"][i], mat.property_set["Center"][i].getSI( Parameter::UNIT_TIME_PS ) );
                if ( QDACC::Math::abs2( mat.property_set["ChirpRate"][i] != 0.0 ) ) Log::L1( " - - Chirp: {}\n", mat.property_set["ChirpRate"][i] );
                if ( QDACC::Math::abs2( mat.property_set["SUPERDelta"][i] != 0.0 ) ) {
                    Log::L1( " - - SUPER Amplitude: {} - {} meV\n", mat.property_set["SUPERDelta"][i], mat.property_set["SUPERDelta"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                    Log::L1( " - - SUPER Frequency: {} - {} meV\n", mat.property_set["SUPERFreq"][i], mat.property_set["SUPERFreq"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                }
                if ( QDACC::Math::abs2( mat.property_set["Phase"][i] != 0.0 ) ) Log::L1( " - - Phase: {}pi\n", mat.property_set["Phase"][i] / QDACC::Math::PI );
                Log::L1( " - - Type: {}{}\n", mat.string_v["Type"][i], mat.string_v["Type"][i] == "gauss" ? std::format( " (Gaussian Amplitude: {})", mat.property_set["GaussAmp"][i] ) : "" );
            }
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using any pulses to exite or drive the system.\n\n" );
    }
    Log::Logger::wrapInBar( "Chirp", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_chirp.size() > 0 ) {
        for ( auto &[name, mat] : input_chirp ) {
            Log::L1( "Chirp: {}\n", name );
            Log::L1( " - Coupled to States:\n" );
            for ( int i = 0; i < mat.string_v["CoupledTo"].size(); i++ ) Log::L1( " - - {} with scaling {}\n", mat.string_v["CoupledTo"][i], mat.property_set["AmpFactor"][i] );
            if ( mat.string.at( "Type" ) != "sine" ) {
                for ( int i = 0; i < mat.property_set["Amplitude"].size(); i++ ) {
                    Log::L1( " - Chirp Point {}:\n", i );
                    Log::L1( " - - Amplitude: {} Hz - {} meV\n", mat.property_set["Amplitude"][i], mat.property_set["Amplitude"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                    Log::L1( " - - Time: {} s - {} ps\n", mat.property_set["Times"][i], mat.property_set["Times"][i].getSI( Parameter::UNIT_TIME_PS ) );
                    Log::L1( " - - Derivative DDT: {}\n", mat.property_set["ddt"][i] );
                }
            } else {
                Log::L1( " - Sine Chirp with {} individual frequencies:\n", mat.property_set.at( "Amplitude" ).size() );
                for ( int i = 0; i < mat.property_set["Amplitude"].size(); i++ ) {
                    Log::L1( " - - Sine {}:\n", i );
                    Log::L1( " - - - Amplitude: {} Hz - {} meV\n", mat.property_set["Amplitude"][i], mat.property_set["Amplitude"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                    Log::L1( " - - - Phase Shift: {} s - {} ps\n", mat.property_set["Times"][i], mat.property_set["Times"][i].getSI( Parameter::UNIT_TIME_PS ) );
                    Log::L1( " - - - Sine Frequency: {} Hz - {} meV\n", mat.property_set["ddt"][i], mat.property_set["ddt"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                }
            }
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using any electronic chirps to shift the system.\n\n" );
    }

    Log::Logger::wrapInBar( "Phonons", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( p_phonon_T >= 0 ) {
        std::vector<std::string> approximations = { "Transformation integral via d/dt chi = -i/hbar*[H,chi] + d*chi/dt onto interaction picture chi(t-tau)",
                                                    "Transformation Matrix U(t,tau)=exp(-i/hbar*H_DQ_L(t)*tau) onto interaction picture chi(t-tau)",
                                                    "No Transformation, only interaction picture chi(t-tau)",
                                                    "Analytical Lindblad formalism",
                                                    "Mixed",
                                                    "Path Integral" };
        Log::L1( "Temperature = {} K\n", p_phonon_T );
        Log::L1( "Ohm Power = {}\n", p_phonon_ohm );
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
        Log::L1( "Transformation approximation used: {} - {}\n", (int)numerics_phonon_approximation_order, approximations.at( (int)numerics_phonon_approximation_order ) );
        // Pathintegral
        if ( numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
            Log::L1( " - Path Integral Settings:\n" );
            Log::L1( " - Backsteps NC: {}\n", p_phonon_nc );
            Log::L1( " - Iterator Stepsize: {}\n", numerics_subiterator_stepsize );
            Log::L1( " - Thresholds: Squared({}), SparsePrune({}), CutoffIterations({}), PropagatorMapping({})\n", numerics_pathintegral_squared_threshold, numerics_pathintegral_sparse_prune_threshold,
                     numerics_pathintegral_dynamiccutoff_iterations_max, numerics_pathintegral_docutoff_propagator );
            Log::L1( " - Used partially summed algorithm?: {}\n", numerics_pathint_partially_summed ? "Yes" : "No" );
            Log::L1( " - Will set the coupling elements to zero for correlation calculations?: {}\n", numerics_pathintegral_set_couplings_zero ? "Yes" : "No" );
            Log::L1( " - Will use QRT for correlation calculations?: {}\n", numerics_pathintegral_use_qrt ? "Yes" : "No" );
        }
        Log::L1( "\n" );
    } else {
        Log::L1( "Not using phonons\n\n" );
    }

    Log::Logger::wrapInBar( "Numerical Parameters" );
    Log::L1( "\n" );
    Log::Logger::wrapInBar( "Time", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Timeborder start: {:.8e} s - {:.2f} ps\n", t_start, t_start * 1E12 );
    Log::L1( "Timeborder end: {:.8e} s - {:.2f} ps\n", t_end, t_end * 1E12 );
    if ( numerics_calculate_till_converged ) {
        Log::L1( "Variable time end at 99.9\% convergence using dm index {} -> {}\n", numerics_groundstate_string, numerics_groundstate );
        if ( numerics_hard_t_max > 0 ) Log::L1( "Or until t > t_hard_max = {}\n", numerics_hard_t_max );
    }
    Log::L1( "Timeborder delta: {:.8e} s - {:.2f} fs \n", t_step, t_step * 1E15 );
    Log::L1( "Subiterator delta: {:.8e} s - {:.2f} fs \n", numerics_subiterator_stepsize, numerics_subiterator_stepsize * 1E15 );
    if ( numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
        Log::L1( "Timeborder delta path integral: {:.8e} s - {:.2f} ps\n", t_step_pathint, t_step_pathint * 1E12 );
    }
    Log::L1( "Output Timeframe: {}", input_conf["DMconfig"].string["interaction_picture"] );
    Log::L1( "\n" );
    Log::Logger::wrapInBar( "G-Function Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( input_correlation.size() > 0 ) {
        // TODO: print grid resolution later
        Log::L1( "Tau-grid resolution is {}\n", numerics_calculate_till_converged ? "to be determined." : std::format( "{}x{}", grid_values.size(), grid_values.size() ) );
        Log::L1( "Interpolator used: {}\n", _get_interpolator_name( numerics_interpolate_method_tau ) );
        Log::L1( "Calculating:\n" );
        for ( auto &[name, all_correlations] : input_correlation )
            for ( auto &mat : all_correlations ) {
                Log::L1( " - {}:\n", name );
                Log::L1( mat.toString( " = " ) );
                /*
                Log::L1( " - {} on mode(s):\n", name );
                for ( auto i = 0; i < mat.string_v["Modes"].size(); i++ ) {
                    const auto &mode = mat.string_v["Modes"][i];
                    if ( name == "Conc" ) {
                        const auto &order = mat.string_v["Order"][i];
                        const auto &method = mat.string_v["Method"][i];
                        Log::L1( " - - {} using matrix evaluation order: '{}'\n", mode, order );
                        Log::L1( " - - {} using method: '{}'\n", mode, method );
                    } else {
                        Log::L1( " - - {}\n", mode );
                    }
                }
                */
            }
        Log::L1( "\n" );
        // Detector Stuff
        if ( input_conf["DetectorTime"].property_set["time_center"].size() > 0 ) {
            Log::L1( "Using {} temporal detection windows with parameters:\n", input_conf["DetectorTime"].property_set["time_center"].size() );
            for ( int i = 0; i < input_conf["DetectorTime"].property_set["time_center"].size(); i++ ) {
                Log::L1( " - Temporal Detection Window {}:\n", i );
                Log::L1( " - - Center: {} - {} ps\n", input_conf["DetectorTime"].property_set["time_center"][i], input_conf["DetectorTime"].property_set["time_center"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Sigma: {} - {} ps\n", input_conf["DetectorTime"].property_set["time_range"][i], input_conf["DetectorTime"].property_set["time_range"][i].getSI( Parameter::UNIT_TIME_PS ) );
                Log::L1( " - - Amplitude: {}\n", input_conf["DetectorTime"].property_set["time_amplitude"][i] );
                Log::L1( " - - Power: {}\n", input_conf["DetectorTime"].property_set["time_power_amplitude"][i] );
                Log::L1( " - - Mode: {}\n", input_conf["DetectorTime"].string_v["time_mode"][i] );
            }
        }
        if ( input_conf["DetectorSpectral"].property_set["spectral_center"].size() > 0 ) {
            Log::L1( "Using {} spectral detection windows with parameters:\n", input_conf["DetectorSpectral"].property_set["spectral_center"].size() );
            for ( int i = 0; i < input_conf["DetectorSpectral"].property_set["spectral_center"].size(); i++ ) {
                Log::L1( " - Spectral Detection Window {}:\n", i );
                Log::L1( " - - Center: {} - {} eV\n", input_conf["DetectorSpectral"].property_set["spectral_center"][i], input_conf["DetectorSpectral"].property_set["spectral_center"][i].getSI( Parameter::UNIT_ENERGY_EV ) );
                Log::L1( " - - Sigma: {} - {} meV\n", input_conf["DetectorSpectral"].property_set["spectral_range"][i], input_conf["DetectorSpectral"].property_set["spectral_range"][i].getSI( Parameter::UNIT_ENERGY_MEV ) );
                Log::L1( " - - Amplitude: {}\n", input_conf["DetectorSpectral"].property_set["spectral_amplitude"][i] );
                Log::L1( " - - Power: {}\n", input_conf["DetectorSpectral"].property_set["spectral_power_amplitude"][i] );
                Log::L1( " - - FT Points: {}\n", input_conf["DetectorSpectral"].property_set["spectral_number_points"][i] );
                Log::L1( " - - Mode: {}\n", input_conf["DetectorSpectral"].string_v["spectral_mode"][i] );
            }
        }
    } else {
        Log::L1( "Not using any G1 or G2 correlation functions.\n\n" );
    }
    Log::L1( "\n" );

    Log::Logger::wrapInBar( "Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Solver used: RK{}{}\n", numerics_rk_order,
             numerics_rk_order != 45 ? "" : std::format( " (Tolerance: {}, Stepdelta: {}, Steplimits: [{},{}])", inputstring_rk45_config, numerics_rk_stepdelta, numerics_rk_stepmin, numerics_rk_stepmax ) );
    if ( numerics_rk_order == 45 and numerics_phonon_nork45 ) Log::L1( "Will NOT use RK45 for the phonon backwards integral!\n" );
    Log::L1( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
    Log::L1( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture ) ? "YES" : "NO" ) );
    Log::L1( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == QDACC::TransformationOrder::Analytical ) ? "Analytic" : "Matrix Exponential" ) );
    Log::L1( "Threads used for primary calculations - {}\nThreads used for Secondary calculations - {}\nThreads used by Eigen: {}\n", numerics_maximum_secondary_threads, numerics_maximum_primary_threads, Eigen::nbThreads() );
    if ( p_phonon_T ) Log::L1( "Cache Phonon Coefficient Matrices? - {}\n", ( numerics_use_saved_coefficients ? "Yes" : "No" ) );
    if ( numerics_interpolate_outputs ) {
        Log::L1( "WARNING: Temporal outputs are interpolated! Interpolator used: {}\n", _get_interpolator_name( numerics_interpolate_method_time ) );
    }
    if ( not numerics_use_function_caching ) Log::L1( "NOT using function caching.\n" );
    if ( not numerics_use_saved_hamiltons ) Log::L1( "NOT using Hamilton caching.\n" );
    Log::L1( "\n" );
    if ( not output_dict.empty() ) {
        Log::L1( "Additional Outputs:\n" );
        std::ranges::for_each( output_dict, []( const auto &el ) { Log::L1( " - {}\n", el ); } );
    }

    for ( int i = 0; i < logfilecounter.size(); i++ ) {
        Log::L1( "Logfile ident number {}: {}\n", i, logfilecounter[i] );
    }
    Log::L1( "\n" );

    Log::Logger::wrapInBar( "Program Log:", Log::BAR_SIZE_FULL, Log::LEVEL_2 );
    Log::L1( "\n" );
}