#include "system/system.h"

using namespace QDACC;

System::System( const std::vector<std::string> &input ) {
    Log::L2( "[System] Creating System Class\n" );
    // Initialize all subclasses with the input vector
    parameters = Parameters( input );
    operatorMatrices = OperatorMatrices( parameters );
    operatorMatrices.output_operators( parameters );
    // Initialize FileOutput with the current system
    FileOutput::init( parameters, operatorMatrices );
    // Initialize / Adjust the remaining system class
    Timer &timer_systeminit = Timers::create( "System Initialization", true /*statistics*/, false /*print*/ ).start();
    Log::L2( "[System] Initialization...\n" );
    if ( !init_system() ) {
        Log::L2( "[System] Initialization failed! Exitting program...\n" );
        Log::Logger::close();
        exit( EXIT_FAILURE );
    }
    Log::L2( "[System] Successful! Elapsed time is {}ms\n", timer_systeminit.getWallTime( Timers::MILLISECONDS ) );
    // Log the operator matrix base
    parameters.log( operatorMatrices.initial_state_vector_ket );
    timer_systeminit.end();
}

bool System::init_system() {
    // Single chirp for single atomic level
    for ( auto &[name, chirpinputs] : parameters.input_chirp ) {
        chirp.emplace_back( chirpinputs, parameters );
        chirp.back().to_file( "chirp_" + name, false /*complex*/, parameters.output_dict.contains( "chirpf" ) /*spectrum*/ );
    }

    // Arbitrary number of pulses onto single atomic level.
    for ( auto &[name, pulseinputs] : parameters.input_pulse ) {
        pulse.emplace_back( pulseinputs, parameters );
        pulse.back().to_file( "pulse_" + name, true /*complex*/, parameters.output_dict.contains( "pulsef" ) /*spectrum*/ );
    }

    // if ( parameters.numerics_use_saved_coefficients )
    //     savedCoefficients.reserve( parameters.numerics_saved_coefficients_max_size );

    // Output Phonon functions if phonons are active
    if ( parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
        initialize_path_integral_functions();
    } else {
        initialize_polaron_frame_functions();
    }

    // Output Phonon J
    if ( parameters.output_dict.contains( "phononJ" ) ) {
        auto &file = FileOutput::add_file( "phonon_spectral" );
        file << std::format( "omega\tJ(omega)\n" );
        for ( double w = parameters.p_phonon_wcutoffdelta; w < 10.0 * parameters.p_phonon_wcutoff; w += parameters.p_phonon_wcutoffdelta ) {
            file << std::format( "{}\t{}\n", w, std::real( dgl_phonons_spectral_density( w ) ) );
        }
    }

    // Time Transformation
    timeTrafoMatrix = MatrixMain( ( Dense( 1.0i * operatorMatrices.H_0 ).exp() ).sparseView() );
    // Check time trafo
    MatrixMain ttrafo = MatrixMain( ( Dense( 1.0i * operatorMatrices.H_0 * 505E-12 ).exp() * operatorMatrices.H_used * Dense( -1.0i * operatorMatrices.H_0 * 505E-12 ).exp() ).sparseView() );
    MatrixMain temp = dgl_timetrafo( operatorMatrices.H_used, 505E-12 );
    if ( double error = std::abs( 1.0 - Dense( temp ).sum() / Dense( ttrafo ).sum() ); error >= 1E-8 ) {
        Log::L1( "[System] Unitary timetransformation error is {:.3f}\%!\n", error * 100.0 );
    }

    // Precalculate Cav Decay, Rad Decay and Dephasing matrices.
    // The DGL then uses the left- and rightside matrices to calculate the Lindbladians
    if ( parameters.p_omega_cavity_loss != 0.0 ) {
        for ( const auto &[transition, data] : operatorMatrices.ph_transitions ) {
            if ( data.direction == 1 ) continue; // Not the right operator
            std::string mode = data.to;
            auto &params = parameters.input_photonic[mode];
            if ( params.property["DecayScaling"] == 0.0 ) continue;
            Log::L2( "[System] Precalculating cavity decay for mode {} with scaling {}\n", mode, params.property["DecayScaling"] );
            const auto factor = 0.5 * parameters.p_omega_cavity_loss * params.property["DecayScaling"];
            cache_cav_decay_left.emplace_back( factor * operatorMatrices.ph_transitions[mode + "b"].hilbert );
            cache_cav_decay_right.emplace_back( operatorMatrices.ph_transitions[mode + "bd"].hilbert );
            cache_cav_decay_rightleft.emplace_back( cache_cav_decay_right.back() * cache_cav_decay_left.back() );
        }
    }
    if ( parameters.p_omega_decay != 0.0 ) {
        for ( const auto &[transition, data] : operatorMatrices.el_transitions ) {
            if ( data.direction == 1 ) continue; // Not the right operator
            const std::string &state = data.to;
            auto &params = parameters.input_electronic[state];
            if ( params.property["DecayScaling"] == 0.0 ) continue;
            Log::L2( "[System] Precalculating radiative decay for state {} with scaling {}\n", state, params.property["DecayScaling"] );
            const std::string &trans_transposed = data.name_transposed;
            cache_rad_decay_left.emplace_back( 0.5 * parameters.p_omega_decay * params.property["DecayScaling"] * operatorMatrices.el_transitions[transition].hilbert );
            cache_rad_decay_right.emplace_back( operatorMatrices.el_transitions[trans_transposed].hilbert );
            cache_rad_decay_rightleft.emplace_back( cache_rad_decay_right.back() * cache_rad_decay_left.back() );
        }
    }
    if ( parameters.p_omega_pure_dephasing != 0.0 ) {
        for ( const auto &[name_a, data_a] : operatorMatrices.el_states ) {
            cache_dephasing_left.emplace_back( -0.5 * parameters.input_electronic[name_a].property["DephasingScaling"] * parameters.p_omega_pure_dephasing * operatorMatrices.el_states[name_a].hilbert );
            MatrixMain right_side = operatorMatrices.zero;
            for ( const auto &[name_b, data_b] : operatorMatrices.el_states ) { // TODO: dephasing über el transitions machen.
                if ( name_a == name_b ) continue;                               // Not the right operator
                if ( parameters.input_electronic[name_b].property["DephasingScaling"] * parameters.input_electronic[name_a].property["DephasingScaling"] == 0.0 ) continue;
                Log::L2( "[System] Precalculating pure dephasing for states {} and {}\n", name_a, name_b );
                right_side += parameters.input_electronic[name_b].property["DephasingScaling"] * operatorMatrices.el_states[name_b].hilbert;
                //-= 0.5 * parameters.input_electronic[name_b].property["DephasingScaling"] * parameters.input_electronic[name_a].property["DephasingScaling"] * parameters.p_omega_pure_dephasing * data_a.hilbert * rho * data_b.hilbert;
            }
            cache_dephasing_right.emplace_back( right_side );
        }
    }
    Log::L2( "[System] Added {} matrices to the cavity decay Lindbladian\n", cache_cav_decay_left.size() );
    Log::L2( "[System] Added {} matrices to the radiative decay Lindbladian\n", cache_rad_decay_left.size() );
    Log::L2( "[System] Added {} matrices to the pure dephasing Lindbladian\n", cache_dephasing_left.size() );
    return true;
}

MatrixMain System::dgl_runge_function( const MatrixMain &rho, const MatrixMain &H, const double t, std::vector<QDACC::SaveState> &past_rhos ) {
    MatrixMain ret = -1.0i * dgl_kommutator( H, rho );

    // Photon Loss
    if ( parameters.p_omega_cavity_loss != 0.0 ) {
        for ( int i = 0; i < cache_cav_decay_left.size(); i++ ) {
            ret += 2.0 * cache_cav_decay_left[i] * rho * cache_cav_decay_right[i] - cache_cav_decay_rightleft[i] * rho - rho * cache_cav_decay_rightleft[i];
        }
    }
    // Radiative Decay
    if ( parameters.p_omega_decay != 0.0 ) {
        for ( int i = 0; i < cache_rad_decay_left.size(); i++ ) {
            ret += 2.0 * cache_rad_decay_left[i] * rho * cache_rad_decay_right[i] - cache_rad_decay_rightleft[i] * rho - rho * cache_rad_decay_rightleft[i];
        }
    }
    // Electronic Dephasing
    if ( parameters.p_omega_pure_dephasing != 0.0 ) {
        for ( int i = 0; i < cache_dephasing_left.size(); i++ ) {
            ret += cache_dephasing_left[i] * rho * cache_dephasing_right[i];
        }
    }

    // Phonons
    if ( parameters.numerics_phonon_approximation_order != QDACC::PhononApproximation::PathIntegral && parameters.p_phonon_T >= 0.0 ) {
        ret += dgl_phonons_pmeq( rho, t, past_rhos );
    }
    return ret;
}

MatrixMain System::dgl_timetrafo( MatrixMain ret, const double t ) {
    if ( not parameters.numerics_use_interactionpicture ) return ret;
    // QDACC::TransformationOrder::Analytical
    if ( parameters.numerics_order_timetrafo == QDACC::TransformationOrder::Analytical ) {
#ifdef USE_SPARSE_MATRIX
        for ( int k = 0; k < ret.outerSize(); ++k ) {
            for ( Sparse::InnerIterator it( ret, k ); it; ++it ) {
                // Convert row/col into respective photon numbers / atomic state
                auto i = it.row();
                auto j = it.col();
                it.valueRef() *= std::exp( t * operatorMatrices.timetrafo_cachematrix( i, j ) );
            }
        }
#else
        for ( int i = 0; i < ret.rows(); i++ )
            for ( int j = 0; j < ret.cols(); j++ ) {
                ret.coeffRef( i, j ) *= std::exp( t * operatorMatrices.timetrafo_cachematrix( i, j ) );
            }
#endif
    }
    // QDACC::TransformationOrder::MatrixExponential
    else if ( parameters.numerics_order_timetrafo == QDACC::TransformationOrder::MatrixExponential ) {
        MatrixMain U = MatrixMain( ( Dense( 1.i * operatorMatrices.H_0 * t ).exp() ).sparseView() );
        ret = ( U * ret * U.adjoint() ).eval(); // aliasing?
    }
    return ret;
}

MatrixMain System::dgl_chirp( const double t ) {
    if ( chirp.empty() ) return operatorMatrices.zero;
    return operatorMatrices.chirp_mat.back() * chirp.back().get( t, !parameters.numerics_use_function_caching );
}

MatrixMain System::dgl_pulse( const double t ) {
    if ( pulse.empty() ) return operatorMatrices.zero;
    MatrixMain ret = operatorMatrices.zero;
    for ( int i = 0; i < pulse.size(); i++ ) {
        auto p = pulse[i].get( t, not parameters.numerics_use_function_caching );
        if ( parameters.numerics_use_rwa )
            ret += operatorMatrices.pulse_mat[2 * i + 1] * p + operatorMatrices.pulse_mat[2 * i] * std::conj( p );
        else
            ret += operatorMatrices.pulse_mat[2 * i + 1] * ( p + std::conj( p ) ) + operatorMatrices.pulse_mat[2 * i] * ( p + std::conj( p ) ); // TODO: verifizieren! trafo selber ausrechnen (in overleaf packen)
    }
    return ret;
}

MatrixMain System::dgl_get_hamilton( const double t ) { return dgl_timetrafo( operatorMatrices.H_used + dgl_pulse( t ) + dgl_chirp( t ), t ); }

void System::calculate_expectation_values( const std::vector<QDACC::SaveState> &rhos, Timer &evalTimer ) {
    // Output expectation Values
    double t_pre = rhos.front().t;
    for ( auto &tup : rhos ) {
        auto rho = tup.mat;
        double t = tup.t;
        double dt = t - t_pre;
        t_pre = t;
        std::string el_out = std::format( "{:.5e}", t );
        std::string ph_out = std::format( "{:.5e}", t );
        std::string el_em = "";
        std::string ph_em = "";
        // Output State Populations and Calculate Emission Probabilities from all radiative transitions:
        for ( auto &[mode, state] : operatorMatrices.ph_states ) {
            double expval = std::real( dgl_expectationvalue<MatrixMain>( rho, state.hilbert, t ) );
            ph_out = std::format( "{:}\t{:.6e}", ph_out, expval );
            if ( parameters.p_omega_cavity_loss > 0.0 ) {
                emission_probabilities[mode] += expval * dt;
                ph_em = std::format( "{:}\t{:.6e}", ph_em, parameters.p_omega_cavity_loss * parameters.input_photonic[mode].property["DecayScaling"] * emission_probabilities[mode] );
            }
        }
        for ( auto &[mode, state] : operatorMatrices.el_states ) {
            double expval = std::real( dgl_expectationvalue<MatrixMain>( rho, state.hilbert, t ) );
            el_out = std::format( "{:}\t{:.6e}", el_out, expval );
            if ( parameters.p_omega_decay > 0.0 and parameters.input_electronic[mode].property["DecayScaling"] != 0.0 ) {
                emission_probabilities[mode] += expval * dt;
                el_em = std::format( "{:}\t{:.6e}", el_em, parameters.p_omega_decay * parameters.input_electronic[mode].property["DecayScaling"] * emission_probabilities[mode] );
            }
        }
        // Output
        if ( not operatorMatrices.ph_states.empty() ) FileOutput::get_file( "photonic" ) << std::format( "{:}{:}\n", ph_out, ph_em );
        if ( not operatorMatrices.el_states.empty() ) FileOutput::get_file( "electronic" ) << std::format( "{:}{:}\n", el_out, el_em );

        if ( parameters.input_conf["DMconfig"].string["output_mode"] != "none" ) {
            FileOutput::get_file( "densitymatrix" ) << std::format( "{:.5e}", t ); //, rho.nonZeros(), rho.rows() * rho.cols() - rho.nonZeros() );
            if ( parameters.input_conf["DMconfig"].string["interaction_picture"] != "int" ) rho = dgl_timetrafo( rho, t );
            if ( parameters.input_conf["DMconfig"].string["output_mode"] == "full" ) {
                for ( int i = 0; i < parameters.maxStates; i++ )
                    for ( int j = 0; j < parameters.maxStates; j++ ) {
                        FileOutput::get_file( "densitymatrix" ) << std::format( "\t{:.10e}", std::real( rho.coeff( i, j ) ) );
                    }
                for ( int i = 0; i < parameters.maxStates; i++ )
                    for ( int j = 0; j < parameters.maxStates; j++ ) {
                        FileOutput::get_file( "densitymatrix" ) << std::format( "\t{:.10e}", std::imag( rho.coeff( i, j ) ) );
                    }
            } else
                for ( int j = 0; j < parameters.maxStates; j++ ) FileOutput::get_file( "densitymatrix" ) << std::format( "\t{:.5e}", std::real( rho.coeff( j, j ) ) );
            FileOutput::get_file( "densitymatrix" ) << "\n";
        }

        // Explicit Photon Matrices
        if ( parameters.output_dict.contains( "photons" ) ) {
            bool int_trafo = parameters.input_conf["DMconfig"].string["interaction_picture"] == "int";
            for ( auto &[mode, state] : operatorMatrices.ph_states ) {
                // Calculate Partial Trace
                const auto base_index = state.base;
                auto current = partial_trace( int_trafo ? rho : dgl_timetrafo( rho, t ), base_index );
                // Output Partial Trace Matrix to File
                FileOutput::get_file( "photons_" + mode ) << std::format( "{:.5e}", t );
                for ( int i = 0; i < current.rows(); i++ )
                    for ( int j = 0; j < current.cols(); j++ ) FileOutput::get_file( "photons_" + mode ) << std::format( "\t{:.5e}", std::real( current.coeff( i, j ) ) );
                // Output Emission Probabilities for diagonal elements
                if ( parameters.p_omega_cavity_loss > 0.0 )
                    for ( int i = 1; i < current.rows(); i++ ) {
                        std::string key = "photons_" + mode + std::to_string( i );
                        emission_probabilities[key] += std::real( current.coeff( i, i ) ) * dt;
                        FileOutput::get_file( "photons_" + mode ) << std::format( "\t{:.5e}", parameters.p_omega_cavity_loss * parameters.input_photonic[mode].property["DecayScaling"] * emission_probabilities[key] );
                    }

                FileOutput::get_file( "photons_" + mode ) << "\n";
            }
        }
        // Output Custom Expectation values
        if ( not operatorMatrices.numerics_custom_expectation_values_operators.empty() ) {
            FileOutput::add_file( "custom_expectation_values" ) << std::format( "{:.5e}", t );
            for ( auto &mat : operatorMatrices.numerics_custom_expectation_values_operators ) {
                FileOutput::get_file( "custom_expectation_values" ) << std::format( "\t{:.5e}", std::real( dgl_expectationvalue<MatrixMain>( rho, mat, t ) ) );
            }
            FileOutput::get_file( "custom_expectation_values" ) << std::format( "\n" );
        }
        evalTimer.iterate();
    }
}

bool System::exit_system( const int failure ) {
    for ( auto &p : pulse ) p.log();
    for ( auto &c : chirp ) c.log();
    Log::L2( "[System-PME] Coefficients: Attempts w/r: {}, Write: {}, Calc: {}, Read: {}, Interpolate/Read: {}.\n", track_getcoefficient_calcattempt, track_getcoefficient_write, track_getcoefficient_calculate,
             track_getcoefficient_read, track_getcoefficient_read_interpolated );
    Log::L2( "[System] Number of approx +/- adjustments: {}\n", globaltries );
    Log::L1( "[System] Maximum RAM used: {} MB\n", getPeakRSS() / 1024 / 1024 );
    FileOutput::close_all();
    return true;
}

bool System::trace_valid( MatrixMain &rho, double t_hit, bool force ) {
    double trace = std::real( get_trace<Scalar>( rho ) );
    parameters.trace.emplace_back( trace );
    if ( trace < 0.99 || trace > 1.01 || force ) {
        if ( force ) std::cout << std::format( "[System] {} Error -> trace check failed at t = {} with trace(rho) = {}\n", QDACC::Message::Prefix::PERROR, t_hit, trace );
        auto &fp_trace = FileOutput::add_file( "trace" );
        for ( int i = 0; i < (int)parameters.trace.size() && parameters.t_step * 1.0 * i < t_hit; i++ ) {
            fp_trace << std::format( "{:.10e} {:.15e}\n", parameters.t_step * 1.0 * ( i + 1 ), parameters.trace.at( i ) );
        }
        FileOutput::close_file( "trace" );
        return false;
    } else {
        return true;
    }
}