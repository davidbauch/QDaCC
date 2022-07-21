#include "system/system.h"

System::System( const std::vector<std::string> &input ) {
    // Set Name of this system.
    name = "Generic Electronic System";
    Log::L2( "[System] Creating System Class for '{}'\n", name );
    // Initialize all subclasses with the input vector
    parameters = Parameters( input );
    operatorMatrices = OperatorMatrices( parameters );
    operatorMatrices.output_operators( parameters );
    // Initialize FileOutput with the current system
    FileOutput::init( parameters, operatorMatrices );
    // Initialize / Adjust the remaining system class
    terminate_message = QDLC::Message::global_normaltermination;
    Timer &timer_systeminit = Timers::create( "System Initialization", true, false );
    Log::L2( "[System] Initialization...\n" );
    timer_systeminit.start();
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
    for ( auto &[mode, p] : parameters.input_chirp ) {
        Chirp::Inputs chirpinputs( parameters.t_start, parameters.t_end, parameters.t_step, p.string["Type"], parameters.numerics_rk_order );
        chirpinputs.add( p.numerical_v["Times"], p.numerical_v["Amplitude"], p.numerical_v["ddt"] );
        chirp.push_back( { chirpinputs } );
    }

    // Arbitrary number of pulses onto single atomic level.
    for ( auto &[name, pulseinputs] : parameters.input_pulse ) {
        pulse.push_back( { pulseinputs, parameters } );
    }
    // pulse_V = Pulse( pulseinputs_V );
    if ( pulse.size() > 0 ) {
        Pulse::fileOutput( parameters.working_directory + "pulse.txt", pulse, parameters.t_start, parameters.t_end, parameters.t_step );
    }
    if ( chirp.size() > 0 ) {
        chirp.back().fileOutput( parameters.working_directory + "chirp.txt" );
    }

    // if ( parameters.numerics_use_saved_coefficients )
    //     savedCoefficients.reserve( parameters.numerics_saved_coefficients_max_size );

    // Output Phonon functions if phonons are active
    if ( parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
        initialize_path_integral_functions();
    } else {
        initialize_polaron_frame_functions();
    }
    // Time Transformation
    timeTrafoMatrix = ( Dense( 1.0i * operatorMatrices.H_0 ).exp() ).sparseView(); //.pruned();
    // Check time trafo
    Sparse ttrafo = ( Dense( 1.0i * operatorMatrices.H_0 * 505E-12 ).exp() * operatorMatrices.H_used * Dense( -1.0i * operatorMatrices.H_0 * 505E-12 ).exp() ).sparseView();
    Sparse temp = dgl_timetrafo( operatorMatrices.H_used, 505E-12 );
    double error = std::abs( 1.0 - Dense( temp ).sum() / Dense( ttrafo ).sum() );
    if ( error >= 1E-8 ) {
        Log::L1( "[System] Unitary timetransformation error is {:.3f}\%!\n", error * 100.0 );
    }
    return true;
}

Sparse System::dgl_runge_function( const Sparse &rho, const Sparse &H, const double t, std::vector<QDLC::SaveState> &past_rhos ) {
    Sparse ret = -1.0i * dgl_kommutator( H, rho );
    Sparse loss = Sparse( rho.rows(), rho.cols() );
    // Photon Loss
    if ( parameters.p_omega_cavity_loss != 0.0 )
        for ( const auto &[transition, data] : operatorMatrices.ph_transitions ) {
            if ( data.direction == 1 )
                continue;
            std::string mode = data.to;
            auto &params = parameters.input_photonic[mode];
            loss += 0.5 * parameters.p_omega_cavity_loss * params.numerical["DecayScaling"] * dgl_lindblad( rho, operatorMatrices.ph_transitions[mode + "b"].hilbert, operatorMatrices.ph_transitions[mode + "bd"].hilbert );
        }
    // Radiative Decay
    if ( parameters.p_omega_decay != 0.0 )
        for ( const auto &[transition, data] : operatorMatrices.el_transitions ) {
            if ( data.direction == 1 )
                continue;
            const std::string &state = data.to;
            auto &params = parameters.input_electronic[state];
            const std::string &trans_transposed = data.name_transposed;
            loss += 0.5 * parameters.p_omega_decay * params.numerical["DecayScaling"] * dgl_lindblad( rho, operatorMatrices.el_transitions[transition].hilbert, operatorMatrices.el_transitions[trans_transposed].hilbert );
        }
    // Electronic Dephasing
    if ( parameters.p_omega_pure_dephasing != 0.0 )
        for ( const auto &[name_a, data_a] : operatorMatrices.el_states ) {
            for ( const auto &[name_b, data_b] : operatorMatrices.el_states ) { // TODO: dephasing über el transitions machen.
                if ( name_a == name_b )
                    continue;
                loss -= 0.5 * parameters.input_electronic[name_b].numerical["DephasingScaling"] * parameters.input_electronic[name_a].numerical["DephasingScaling"] * parameters.p_omega_pure_dephasing * data_a.hilbert * rho * data_b.hilbert;
            }
        }

    if ( parameters.numerics_phonon_approximation_order != PHONON_PATH_INTEGRAL && parameters.p_phonon_T >= 0.0 ) {
        Sparse phonons = dgl_phonons_pmeq( rho, t, past_rhos );
        Log::L3( "[System] Lindblad Equation:\n{}\nLosses:\n{}\nPhonon Part:\n{}\n", Dense( ret ).format( operatorMatrices.output_format ), Dense( loss ).format( operatorMatrices.output_format ), Dense( phonons ).format( operatorMatrices.output_format ) );
        ret += phonons; // dgl_phonons_pmeq( rho, t, past_rhos );
    }

    return ret + loss;
}

Sparse System::dgl_timetrafo( Sparse ret, const double t ) {
    if ( parameters.numerics_use_interactionpicture ) {
        // TIMETRANSFORMATION_ANALYTICAL
        if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) {
            // std::vector<Eigen::Triplet<Scalar>> ret_v;
            for ( int k = 0; k < ret.outerSize(); ++k ) {
                for ( Sparse::InnerIterator it( ret, k ); it; ++it ) {
                    // Convert row/col into respective photon numbers / atomic state
                    auto i = it.row();
                    auto j = it.col();
                    // ret.coeffRef( i, j ) *= std::exp( t * operatorMatrices.timetrafo_cachematrix( i, j ) ); // it.value() = ?
                    it.valueRef() *= std::exp( t * operatorMatrices.timetrafo_cachematrix( i, j ) ); // it.value() = ?
                }
            }
        }
        // TIMETRANSFORMATION_MATRIXEXPONENTIAL
        else if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_MATRIXEXPONENTIAL ) {
            Sparse U = ( Dense( 1.0i * operatorMatrices.H_0 * t ).exp() ).sparseView();
            return U * ret * U.adjoint();
        }
    }
    return ret;
}

Sparse System::dgl_chirp( const double t ) {
    if ( chirp.empty() )
        return Sparse( parameters.maxStates, parameters.maxStates );
    return operatorMatrices.chirp_mat.back() * chirp.back().get( t, !parameters.numerics_use_function_caching );
}

Sparse System::dgl_pulse( const double t ) {
    Sparse ret = Sparse( parameters.maxStates, parameters.maxStates );
    if ( pulse.empty() ) {
        return ret;
    }
    for ( int i = 0; i < pulse.size(); i++ ) {
        auto p = pulse[i].get( t, not parameters.numerics_use_function_caching );
        if ( parameters.numerics_use_rwa )
            ret += operatorMatrices.pulse_mat[2 * i + 1] * p + operatorMatrices.pulse_mat[2 * i] * std::conj( p );
        else
            ret += operatorMatrices.pulse_mat[2 * i + 1] * ( p + std::conj( p ) ) + operatorMatrices.pulse_mat[2 * i] * ( p + std::conj( p ) );
    }
    return ret;
}

void System::calculate_expectation_values( const std::vector<QDLC::SaveState> &rhos, Timer &evalTimer ) {
    // Output expectation Values
    double t_pre = rhos.front().t;
    for ( auto &tup : rhos ) {
        auto rho = tup.mat;
        double t = tup.t;
        double dt = t - t_pre;
        t_pre = t;
        std::string el_out = fmt::format( "{:.5e}", t );
        std::string ph_out = fmt::format( "{:.5e}", t );
        std::string el_em = "";
        std::string ph_em = "";
        // Output State Populations and Calculate Emission Probabilities from all radiative transitions:
        for ( auto &[mode, state] : operatorMatrices.ph_states ) {
            double expval = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, state.hilbert, t ) );
            ph_out = fmt::format( "{:}\t{:.6e}", ph_out, expval );
            if ( parameters.p_omega_cavity_loss > 0.0 ) {
                emission_probabilities[mode] += expval * dt;
                ph_em = fmt::format( "{:}\t{:.6e}", ph_em, parameters.p_omega_cavity_loss * parameters.input_photonic[mode].numerical["DecayScaling"] * emission_probabilities[mode] );
            }
        }
        for ( auto &[mode, state] : operatorMatrices.el_states ) {
            double expval = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, state.hilbert, t ) );
            el_out = fmt::format( "{:}\t{:.6e}", el_out, expval );
            if ( parameters.p_omega_decay > 0.0 and parameters.input_electronic[mode].numerical["DecayScaling"] != 0.0 ) {
                emission_probabilities[mode] += expval * dt;
                el_em = fmt::format( "{:}\t{:.6e}", el_em, parameters.p_omega_decay * parameters.input_electronic[mode].numerical["DecayScaling"] * emission_probabilities[mode] );
            }
        }
        // Output
        if ( not operatorMatrices.ph_states.empty() )
            FileOutput::get_file( "photonic" ) << fmt::format( "{:}\t{:}\n", ph_out, ph_em );
        if ( not operatorMatrices.el_states.empty() )
            FileOutput::get_file( "electronic" ) << fmt::format( "{:}\t{:}\n", el_out, el_em );

        if ( parameters.input_conf["DMconfig"].string["output_mode"] != "none" ) {
            FileOutput::get_file( "densitymatrix" ) << fmt::format( "{:.5e}\t", t ); //, rho.nonZeros(), rho.rows() * rho.cols() - rho.nonZeros() );
            if ( parameters.input_conf["DMconfig"].string["interaction_picture"] != "int" )
                rho = dgl_timetrafo( rho, t );
            if ( parameters.input_conf["DMconfig"].string["output_mode"] == "full" ) {
                for ( int i = 0; i < parameters.maxStates; i++ )
                    for ( int j = 0; j < parameters.maxStates; j++ ) {
                        FileOutput::get_file( "densitymatrix" ) << fmt::format( "{:.10e}\t", std::real( rho.coeff( i, j ) ) );
                    }
                for ( int i = 0; i < parameters.maxStates; i++ )
                    for ( int j = 0; j < parameters.maxStates; j++ ) {
                        FileOutput::get_file( "densitymatrix" ) << fmt::format( "{:.10e}\t", std::imag( rho.coeff( i, j ) ) );
                    }
            } else
                for ( int j = 0; j < parameters.maxStates; j++ )
                    FileOutput::get_file( "densitymatrix" ) << fmt::format( "{:.5e}\t", std::real( rho.coeff( j, j ) ) );
            FileOutput::get_file( "densitymatrix" ) << "\n";
        }
        evalTimer.iterate();
    }
}

Sparse System::dgl_get_hamilton( const double t ) {
    return dgl_timetrafo( operatorMatrices.H_used + dgl_pulse( t ) + dgl_chirp( t ), t );
}

bool System::exit_system( const int failure ) {
    for ( auto &p : pulse )
        p.log();
    for ( auto &c : chirp )
        c.log();
    Log::L2( "[System-PME] Coefficients: Attempts w/r: {}, Write: {}, Calc: {}, Read: {}, Interpolate/Read: {}.\n", track_getcoefficient_calcattempt, track_getcoefficient_write, track_getcoefficient_calculate, track_getcoefficient_read, track_getcoefficient_read_interpolated );
    Log::L2( "[System] Number of approx +/- adjustments: {}\n", globaltries );
    Log::L1( "[System] Maximum RAM used: {} MB\n", getPeakRSS() / 1024 / 1024 );
    FileOutput::close_all();
    return true;
}

bool System::trace_valid( Sparse &rho, double t_hit, bool force ) {
    double trace = std::real( get_trace<Scalar>( rho ) );
    parameters.trace.emplace_back( trace );
    if ( trace < 0.99 || trace > 1.01 || force ) {
        if ( force )
            fmt::print( "[System] {} {} -> trace check failed at t = {} with trace(rho) = {}\n", QDLC::Message::Prefix::PERROR, QDLC::Message::global_error_divergent, t_hit, trace );
        terminate_message = QDLC::Message::global_error_divergent;
        auto &fp_trace = FileOutput::add_file( "trace" );
        for ( int i = 0; i < (int)parameters.trace.size() && parameters.t_step * 1.0 * i < t_hit; i++ ) {
            fp_trace << fmt::format( "{:.10e} {:.15e}\n", parameters.t_step * 1.0 * ( i + 1 ), parameters.trace.at( i ) );
        }
        FileOutput::close_file( "trace" );
        return false;
    } else {
        return true;
    }
}