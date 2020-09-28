#include "system/system.h"

System::System( const std::vector<std::string> &input ) : Operatable() {
    // Set Name of this system.
    name = "Biexciton (4NS)";
    logs.level2( "Creating System Class for '{}'\n", name );
    // Initialize all subclasses with the input vector
    parameters = Parameters( input );
    operatorMatrices = OperatorMatrices( parameters );
    operatorMatrices.outputOperators( parameters );
    // Log the operator matrix base
    parameters.log( operatorMatrices.base );
    // Create all possible file outputs
    fileoutput = FileOutput( {"densitymatrix.txt", "atomicinversion.txt", "photonpopulation.txt"}, parameters, operatorMatrices );
    // Initialize / Adjust the remaining system class
    terminate_message = global_message_normaltermination;
    Timer &timer_systeminit = createTimer( "System Initialization", true, false );
    logs.level2( "System initialization... " );
    timer_systeminit.start();
    if ( !init_system() ) {
        logs.level2( "System initialization failed! Exitting program...\n" );
        logs.close();
        exit( EXIT_FAILURE );
    }
    logs.level2( "successful. Elapsed time is {}ms\n", timer_systeminit.getWallTime( TIMER_MILLISECONDS ) );
    timer_systeminit.end();
}

bool System::init_system() {
    // Single chirp for single atomic level
    Chirp::Inputs chirpinputs( parameters.t_start, parameters.t_end, parameters.t_step, parameters.chirp_type, parameters.numerics_order_highest );
    chirpinputs.add( parameters.chirp_t, parameters.chirp_y, parameters.chirp_ddt );
    chirp = Chirp( chirpinputs );
    if ( parameters.chirp_total != 0 )
        chirp.fileOutput( parameters.subfolder + "chirp.txt" );

    // Arbitrary number of pulses onto single atomic level.
    Pulse::Inputs pulseinputs_H( parameters.t_start, parameters.t_end, parameters.t_step, parameters.numerics_order_highest );
    pulseinputs_H.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "H" );
    pulseinputs_H.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "+", 1.0 / std::sqrt( 2.0 ) );
    pulseinputs_H.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "-", 1.0 / std::sqrt( 2.0 ) );
    pulse_H = Pulse( pulseinputs_H );
    Pulse::Inputs pulseinputs_V( parameters.t_start, parameters.t_end, parameters.t_step, parameters.numerics_order_highest );
    pulseinputs_V.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "V" );
    pulseinputs_V.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "+", 1i / std::sqrt( 2.0 ) );
    pulseinputs_V.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_omega_chirp, parameters.pulse_type, parameters.pulse_pol, "-", -1i / std::sqrt( 2.0 ) );
    pulse_V = Pulse( pulseinputs_V );
    if ( parameters.pulse_amp.at( 0 ) != 0 ) {
        Pulse::fileOutput( parameters.subfolder + "pulse.txt", {pulse_H, pulse_V} );
    }

    if ( parameters.numerics_use_saved_coefficients )
        savedCoefficients.reserve( parameters.numerics_saved_coefficients_max_size );

    // Output Phonon functions if phonons are active
    if ( parameters.p_phonon_T != 0 ) {
        phi_vector.reserve( std::ceil( parameters.p_phonon_tcutoff / parameters.t_step * 3.1 ) );
        phi_vector.emplace_back( dgl_phonons_phi( parameters.t_step * 0.01 ) );
        for ( double tau = parameters.t_step; tau < parameters.p_phonon_tcutoff * 3.1; tau += parameters.t_step ) {
            phi_vector.emplace_back( dgl_phonons_phi( tau ) );
        }
        FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
        fmt::print( fp_phonons, "t\treal(phi(t))\timag(phi(t))\treal(g_u(t))\timag(g_u(t))\treal(g_g(t))\timag(g_g(t))\n" );
        for ( double t = parameters.t_start; t < 3.0 * parameters.p_phonon_tcutoff; t += parameters.t_step ) {
            auto greenu = dgl_phonons_greenf( t, 'u' );
            auto greeng = dgl_phonons_greenf( t, 'g' );
            fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", t, std::real( phi_vector.at( std::floor( t / parameters.t_step ) ) ), std::imag( phi_vector.at( std::floor( t / parameters.t_step ) ) ), std::real( greenu ), std::imag( greenu ), std::real( greeng ), std::imag( greeng ) );
        }
        std::fclose( fp_phonons );
        if ( parameters.output_coefficients ) {
            fp_phonons = std::fopen( ( parameters.subfolder + "phonons_lb.txt" ).c_str(), "w" );
            fmt::print( fp_phonons, "t\tL_a_+\tL_a_-\tL_c_+\tL_c_-\n" );
            for ( double t = parameters.t_start; t < parameters.t_end; t += parameters.t_step ) {
                fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\n", t, dgl_phonons_lindblad_coefficients( t, 'L', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'L', -1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', -1.0 ) );
            }
            std::fclose( fp_phonons );
        }
    }

    // Time Transformation
    timeTrafoMatrix = ( Dense( 1i * operatorMatrices.H_0 ).exp() ).sparseView().pruned();
    // Check time trafo
    Sparse ttrafo = ( Dense( 1i * operatorMatrices.H_0 * 125E-12 ).exp() * operatorMatrices.H_used * Dense( -1i * operatorMatrices.H_0 * 125E-12 ).exp() ).sparseView();
    Sparse temp = dgl_timetrafo( operatorMatrices.H_used, 125E-12 ) - ttrafo;
    if ( std::abs( Dense( temp ).sum() / parameters.p_omega_atomic_G_H ) >= 0.0001 ) {
        logs( "Unitary timetransformation error is bigger than 0.01% of atomic transition energy!\n\n" );
    }
    return true;
}

Sparse System::dgl_rungeFunction( const Sparse &rho, const Sparse &H, const double t, std::vector<SaveState> &past_rhos ) {
    Sparse ret = -1i * dgl_kommutator( H, rho );
    // Photon losses
    if ( parameters.p_omega_cavity_loss != 0.0 ) {
        ret += parameters.p_omega_cavity_loss / 2.0 * dgl_lindblad( rho, operatorMatrices.photon_annihilate_H, operatorMatrices.photon_create_H ); //BUGFIX: kappa/2.0
        ret += parameters.p_omega_cavity_loss / 2.0 * dgl_lindblad( rho, operatorMatrices.photon_annihilate_V, operatorMatrices.photon_create_V ); //BUGFIX: kappa/2.0
    }
    // Pure Dephasing
    if ( parameters.p_omega_pure_dephasing != 0.0 ) {
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_ground * rho * operatorMatrices.atom_state_H + operatorMatrices.atom_state_H * rho * operatorMatrices.atom_state_ground );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_ground * rho * operatorMatrices.atom_state_V + operatorMatrices.atom_state_V * rho * operatorMatrices.atom_state_ground );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_H * rho * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_biexciton * rho * operatorMatrices.atom_state_H );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_V * rho * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_biexciton * rho * operatorMatrices.atom_state_V );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_H * rho * operatorMatrices.atom_state_V + operatorMatrices.atom_state_V * rho * operatorMatrices.atom_state_H );
        ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_ground * rho * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_biexciton * rho * operatorMatrices.atom_state_ground );
    }
    // Radiative decay
    if ( parameters.p_omega_decay != 0.0 ) {
        ret += parameters.p_omega_decay / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H, operatorMatrices.atom_sigmaplus_G_H );
        ret += parameters.p_omega_decay / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V, operatorMatrices.atom_sigmaplus_G_V );
        ret += parameters.p_omega_decay / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B, operatorMatrices.atom_sigmaplus_H_B );
        ret += parameters.p_omega_decay / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B, operatorMatrices.atom_sigmaplus_V_B );
    }
    if ( parameters.p_phonon_T != 0.0 ) {
        ret += dgl_phonons( rho, t, past_rhos );
    }
    return ret.pruned();
}

Sparse System::dgl_timetrafo( const Sparse &A, const double t ) {
    Sparse ret = A;
    if ( parameters.numerics_use_interactionpicture == 1 ) {
        // TIMETRANSFORMATION_ANALYTICAL
        if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) {
            std::vector<Eigen::Triplet<Scalar>> ret_v;
            for ( int k = 0; k < A.outerSize(); ++k ) {
                for ( Sparse::InnerIterator it( A, k ); it; ++it ) {
                    // Convert row/col into respective photon numbers / atomic state
                    int i = it.row() % 4;
                    int photon_H_i = std::floor( it.row() / ( 4.0 * ( parameters.p_max_photon_number + 1 ) ) );
                    int photon_V_i = ( (int)std::floor( it.row() / 4.0 ) ) % ( parameters.p_max_photon_number + 1 );
                    int j = it.col() % 4;
                    int photon_H_j = std::floor( it.col() / ( 4.0 * ( parameters.p_max_photon_number + 1 ) ) );
                    int photon_V_j = ( (int)std::floor( it.col() / 4.0 ) ) % ( parameters.p_max_photon_number + 1 );
                    // Add new Value to triplet list
                    Scalar factor = std::exp( 1i * t * ( parameters.p_omega_atomic_G_H * ( delta( i, 1 ) - delta( j, 1 ) ) + parameters.p_omega_atomic_G_V * ( delta( i, 2 ) - delta( j, 2 ) ) + parameters.p_omega_atomic_B * ( delta( i, 3 ) - delta( j, 3 ) ) + parameters.p_omega_cavity_H * ( photon_H_i - photon_H_j ) + parameters.p_omega_cavity_V * ( photon_V_i - photon_V_j ) ) );
                    //Scalar factor = std::exp( 1i * t * ( parameters.p_omega_atomic_G_H / 2.0 * ( delta( i, 1 ) - delta( i, 0 ) - delta( j, 1 ) + delta( j, 0 ) ) + parameters.p_omega_atomic_G_V / 2.0 * ( delta( i, 2 ) - delta( i, 0 ) - delta( j, 2 ) + delta( j, 0 ) ) + parameters.p_omega_atomic_H_B / 2.0 * ( delta( i, 3 ) - delta( i, 1 ) - delta( j, 3 ) + delta( j, 1 ) ) + parameters.p_omega_atomic_V_B / 2.0 * ( delta( i, 3 ) - delta( i, 2 ) - delta( j, 3 ) + delta( j, 2 ) ) + parameters.p_omega_cavity_H * ( photon_H_i - photon_H_j ) + parameters.p_omega_cavity_V * ( photon_V_i - photon_V_j ) ) );
                    ret_v.emplace_back( it.row(), it.col(), it.value() * factor );
                }
            }
            // Generate new Matrix from triplet list
            ret.setFromTriplets( ret_v.begin(), ret_v.end() );
        }
        // TIMETRANSFORMATION_MATRIXEXPONENTIAL
        else if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_MATRIXEXPONENTIAL ) {
            Sparse U = ( Dense( 1i * operatorMatrices.H_0 * t ).exp() ).sparseView();
            ret = U * A * U.adjoint();
        }
    }
    return ret;
}

Sparse System::dgl_chirp( const double t ) {
    if ( parameters.chirp_total == 0 )
        return Sparse( parameters.maxStates, parameters.maxStates );
    return ( 2.0*operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_H + operatorMatrices.atom_state_V ) * chirp.get( t ); //Experimental; corrected chirp
}

Sparse System::dgl_pulse( const double t ) {
    if ( parameters.pulse_amp.at( 0 ) == 0 ) {
        return Sparse( parameters.maxStates, parameters.maxStates );
    }
    return ( ( operatorMatrices.atom_sigmaminus_G_H + operatorMatrices.atom_sigmaminus_H_B ) * std::conj( pulse_H.get( t ) ) + ( operatorMatrices.atom_sigmaminus_G_V + operatorMatrices.atom_sigmaminus_V_B ) * std::conj( pulse_V.get( t ) ) + ( operatorMatrices.atom_sigmaplus_G_H + operatorMatrices.atom_sigmaplus_H_B ) * pulse_H.get( t ) + ( operatorMatrices.atom_sigmaplus_G_V + operatorMatrices.atom_sigmaplus_V_B ) * pulse_V.get( t ) );
    //return 0.5 * ( operatorMatrices.atom_sigmaplus_G_H * pulse_H.get( t ) + operatorMatrices.atom_sigmaminus_G_H * std::conj( pulse_H.get( t ) ) + operatorMatrices.atom_sigmaplus_G_V * pulse_V.get( t ) + operatorMatrices.atom_sigmaminus_G_V * std::conj( pulse_V.get( t ) ) );
}

Sparse System::dgl_phonons_rungefunc( const Sparse &chi, const double t ) {
    double chirpcorrection = chirp.get( t ) + t * ( chirp.get( t ) - chirp.derivative( t ) );
    // FIX: für w_b-w_xi muss da -chirpcorrection statt +chirpcorrection!
    //FIXME: E_B(tau)!! -> 2de3lta - delta -> +delta statt -delta für B-X
    Sparse explicit_time = 1i * ( parameters.p_omega_atomic_G_H + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_G_H + 1i * ( parameters.p_omega_atomic_H_B + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_H_B + 1i * ( ( parameters.p_omega_atomic_H_B + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_H_B + ( parameters.p_omega_atomic_G_H + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_G_H - parameters.p_omega_cavity_H * operatorMatrices.projector_atom_sigmaplus_G_H - parameters.p_omega_cavity_H * operatorMatrices.projector_atom_sigmaplus_H_B ) * operatorMatrices.projector_photon_annihilate_H + pulse_H.derivative( t ) * ( operatorMatrices.projector_atom_sigmaplus_G_H + operatorMatrices.projector_atom_sigmaplus_H_B );
    explicit_time += 1i * ( parameters.p_omega_atomic_G_V + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_G_V + 1i * ( parameters.p_omega_atomic_V_B + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_V_B + 1i * ( ( parameters.p_omega_atomic_V_B + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_V_B + ( parameters.p_omega_atomic_G_V + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_G_V - parameters.p_omega_cavity_V * operatorMatrices.projector_atom_sigmaplus_G_V - parameters.p_omega_cavity_V * operatorMatrices.projector_atom_sigmaplus_V_B ) * operatorMatrices.projector_photon_annihilate_V + pulse_V.derivative( t ) * ( operatorMatrices.projector_atom_sigmaplus_G_V + operatorMatrices.projector_atom_sigmaplus_V_B );

    Sparse hamilton = dgl_getHamilton( t );

    return -1i * dgl_kommutator( hamilton, chi ) + explicit_time.cwiseProduct( project_matrix_sparse( chi ) );
}

Sparse System::dgl_phonons_chiToX( const Sparse &chi, const char mode ) {
    Sparse adjoint = chi.adjoint();
    if ( mode == 'g' ) {
        return chi + adjoint;
    }
    return 1i * ( chi - adjoint );
}

Scalar System::dgl_phonons_greenf( double t, const char mode ) {
    int i = std::floor( t / parameters.t_step );
    auto phi = phi_vector.at( i ); //dgl_phonons_phi( t );
    if ( mode == 'g' ) {
        return parameters.p_phonon_b * parameters.p_phonon_b * ( std::cosh( phi ) - 1.0 );
    }
    return parameters.p_phonon_b * parameters.p_phonon_b * std::sinh( phi );
}

Scalar System::dgl_phonons_phi( const double t ) {
    Scalar integral = 0;
    double stepsize = 0.01 * parameters.p_phonon_wcutoff;
    for ( double w = stepsize; w < 10 * parameters.p_phonon_wcutoff; w += stepsize ) {
        integral += stepsize * ( parameters.p_phonon_alpha * w * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff ) * ( std::cos( w * t ) / std::tanh( 1.0545718E-34 * w / 2.0 / 1.3806488E-23 / parameters.p_phonon_T ) - 1i * std::sin( w * t ) ) );
    }
    return integral;
}

double System::dgl_phonons_lindblad_coefficients( double t, double omega_atomic, const char mode, const char level, const double sign ) {
    if ( t == 0.0 )
        t = parameters.t_step * 0.01;
    double ret = 0;
    double step = parameters.t_step;
    if ( mode == 'L' ) {
        double bpulsesquared, delta, nu;
        if ( level == 'H' ) {
            bpulsesquared = std::pow( std::abs( parameters.p_phonon_b * pulse_H.get( t ) ), 2.0 ); //TODO: chirp correction for omega_atomic!!!!
            delta = parameters.pulse_omega.at( 0 ) - omega_atomic;                                 //FIXME : different pulse frequencies
        } else {
            bpulsesquared = std::pow( std::abs( parameters.p_phonon_b * pulse_V.get( t ) ), 2.0 );
            delta = parameters.pulse_omega.at( 0 ) - omega_atomic; //FIXME : different pulse frequencies
        }
        nu = std::sqrt( bpulsesquared + delta * delta );
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += step ) {
            Scalar f = ( delta * delta * std::cos( nu * tau ) + bpulsesquared ) / std::pow( nu, 2.0 );
            ret += std::real( ( std::cosh( phi_vector.at( i ) ) - 1.0 ) * f + std::sinh( phi_vector.at( i ) ) * std::cos( nu * tau ) ) - sign * std::imag( ( std::exp( phi_vector.at( i ) ) - 1.0 ) * delta * std::sin( nu * tau ) / nu );
            i++;
        }
        ret *= 2.0 * bpulsesquared * step;
    } else if ( mode == 'C' ) {
        double delta;
        if ( mode == 'H' ) {
            delta = parameters.p_omega_cavity_H - omega_atomic;
        } else {
            delta = parameters.p_omega_cavity_V - omega_atomic;
        }
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += step ) {
            ret += std::real( std::exp( 1i * sign * delta * tau ) * ( std::exp( phi_vector.at( i ) ) - 1.0 ) );
            i++;
        }
        ret *= parameters.p_phonon_b * parameters.p_phonon_b * parameters.p_omega_coupling * parameters.p_omega_coupling * step;
    }
    return ret;
}

int System::dgl_get_coefficient_index( const double t, const double tau ) {
    if ( parameters.numerics_use_saved_coefficients ) {
        // Look for already calculated coefficient. if not found, calculate and save new coefficient
        if ( savedCoefficients.size() > 0 && t <= savedCoefficients.back().t && tau <= savedCoefficients.back().tau ) {
            track_getcoefficient_calcattempt++;
            double mult = ( parameters.numerics_phonon_approximation_markov1 ? std::floor( parameters.p_phonon_tcutoff / parameters.t_step ) : 1.0 );
            int add = ( parameters.numerics_phonon_approximation_markov1 ? (int)std::floor( ( parameters.p_phonon_tcutoff - tau ) / parameters.t_step ) : 0 );
            int approx = std::max( 0, std::min( (int)savedCoefficients.size() - 1, (int)( ( std::floor( t / parameters.t_step * 2.0 ) - 1.0 ) * mult ) ) - add ); //TODO: *4.0 wenn RK5, *parameters.rungekutta_timestep_multiplier
            int tries = 0;
            while ( approx < (int)savedCoefficients.size() - 1 && t > savedCoefficients.at( approx ).t ) {
                approx++;
                tries++;
            }
            while ( approx > 0 && t < savedCoefficients.at( approx ).t ) {
                approx--;
                tries++;
            }
            while ( approx < (int)savedCoefficients.size() - 1 && tau > savedCoefficients.at( approx ).tau ) {
                approx++;
                tries++;
            }
            while ( approx > 0 && tau < savedCoefficients.at( approx ).tau ) {
                approx--;
                tries++;
            }
            //logs.level2("Approx = {}, size = {}\n",approx,savedCoefficients.size());
            if ( approx < (int)savedCoefficients.size() ) {
                if ( t != savedCoefficients.at( approx ).t || tau != savedCoefficients.at( approx ).tau ) {
                    // This situation may occur during multithreading.
                    track_getcoefficient_read_but_unequal++;
                    //logs.level2( "Coefficient time mismatch! t = {} but coefficient.t = {}, tau = {} but coefficient.tau = {}\n", t, savedCoefficients.at( approx ).t, tau, savedCoefficients.at( approx ).tau );
                } else {
                    track_getcoefficient_read++;
                    globaltries += tries;
                    //logs.level2( "Coefficient time match! t = {} but coefficient.t = {}, tau = {} but coefficient.tau = {}\n", t, savedCoefficients.at( approx ).t, tau, savedCoefficients.at( approx ).tau );
                    return approx;
                }
            }
        }
    }
    return -1;
}

void System::dgl_save_coefficient( const Sparse &coefficient1, const Sparse &coefficient2, const double t, const double tau ) {
    track_getcoefficient_calculate++;
    if ( parameters.numerics_use_saved_coefficients && savedCoefficients.size() < parameters.numerics_saved_coefficients_max_size ) {
        track_getcoefficient_write++;
#pragma omp critical
        savedCoefficients.emplace_back( SaveStateTau( coefficient1, coefficient2, t, tau ) );
        // If use cutoff and vector saved more than 5 timesteps worth of matrices, delete other matrices
        if ( parameters.numerics_saved_coefficients_cutoff > 0 && savedCoefficients.size() > parameters.numerics_saved_coefficients_cutoff ) {
#pragma omp critical
            savedCoefficients.erase( savedCoefficients.begin() );
        }
        //logs.level2("Saved coefficient for t = {}, tau = {}\n",t,tau);
    }
}

Sparse System::dgl_phonons_chi( const double t ) {
    return dgl_timetrafo( parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H + ( operatorMatrices.atom_sigmaplus_G_H + operatorMatrices.atom_sigmaplus_H_B ) * pulse_H.get( t ) + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V + ( operatorMatrices.atom_sigmaplus_G_V + operatorMatrices.atom_sigmaplus_V_B ) * pulse_V.get( t ), t );
}

Sparse System::dgl_phonons_calculate_transformation( Sparse &chi_tau, double t, double tau ) {
    if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_BACKWARDS_INTEGRAL ) {
        return Solver::calculate_definite_integral( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), -parameters.t_step ).mat;
    } else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_TRANSFORMATION_MATRIX ) {
        Sparse U = ( Dense( -1i * dgl_getHamilton( t ) * tau ).exp() ).sparseView();
        return ( U * chi_tau * U.adjoint() ).eval();
    } else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_MIXED ) {
        double error = std::abs( pulse_H.get( t ) + pulse_V.get( t ) );
        if ( ( pulse_H.maximum > 0 && error > pulse_H.maximum * 0.1 ) || ( pulse_V.maximum > 0 && error > pulse_V.maximum * 0.1 ) || chirp.derivative(t) != 0) {
            return Solver::calculate_definite_integral( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), -parameters.t_step ).mat;
        }
        return chi_tau;
    } else {
        return chi_tau;
    }
}

Sparse System::dgl_phonons( const Sparse &rho, const double t, const std::vector<SaveState> &past_rhos ) {
    Sparse ret = Sparse( parameters.maxStates, parameters.maxStates );
    Sparse chi = dgl_phonons_chi( t );
    std::vector<Sparse> threadmap_1, threadmap_2;

    // Most precise approximation used. Caclulate polaron fram Chi by integrating backwards from t to t-tau.
    if ( parameters.numerics_phonon_approximation_order != PHONON_APPROXIMATION_LINDBLAD_FULL ) {
        Sparse XUT = dgl_phonons_chiToX( chi, 'u' );
        Sparse XGT = dgl_phonons_chiToX( chi, 'g' );
        int _taumax = (int)std::min( parameters.p_phonon_tcutoff / parameters.t_step, t / parameters.t_step );
        // Use markov approximation
        if ( parameters.numerics_phonon_approximation_markov1 ) {
            // Temporary variables
            Sparse chi_tau_back_u, chi_tau_back_g, integrant;
            // Look if we already calculated coefficient sum for this specific t-value
            int index = dgl_get_coefficient_index( t, 0 );
            if ( index != -1 ) {
                // Index was found, chi(t-tau) sum used from saved vector
                chi_tau_back_u = savedCoefficients.at( index ).mat1;
                chi_tau_back_g = savedCoefficients.at( index ).mat2;
            } else {
                // Index was not found, (re)calculate chi(t-tau) sum
                // Initialize temporary matrices to zero for threads to write to
                init_sparsevector( threadmap_1, parameters.maxStates, parameters.numerics_phonons_maximum_threads );
                init_sparsevector( threadmap_2, parameters.maxStates, parameters.numerics_phonons_maximum_threads );
                // Calculate backwards integral and sum it into threadmaps. Threadmaps will later be summed into one coefficient matrix.
#pragma omp parallel for ordered schedule( dynamic ) shared( savedCoefficients ) num_threads( parameters.numerics_phonons_maximum_threads )
                for ( int _tau = 0; _tau < _taumax; _tau++ ) {
                    Sparse chi_tau_back_u, chi_tau_back_g, chi_tau, chi_tau_back;
                    double tau = ( 1.0 * _tau ) * parameters.t_step;
                    chi_tau = dgl_phonons_chi( t - tau );
                    chi_tau_back = dgl_phonons_calculate_transformation( chi_tau, t, tau );
                    chi_tau_back_u = dgl_phonons_chiToX( chi_tau_back, 'u' );
                    chi_tau_back_g = dgl_phonons_chiToX( chi_tau_back, 'g' );
                    auto thread = omp_get_thread_num();
                    threadmap_1.at( thread ) += dgl_phonons_greenf( tau, 'u' ) * chi_tau_back_u;
                    threadmap_2.at( thread ) += dgl_phonons_greenf( tau, 'g' ) * chi_tau_back_g;
                }
                // Sum all contributions from threadmaps into one coefficient
                chi_tau_back_u = std::accumulate( threadmap_1.begin(), threadmap_1.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
                chi_tau_back_g = std::accumulate( threadmap_2.begin(), threadmap_2.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
                // Save coefficients
                dgl_save_coefficient( chi_tau_back_u, chi_tau_back_g, t, 0 );
            }
            // Calculate phonon contributions from (saved/calculated) coefficients and rho(t)
            integrant = dgl_kommutator( XUT, ( chi_tau_back_u * rho ).eval() );
            integrant += dgl_kommutator( XGT, ( chi_tau_back_g * rho ).eval() );
            Sparse adjoint = integrant.adjoint();
            ret -= ( integrant + adjoint ) * parameters.t_step;
        } else {
            init_sparsevector( threadmap_1, parameters.maxStates, parameters.numerics_phonons_maximum_threads );
            // Dont use markov approximation; Full integral
#pragma omp parallel for ordered schedule( dynamic ) shared( savedCoefficients ) num_threads( parameters.numerics_phonons_maximum_threads )
            for ( int _tau = 0; _tau < _taumax; _tau++ ) {
                Sparse chi_tau_back_u, chi_tau_back_g, chi_tau, chi_tau_back, integrant;
                int rho_index = std::max( 0, (int)past_rhos.size() - 1 - _tau );
                double tau = ( 1.0 * _tau ) * parameters.t_step;
                // Check if chi(t-tau) has to be recalculated, or can just be retaken from saved matrices, since chi(t-tau) is only dependant on time
                int index = dgl_get_coefficient_index( t, tau );
                if ( index != -1 ) {
                    // Index was found, chi(t-tau) used from saved vector
                    chi_tau_back_u = savedCoefficients.at( index ).mat1;
                    chi_tau_back_g = savedCoefficients.at( index ).mat2;
                } else {
                    // Index not found, or no saving is used. Recalculate chi(t-tau).
                    chi_tau = dgl_phonons_chi( t - tau );
                    chi_tau_back = chi_tau_back = dgl_phonons_calculate_transformation( chi_tau, t, tau );
                    chi_tau_back_u = dgl_phonons_chiToX( chi_tau_back, 'u' );
                    chi_tau_back_g = dgl_phonons_chiToX( chi_tau_back, 'g' );
                    // If saving is used, save current chi(t-tau). If markov approximation is used, only save final contributions
                    dgl_save_coefficient( chi_tau_back_u, chi_tau_back_g, t, tau );
                }
                integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( XUT, ( chi_tau_back_u * past_rhos.at( rho_index ).mat ).eval() );
                integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( XGT, ( chi_tau_back_g * past_rhos.at( rho_index ).mat ).eval() );
                Sparse adjoint = integrant.adjoint();
                auto thread = omp_get_thread_num();
                threadmap_1.at( thread ) += ( integrant + adjoint ) * parameters.t_step;
            }

            ret -= std::accumulate( threadmap_1.begin(), threadmap_1.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
        }
    } else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_LINDBLAD_FULL ) {
        // H
        double chirpcorrection = chirp.get( t );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H + chirpcorrection, 'L', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H, operatorMatrices.atom_sigmaplus_G_H );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H + chirpcorrection, 'L', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_H, operatorMatrices.atom_sigmaminus_G_H );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H + chirpcorrection, 'C', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H * operatorMatrices.photon_create_H, operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H + chirpcorrection, 'C', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H, operatorMatrices.atom_sigmaminus_G_H * operatorMatrices.photon_create_H );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B + chirpcorrection, 'L', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B, operatorMatrices.atom_sigmaplus_H_B );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B + chirpcorrection, 'L', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_H_B, operatorMatrices.atom_sigmaminus_H_B );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B + chirpcorrection, 'C', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B * operatorMatrices.photon_create_H, operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B + chirpcorrection, 'C', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H, operatorMatrices.atom_sigmaminus_H_B * operatorMatrices.photon_create_H );
        // V
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V + chirpcorrection, 'L', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V, operatorMatrices.atom_sigmaplus_G_V );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V + chirpcorrection, 'L', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_V, operatorMatrices.atom_sigmaminus_G_V );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V + chirpcorrection, 'C', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V * operatorMatrices.photon_create_V, operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V + chirpcorrection, 'C', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V, operatorMatrices.atom_sigmaminus_G_V * operatorMatrices.photon_create_V );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B + chirpcorrection, 'L', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B, operatorMatrices.atom_sigmaplus_V_B );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B + chirpcorrection, 'L', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_V_B, operatorMatrices.atom_sigmaminus_V_B );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B + chirpcorrection, 'C', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B * operatorMatrices.photon_create_V, operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V );
        ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B + chirpcorrection, 'C', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V, operatorMatrices.atom_sigmaminus_V_B * operatorMatrices.photon_create_V );
    }
    return ret;
}

Scalar System::dgl_raman_population_increment( const std::vector<SaveState> &past_rhos, const char mode, const Scalar before, const double t ) {
    Scalar ret = 0;
    double chirpcorrection = chirp.get( t );
    double w1 = ( mode == 'h' ? parameters.p_omega_atomic_G_H + chirpcorrection : parameters.p_omega_atomic_G_V + chirpcorrection );
    double w2 = ( mode == 'h' ? parameters.p_omega_atomic_H_B + chirpcorrection : parameters.p_omega_atomic_V_B + chirpcorrection );
    double wc = ( mode == 'h' ? parameters.p_omega_cavity_H : parameters.p_omega_cavity_V );
    double sigma1 = parameters.p_omega_pure_dephasing + parameters.p_omega_decay;
    double sigma2 = parameters.p_omega_pure_dephasing + 3. * parameters.p_omega_decay;
    Sparse op = operatorMatrices.atom_sigmaminus_G_B * ( mode == 'h' ? operatorMatrices.photon_create_H : operatorMatrices.photon_create_V );
    double A = std::exp( -parameters.p_omega_cavity_loss * parameters.t_step );
    Scalar B, R;
#pragma omp parallel for ordered schedule( dynamic ) shared( past_rhos ) num_threads( parameters.numerics_phonons_maximum_threads )
    for ( long unsigned int i = 0; i < past_rhos.size(); i++ ) {
        //for ( SaveState savestate : past_rhos ) {
        double tdd = past_rhos.at( i ).t;
        B = std::exp( -1i * ( w2 - wc - 0.5i * ( parameters.p_omega_cavity_loss + sigma2 ) ) * ( t - tdd ) ) - std::exp( -1i * ( w1 - wc - 0.5i * ( parameters.p_omega_cavity_loss + sigma1 ) ) * ( t - tdd ) );
        R = dgl_expectationvalue<Sparse, Scalar>( past_rhos.at( i ).mat, op, tdd ) * ( mode == 'h' ? std::conj( pulse_H.get( tdd ) ) : std::conj( pulse_V.get( tdd ) ) );
        //fmt::print("t = {}, tau = {}, A = {}, B = {}, R = {}\n",t,tdd,A,B,R);
#pragma omp critical
        ret += B * R;
    }
    return A * before + parameters.t_step * parameters.t_step * ret;
}

void System::expectationValues( const Sparse &rho, const double t, const std::vector<SaveState> &past_rhos ) {
    // Calculate reused expectation values
    double exp_photon_H = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.photon_n_H, t ) );
    double exp_photon_V = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.photon_n_V, t ) );
    // Calculating photon emission probability through photon expectation values
    photonemissionprob_integral_H += exp_photon_H * parameters.t_step * parameters.p_omega_cavity_loss;
    photonemissionprob_integral_V += exp_photon_V * parameters.t_step * parameters.p_omega_cavity_loss;
    // Calculating raman photon expectation values and emission probability
    if ( parameters.numerics_output_raman_population ) {
        ramanphotonpopulation_integral_H = dgl_raman_population_increment( past_rhos, 'h', ramanphotonpopulation_integral_H, t );
        ramanphotonpopulation_integral_V = dgl_raman_population_increment( past_rhos, 'v', ramanphotonpopulation_integral_V, t );
        ramanphotonemissionprob_integral_H += std::real( 2.0 * parameters.p_omega_coupling * ramanphotonpopulation_integral_H ) * parameters.t_step * parameters.p_omega_cavity_loss;
        ramanphotonemissionprob_integral_V += std::real( 2.0 * parameters.p_omega_coupling * ramanphotonpopulation_integral_V ) * parameters.t_step * parameters.p_omega_cavity_loss;
    }
    // Output expectation values
    double ev_operators_electronic_H = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.atom_state_H, t ) );
    double ev_operators_electronic_V = std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.atom_state_V, t ) );
    fmt::print( fileoutput.fp_atomicinversion, "{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\n", t, std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.atom_state_ground, t ) ), ev_operators_electronic_H, ev_operators_electronic_V, std::real( dgl_expectationvalue<Sparse, Scalar>( rho, operatorMatrices.atom_state_biexciton, t ) ) );
    fmt::print( fileoutput.fp_photonpopulation, "{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}", t, exp_photon_H, exp_photon_V, photonemissionprob_integral_H, photonemissionprob_integral_V );
    if ( parameters.numerics_output_raman_population ) {
        fmt::print( fileoutput.fp_photonpopulation, "\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}", std::real(ramanphotonpopulation_integral_H), std::real(ramanphotonpopulation_integral_V), ramanphotonemissionprob_integral_H, ramanphotonemissionprob_integral_V );
    }

    if ( parameters.numerics_output_electronic_emission ) {
        electronic_emissionprob_integral_H += ev_operators_electronic_H * parameters.t_step * parameters.p_omega_decay;
        electronic_emissionprob_integral_V += ev_operators_electronic_V * parameters.t_step * parameters.p_omega_decay;
        fmt::print( fileoutput.fp_photonpopulation, "\t{:.5e}\t{:.5e}", electronic_emissionprob_integral_H, electronic_emissionprob_integral_V );
    }
    fmt::print( fileoutput.fp_photonpopulation, "\n" );

    if ( !parameters.output_no_dm ) {
        fmt::print( fileoutput.fp_densitymatrix, "{:.5e}\t", t );
        if ( parameters.output_full_dm ) {
            for ( int i = 0; i < parameters.maxStates; i++ )
                for ( int j = 0; j < parameters.maxStates; j++ ) {
                    fmt::print( fileoutput.fp_densitymatrix, "{:.5e}\t", std::real( rho.coeff( i, j ) ) );
                }
        } else
            for ( int j = 0; j < parameters.maxStates; j++ )
                fmt::print( fileoutput.fp_densitymatrix, "{:.5e}\t", std::real( rho.coeff( j, j ) ) );
        fmt::print( fileoutput.fp_densitymatrix, "\n" );
    }
}

Sparse System::dgl_getHamilton( const double t ) {
    return dgl_timetrafo( parameters.p_phonon_b * ( operatorMatrices.H_used + dgl_pulse( t ) ) + dgl_chirp( t ), t );
}

bool System::command( unsigned int index ) {
    // The only reason for classes using this system class to set the main programs maximum threads to 1 is, if the usual T-direction is already done.
    if ( index == Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM ) {
        parameters.numerics_phonons_maximum_threads = 1;
        logs.level2( "Set maximum number of Threads for primary calculations to {}\n", parameters.numerics_phonons_maximum_threads );
        if ( parameters.numerics_use_saved_coefficients ) {
            // Sort after t
            logs.level2( "Sorting saved coefficients by t... " );
            std::sort( savedCoefficients.begin(), savedCoefficients.end(), Save_State_sort_t );
            // Sort after tau
            logs.level2( "Done, sorting saved coefficients chunkwise by tau... " );
            int current_start = 0;
            int i;
            for ( i = 0; i < (int)savedCoefficients.size(); i++ ) {
                if ( savedCoefficients.at( i ).t > savedCoefficients.at( current_start ).t || i >= (int)savedCoefficients.size() ) {
                    std::sort( savedCoefficients.begin() + current_start, savedCoefficients.begin() + i, Save_State_sort_tau );
                    current_start = i;
                }
            }
            // Sort last entries
            std::sort( savedCoefficients.begin() + current_start, savedCoefficients.begin() + i, Save_State_sort_tau );
            // Just in case, test coefficients:
            //for ( i = 0; i < (int)savedCoefficients.size() - 1; i++ ) {
            //    if ( ( savedCoefficients.at( i ).t == savedCoefficients.at( i + 1 ).t && ( savedCoefficients.at( i ).tau >= savedCoefficients.at( i + 1 ).tau ) ) || savedCoefficients.at( i ).t > savedCoefficients.at( i + 1 ).t ) {
            //        logs.level2( "Coefficient mismatch after sorting! i = {}, t = {}, tau = {}, i+1 = {}, t = {}, tau = {}\n", i, savedCoefficients.at( i ).t, savedCoefficients.at( i ).tau, i + 1, savedCoefficients.at( i + 1 ).t, savedCoefficients.at( i + 1 ).tau );
            //    }
            //logs.level2( "T = {}, tau 0 {}\n", savedCoefficients.at( i ).t, savedCoefficients.at( i ).tau );
            //}
            logs.level2( "Done!\n" );
        }
    }
    return true;
}

bool System::exit_system( const int failure ) {
    pulse_H.log();
    pulse_V.log();
    chirp.log();
    logs.level2( "Coefficients: Attempts w/r: {}, Write: {}, Calc: {}, Read: {}, Read-But-Not-Equal: {}. Done!\n", track_getcoefficient_calcattempt, track_getcoefficient_write, track_getcoefficient_calculate, track_getcoefficient_read, track_getcoefficient_read_but_unequal );
    logs.level2( "Number of approx+/- adjustments: {}\n", globaltries );
    logs( "Maximum RAM used: {} MB\n", getPeakRSS() / 1024 / 1024 );
    fileoutput.close();
    return true;
}

bool System::traceValid( MatType &rho, double t_hit, bool force ) {
    double trace = std::real( getTrace<Scalar>( rho ) );
    parameters.trace.emplace_back( trace );
    if ( trace < 0.99 || trace > 1.01 || force ) {
        if ( force )
            fmt::print( "{} {} -> trace check failed at t = {} with trace(rho) = {}\n", PREFIX_ERROR, global_message_error_divergent, t_hit, trace );
        terminate_message = global_message_error_divergent;
        parameters.numerics_calculate_spectrum = 0;
        FILE *fp_trace = std::fopen( ( parameters.subfolder + "trace.txt" ).c_str(), "w" );
        for ( int i = 0; i < (int)parameters.trace.size() && parameters.t_step * 1.0 * i < t_hit; i++ ) {
            fmt::print( fp_trace, "{:.10e} {:.15e}\n", parameters.t_step * 1.0 * ( i + 1 ), parameters.trace.at( i ) );
        }
        std::fclose( fp_trace );
        return false;
    } else {
        return true;
    }
}