#include "system/system.h"

// Phonon Contributions calculated via the Polaron Frame Approach
Scalar System::dgl_phonons_phi( const double t ) {
    Scalar integral = 0;
    double stepsize = 0.01 * parameters.p_phonon_wcutoff;
    //double eV7 = convertParam<double>( "7.0eV" );
    //double eV35 = -convertParam<double>( "3.5eV" );
    //double v_c = 5110.0;
    //double a_e = 5E-9;       //3E-9;
    //double a_h = 0.87 * a_e; //a_e / 1.15;
    //double rho = 5370.0;
    for ( double w = stepsize; w < 10 * parameters.p_phonon_wcutoff; w += stepsize ) {
        double J = parameters.p_phonon_alpha * w * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff );
        //double J = w * parameters.hbar * std::pow( eV7 * std::exp( -w * w * a_e * a_e / ( 4. * v_c * v_c ) ) - eV35 * std::exp( -w * w * a_h * a_h / ( 4. * v_c * v_c ) ), 2. ) / ( 4. * 3.1415 * 3.1415 * rho * std::pow( v_c, 5. ) );
        integral += stepsize * ( J * ( std::cos( w * t ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) - 1i * std::sin( w * t ) ) );
    }
    return integral;
}

void System::initialize_polaron_frame_functions() {
    if ( parameters.p_phonon_T >= 0 ) {
        // Initialize Phi(tau)
        phi_vector[0.0] = dgl_phonons_phi( parameters.t_step * 0.001 );
        for ( double tau = parameters.t_step; tau < parameters.p_phonon_tcutoff * 3.1; tau += parameters.t_step ) {
            phi_vector[tau] = dgl_phonons_phi( tau );
        }

        // Output Phonon Functions
        FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
        fmt::print( fp_phonons, "t\treal(phi(t))\timag(phi(t))\treal(g_u(t))\timag(g_u(t))\treal(g_g(t))\timag(g_g(t))\n" );
        for ( double t = parameters.t_start; t < 3.0 * parameters.p_phonon_tcutoff; t += parameters.t_step ) {
            auto greenu = dgl_phonons_greenf( t, 'u' );
            auto greeng = dgl_phonons_greenf( t, 'g' );
            fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", t, std::real( phi_vector[t] ), std::imag( phi_vector[t] ), std::real( greenu ), std::imag( greenu ), std::real( greeng ), std::imag( greeng ) );
        }
        std::fclose( fp_phonons );
        if ( parameters.output_coefficients ) {
            //fp_phonons = std::fopen( ( parameters.subfolder + "phonons_lb.txt" ).c_str(), "w" );
            //fmt::print( fp_phonons, "t\tL_a_+\tL_a_-\tL_c_+\tL_c_-\n" );
            //for ( double t = parameters.t_start; t < parameters.t_end; t += parameters.t_step ) {
            //    fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\n", t, dgl_phonons_lindblad_coefficients( t, 'L', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'L', -1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', -1.0 ) );
            //}
            //std::fclose( fp_phonons );
        }
    }
}

Sparse System::dgl_phonons_rungefunc( const Sparse &chi, const double t ) {
    //TODO Maybe? Cache explicit times?
    double chirpcorrection = chirp.size() > 0 ? ( chirp.back().get( t ) + t * ( chirp.back().get( t ) - parameters.scaleVariable( chirp.back().derivative( t ), parameters.scale_value ) ) ) : 0;
    Sparse explicit_time = Sparse( chi.rows(), chi.cols() );
    for ( auto &[mode, param] : operatorMatrices.el_transitions ) {
        if ( param.direction == -1 )
            continue;
        explicit_time += ( param.energy + chirpcorrection ) * param.projector;
    }
    for ( auto &[mode, param] : parameters.input_photonic ) { //FIXME: ohne cav, keine übergäng!!
        for ( auto transition : param.string_v["CoupledTo"] ) {
            std::reverse( transition.begin(), transition.end() );
            explicit_time += 1.0i * ( operatorMatrices.el_transitions[transition].energy + chirpcorrection - operatorMatrices.ph_states[mode].energy ) * operatorMatrices.el_transitions[transition].projector * operatorMatrices.ph_transitions[mode + "b"].projector;
        }
    }
    int p = 0;
    int pp = 0;
    for ( auto &[mode, param] : parameters.input_pulse ) {
        for ( auto transition : param.string_v["CoupledTo"] ) {
            std::reverse( transition.begin(), transition.end() );
            explicit_time += 1.0i * ( pulse[p].get( t ) + 1.0i * parameters.scaleVariable( pulse[p].derivative( t ), parameters.scale_value ) ) * operatorMatrices.polaron_pulse_factors_explicit_time[pp++];
        }
        p++;
    }

    //Sparse explicit_time = 1i * ( parameters.p_omega_atomic_G_H + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_G_H
    // + 1i * ( parameters.p_omega_atomic_H_B + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_H_B
    // + 1i * (
    //  ( parameters.p_omega_atomic_H_B + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_H_B
    //      + ( parameters.p_omega_atomic_G_H + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_G_H
    //      - parameters.p_omega_cavity_H * operatorMatrices.projector_atom_sigmaplus_G_H
    //      - parameters.p_omega_cavity_H * operatorMatrices.projector_atom_sigmaplus_H_B
    //  ) * operatorMatrices.projector_photon_annihilate_H
    // + parameters.scaleVariable( pulse_H.derivative( t ), parameters.scale_value ) * ( operatorMatrices.projector_atom_sigmaplus_G_H + operatorMatrices.projector_atom_sigmaplus_H_B );
    //explicit_time += 1i * ( parameters.p_omega_atomic_G_V + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_G_V + 1i * ( parameters.p_omega_atomic_V_B + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_V_B + 1i * ( ( parameters.p_omega_atomic_V_B + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_V_B + ( parameters.p_omega_atomic_G_V + chirpcorrection ) * operatorMatrices.projector_atom_sigmaplus_G_V - parameters.p_omega_cavity_V * operatorMatrices.projector_atom_sigmaplus_G_V - parameters.p_omega_cavity_V * operatorMatrices.projector_atom_sigmaplus_V_B ) * operatorMatrices.projector_photon_annihilate_V + parameters.scaleVariable( pulse_V.derivative( t ), parameters.scale_value ) * ( operatorMatrices.projector_atom_sigmaplus_G_V + operatorMatrices.projector_atom_sigmaplus_V_B );
    explicit_time = parameters.scaleVariable( explicit_time, 1.0 / parameters.scale_value );

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
    if ( phi_vector.count( t ) == 0 )
        phi_vector[t] = dgl_phonons_phi( t ); //std::cout << "Time " << t << " is not in phi vecotr\n";
    auto phi = phi_vector[t];
    if ( mode == 'g' ) {
        return parameters.p_phonon_b * parameters.p_phonon_b * ( std::cosh( phi ) - 1.0 );
    }
    return parameters.p_phonon_b * parameters.p_phonon_b * std::sinh( phi );
}

double System::dgl_phonons_lindblad_coefficients( double t, double energy, double coupling, Scalar pulse, const char mode, const double sign ) {
    if ( t == 0.0 )
        t = parameters.t_step * 0.01;
    if ( phi_vector.count( t ) == 0 ) //TODO: get function for phi?
        phi_vector[t] = dgl_phonons_phi( t );
    double ret = 0;
    if ( mode == 'L' ) {
        double bpulsesquared = std::pow( std::abs( parameters.p_phonon_b * pulse ), 2.0 );

        double nu = std::sqrt( bpulsesquared + energy * energy );
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += parameters.t_step ) {
            Scalar f = ( energy * energy * std::cos( nu * tau ) + bpulsesquared ) / std::pow( nu, 2.0 );
            ret += std::real( ( std::cosh( phi_vector[t] ) - 1.0 ) * f + std::sinh( phi_vector[t] ) * std::cos( nu * tau ) ) - sign * std::imag( ( std::exp( phi_vector[t] ) - 1.0 ) * energy * std::sin( nu * tau ) / nu );
            i++;
        }
        ret *= 2.0 * bpulsesquared * parameters.t_step;
    } else if ( mode == 'C' ) {
        int i = 0;
        for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += parameters.t_step ) {
            ret += std::real( std::exp( 1i * sign * energy * tau ) * ( std::exp( phi_vector[t] ) - 1.0 ) );
            i++;
        }
        ret *= parameters.p_phonon_b * parameters.p_phonon_b * coupling * coupling * parameters.t_step;
    }
    return ret;
}

Sparse System::dgl_phonons_chi( const double t ) {
    // Electron-Cavity
    Sparse ret = operatorMatrices.polaron_factors[0];
    // Electron-Pulse
    for ( int p = 1; p < parameters.input_pulse.size() + 1; p++ ) {
        ret += operatorMatrices.polaron_factors[p] * pulse[p - 1].get( t );
    }
    return dgl_timetrafo( ret, t );
    //return dgl_timetrafo( parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H + ( operatorMatrices.atom_sigmaplus_G_H + operatorMatrices.atom_sigmaplus_H_B ) * pulse_H.get( t ) + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V + ( operatorMatrices.atom_sigmaplus_G_V + operatorMatrices.atom_sigmaplus_V_B ) * pulse_V.get( t ), t );
    //return dgl_timetrafo( parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V + parameters.p_omega_coupling * operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V, t );
}

Sparse System::dgl_phonons_calculate_transformation( Sparse &chi_tau, double t, double tau ) {
    if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_BACKWARDS_INTEGRAL ) {
        return Solver::calculate_definite_integral( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.t_step, parameters.numerics_rk_tol, parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() ).mat;
    } else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_TRANSFORMATION_MATRIX ) {
        Sparse U = ( Dense( -1.0i * dgl_getHamilton( t ) * tau ).exp() ).sparseView();
        return ( U * chi_tau * U.adjoint() ).eval();
    } else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_MIXED ) {
        double threshold = 0;
        for ( auto &p : pulse )
            threshold += std::abs( p.get( t ) / p.maximum ); //TODO: + photon zahl, wenn photon number > 0.1 oder so dann auch. für große kopplungen gibts sonst starke abweichungen. vil. number*g draufaddieren.
        if ( threshold > 1E-4 || ( chirp.size() > 0 && chirp.back().derivative( t ) != 0 ) ) { //TODO: threshold als parameter
            return Solver::calculate_definite_integral( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.t_step, parameters.numerics_rk_tol, parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() ).mat;
        }
    }
    return chi_tau;
}

Sparse System::dgl_phonons_pmeq( const Sparse &rho, const double t, const std::vector<SaveState> &past_rhos ) {
    Sparse ret = Sparse( parameters.maxStates, parameters.maxStates );
    Sparse chi = dgl_phonons_chi( t );
    std::vector<Sparse> threadmap_1, threadmap_2;

    // Most precise approximation used. Caclulate polaron fram Chi by integrating backwards from t to t-tau.
    if ( parameters.numerics_phonon_approximation_order != PHONON_APPROXIMATION_LINDBLAD_FULL ) {
        Sparse XUT = dgl_phonons_chiToX( chi, 'u' );
        Sparse XGT = dgl_phonons_chiToX( chi, 'g' );
        int _taumax = (int)std::min( parameters.p_phonon_tcutoff / parameters.t_step, t / parameters.t_step );
        // Index vector for thread ordering
        //std::vector<int> thread_index;
        //thread_index.reserve( _taumax );
        //int thread_increment = std::ceil( _taumax / parameters.numerics_maximum_threads );
        //for ( int thread = 0; thread < parameters.numerics_maximum_threads; thread++ ) {
        //    for ( int cur = 0; cur < parameters.numerics_maximum_threads; cur++ ) {
        //        int vec_index = thread + cur * parameters.numerics_maximum_threads;
        //        if ( vec_index < _taumax ) {
        //            thread_index.emplace_back( vec_index );
        //        } else {
        //            break;
        //        }
        //    }
        //}
        // debug
        //Log::L2("Vec index:\n");
        //for ( auto &a : thread_index) {
        //    Log::L2("{}, ",a);
        //}
        //Log::L2("\nVec index done!\n");
        // Use markov approximation
        if ( parameters.numerics_phonon_approximation_markov1 ) {
            // Temporary variables
            Sparse chi_tau_back_u, chi_tau_back_g, integrant;
            // Look if we already calculated coefficient sum for this specific t-value
            //int index = dgl_get_coefficient_index( t, 0 );
            std::pair<double, double> time_pair = std::make_pair( t, 0.0 );
            // FIXME: multithreading mit der map is RIP wenn RK 45
            if ( parameters.numerics_use_saved_coefficients && savedCoefficients.count( time_pair ) > 0 ) {
                // Index was found, chi(t-tau) sum used from saved vector
                auto &coeff = savedCoefficients[time_pair];
                chi_tau_back_u = coeff.mat1; //savedCoefficients.at( index ).mat1;
                chi_tau_back_g = coeff.mat2; //savedCoefficients.at( index ).mat2;
            } else {
                // Index was not found, (re)calculate chi(t-tau) sum
                // Initialize temporary matrices to zero for threads to write to
                init_sparsevector( threadmap_1, parameters.maxStates, parameters.numerics_phonons_maximum_threads );
                init_sparsevector( threadmap_2, parameters.maxStates, parameters.numerics_phonons_maximum_threads );
                // Calculate backwards integral and sum it into threadmaps. Threadmaps will later be summed into one coefficient matrix.
#pragma omp parallel for ordered schedule( dynamic ) shared( savedCoefficients ) num_threads( parameters.numerics_phonons_maximum_threads )
                for ( int cur_thread = 0; cur_thread < parameters.numerics_phonons_maximum_threads; cur_thread++ ) {
                    for ( int cur = 0; cur < _taumax; cur += parameters.numerics_phonons_maximum_threads ) {
                        int _tau = cur_thread + cur;
                        //for ( int _tau = 0; _tau < _taumax; _tau++ ) {
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
                }
                // Sum all contributions from threadmaps into one coefficient
                chi_tau_back_u = std::accumulate( threadmap_1.begin(), threadmap_1.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
                chi_tau_back_g = std::accumulate( threadmap_2.begin(), threadmap_2.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
                // Save coefficients
                //dgl_save_coefficient( chi_tau_back_u, chi_tau_back_g, t, 0 );
                if ( parameters.numerics_use_saved_coefficients ) {
#pragma omp critical
                    savedCoefficients[time_pair] = SaveStateTau( chi_tau_back_u, chi_tau_back_g, t, 0 );
                }
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
            for ( int cur_thread = 0; cur_thread < parameters.numerics_phonons_maximum_threads; cur_thread++ ) {
                for ( int cur = 0; cur < _taumax; cur += parameters.numerics_phonons_maximum_threads ) {
                    int _tau = cur_thread + cur;
                    //for ( int _tau = 0; _tau < _taumax; _tau++ ) {
                    Sparse chi_tau_back_u, chi_tau_back_g, chi_tau, chi_tau_back, integrant;
                    int rho_index = std::max( 0, (int)past_rhos.size() - 1 - _tau );
                    double tau = ( 1.0 * _tau ) * parameters.t_step;
                    // Check if chi(t-tau) has to be recalculated, or can just be retaken from saved matrices, since chi(t-tau) is only dependant on time
                    //int index = dgl_get_coefficient_index( t, tau );
                    std::pair<double, double> time_pair = std::make_pair( t, tau );
                    if ( parameters.numerics_use_saved_coefficients && savedCoefficients.count( time_pair ) > 0 ) {
// Index was found, chi(t-tau) used from saved vector
#pragma omp critical
                        {
                            auto &coeff = savedCoefficients[time_pair];
                            chi_tau_back_u = coeff.mat1;
                            chi_tau_back_g = coeff.mat2;
                        }
                    } else {
                        // Index not found, or no saving is used. Recalculate chi(t-tau).
                        chi_tau = dgl_phonons_chi( t - tau );
                        chi_tau_back = chi_tau_back = dgl_phonons_calculate_transformation( chi_tau, t, tau );
                        chi_tau_back_u = dgl_phonons_chiToX( chi_tau_back, 'u' );
                        chi_tau_back_g = dgl_phonons_chiToX( chi_tau_back, 'g' );
                        // If saving is used, save current chi(t-tau). If markov approximation is used, only save final contributions
                        //dgl_save_coefficient( chi_tau_back_u, chi_tau_back_g, t, tau );
                        if ( parameters.numerics_use_saved_coefficients ) {
#pragma omp critical
                            savedCoefficients[time_pair] = SaveStateTau( chi_tau_back_u, chi_tau_back_g, t, tau );
                        }
                    }
                    integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( XUT, ( chi_tau_back_u * past_rhos.at( rho_index ).mat ).eval() );
                    integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( XGT, ( chi_tau_back_g * past_rhos.at( rho_index ).mat ).eval() );
                    Sparse adjoint = integrant.adjoint();
                    auto thread = omp_get_thread_num();
                    threadmap_1.at( thread ) += ( integrant + adjoint ) * parameters.t_step;
                }
            }
            ret -= std::accumulate( threadmap_1.begin(), threadmap_1.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
        }
    } else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_LINDBLAD_FULL ) {
        double chirpcorrection = chirp.size() > 0 ? chirp.back().get( t ) : 0;
        for ( auto &[name, mat] : parameters.input_pulse ) {
            for ( int p = 0; p < mat.string_v["CoupledTo"].size(); p++ ) {
                auto &mode = mat.string_v["CoupledTo"][p];
                auto &transition = operatorMatrices.el_transitions[mode].hilbert;
                auto mode_transposed = mode;
                std::reverse( mode_transposed.begin(), mode_transposed.end() );
                auto &transition_transposed = operatorMatrices.el_transitions[mode_transposed].hilbert;
                //std::cout << p << " m1 = " << mode << ", m2 = " << mode_transposed << std::endl;
                auto delta_E = operatorMatrices.el_transitions[mode].energy - mat.numerical_v["Frequency"][p] + chirpcorrection;
                //std::cout << "Pulsed Transition " << mode << " -> coeeff = " << dgl_phonons_lindblad_coefficients( t, delta_E, 0.0, pulse[p].get( t ), 'L', operatorMatrices.el_transitions[mode].direction ) << std::endl;
                ret += dgl_phonons_lindblad_coefficients( t, delta_E, 0.0, pulse[p].get( t ), 'L', -1 ) * dgl_lindblad( rho, transition, transition_transposed );
                ret += dgl_phonons_lindblad_coefficients( t, delta_E, 0.0, pulse[p].get( t ), 'L', 1 ) * dgl_lindblad( rho, transition_transposed, transition );
            }
        }
        for ( auto &[name, mat] : parameters.input_photonic ) {
            for ( auto &mode : mat.string_v["CoupledTo"] ) {
                int c = 0;
                auto &transition = operatorMatrices.el_transitions[mode].hilbert;
                auto mode_transposed = mode;
                std::reverse( mode_transposed.begin(), mode_transposed.end() );
                auto &transition_transposed = operatorMatrices.el_transitions[mode_transposed].hilbert;
                auto &optical_transition = operatorMatrices.ph_transitions[name + "b"].hilbert;
                auto &optical_transition_transposed = operatorMatrices.ph_transitions[name + "bd"].hilbert;
                auto delta_E = operatorMatrices.el_transitions[mode].energy - mat.numerical["Frequency"] + chirpcorrection;
                //std::cout << "Transition " << mode << "-" << name + "bd (" << dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c++] * parameters.p_omega_coupling, 0.0, 'C', -1 ) << "),   " << mode_transposed << "-" << name + "b (" << dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c++] * parameters.p_omega_coupling, 0.0, 'C', 1 ) << ")" << std::endl;
                ret += dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c++] * parameters.p_omega_coupling, 0.0, 'C', -1 ) * dgl_lindblad( rho, transition * optical_transition_transposed, transition_transposed * optical_transition );
                ret += dgl_phonons_lindblad_coefficients( t, delta_E, mat.numerical_v["CouplingScaling"][c++] * parameters.p_omega_coupling, 0.0, 'C', 1 ) * dgl_lindblad( rho, transition_transposed * optical_transition, transition * optical_transition_transposed );
            }
        }
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H + chirpcorrection, 'L', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H, operatorMatrices.atom_sigmaplus_G_H );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H + chirpcorrection, 'L', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_H, operatorMatrices.atom_sigmaminus_G_H );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H + chirpcorrection, 'C', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H * operatorMatrices.photon_create_H, operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H + chirpcorrection, 'C', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H, operatorMatrices.atom_sigmaminus_G_H * operatorMatrices.photon_create_H );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B + chirpcorrection, 'L', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B, operatorMatrices.atom_sigmaplus_H_B );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B + chirpcorrection, 'L', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_H_B, operatorMatrices.atom_sigmaminus_H_B );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B + chirpcorrection, 'C', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B * operatorMatrices.photon_create_H, operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B + chirpcorrection, 'C', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H, operatorMatrices.atom_sigmaminus_H_B * operatorMatrices.photon_create_H );
        //
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V + chirpcorrection, 'L', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V, operatorMatrices.atom_sigmaplus_G_V );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V + chirpcorrection, 'L', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_V, operatorMatrices.atom_sigmaminus_G_V );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V + chirpcorrection, 'C', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V * operatorMatrices.photon_create_V, operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V + chirpcorrection, 'C', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V, operatorMatrices.atom_sigmaminus_G_V * operatorMatrices.photon_create_V );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B + chirpcorrection, 'L', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B, operatorMatrices.atom_sigmaplus_V_B );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B + chirpcorrection, 'L', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_V_B, operatorMatrices.atom_sigmaminus_V_B );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B + chirpcorrection, 'C', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B * operatorMatrices.photon_create_V, operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V );
        //ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B + chirpcorrection, 'C', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V, operatorMatrices.atom_sigmaminus_V_B * operatorMatrices.photon_create_V );
    }
    return ret;
}