#include "System_Exciton4NS_Parameters.cpp"
#include "System_Exciton4NS_OperatorMatrices.cpp"
#include "System_Exciton4NS_FileOutput.cpp"
#include "../system.cpp"
#include "../../chirp.h"
#include "../../pulse.h"
#include "../../solver.h"

class System : public System_Parent {
   public:
    // System Components
    Chirp chirp_H;
    Chirp chirp_V;
    Pulse pulse_H;
    Pulse pulse_V;
    FileOutput fileoutput;
    std::vector<std::complex<double>> phi_vector;

    DenseMat timeTrafoMatrix;

    System(){};
    System( const std::vector<std::string> &input ) {
        // Set Name of this system. 2 Level atomic system coupled to single mode optical field
        name = "Biexciton (4NS)";
        logs.level2( "Creating System Class for '{}'\n", name );
        parameters = Parameters( input );
        parameters.log();
        operatorMatrices = OperatorMatrices( parameters );
        operatorMatrices.outputOperators( parameters );
        fileoutput = FileOutput( {"densitymatrix.txt", "atomicinversion.txt", "photonpopulation.txt"}, parameters, operatorMatrices );
        init();
    }

    bool init_system() {
        // Single chirp for single atomic level
        Chirp::Inputs chirpinputs( parameters.t_start, parameters.t_end, parameters.t_step, parameters.chirp_type, parameters.numerics_order_highest );
        chirpinputs.add( parameters.chirp_t, parameters.chirp_y, parameters.chirp_ddt );
        chirp_H = Chirp( chirpinputs );
        chirp_V = Chirp( chirpinputs );
        if ( parameters.chirp_total != 0 )
            chirp_H.fileOutput( parameters.subfolder + "chirp.txt" );

        // Arbitrary number of pulses onto single atomic level. As of now, only a single pulse input is supported via input.
        Pulse::Inputs pulseinputs( parameters.t_start, parameters.t_end, parameters.t_step, parameters.numerics_order_highest );
        pulseinputs.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_type );
        pulse_H = Pulse( pulseinputs );
        pulse_V = Pulse( pulseinputs );
        if ( parameters.pulse_amp.at( 0 ) != 0 )
            pulse_H.fileOutput( parameters.subfolder + "pulse.txt" );

        // Output Phonon functions if phonons are active
        if ( parameters.p_phonon_T != 0 ) {
            phi_vector.reserve( std::ceil( parameters.p_phonon_tcutoff / getTimeStep() * 3.1 ) );
            phi_vector.emplace_back( dgl_phonons_phi( getTimeStep() * 0.01 ) );
            for ( double tau = getTimeStep(); tau < parameters.p_phonon_tcutoff * 3.1; tau += getTimeStep() ) {
                phi_vector.emplace_back( dgl_phonons_phi( tau ) );
            }
            FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
            fmt::print( fp_phonons, "t\treal(phi(t))\timag(phi(t))\treal(g_u(t))\timag(g_u(t))\treal(g_g(t))\timag(g_g(t))\n" );
            for ( double t = getTimeborderStart(); t < 3.0 * parameters.p_phonon_tcutoff; t += getTimeStep() ) {
                auto greenu = dgl_phonons_greenf( t, 'u' );
                auto greeng = dgl_phonons_greenf( t, 'g' );
                fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", t, std::real( phi_vector.at( std::floor( t / getTimeStep() ) ) ), std::imag( phi_vector.at( std::floor( t / getTimeStep() ) ) ), std::real( greenu ), std::imag( greenu ), std::real( greeng ), std::imag( greeng ) );
            }
            std::fclose( fp_phonons );
            if ( parameters.output_coefficients ) {
                fp_phonons = std::fopen( ( parameters.subfolder + "phonons_lb.txt" ).c_str(), "w" );
                fmt::print( fp_phonons, "t\tL_a_+\tL_a_-\tL_c_+\tL_c_-\n" );
                for ( double t = getTimeborderStart(); t < getTimeborderEnd(); t += getTimeStep() ) {
                    fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\n", t, dgl_phonons_lindblad_coefficients( t, 'L', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'L', -1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', -1.0 ) );
                }
                std::fclose( fp_phonons );
            }
        }

        // Time Transformation
        timeTrafoMatrix = ( 1i * operatorMatrices.H_0 ).exp();
        // Check time trafo
        DenseMat ttrafo = ( 1i * operatorMatrices.H_0 * 125E-12 ).exp() * operatorMatrices.H_used * ( -1i * operatorMatrices.H_0 * 125E-12 ).exp();
        DenseMat temp = dgl_timetrafo( operatorMatrices.H_used, 125E-12 ) - ttrafo;
        if ( std::abs( temp.sum() / parameters.p_omega_atomic_G_H ) >= 0.0001 ) {
            logs( "Unitary timetransformation error is bigger than 0.01% of atomic transition energy!\n\n" );
        }
        return true;
    }

    DenseMat dgl_rungeFunction( const DenseMat &rho, const DenseMat &H, const double t, std::vector<SaveState> &past_rhos ) {
        DenseMat ret = -1i * dgl_kommutator( H, rho );
        // Photone losses
        if ( parameters.p_omega_cavity_loss != 0.0 ) {
            //ret += parameters.p_omega_cavity_loss / 2.0 * ( 2 * operatorMatrices.photon_annihilate_H * rho * operatorMatrices.photon_create_H - dgl_antikommutator( operatorMatrices.photon_n_H, rho ) );
            ret += parameters.p_omega_cavity_loss * dgl_lindblad( rho, operatorMatrices.photon_annihilate_H, operatorMatrices.photon_create_H );
            ret += parameters.p_omega_cavity_loss * dgl_lindblad( rho, operatorMatrices.photon_annihilate_V, operatorMatrices.photon_create_V );
        }
        // Pure Dephasing
        if ( parameters.p_omega_pure_dephasing != 0.0 ) {
            ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_ground * rho * operatorMatrices.atom_state_H + operatorMatrices.atom_state_H * rho * operatorMatrices.atom_state_ground );
            ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_ground * rho * operatorMatrices.atom_state_V + operatorMatrices.atom_state_V * rho * operatorMatrices.atom_state_ground );
            ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_H * rho * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_biexciton * rho * operatorMatrices.atom_state_H );
            ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_state_V * rho * operatorMatrices.atom_state_biexciton + operatorMatrices.atom_state_biexciton * rho * operatorMatrices.atom_state_V );
        }
        // Radiative decay
        if ( parameters.p_omega_decay != 0.0 ) {
            ret += parameters.p_omega_pure_dephasing / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H, operatorMatrices.atom_sigmaplus_G_H );
            ret += parameters.p_omega_pure_dephasing / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V, operatorMatrices.atom_sigmaplus_G_V );
            ret += parameters.p_omega_pure_dephasing / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B, operatorMatrices.atom_sigmaplus_H_B );
            ret += parameters.p_omega_pure_dephasing / 2.0 * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B, operatorMatrices.atom_sigmaplus_V_B );
            //ret += parameters.p_phonon_b * parameters.p_phonon_b * parameters.p_omega_decay * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus, operatorMatrices.atom_sigmaplus );
        }
        if ( parameters.p_phonon_T != 0.0 ) {
            ret += dgl_phonons( rho, t, past_rhos );
        }
        return ret;
    }

    DenseMat dgl_timetrafo( const DenseMat &A, const double t ) {
        DenseMat ret = A;
        if ( parameters.numerics_use_interactionpicture == 1 ) {
            // TIMETRANSFORMATION_PRECALCULATED
            if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_PRECALCULATED ) {
                int i, j, pn, pm;
                for ( int n = 0; n < A.rows(); n++ ) {
                    i = n % 2;
                    pn = (int)( n / 2 );
                    for ( int m = 0; m < A.cols(); m++ ) {
                        if ( A( n, m ) != 0.0 ) {
                            j = m % 2;
                            pm = (int)( m / 2 );
                            // FIXME: calculate analytical transformatiuon!
                            //ret( n, m ) = A( n, m ) * std::exp( 1i * t * ( ( parameters.p_omega_atomic_G_H ) / 2. * ( delta( i, 1 ) - delta( i, 0 ) - delta( j, 1 ) + delta( j, 0 ) ) + parameters.p_omega_cavity_H * ( pn - pm ) ) );
                        }
                    }
                }
            }
            // TIMETRANSFORMATION_MATRIXEXPONENTIAL
            else if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_MATRIXEXPONENTIAL ) {
                DenseMat U = ( 1i * operatorMatrices.H_0 * t ).exp(); // auto slower than direct type, because of doubled evaluation
                ret = U * A * U.conjugate();
            } else if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_MATRIXEXPONENTIAL_PRECALCULATED ) {
                DenseMat U = timeTrafoMatrix.array().pow( t ); // auto slower than direct type, because of doubled evaluation
                ret = U * A * U.conjugate();
            }
        }
        return ret;
    }

    DenseMat dgl_chirp( const double t ) {
        if ( parameters.chirp_total == 0 )
            return DenseMat::Zero( parameters.maxStates, parameters.maxStates );
        return 0.5 * operatorMatrices.atom_inversion_G_H * chirp_H.get( t );
    }
    DenseMat dgl_pulse( const double t ) {
        if ( parameters.pulse_amp.at( 0 ) == 0 )
            return DenseMat::Zero( parameters.maxStates, parameters.maxStates );
        return 0.5 * ( operatorMatrices.atom_sigmaplus_G_H * pulse_H.get( t ) + operatorMatrices.atom_sigmaminus_G_H * std::conj( pulse_H.get( t ) ) );
    }

    DenseMat dgl_phonons_rungefunc( const DenseMat &chi, const double t ) {
        // H
        DenseMat explicit_time = 1i * parameters.p_omega_atomic_G_H * parameters.p_omega_coupling * ( project_matrix( operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H ).cwiseProduct( chi ) );
        explicit_time += 1i * parameters.p_omega_atomic_H_B * parameters.p_omega_coupling * ( project_matrix( operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H ).cwiseProduct( chi ) );
        explicit_time += 1i * parameters.p_omega_atomic_G_H * pulse_H.get( t ) * ( project_matrix( operatorMatrices.atom_sigmaplus_G_H ).cwiseProduct( chi ) );
        explicit_time -= 1i * parameters.p_omega_coupling * parameters.p_omega_cavity_H * ( project_matrix( operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H ).cwiseProduct( chi ) );
        explicit_time -= 1i * parameters.p_omega_coupling * parameters.p_omega_cavity_H * ( project_matrix( operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H ).cwiseProduct( chi ) );
        explicit_time += 1i * ( pulse_H.get( t ) - pulse_H.get( std::max( 0.0, t - parameters.t_step ) ) ) / parameters.t_step * ( project_matrix( operatorMatrices.atom_sigmaplus_G_H ).cwiseProduct( chi ) );
        explicit_time += 1i * ( pulse_H.get( t ) - pulse_H.get( std::max( 0.0, t - parameters.t_step ) ) ) / parameters.t_step * ( project_matrix( operatorMatrices.atom_sigmaplus_H_B ).cwiseProduct( chi ) );
        // V
        explicit_time += 1i * parameters.p_omega_atomic_G_V * parameters.p_omega_coupling * ( project_matrix( operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V ).cwiseProduct( chi ) );
        explicit_time += 1i * parameters.p_omega_atomic_V_B * parameters.p_omega_coupling * ( project_matrix( operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V ).cwiseProduct( chi ) );
        explicit_time += 1i * parameters.p_omega_atomic_G_V * pulse_V.get( t ) * ( project_matrix( operatorMatrices.atom_sigmaplus_G_V ).cwiseProduct( chi ) );
        explicit_time -= 1i * parameters.p_omega_coupling * parameters.p_omega_cavity_V * ( project_matrix( operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V ).cwiseProduct( chi ) );
        explicit_time -= 1i * parameters.p_omega_coupling * parameters.p_omega_cavity_V * ( project_matrix( operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V ).cwiseProduct( chi ) );
        explicit_time += 1i * ( pulse_V.get( t ) - pulse_V.get( std::max( 0.0, t - parameters.t_step ) ) ) / parameters.t_step * ( project_matrix( operatorMatrices.atom_sigmaplus_G_V ).cwiseProduct( chi ) );
        explicit_time += 1i * ( pulse_V.get( t ) - pulse_V.get( std::max( 0.0, t - parameters.t_step ) ) ) / parameters.t_step * ( project_matrix( operatorMatrices.atom_sigmaplus_V_B ).cwiseProduct( chi ) );
        DenseMat hamilton = dgl_getHamilton( t );
        return -1i * dgl_kommutator( hamilton, chi ) + explicit_time;
    }

    DenseMat dgl_phonons_chiToX( const DenseMat &chi, const char mode = 'u' ) {
        if ( mode == 'g' ) {
            return chi + chi.adjoint().eval();
        }
        return 1i * ( chi - chi.adjoint().eval() );
    }

    std::complex<double> dgl_phonons_greenf( double t, const char mode = 'u' ) {
        int i = std::floor( t / getTimeStep() );
        auto phi = phi_vector.at( i ); //dgl_phonons_phi( t );
        if ( mode == 'g' ) {
            return parameters.p_phonon_b * parameters.p_phonon_b * ( std::cosh( phi ) - 1.0 );
        }
        return parameters.p_phonon_b * parameters.p_phonon_b * std::sinh( phi );
    }

    std::complex<double> dgl_phonons_phi( const double t ) {
        std::complex<double> integral = 0;
        double stepsize = 0.01 * parameters.p_phonon_wcutoff;
        for ( double w = stepsize; w < 10 * parameters.p_phonon_wcutoff; w += stepsize ) {
            integral += stepsize * ( parameters.p_phonon_alpha * w * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff ) * ( std::cos( w * t ) / std::tanh( 1.0545718E-34 * w / 2.0 / 1.3806488E-23 / parameters.p_phonon_T ) - 1i * std::sin( w * t ) ) );
        }
        return integral;
    }

    double dgl_phonons_lindblad_coefficients( double t, double omega_atomic, const char mode = 'L', const char level = 'H', const double sign = 1.0 ) {
        if ( t == 0.0 )
            t = parameters.t_step * 0.01;
        double ret = 0;
        double step = parameters.t_step;
        if ( mode == 'L' ) {
            double bpulsesquared, delta, nu;
            if ( level == 'H' ) {
                bpulsesquared = std::pow( std::abs( parameters.p_phonon_b * pulse_H.get( t ) ), 2.0 );
                delta = parameters.pulse_omega.at( 0 ) - omega_atomic; //FIXME : different pulse frequencies
            } else {
                bpulsesquared = std::pow( std::abs( parameters.p_phonon_b * pulse_V.get( t ) ), 2.0 );
                delta = parameters.pulse_omega.at( 0 ) - omega_atomic; //FIXME : different pulse frequencies
            }
            nu = std::sqrt( bpulsesquared + delta * delta );
            int i = 0;
            for ( double tau = 0; tau < parameters.p_phonon_tcutoff; tau += step ) {
                std::complex<double> f = ( delta * delta * std::cos( nu * tau ) + bpulsesquared ) / std::pow( nu, 2.0 );
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

    DenseMat dgl_phonons( const DenseMat &rho, const double t, const std::vector<SaveState> &past_rhos ) {
        // Determine Chi
        DenseMat ret = DenseMat::Zero( parameters.maxStates, parameters.maxStates );
        DenseMat chi = dgl_timetrafo( operatorMatrices.atom_sigmaplus_G_H * parameters.p_omega_coupling * operatorMatrices.photon_annihilate_H + operatorMatrices.atom_sigmaplus_H_B * parameters.p_omega_coupling * operatorMatrices.photon_annihilate_H + ( operatorMatrices.atom_sigmaplus_G_H + operatorMatrices.atom_sigmaplus_H_B ) * pulse_H.get( t ), t );
        DenseMat integrant = ret;
        DenseMat rho_used = rho;

        if ( parameters.numerics_phonon_approximation_2 == PHONON_APPROXIMATION_BACKWARDS_INTEGRAL ) {
            // Calculate Chi(t) backwards to Chi(t-tau)
            for ( double tau = 0; tau < std::min( parameters.p_phonon_tcutoff, t ); tau += parameters.t_step ) { //for ( auto chit : chis ) {
                if ( !parameters.numerics_phonon_approximation_1 ) {
                    int tau_index = tau / parameters.t_step;
                    rho_used = past_rhos.at( std::max( 0, (int)past_rhos.size() - 1 - tau_index ) ).mat;
                }
                DenseMat chi_tau = dgl_timetrafo( operatorMatrices.atom_sigmaplus_G_H * parameters.p_omega_coupling * operatorMatrices.photon_annihilate_H + operatorMatrices.atom_sigmaplus_H_B * parameters.p_omega_coupling * operatorMatrices.photon_annihilate_H + ( operatorMatrices.atom_sigmaplus_G_H + operatorMatrices.atom_sigmaplus_H_B ) * pulse_H.get( t - tau ), t - tau );
                DenseMat chi_tau_back = ODESolver::calculate_definite_integral( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), -parameters.t_step ).mat;
                integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( dgl_phonons_chiToX( chi, 'u' ), ( dgl_phonons_chiToX( chi_tau_back, 'u' ) * rho_used ).eval() );
                integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( dgl_phonons_chiToX( chi, 'g' ), ( dgl_phonons_chiToX( chi_tau_back, 'g' ) * rho_used ).eval() );
                ret -= ( integrant + integrant.adjoint().eval() ) * parameters.t_step;
            }
        }
        if ( parameters.numerics_phonon_approximation_2 == PHONON_APPROXIMATION_TRANSFORMATION_MATRIX ) {
            for ( double tau = 0; tau < std::min( parameters.p_phonon_tcutoff, t ); tau += parameters.t_step ) {
                if ( !parameters.numerics_phonon_approximation_1 ) {
                    int tau_index = tau / parameters.t_step;
                    rho_used = past_rhos.at( std::max( 0, (int)past_rhos.size() - 1 - tau_index ) ).mat;
                }
                DenseMat chi_tau = dgl_timetrafo( operatorMatrices.atom_sigmaplus_G_H * parameters.p_omega_coupling * operatorMatrices.photon_annihilate_H + operatorMatrices.atom_sigmaplus_H_B * parameters.p_omega_coupling * operatorMatrices.photon_annihilate_H + ( operatorMatrices.atom_sigmaplus_G_H + operatorMatrices.atom_sigmaplus_H_B ) * pulse_H.get( t - tau ), t - tau );
                DenseMat U = ( -1i * dgl_getHamilton( t ) * tau ).exp();
                integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( dgl_phonons_chiToX( chi, 'u' ), ( dgl_phonons_chiToX( U * chi_tau * U.adjoint().eval(), 'u' ) * rho_used ).eval() );
                integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( dgl_phonons_chiToX( chi, 'g' ), ( dgl_phonons_chiToX( U * chi_tau * U.adjoint().eval(), 'g' ) * rho_used ).eval() );
                ret -= ( integrant + integrant.adjoint().eval() ) * parameters.t_step;
            }
        } else if ( parameters.numerics_phonon_approximation_2 == PHONON_APPROXIMATION_TIMETRANSFORMATION ) {
            for ( double tau = 0; tau < std::min( parameters.p_phonon_tcutoff, t ); tau += parameters.t_step ) {
                if ( !parameters.numerics_phonon_approximation_1 ) {
                    int tau_index = tau / parameters.t_step;
                    rho_used = past_rhos.at( std::max( 0, (int)past_rhos.size() - 1 - tau_index ) ).mat;
                }
                DenseMat chi_tau = dgl_timetrafo( operatorMatrices.atom_sigmaplus_G_H * parameters.p_omega_coupling * operatorMatrices.photon_annihilate_H + operatorMatrices.atom_sigmaplus_H_B * parameters.p_omega_coupling * operatorMatrices.photon_annihilate_H + ( operatorMatrices.atom_sigmaplus_G_H + operatorMatrices.atom_sigmaplus_H_B ) * pulse_H.get( t - tau ), t - tau );
                integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( dgl_phonons_chiToX( chi, 'u' ), ( dgl_phonons_chiToX( chi_tau, 'u' ) * rho_used ).eval() );
                integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( dgl_phonons_chiToX( chi, 'g' ), ( dgl_phonons_chiToX( chi_tau, 'g' ) * rho_used ).eval() );
                ret -= ( integrant + integrant.adjoint().eval() ) * parameters.t_step;
            }
        } else if ( parameters.numerics_phonon_approximation_2 == PHONON_APPROXIMATION_LINDBLAD_FULL ) {
            // H
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H, 'L', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H, operatorMatrices.atom_sigmaplus_G_H );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H, 'L', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_H, operatorMatrices.atom_sigmaminus_G_H );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H, 'C', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_H * operatorMatrices.photon_create_H, operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_H, 'C', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_H * operatorMatrices.photon_annihilate_H, operatorMatrices.atom_sigmaminus_G_H * operatorMatrices.photon_create_H );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B, 'L', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B, operatorMatrices.atom_sigmaplus_H_B );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B, 'L', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_H_B, operatorMatrices.atom_sigmaminus_H_B );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B, 'C', 'H', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_H_B * operatorMatrices.photon_create_H, operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_H_B, 'C', 'H', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_H_B * operatorMatrices.photon_annihilate_H, operatorMatrices.atom_sigmaminus_H_B * operatorMatrices.photon_create_H );
            // V
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V, 'L', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V, operatorMatrices.atom_sigmaplus_G_V );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V, 'L', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_V, operatorMatrices.atom_sigmaminus_G_V );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V, 'C', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_G_V * operatorMatrices.photon_create_V, operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_G_V, 'C', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_G_V * operatorMatrices.photon_annihilate_V, operatorMatrices.atom_sigmaminus_G_V * operatorMatrices.photon_create_V );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B, 'L', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B, operatorMatrices.atom_sigmaplus_V_B );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B, 'L', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_V_B, operatorMatrices.atom_sigmaminus_V_B );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B, 'C', 'V', -1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus_V_B * operatorMatrices.photon_create_V, operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V );
            ret += dgl_phonons_lindblad_coefficients( t, parameters.p_omega_atomic_V_B, 'C', 'V', 1.0 ) * dgl_lindblad( rho, operatorMatrices.atom_sigmaplus_V_B * operatorMatrices.photon_annihilate_V, operatorMatrices.atom_sigmaminus_V_B * operatorMatrices.photon_create_V );
        }
        return ret;
    }

    void expectationValues( const DenseMat &rho, const double t ) {
        std::fprintf( fileoutput.fp_atomicinversion, "%.15e\t%.15e\n", t, std::real( dgl_expectationvalue( rho, operatorMatrices.atom_inversion_G_B, t ) ) );
        std::fprintf( fileoutput.fp_photonpopulation, "%.15e\t%.15e\n", t, std::real( dgl_expectationvalue( rho, operatorMatrices.photon_n_H, t ) ) );
        std::fprintf( fileoutput.fp_densitymatrix, "%e\t", t );
        if ( parameters.output_full_dm ) {
            for ( int i = 0; i < parameters.maxStates; i++ )
                for ( int j = 0; j < parameters.maxStates; j++ ) {
                    std::fprintf( fileoutput.fp_densitymatrix, "%e\t", std::real( rho( i, j ) ) );
                }
        } else
            for ( int j = 0; j < parameters.maxStates; j++ )
                std::fprintf( fileoutput.fp_densitymatrix, "%e\t", std::real( rho( j, j ) ) );
        std::fprintf( fileoutput.fp_densitymatrix, "\n" );
    }

    DenseMat dgl_getHamilton( const double t ) {
        return dgl_timetrafo( parameters.p_phonon_b * ( operatorMatrices.H_used + dgl_chirp( t ) + dgl_pulse( t ) ), t );
    }

    bool exit_system( const int failure = 0 ) {
        pulse_H.log();
        chirp_H.log();
        fileoutput.close();
        return true;
    }
};