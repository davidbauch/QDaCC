#include "System_Exziton2NS_Parameters.cpp"
#include "System_Exziton2NS_OperatorMatrices.cpp"
#include "System_Exziton2NS_FileOutput.cpp"
#include "../system.cpp"
#include "../../chirp.h"
#include "../../pulse.h"
#include "../../solver.h"

class System : public System_Parent {
   public:
    // System Components
    Chirp chirp;
    Pulse pulse;
    FileOutput fileoutput;

    MatrixXcd timeTrafoMatrix;

    System(){};
    System( const std::vector<std::string> &input ) {
        // Set Name of this system. 2 Level atomic system coupled to single mode optical field
        name = "Exziton (2NS)";
        logs.level2( "Creating System Class for '{}'\n", name );
        parameters = Parameters( input );
        parameters.log();
        operatorMatrices = OperatorMatrices( parameters );
        operatorMatrices.outputOperators( parameters );
        fileoutput = FileOutput( {"densitymatrix.txt", "atomicinversion.txt", "photonpopulation.txt"}, parameters );
        init();
    }

    bool init_system() {
        logs.inBar( "Program Log" );
        // Single chirp for single atomic level
        Chirp::Inputs chirpinputs( parameters.t_start, parameters.t_end, parameters.t_step, parameters.chirp_type, parameters.numerics_order_highest );
        chirpinputs.add( parameters.chirp_t, parameters.chirp_y, parameters.chirp_ddt );
        chirp = Chirp( chirpinputs );
        if ( parameters.chirp_total != 0 )
            chirp.fileOutput( parameters.subfolder + "chirp.txt" );

        // Arbitrary number of pulses onto single atomic level. As of now, only a single pulse input is supported via input.
        Pulse::Inputs pulseinputs( parameters.t_start, parameters.t_end, parameters.t_step, parameters.numerics_order_highest );
        pulseinputs.add( parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_type );
        pulse = Pulse( pulseinputs );
        if ( parameters.pulse_amp.at( 0 ) != 0 )
            pulse.fileOutput( parameters.subfolder + "pulse.txt" );

        // Output Phonon functions of phonons are active
        if ( parameters.p_phonon_T != 0 ) {
            FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
            fmt::print( fp_phonons, "t\tphi(t)\tg_u(t)\tg_g(t)\n" );
            for ( double t = getTimeborderStart() + getTimeStep(); t < 3 * parameters.p_photon_tcutoff; t += getTimeStep() ) {
                fmt::print( fp_phonons, "{}\t{}\t{}\t{}\n", t, std::abs( dgl_phonons_phi( t ) ), std::abs( dgl_phonons_greenf( t, 'u' ) ), std::abs( dgl_phonons_greenf( t, 'g' ) ) );
            }
            std::fclose( fp_phonons );
        }

        // Time Transformation
        timeTrafoMatrix = ( 1i * operatorMatrices.H_0 ).exp();
        // Check time trafo
        MatrixXcd ttrafo = ( 1i * operatorMatrices.H_0 * 125E-12 ).exp() * operatorMatrices.H_used * ( -1i * operatorMatrices.H_0 * 125E-12 ).exp();
        MatrixXcd temp = dgl_timetrafo( operatorMatrices.H_used, 125E-12 ) - ttrafo;
        if ( std::abs( temp.sum() / parameters.p_omega_atomic ) >= 0.0001 ) {
            logs( "Unitary timetransformation error is bigger than 0.01% of atomic transition energy!\n\n" );
        }
        return true;
    }

    MatrixXcd dgl_rungeFunction( const MatrixXcd &rho, const MatrixXcd &H, const double t ) {
        MatrixXcd ret = -1i * dgl_kommutator( H, rho );
        // Photone losses
        if ( parameters.p_omega_cavity_loss != 0.0 ) {
            ret += parameters.p_omega_cavity_loss / 2.0 * ( 2 * operatorMatrices.photon_annihilate * rho * operatorMatrices.photon_create - dgl_antikommutator( operatorMatrices.photon_n, rho ) ); /* ret + p_omega_cavity_loss/2*(2*b*rho*b^+ - [b^+b,rho]) */
            //ret += parameters.p_omega_cavity_loss/2.0*dgl_lindblad(rho,operatorMatrices.photon_annihilate,operatorMatrices.photon_create);
        }
        // Pure Dephasing
        if ( parameters.p_omega_pure_dephasing != 0.0 ) {
            ret -= parameters.p_omega_pure_dephasing / 2.0 * ( operatorMatrices.atom_ground * rho * operatorMatrices.atom_exited + operatorMatrices.atom_exited * rho * operatorMatrices.atom_ground ); /* -p_omega_pure_dephasing/2*(|g><g|rho|e><e| + |e><e|rho|g><g|) */
            //ret += parameters.p_omega_pure_dephasing/2.0*dgl_lindblad(rho,operatorMatrices.atom_exited,operatorMatrices.atom_ground);
        }
        if ( parameters.p_omega_decay != 0.0 ) {
            ret += parameters.p_omega_decay * dgl_lindblad( rho, operatorMatrices.atom_sigmaminus, operatorMatrices.atom_sigmaplus );
        }
        if ( parameters.p_phonon_T != 0.0 ) {
            ret -= dgl_phonons( rho, t );
        }
        return ret;
    }

    MatrixXcd dgl_timetrafo( const MatrixXcd &A, const double t ) {
        MatrixXcd ret = A;
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
                            ret( n, m ) = A( n, m ) * std::exp( 1i * t * ( ( parameters.p_omega_atomic ) / 2. * ( delta( i, 1 ) - delta( i, 0 ) - delta( j, 1 ) + delta( j, 0 ) ) + parameters.p_omega_cavity * ( pn - pm ) ) );
                        }
                    }
                }
            }
            // TIMETRANSFORMATION_MATRIXEXPONENTIAL
            else if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_MATRIXEXPONENTIAL ) {
                MatrixXcd U = ( 1i * operatorMatrices.H_0 * t ).exp(); // auto slower than direct type, because of doubled evaluation
                ret = U * A * U.conjugate();
            } else if ( parameters.numerics_order_timetrafo == TIMETRANSFORMATION_MATRIXEXPONENTIAL_PRECALCULATED ) {
                MatrixXcd U = timeTrafoMatrix.array().pow( t ); // auto slower than direct type, because of doubled evaluation
                ret = U * A * U.conjugate();
            }
        }
        return ret;
    }

    MatrixXcd dgl_chirp( const double t ) {
        if ( parameters.chirp_total == 0 )
            return MatrixXcd::Zero( parameters.maxStates, parameters.maxStates );
        return 0.5 * operatorMatrices.atom_inversion * chirp.get( t );
    }
    MatrixXcd dgl_pulse( const double t ) {
        if ( parameters.pulse_amp.at( 0 ) == 0 )
            return MatrixXcd::Zero( parameters.maxStates, parameters.maxStates );
        return 0.5 * ( operatorMatrices.atom_sigmaplus * pulse.get( t ) + operatorMatrices.atom_sigmaminus * std::conj( pulse.get( t ) ) );
    }

    MatrixXcd dgl_phonons_rungefunc( const MatrixXcd &chi, const double t ) {
        MatrixXcd explicit_time = 1i * parameters.p_omega_atomic * operatorMatrices.atom_sigmaplus * parameters.p_omega_coupling * operatorMatrices.photon_annihilate + 1i * parameters.p_omega_atomic * operatorMatrices.atom_sigmaplus * pulse.get( t ) - 1i * operatorMatrices.atom_sigmaplus * parameters.p_omega_cavity * parameters.p_omega_coupling * operatorMatrices.photon_annihilate + operatorMatrices.atom_sigmaplus * ( pulse.get( t ) - pulse.get( std::max( 0.0, t - parameters.t_step ) ) ) / parameters.t_step;
        MatrixXcd hamilton = parameters.p_phonon_b * dgl_getHamilton( t );
        return -1i * dgl_kommutator( hamilton, chi ) + explicit_time;
    }

    MatrixXcd dgl_phonons_chiToX( const MatrixXcd &chi, const char mode = 'u' ) {
        if ( mode == 'g' ) {
            return chi + chi.adjoint().eval();
        }
        return 1i * ( chi - chi.adjoint().eval() );
    }

    std::complex<double> dgl_phonons_greenf( const double t, const char mode = 'u' ) {
        auto phi = dgl_phonons_phi( t );
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

    MatrixXcd dgl_phonons( const MatrixXcd &rho, const double t ) {
        // Determine Chi
        MatrixXcd ret = MatrixXcd::Zero( parameters.maxStates, parameters.maxStates );
        MatrixXcd chi = operatorMatrices.atom_sigmaplus * parameters.p_omega_coupling * operatorMatrices.photon_annihilate + operatorMatrices.atom_sigmaplus * pulse.get( t );
        MatrixXcd integrant;

        if ( parameters.numerics_phonon_approximation == PHONON_APPROXIMATION_FULL ) {
            // Calculate Chi(t) backwards to Chi(t-tau)
            MatrixXcd chit = dgl_timetrafo( chi, t );
            std::vector<ODESolver::SaveState> chis = ODESolver::calculate_definite_integral_vec( chit, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( 0.0, t - parameters.p_photon_tcutoff ), -parameters.t_step );
            //fmt::print("Times for chi ({}): \n",t);
            //for ( int i = 1; i < (int)( chis.size() ); i++ ) {
            //    fmt::print("\ti = {} -> t = {}\n",i,chis.at(i).t);
            //}
            for ( int i = 1; i < (int)( chis.size() ); i++ ) {
                double tau = parameters.t_step * (double)i;
                integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( dgl_phonons_chiToX( chis.at( 0 ).mat, 'u' ), ( dgl_phonons_chiToX( chis.at( i ).mat, 'u' ) * rho ).eval() );
                integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( dgl_phonons_chiToX( chis.at( 0 ).mat, 'g' ), ( dgl_phonons_chiToX( chis.at( i ).mat, 'g' ) * rho ).eval() );
                //integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( dgl_phonons_chiToX( dgl_timetrafo( chi, t ), 'u' ), ( dgl_phonons_chiToX( dgl_timetrafo( chi_tau, t - tau ), 'u' ) * rho ).eval() );
                //integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( dgl_phonons_chiToX( dgl_timetrafo( chi, t ), 'g' ), ( dgl_phonons_chiToX( dgl_timetrafo( chi_tau, t - tau ), 'g' ) * rho ).eval() );
                ret += ( integrant + integrant.adjoint().eval() ) * parameters.t_step;
            }
        } else if ( parameters.numerics_phonon_approximation == PHONON_APPROXIMATION_MARKOV_1 ) {
            for ( int i = 1; i < (int)( parameters.p_photon_tcutoff / parameters.t_step ); i++ ) {
                double tau = parameters.t_step * (double)i;
                MatrixXcd chi_tau = operatorMatrices.atom_sigmaplus * parameters.p_omega_coupling * operatorMatrices.photon_annihilate + operatorMatrices.atom_sigmaplus * pulse.get( t - tau );
                //integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( dgl_phonons_chiToX( chis.at( 0 ).mat, 'u' ), (dgl_phonons_chiToX( chis.at( i ).mat, 'u' )*rho).eval() );
                //integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( dgl_phonons_chiToX( chis.at( 0 ).mat, 'g' ), (dgl_phonons_chiToX( chis.at( i ).mat, 'g' )*rho).eval() );
                integrant = dgl_phonons_greenf( tau, 'u' ) * dgl_kommutator( dgl_phonons_chiToX( dgl_timetrafo( chi, t ), 'u' ), ( dgl_phonons_chiToX( dgl_timetrafo( chi_tau, t - tau ), 'u' ) * rho ).eval() );
                integrant += dgl_phonons_greenf( tau, 'g' ) * dgl_kommutator( dgl_phonons_chiToX( dgl_timetrafo( chi, t ), 'g' ), ( dgl_phonons_chiToX( dgl_timetrafo( chi_tau, t - tau ), 'g' ) * rho ).eval() );
                ret += ( integrant + integrant.adjoint().eval() ) * parameters.t_step;
            }
        } else if ( parameters.numerics_phonon_approximation == PHONON_APPROXIMATION_MARKOV_2 ) {
        }
        return ret;
    }

    void expectationValues( const MatrixXcd &rho, const double t ) {
        std::fprintf( fileoutput.fp_atomicinversion, "%.15e\t%.15e\n", t, std::real( dgl_expectationvalue( rho, operatorMatrices.atom_inversion, t ) ) );
        std::fprintf( fileoutput.fp_photonpopulation, "%.15e\t%.15e\n", t, std::real( dgl_expectationvalue( rho, operatorMatrices.photon_n, t ) ) );
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

    MatrixXcd dgl_getHamilton( const double t ) {
        return dgl_timetrafo( operatorMatrices.H_used + dgl_chirp( t ) + dgl_pulse( t ), t );
    }

    bool exit_system( const int failure = 0 ) {
        pulse.log();
        chirp.log();
        fileoutput.close();
        return true;
    }
};