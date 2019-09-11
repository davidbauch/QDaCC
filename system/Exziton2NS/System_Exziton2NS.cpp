#include "System_Exziton2NS_Parameters.cpp"
#include "System_Exziton2NS_OperatorMatrices.cpp"
#include "System_Exziton2NS_FileOutput.cpp"
#include "../system.cpp"
#include "../../chirp.h" 
#include "../../pulse.h"

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
        // Single chirp for single atomic level
        Chirp::Inputs chirpinputs(parameters.t_start, parameters.t_end, parameters.t_step, parameters.chirp_type, parameters.numerics_order_highest);
        chirpinputs.add(parameters.chirp_t, parameters.chirp_y, parameters.chirp_ddt);
        chirp = Chirp(chirpinputs);
        if ( parameters.chirp_total != 0 )
            chirp.fileOutput( parameters.subfolder + "chirp.txt" );

        // Arbitrary number of pulses onto single atomic level. As of now, only a single pulse input is supported via input.
        Pulse::Inputs pulseinputs(parameters.t_start, parameters.t_end, parameters.t_step, parameters.numerics_order_highest);
        pulseinputs.add(parameters.pulse_center, parameters.pulse_amp, parameters.pulse_sigma, parameters.pulse_omega, parameters.pulse_type);
        pulse = Pulse( pulseinputs );
        if ( parameters.pulse_amp.at(0) != 0 )
            pulse.fileOutput( parameters.subfolder + "pulse.txt" );
        
        // Time Transformation
        timeTrafoMatrix = ( 1i * operatorMatrices.H_0 ).exp();
        // Check time trafo
        MatrixXcd ttrafo = ( 1i * operatorMatrices.H_0 * 125E-12 ).exp() * operatorMatrices.H_used * ( -1i * operatorMatrices.H_0 * 125E-12 ).exp();
        MatrixXcd temp = dgl_timetrafo( operatorMatrices.H_used, 125E-12 ) - ttrafo;
        if ( std::abs( temp.sum()/parameters.p_omega_atomic ) >= 0.0001 ) {
            logs( "Unitary timetransformation error is bigger than 0.01% of atomic transition energy!\n\n" );
        }
        return true;
    }

    double getTimeborderStart() const {
        return parameters.t_start;
    }

    double getTimeborderEnd() const {
        return parameters.t_end;
    }

    double getTimeStep() const {
        return parameters.t_step;
    }

    MatrixXcd getRho0() const {
        return operatorMatrices.rho;
    }

    MatrixXcd dgl_rungeFunction( const MatrixXcd &rho, const MatrixXcd &H, const double t ) const {
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
        return ret;
    }

    MatrixXcd dgl_timetrafo( const MatrixXcd &A, const double t ) const {
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

    MatrixXcd dgl_chirp( const double t ) const {
        if ( parameters.chirp_total == 0 )
            return MatrixXcd::Zero( parameters.maxStates, parameters.maxStates );
        return 0.5 * operatorMatrices.atom_inversion * chirp.get( t );
    }
    MatrixXcd dgl_pulse( const double t ) const {
        if ( parameters.pulse_amp.at(0) == 0 )
            return MatrixXcd::Zero( parameters.maxStates, parameters.maxStates );
        return 0.5 * ( operatorMatrices.atom_sigmaplus * pulse.get( t ) + operatorMatrices.atom_sigmaminus * std::conj( pulse.get( t ) ) );
    }

    MatrixXcd dgl_phonons( const double t ) const {
        return MatrixXcd::Zero( 1, 1 );
    }

    void expectationValues( const MatrixXcd &rho, const double t ) const {
        std::fprintf( fileoutput.fp_atomicinversion, "%.15e\t%.15e\n", t, std::real( dgl_expectationvalue( rho, operatorMatrices.atom_inversion, t ) ) );
        std::fprintf( fileoutput.fp_photonpopulation, "%.15e\t%.15e\n", t, std::real( dgl_expectationvalue( rho, operatorMatrices.photon_n, t ) ) );
        std::fprintf( fileoutput.fp_densitymatrix, "%e\t", t );
        for ( int j = 0; j < parameters.maxStates; j++ )
            std::fprintf( fileoutput.fp_densitymatrix, "%e\t", std::real( rho( j, j ) ) );
        std::fprintf( fileoutput.fp_densitymatrix, "\n" );
    }

    MatrixXcd dgl_getHamilton( double t ) const {
        return dgl_timetrafo( operatorMatrices.H_used + dgl_chirp( t ) + dgl_pulse( t ), t );
    }

    bool exit_system( const int failure = 0 ) {
        fileoutput.close();
        return true;
    }
};