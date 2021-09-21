#pragma once

#include "global.h"
#include "system/parameters.h"
#include "system/fileoutput.h"
#include "system/operatormatrices.h"
#include "chirp.h"
#include "pulse.h"
#include "solver/solver.h"

class System {
   public:
    std::string name;
    std::string terminate_message;

    std::vector<std::string> arguments;

    // ##### System Components #####
    std::vector<Chirp> chirp;
    //Pulse pulse_H; //REMOVE
    //Pulse pulse_V; //REMOVE
    std::vector<Pulse> pulse;
    FileOutput fileoutput;
    Parameters parameters;
    OperatorMatrices operatorMatrices;

    // ##### Cache Components #####
    // Runtime efficient caching vector
    std::map<double, Scalar> phi_vector;                                 // Vector of saved phonon-phi function
    std::vector<Scalar> phi_vector_int;                                  // Same, but with integers for path integral
    std::map<std::pair<double, double>, SaveStateTau> savedCoefficients; // Vector of saved coefficients for e.g. phonon terms.

    // ##### Helper Variables #####
    std::map<std::string, double> emission_probabilities;
    // Helpervariables for photon emission probability
    double photonemissionprob_integral_H = 0;
    double photonemissionprob_integral_V = 0;
    double electronic_emissionprob_integral_H = 0;
    double electronic_emissionprob_integral_V = 0;
    // Helpervariables for raman emission process
    double ramanphotonemissionprob_integral_H = 0;
    double ramanphotonemissionprob_integral_V = 0;
    Scalar ramanphotonpopulation_integral_H = 0;
    Scalar ramanphotonpopulation_integral_V = 0;

    // Coefficient Tracking variables
    int track_getcoefficient_read = 0;             // Read coefficient from vector
    int track_getcoefficient_write = 0;            // Wrote coefficient to vector
    int track_getcoefficient_calculate = 0;        // Sucessfully calculated coefficient
    int track_getcoefficient_read_but_unequal = 0; // Tried to read coefficient, but didnt find any. Recalculate
    int track_getcoefficient_calcattempt = 0;      // Attempt to calculate coefficient
    int globaltries = 0;

    // Time Trafo Matrix Caching
    Sparse timeTrafoMatrix;

    // ##### Content #####
    // Constructor
    System();
    System( const std::vector<std::string> &input );

    // Initializes the System class. Will be called inside the system constructor
    Sparse dgl_timetrafo( Sparse ret, const double t );

    // Calculates the differential equation for different input times for evaluation with the any-order Runge-Kutta solver
    Sparse dgl_rungeFunction( const Sparse &rho, const Sparse &H, const double t, std::vector<SaveState> &past_rhos );

    // Initializes all system parameters, cache variables and precalculates all functions that allow for caching
    bool init_system();

    // Finalizes all system parameters, empties cached variables and saves all outputs.
    bool exit_system( const int failure = 0 );

    // Helperfunction for the solver class. This method allows the solver to change system parameters
    bool command( unsigned int index = 0 );

    // Calculates the chirped Hamilton operator
    Sparse dgl_chirp( const double t );

    // Calculates the pulse Hamilton operator
    Sparse dgl_pulse( const double t );

    // Integrates the Raman photon population. Very runtime costly. TODO: Multithread/Optimize integral.
    Scalar dgl_raman_population_increment( const std::vector<SaveState> &past_rhos, const char mode, const Scalar before, const double t );

    // Calculates and outputs expectation values for all available observables
    void expectationValues( const std::vector<SaveState> &rhos, Timer &evalTimer );

    // Calculates or returns the cached(if allowed) Hamiltonian for current time t.
    // This function is important because it allows for e.g. path integral to use different definitions of H
    // where no phonon <B> is incorporated.
    Sparse dgl_getHamilton( const double t );

    // Validates trace is still contained, if not, outputs trace in file and returns false.
    // Can also be forced by setting force = true
    bool traceValid( Sparse &rho, double t_hit, bool force = false );

    // ##### Phonon Functions and Helperfunctions #####

    // Calculates the differential for the phonon matrix. Uses the dgl_getHamilton function.
    Sparse dgl_phonons_rungefunc( const Sparse &chi, const double t );

    // Calculates the phonon X Matrix from a given Chi
    Sparse dgl_phonons_chiToX( const Sparse &chi, const char mode = 'u' );

    // Either calculates or returns the polaron green function from cache
    Scalar dgl_phonons_greenf( double t, const char mode = 'u' );

    // Claculates the Phi function for the polaron green function. This function will be cached
    Scalar dgl_phonons_phi( const double t );

    // Calculates the Lindbladian coefficients for the analytical phonon contributions using the polaron frame
    double dgl_phonons_lindblad_coefficients( double t, double energy, double coupling, Scalar pulse, const char mode = 'L', const double sign = 1.0 );

    // Initializes the polaron functions
    void initialize_polaron_frame_functions();

    // Calculates Chi(t,0)
    Sparse dgl_phonons_chi( const double t );

    // Calculates Chi(t,tau > 0) after Chi(t,0) was calculated with dgl_phonons_chi.
    // Uses different approaches
    Sparse dgl_phonons_calculate_transformation( Sparse &chi_tau, double t, double tau );

    // Calculates L_phonons(t)
    Sparse dgl_phonons_pmeq( const Sparse &rho, const double t, const std::vector<SaveState> &past_rhos );

    // Path integral phonon functions
    void initialize_path_integral_functions();
    Scalar dgl_phonons_kernel( const double t, const double t_step );
    Scalar dgl_phonon_S_function( const int t_delta, const int i_n, const int j_n, const int i_nd, const int j_nd );
    // ##### Template Functions #####

    // Calculates the expectation values for a given operator
    template <typename T, typename R>
    inline R dgl_expectationvalue( const T &rho, const T &op, const double t ) {
        return getTrace<R>( ( rho * dgl_timetrafo( op, t ) ).eval() );
    }

    // Calculates the transformed density matrix for the first order correlation function
    template <typename T>
    inline T dgl_calc_rhotau( const T &rho, const T &op, const double t ) {
        return dgl_timetrafo( op, t ) * rho;
    }

    // Calculates the transformed density matrix for the second order correlation function
    template <typename T>
    inline T dgl_calc_rhotau_2( const T &rho, const T &op1, const T &op2, const double t ) {
        return dgl_timetrafo( op1, t ) * rho * dgl_timetrafo( op2, t );
    }

    // Functions to determine numerical matrix trace from input types supporting .trace()
    // @param &mat: Input dense matrix
    // @return Matrix trace
    template <typename T>
    inline T getTrace( const Dense &mat ) const {
        return mat.trace();
    }
    // Functions to determine numerical matrix trace from sparse input types not supporting .trace()
    // @param &mat: Input sparse matrix
    // @return Matrix trace of type T
    template <typename T>
    inline T getTrace( const Sparse &mat ) const {
        //return getTrace<T>( Dense( mat ) );
        T ret = (T)0.0;
        for ( int k = 0; k < mat.outerSize(); ++k ) {
            for ( Sparse::InnerIterator it( mat, k ); it; ++it ) {
                if ( it.row() == it.col() ) {
                    ret += it.value();
                }
            }
        }
        return ret;
    }

    // Function to calculate the partial trace for a specific base index
    Dense partialTrace( const Sparse &mat, int i ) {
        Dense ret = Dense::Zero( operatorMatrices.base_selfhilbert.at( i ).rows(), operatorMatrices.base_selfhilbert.at( i ).cols() );
        for ( int k = 0; k < mat.outerSize(); ++k ) {
            for ( Sparse::InnerIterator it( mat, k ); it; ++it ) {
                int l = it.row();
                int j = it.col();
                int hi = std::real( operatorMatrices.base_hilbert_index[i]( l, j ) ) - 1;
                int hj = std::imag( operatorMatrices.base_hilbert_index[i]( l, j ) ) - 1;
                if ( hi >= 0 and hj >= 0 )
                    ret( hi, hj ) += it.value();
            }
        }
        return ret;
    }

    // Calculates the Kommutator of two matrices of identical type
    // @param &A,&B: Input matrices
    // @return Kommutator of A and B where [A,B] = AB-BA
    template <typename T>
    inline T dgl_kommutator( const T &A, const T &B ) {
        return A * B - B * A;
    }

    // Calculates the Anti-Kommutator of two matrices of identical type
    // @param &A,&B: Input matrices
    // @return Anti-Kommutator of A and B where [A,B]^- = BA-AB
    template <typename T>
    inline T dgl_antikommutator( const T &A, const T &B ) {
        return A * B + B * A;
    }

    // Calculates the Lindbladian of two input matrices where L = 2*A*rho*B - B*A*rho - rho*B*A
    // @param &rho: Input density matrix
    // @param &op: Input matrix O
    // @param &opd: Input matrix O^dagger
    // @return Returns the Lindbladian of rho, O and O^dagger
    template <typename T, typename T2, typename T3>
    inline T dgl_lindblad( const T &rho, const T2 &op, const T3 &opd ) {
        return 2.0 * op * rho * opd - opd * op * rho - rho * opd * op;
    }
};