#pragma once

#include "global.h"
#include "system/parameters.h"
#include "system/fileoutput.h"
#include "system/operatormatrices.h"
#include "system/evaluable/chirp.h"
#include "system/evaluable/pulse.h"
#include "solver/solver.h"

namespace QDACC {

class System {
   public:
    // Vector of input arguments. The argv elements will get pushed into this vector
    std::vector<std::string> arguments;

    // Shift in Energy for any state or cavity mode
    std::vector<Chirp> chirp;
    // Pulse for any transition or cavity mode
    std::vector<Pulse> pulse;
    // FileOutput Object, handles all fileoutputs like opening/closing/header etc.
    FileOutput fileoutput;
    // Parameters object, holds all system and numerical parameters
    Parameters parameters;
    // OperatorMatrices object, holds all state- and transition matrices as well as some helpervariables
    OperatorMatrices operatorMatrices;

    // Runtime efficient caching vector
    // Map of saved phonon-phi function
    std::map<double, Scalar> phi_vector;
    // Vector of saved phonon-phi function
    std::vector<Scalar> phi_vector_int;
    // Vector of saved coefficients for e.g. phonon terms.
    std::map<double, std::map<double, QDACC::SaveStateTau>> savedCoefficients;
    // ##### Helper Variables #####
    std::map<std::string, double> emission_probabilities;

    // Coefficient Tracking variables
    int track_getcoefficient_read = 0;              // Read coefficient from vector
    int track_getcoefficient_read_interpolated = 0; // Read interpolated coefficient from vector
    int track_getcoefficient_write = 0;             // Wrote coefficient to vector
    int track_getcoefficient_calculate = 0;         // Sucessfully calculated coefficient
    int track_getcoefficient_calcattempt = 0;       // Attempt to calculate coefficient
    int globaltries = 0;

    // Time Trafo Matrix Caching
    Sparse timeTrafoMatrix;

    /**
     * @brief Constructs System
     *
     */
    System();
    System( const std::vector<std::string> &input );

    /**
     * @brief Transforms a Matrix into the interaction picture if enabled.
     *
     * @return Sparse input matrix in the interaction picture
     */
    Sparse dgl_timetrafo( Sparse ret, const double t );

    /**
     * @brief Calculates the differential equation for different input times for evaluation with the any-order Runge-Kutta solver
     *
     * @return Sparse Matrix Runge Function
     */
    Sparse dgl_runge_function( const Sparse &rho, const Sparse &H, const double t, std::vector<QDACC::SaveState> &past_rhos );

    /**
     * @brief Initializes all system parameters, cache variables and precalculates all functions that allow for caching
     *
     * @return true If Initialization was succesfull
     * @return false Otherwise
     */
    bool init_system();

    /**
     * @brief Finalizes all system parameters, empties cached variables and saves all outputs.
     *
     * @return true If Finalization was succesfull
     * @return false Otherwise
     */
    bool exit_system( const int failure = 0 );

    /**
     * @brief Calculates the chirped Hamilton operator
     *
     * @return Sparse Matrix Chirp Contribution to add onto the Hamilton Operator
     */
    Sparse dgl_chirp( const double t );

    /**
     * @brief Calculates the pulse Hamilton operator
     *
     * @return Sparse Matrix Pulse Contribution to add onto the Hamilton Operator
     */
    Sparse dgl_pulse( const double t );

    /**
     * @brief Calculates and outputs expectation values for all available observables
     *
     */
    void calculate_expectation_values( const std::vector<QDACC::SaveState> &rhos, Timer &evalTimer );

    /**
     * @brief Calculates or returns the cached(if allowed) Hamiltonian for current time t.
     *
     * @return Sparse Matrix Hamilton Operator
     */
    Sparse dgl_get_hamilton( const double t );

    /**
     * @brief Validates trace is still contained, if not, outputs trace in file and returns false. Can also be forced by setting force = true.
     *
     * @return true Trace is valid
     * @return false Trace is not valid
     */
    bool trace_valid( Sparse &rho, double t_hit, bool force = false );

    // ##### Phonon Functions and Helperfunctions #####

    /**
     * @brief The Polaron Transformed Operators are calculated by evaluating the Time Transformed X_g and X_u. This Runge Function is used to calculate this time transformation using the von-Neumann euqation instead of
     * numerically evaluating the tim transformation by e.g. using the numerical matrix exponential. The noteworthy difference here is to also include the explicit time dependency of the transformed operators.
     *
     * @param chi Current Chi(t)
     * @param t Current Time
     * @return Sparse: Applied Runge function
     */
    Sparse dgl_phonons_rungefunc( const Sparse &chi, const double t );

    /**
     * @brief Calculates the Polaron Green Function. This function is currently not cached.
     *
     * @param t Current Time
     * @param mode Either "g" or "u"
     * @return Scalar: G_u/g(t)
     */
    Scalar dgl_phonons_greenf( double tau, const char mode = 'u' );

    /**
     * @brief Same but with phonon coupling scaling, returning a Sparse Matrix to do cwiseMultiplication with. This function will cache the matrices to avoid expensive recalculations.
     *
     * @return Sparse& Matrix Reference with Green Function per element.
     */
    Sparse &dgl_phonons_greenf_matrix( double tau, const char mode = 'u' );

    /**
     * @brief Calculates the Phonon Correlation Kerbal Phi(tau) = int_0^inf dw J(w) / w^2 * [coth(hbar w/(2 k_b T)) * cos(w t) - i*sin(w t)]
     *
     * @return Scalar value Phi(t)
     */
    Scalar dgl_phonons_phi( const double tau );

    /**
     * @brief Calculates the Phonon Spectral Density
     *
     * @param w Input frequency
     * @return J(w)
     */
    double dgl_phonons_spectral_density( const double w );

    /**
     * @brief Calculates to Lindblad Coefficients for the analytical solution of the Polaron Contributions.
     *
     * @param energy Transition Energy
     * @param coupling Cavity Coupling
     * @param pulse Current Pulse values
     * @param mode Mode
     * @param sign +1 or -1
     * @return double: L(t)
     */
    double dgl_phonons_lindblad_coefficients( const double energy, const double coupling, const Scalar pulse, const char mode = 'L', const double scaling = 1.0, const double sign = 1.0 );

    Sparse dgl_phonons_lindblad_contribution( const double t, const Sparse &rho );
    Sparse dgl_phonons_integrated_contribution( const double t, const Sparse &rho, const std::vector<QDACC::SaveState> &past_rhos );

    /**
     * @brief Initializes the Polaron Frame Functions by precalculating the Phi(tau) function and the corresponding Green functions
     *
     */
    void initialize_polaron_frame_functions();

    /**
     * @brief Calculates Chi(t,0)
     *
     * @return Sparse Matrix
     */
    Sparse dgl_phonons_chi( const double t );

    /**
     * @brief Evaluates the Transformation U(t,tau) Chi(t) U(t,tau)^+ = Chi(t,tau) from tau' = 0 to tau' = tau
     *
     * @param chi_tau Current Chi(t)
     * @param t Current Time
     * @param tau Current Delay
     * @return Sparse: Transformed Chi(t) = Chi(t,tau)
     */

    Sparse dgl_phonons_calculate_transformation( double t, double tau );

    /**
     * @brief Calculates the Polaron Master Equation Contribution L_phonons(t)
     *
     * @return Sparse Matrix to add onto the usual von-Neumann equatation / the runge function
     */
    Sparse dgl_phonons_pmeq( const Sparse &rho, const double t, const std::vector<QDACC::SaveState> &past_rhos );

    /**
     * @brief Initializes the Path Integral Functions by precalculating the Kernel functions
     *
     */
    void initialize_path_integral_functions();

    /**
     * @brief Calculates the Phonon Kernel Functions for the Path Integral Method
     *
     * @return Scalar valued K(t)
     */
    Scalar dgl_phonons_kernel( const double t, const double t_step );

    /**
     * @brief Calculates the phonon Correlation function for the Path Integral Method
     *
     * @return Scalar valued S(i,j,id,jd)
     */
    Scalar dgl_phonon_memory_function( const int t_delta, const int i_n, const int j_n, const int i_nd, const int j_nd );

    /**
     * @brief Calculates the expectation values for a given operator
     *
     * @return Expectation Value <op>
     */
    template <typename T, typename R>
    inline R dgl_expectationvalue( const T &rho, const T &op, const double t ) {
        return get_trace<R>( ( rho * dgl_timetrafo( op, t ) ).eval() );
    }

    /**
     * @brief Calculates the transformed two-time density matrix for the first order correlation function
     *
     * @return T op*rho
     */
    template <typename T>
    inline T dgl_calc_rhotau( const T &rho, const T &op, const double t ) {
        return op * rho; // dgl_timetrafo( op, t ) * rho;
    }

    /**
     * @brief Calculates the transformed two-time density matrix for the second order correlation function.
     *
     * @return op1*rho*op2
     */
    template <typename T>
    inline T dgl_calc_rhotau_2( const T &rho, const T &op1, const T &op2, const double t ) {
        // return op1 * rho * op2; // dgl_timetrafo( op1, t ) * rho * dgl_timetrafo( op2, t );
        return dgl_timetrafo( op1, t ) * rho * dgl_timetrafo( op2, t );
    }

    /**
     * @brief Wrapper to calculate the matrix trace from Dense input types
     *
     * @return Trace(mat)
     */
    template <typename T>
    inline T get_trace( const Dense &mat ) const {
        return mat.trace();
    }

    /**
     * @brief Wrapper to calculate the matrix trace from Sparse input types not supporting .trace()
     *
     * @return Trace(mat)
     */
    template <typename T>
    inline T get_trace( const Sparse &mat ) const {
        // return get_trace<T>( Dense( mat ) );
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

    /**
     * @brief Function to calculate the partial trace for a specific base index. The base index is an internal value determined automatically on generation of the matrices.
     * The base index from a base state vector |i,j,k,l,...> can be one of the i,j,k,l indices to partially trace over.
     * @return Dense
     */
    Dense partial_trace( const Sparse &mat, int i ) {
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

    /**
     * @brief Alias for the commutator return type, which is specified explicitely using Eigen's reference objects.
     *
     * @tparam T1
     * @tparam T2
     */
    template <typename T1, typename T2>
    using commutator_return_type = Eigen::CwiseBinaryOp<Eigen::internal::scalar_difference_op<typename T1::Scalar, typename T1::Scalar>, const Eigen::Product<T1, T2, 2>, const Eigen::Product<T2, T1, 2>>;
    template <typename T1, typename T2>
    using anticommutator_return_type = Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<typename T1::Scalar, typename T1::Scalar>, const Eigen::Product<T1, T2, 2>, const Eigen::Product<T2, T1, 2>>;

    /**
     * @brief Calculates the Kommutator [A,B] of two matrices A and B of identical type
     *
     * @return Eigen Reference for A*B-B*A
     */
    template <typename T1, typename T2>
    inline commutator_return_type<T1, T2> dgl_kommutator( const T1 &A, const T2 &B ) {
        return A * B - B * A;
    }

    /**
     * @brief Calculates the Anti-Kommutator [A,B]_- of two matrices A and B of identical type
     *
     * @return Eigen Reference for A*B+B*A
     */
    template <typename T1, typename T2>
    inline anticommutator_return_type<T1, T2> dgl_antikommutator( const T1 &A, const T2 &B ) {
        return A * B + B * A;
    }

    /**
     * @brief Calculates the Lindbladian of two input matrices where L = 2*op*rho*opd - opd*op*rho - rho*opd*op
     *
     * @return T
     */
    template <typename T, typename T2, typename T3>
    inline T dgl_lindblad( const T &rho, const T2 &op, const T3 &opd ) {
        return 2.0 * op * rho * opd - opd * op * rho - rho * opd * op;
    }
};

} // namespace QDACC