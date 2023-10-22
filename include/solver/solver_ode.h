#pragma once

#include "global.h"
#include "system/system.h"
#include "solver/solver_tensor_map.h"
#include "types/cache_matrix.h"

namespace QDACC::Numerics {

// Description: ODESolver class provides both Runge-Kutta functions of different orders and functions for different numerical operations
class ODESolver {
   private:
    // Hamilton Cache Tracker
    // TODO: Move to system class OR move phonon tracking variables to solver class.
    int track_gethamilton_read, track_gethamilton_write, track_gethamilton_calc, track_gethamilton_calcattempt;

    // Matrix Cache Map
    std::map<std::string, MultidimensionalCacheMatrix> cache;

    template <typename T>
    using t_output_type = std::map<std::string, std::map<std::string, std::vector<T>>>;
    // Vector Cache Map which are output to file
    t_output_type<Scalar> to_output;
    // Matrix Cache Map which are output to file
    t_output_type<Dense> to_output_m;

    // Cached Entries
    // Vector for saved matrix-time tuples for densitymatrix
    std::vector<QDACC::SaveState> savedStates;
    // Vector for saved matrix-time tuples for hamilton operators
    std::map<double, Sparse> savedHamiltons;
    // Propagators for the path integral. Used for their corresponding correlation functions.
    std::map<double, std::vector<std::vector<Sparse>>> pathint_propagator;
    // Maps t_t onto i for accessing the savedSate Vector via doubles.
    std::map<double, size_t> rho_index_map;

    // Saves the RK45 errors for accepted timesteps
    std::vector<std::tuple<double, double>> rk_error_accepted;
    // Saves the RK45 errors for attempted and accepted timesteps
    std::vector<std::tuple<double, double, double, int>> rk_error;

    // Path Integral Helper Variables
    std::vector<int> pathint_tensor_dimensions;

    // Detector matrix
    std::map<std::string, MultidimensionalCacheMatrix> detector_temporal_mask;
    using t_detector_frequency_mask = std::vector<std::tuple<double, double, double>>;
    std::map<std::string, t_detector_frequency_mask> detector_frequency_mask;
    // Dense detector_frequency_mask_cache;

    /**
     * @brief Saves a tuple of a complex (density-)matrix and time, ensuring times and matrices don't get mixed up
     *
     */
    void saveState( const Sparse &mat, const double t, std::vector<QDACC::SaveState> &savedStates );

    /**
     * @brief Saves a tuple of a complex (Hamilton-)matrix and time, ensuring times and matrices don't get mixed up
     *
     */
    void save_hamilton( const Sparse &mat, const double t );

    /**
     * @brief Iterates Runge-Kutta of order 4 at time t onto rho using the systems hamilton operator.
     *
     */
    Sparse iterateRungeKutta4( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDACC::SaveState> &savedStates );

    /**
     * @brief Iterates Runge-Kutta of order 5 at time t onto rho using the systems hamilton operator.
     *
     */
    Sparse iterateRungeKutta5( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDACC::SaveState> &savedStates );

    /**
     * @brief Iterates Runge-Kutta of order 5 at time t onto rho using the systems hamilton operator. Also returns the error.
     *
     */
    std::pair<Sparse, double> iterateRungeKutta45( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDACC::SaveState> &savedStates );

    /**
     * @brief Get the Time at index i
     *
     */
    double &get_time_at( int i );

    /**
     * @brief Get the Density Matrix at index i
     *
     */
    Sparse &get_rho_at( int i );

    /**
     * @brief Gatheres the Hamiltonoperator for time t using a systems get_hamilton() function. This workaround enables the saving of already calculated matrices for dublicate uses.
     *
     */
    Sparse getHamilton( System &s, const double t );

    /**
     * @brief Calculates the G1(tau) function. Calculates <b^+(t) * b(t+tau)> via quantum regression theorem. Logs and outputs progress.
     *
     */
    void calculate_g1( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &purpose = "unknown" );
    void calculate_g1( System &s, const std::vector<std::string> &s_op_i, const std::string &s_op_j, const std::vector<std::string> &purposes );

    /**
     * @brief Calculates the G2 Correlation Function
     * The G2 function will be evaluated in time only once, but the expectation value loop can be evaulated multiple times.
     *
     */
    void calculate_g2( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &purpose = "unknown" );
    void calculate_g2( System &s, const std::string &s_op_i, const std::vector<std::string> &s_op_j, const std::vector<std::string> &s_op_k, const std::string &s_op_l, const std::vector<std::string> &purposes );

    void calculate_g3( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &s_op_m, const std::string &s_op_n, const std::string &purpose = "unknown" );
    void calculate_g3( System &s, const std::string &s_op_i, const std::vector<std::string> &s_op_j, const std::vector<std::string> &s_op_k, const std::vector<std::string> &s_op_l, const std::vector<std::string> &s_op_m, const std::string &s_op_n, const std::vector<std::string> &purposes );

   public:
    ODESolver(){};
    ODESolver( System &s );

    /**
     * @brief Iterates Runge-Kutta with given order depending on the systems settings.
     *
     */
    Sparse iterate( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDACC::SaveState> &savedStates );

    /**
     * @brief Calculates the normal t-direction via solving the von-Neumann equation for rho. May save some of the density matrices for later uses. Logs the calculation and outputs progress.
     *
     */
    bool calculate_t_direction( System &s );

    /**
     * @brief Takes a superposition of cavity modes or annihilator-type operator strings (e.g. 'h+GH+HZ') and returns a tuple of [creator, annihilator] strings
     *
     */
    std::tuple<std::string, std::string> get_operator_strings( System &s, const std::string &operators );

    /**
     * @brief Returns an intentifier string for a combination of G-function order and a superposition of operators
     *
     */
    std::string get_operators_purpose( const std::vector<std::string> &operators );

    /**
     * @brief Returns a tuple of [creator, annihilator] Sparse Operator Matrices from a given Operator String
     *
     */
    Sparse get_operators_matrix( System &s, const std::string &s_op );
    std::tuple<Sparse, Sparse> get_operators_matrices( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator );

    /**
     * @brief Calculates the Spectrum for a given operator combination. This function requires the calculation of 1 G1 Function
     *
     */
    bool calculate_spectrum( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator, double frequency_center, double frequency_range, int resolution, int order, bool normalize, std::string s_g = "" );

    /**
     * @brief Calculates the Indistinguishability and Visibility for a given operator combination. This function requires the calculation of 1 G1 and 1 G2 Function.
     *
     */
    bool calculate_indistinguishability( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator );

    /**
     * @brief Calculates the Concurrence for a given operator combination. This function requires the calculation of 6 G2 Functions and will produce 4 values for the concurrence.
     *
     */
    bool calculate_concurrence( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2, const int matrix_priority_evaluation, const double spec_center, const double spec_range, const double spec_res );

    /**
     * @brief Calculates the Wigner Function for a given operator combination
     *
     */
    bool calculate_wigner( System &s, const std::string &s_mode, const double x, const double y, const int resolution, const int skip );

    /**
     * @brief Chains all calculations requiring G1 and G2 functions as well as Wigner functions and Raman populations.
     *
     */
    bool calculate_advanced_photon_statistics( System &s );

    /**
     * @brief Calculates a single Path Integral Propagator
     *
     */
    Sparse calculate_propagator_single( System &s, size_t tensor_dim, double t0, double t_step, int i, int j, std::vector<QDACC::SaveState> &output, const Sparse &one );

    /**
     * @brief Calculates the Path Integral Propagator for all index combinations. This function also saves the resuling propagators and returns a reference to the cached object.
     *
     */
    std::vector<std::vector<Sparse>> &calculate_propagator_vector( System &s, size_t tensor_dim, double t0, double t_step, std::vector<QDACC::SaveState> &output );

    /**
     * @brief Iterates rho0 from t_start to t_end using RK4, RK5 or RK45. If the order is 45, this function calls calculate_runge_kutta_45.
     *
     */
    bool calculate_runge_kutta( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDACC::SaveState> &output, bool do_output = true );

    /**
     * @brief Iterates rho0 from t_start to t_end using RK45. This function does NOT interpolate the result but instead returns all calculated densitymatrices as vector.
     *
     */
    bool calculate_runge_kutta_45( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDACC::SaveState> &output, bool do_output = true );

    /**
     * @brief Calculates rho0 from t_start to t_end using the Path Integral Method. If correlation functions are evaluated, this function calls calculate_path_integral_correlation multiple times per iteration.
     *
     */
    bool calculate_path_integral( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDACC::SaveState> &output, bool do_output = true );

    /**
     * @brief Iterates rho0 from t_start to t_end using the Path Integral Method.
     */
    Tensor iterate_path_integral( System &s, Tensor &adm_tensor, std::vector<std::vector<Sparse>> &propagator, const int max_index );
    Tensor iterate_path_integral_gpu( System &s, Tensor &adm_tensor, std::vector<std::vector<Sparse>> &propagator, const int max_index );

    /**
     * @brief Calculates rho0 from t_start to t_end incorporating the corresponding correlation modification.
     *
     */
    bool calculate_path_integral_correlation( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, System &s, std::vector<QDACC::SaveState> &output, const Sparse &op_l, const Sparse &op_i, int adm_multithreading_cores );

    /**
     * @brief Creates a table file containing the available propagation paths in the .dot format
     *
     */
    bool visualize_path( Sparse &rho0, System &s );

    /**
     * @brief Integrates the Raman photon population. Very runtime costly. TODO: Multithread/Optimize integral.
     * FIXME: Currently broken lol.
     */
    bool calculate_raman_population( System &s, const std::string &electronic_transition1, const std::string &electronic_transition2, const std::string &optical_transition, const std::string &pulse_mode );

    /**
     * @brief Outputs numerical data like the RK45 errors.
     *
     */
    bool output_numerical_data( System &s );

    /**
     * @brief Applies the detector mapping function to a G1 or G2 correlation function.
     *
     */
    void apply_detector_function( System &s, MultidimensionalCacheMatrix &mat, const std::string &mat_mode = "G" );

    /**
     * @brief Initializes the detector functions for G1 and G2 correlation functions.
     * This function requireds a valid MultidimensionalCacheMatrix object for the dimensions, so it can
     * only be called after the MultidimensionalCacheMatrix has been initialized.
     */
    void initialize_detector_functions( System &s, MultidimensionalCacheMatrix &mat );

    /**
     * @brief Calculates the Eigenvalues of the Hamilton operators H_0, H_I and H_0+H_I and outputs them to "hamilton_eigenvalues.txt"
     *
     */
    void calculate_hamilton_eigenvalues( System &s, const int power = 1 );

    /**
    * @brief Adds value to the output vector for name using key as the individual dataset name.
    * If to_output[name] already contains key, then key is appended with a number to avoid overwriting.
    * @param name: [std::string] Name of the dataset
    * @param key: [std::string] Name of the individual dataset
    * @param value: [T] Value to be added to the dataset
    * @param where: [t_output_type<T>] Reference to the output vector
    * @param overwrite: [bool] If true, overwrites the dataset if it already exists anyways
    */
    template <typename T>
    std::vector<T> &add_to_output( const std::string &name, std::string key, const std::vector<T> &value, t_output_type<T> &where, bool overwrite = false ) {
        int count;
        if ( not overwrite and where.contains( name ) and ( count = where[name].count( key ) ) > 0 ) {
            key += "_" + std::to_string( count );
        }
        where[name][key] = value;
        return where[name][key];
    }
};

} // namespace QDACC::Numerics