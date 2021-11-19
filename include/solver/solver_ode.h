#pragma once

#include "global.h"
#include "system/system.h"
#include "misc/interpolant.h"
#include "solver/solver_tensor_map.h"

namespace QDLC {

namespace Numerics {

// Description: ODESolver class provides both Runge-Kutta functions of different orders and functions for different numerical operations
class ODESolver {
   private:
    int track_gethamilton_read, track_gethamilton_write, track_gethamilton_calc, track_gethamilton_calcattempt;
    int dim; //TODO: //Remove
    // Cache Matrices that allow for variable time steps.
    std::map<std::string, Dense> cache;
    std::map<std::string, std::map<std::string, std::vector<Scalar>>> to_output;
    std::map<std::string, std::map<std::string, std::vector<Dense>>> to_output_m;

    // Cached Entries
    std::vector<QDLC::SaveState> savedStates;                              // Vector for saved matrix-time tuples for densitymatrix
    std::map<double, Sparse> savedHamiltons;                               // Vector for saved matrix-time tuples for hamilton operators
    std::map<double, std::vector<std::vector<Sparse>>> pathint_propagator; // Propagators for the path integral. Used for their corresponding correlation functions.
    std::map<double, size_t> rho_index_map;                                    // Maps t_t onto i for accessing the savedSate Vector via doubles.

    // Path Integral Helper Variables
    std::vector<int> pathint_tensor_dimensions;

    // Description: Saves a tuple of a complex (density-)matrix and time, ensuring times and matrices don't get mixed up
    // Type: ODESolver private function
    // @param mat: [&Sparse] Matrix to save
    // @param t: [double] Corresponding time
    // @return: [void]
    void saveState( const Sparse &mat, const double t, std::vector<QDLC::SaveState> &savedStates );

    // Description: Saves a tuple of a complex (Hamilton-)matrix and time, ensuring times and matrices don't get mixed up
    // Type: ODESolver private function
    // @param mat: [&Sparse] Matrix to save
    // @param t: [double] Corresponding time
    // @return [void]
    void saveHamilton( const Sparse &mat, const double t );

    // Description: Iterates Runge-Kutta of order 4 at time t onto rho using the systems hamilton operator.
    // Type: ODESolver private function
    // @param rho: [&Sparse] Input (density-) matrix
    // @param s: [&System] Class providing set of system functions
    // @param t: [double] Time to iterate at
    // @return: [Sparse] rho at time t+t_step

    Sparse iterateRungeKutta4( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDLC::SaveState> &savedStates );

    // Description: Iterates Runge-Kutta of order 5 at time t onto rho using the systems hamilton operator.
    // Type: ODESolver private function
    // @param rho: [&Sparse] Input (density-) matrix
    // @param s: [&System] Class providing set of system functions
    // @param t: [double] Time to iterate at
    // @return: [Sparse] rho at time t+t_step
    Sparse iterateRungeKutta5( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDLC::SaveState> &savedStates );

    std::pair<Sparse, double> iterateRungeKutta45( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDLC::SaveState> &savedStates );

    // Desciption: Function to calculate the number of iterations used for tau direction calculations
    // Type: ODESolver private function
    // @param s: [&System] Class providing set of system functions
    // @return [int] Number of tau-direction iterations
    int getIterationNumberTau( System &s );

    int getIterationNumberSpectrum( System &s );
    double &getTimeAt( int i );
    Sparse &getRhoAt( int i );

    // Description: Gatheres a Hamiltonoperator using a systems get_hamilton() function. This workaround enables the saving of already calculated matrices for dublicate uses. Seems to have almost no influence on runtime.
    // Type: ODESolver private function
    // @param s: [&System] Class providing set of system functions
    // @param t: [double] Time to return Hamiltonoperator at
    // @return: [Sparse] hamilton matrix of type Sparse
    Sparse getHamilton( System &s, const double t, bool use_saved_hamiltons = false );

    // Description: Resets all temporary variables. Currently: out, akf_mat, track variables
    // Type: ODESolver private function
    // @param s: [&System] Class providing set of system functions
    // @return: [int] dimensions of temporary variables
    int reset( System &s );

    // @return: [tuple] Used Operator Matrices in order of input parameters
    std::tuple<Sparse, Sparse> calculate_g1( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator, std::string purpose = "unknown" ); //std::vector<std::vector<QDLC::SaveScalar>>
    std::tuple<Sparse, Sparse, Sparse, Sparse> calculate_g2( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2, std::string purpose = "unknown" );

    // Description: Stretches Data evaluated at various timesteps onto an equidistant grid, such that g1 and g2 can work with said grid.
    bool scale_grid( System &s, Dense &cache, std::vector<std::vector<QDLC::SaveScalar>> &cache_noneq );

   public:
    ODESolver(){};
    ODESolver( System &s );

    // Description: Iterates Runge-Kutta with given order depending on the systems settings.
    // Type: ODESolver public function
    // @param rho: [&Sparse] Input (density-) matrix
    // @param s: [&System] Class providing set of system functions
    // @param t: [double] Time to iterate at
    // @param dir: [int] Time direction. Solver order can differ in both t and tau direction, depending on the systems settings. Default value is DIR_T
    // @return: [Sparse] rho at time t+t_step
    Sparse iterate( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDLC::SaveState> &savedStates, const int dir = DIR_T );

    // Description: Calculates the normal t-direction via solving the von-Neumann equation for rho. May save some of the density matrices for later uses. Logs the calculation and outputs progress.
    // Type: ODESolver public function
    // @param s: [&System] Class providing set of system functions
    // @return: [bool] True if calculations are sucessfull, else false
    bool calculate_t_direction( System &s );

    //bool calculate_g2_0( System &s, const Sparse &op_creator, const Sparse &op_annihilator, std::string fileOutputName ); //moved to advancedPhotonStatistics
    std::tuple<std::string, std::string> get_operator_strings( const std::string &operators );
    std::string get_operators_purpose( const std::vector<std::string> &operators, int order );
    std::tuple<Sparse, Sparse> get_operators_matrices( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator );
    bool calculate_spectrum( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator, double frequency_center, double frequency_range, int resolution );
    bool calculate_indistinguishability( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator );
    bool calculate_concurrence( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2 );
    bool calculate_wigner( System &s, const std::string &s_mode, const double x, const double y, const int resolution, const int skip );
    bool calculate_advanced_photon_statistics( System &s );

    Sparse calculate_propagator_single( System &s, size_t tensor_dim, double t0, double t_step, int i, int j, std::vector<QDLC::SaveState> &output, const Sparse &one );
    std::vector<std::vector<Sparse>> &calculate_propagator_vector( System &s, size_t tensor_dim, double t0, double t_step, std::vector<QDLC::SaveState> &output );

    bool calculate_runge_kutta( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output = true );
    bool calculate_runge_kutta_45( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output = true, bool interpolate = true, double tolerance = 1E-4 );
    bool calculate_path_integral( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output = true );
    bool calculate_path_integral_correlation( Tensor adms, Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, size_t total_progressbar_iterations, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output, const std::vector<Sparse> &matrices, const std::vector<std::vector<int>> &adm_multithreaded_indices, int adm_multithreading_cores );
};

} // namespace Numerics

} // namespace QDLC