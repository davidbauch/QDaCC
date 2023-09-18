#pragma once
// Dependencies
#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "misc/timer.h"
#include "system/parameter.h"

#define GLOBAL_PROGRAM_VERSION "4.0.0"
#define GLOBAL_PROGRAM_LASTCHANGE ""

namespace QDLC {

class Parameters {
   public:
    // Numerical Parameters
    Parameter hbar, kb;

    // Working Directory
    std::string working_directory;

    // Numerical Approximations
    bool numerics_use_interactionpicture, numerics_use_rwa;

    // If True, the program will extend itself until the system converges into the defined groundstate or t > numerics_hard_t_max.
    bool numerics_calculate_till_converged;
    double numerics_hard_t_max;

    // Groundstate Parameters
    int numerics_groundstate;
    std::string numerics_groundstate_string;

    // Order of the Timetransformation. Can be Analytical or Numerical
    QDLC::TransformationOrder numerics_order_timetrafo;

    // Advanced Output Settings for Debugging
    int output_advanced_log, output_handlerstrings, output_coefficients;

    // Number of Gridpoints for Tau-Calculations
    int grid_resolution;
    // Grid Helpervariables to allow for the Precalculation of the entire Grid
    // Calculate T and Correlation Grid for these timesteps. maps Index -> t_t, t_step
    std::vector<double> grid_steps, grid_values;
    // maps t_t -> Index
    std::map<double, size_t> grid_value_indices;
    // Number of Iterations to skip for the grid. E.g. :....:....:....: (every fith value is used for grid)
    int iterations_t_skip;

    // Assumed maximum iterations for the main t-direction calculation
    int iterations_t_max;

    // Vector to save trace into.
    std::vector<double> trace;

    // Output Switches
    // Interpolate Outputs
    bool numerics_interpolate_outputs;
    std::set<std::string> output_dict;

    // Maximum Threads for the Program to use in general
    int numerics_maximum_primary_threads;
    // Maximum Threads for the Phonon- and other Subroutines to use
    int numerics_maximum_secondary_threads;

    // Phonon Parameters
    // Enables the first Markov Approximation
    bool numerics_phonon_approximation_markov1;
    // Defines the Approximation Method used for the Polaron Timetransformation
    QDLC::PhononApproximation numerics_phonon_approximation_order;
    // Quantum Dot / Phonon Parameters
    Parameter p_phonon_b, p_phonon_ohm, p_phonon_alpha, p_phonon_wcutoff, p_phonon_T, p_phonon_tcutoff, p_phonon_wcutoffdelta, p_phonon_pure_dephasing;
    Parameter p_phonon_qd_de, p_phonon_qd_dh, p_phonon_qd_rho, p_phonon_qd_cs, p_phonon_qd_ratio, p_phonon_qd_ae;
    // Additional Scalings for phonon adjustment.
    double p_phonon_adjust_rad, p_phonon_adjust_dep, p_phonon_adjust_b;
    // When false, RK45 for the Polaron Backwards integral is disabled. Default should be false.
    // TODO: currently, enabling this would lead to wrong results because the PME subroutine expects equidistant matrices. Need to interpolate then.
    bool numerics_phonon_nork45;

    // Path integral Parameters
    // Tensor Rank
    int p_phonon_nc;
    // Enables the use of the partially summed PI algorithm, greatly increasing its efficiency. This should always be true and only set false for debugging purposes.
    bool numerics_pathint_partially_summed;
    // Squared threshold for elements to be considered zero
    double numerics_pathintegral_squared_threshold;
    // Threshold for the .prune() method
    double numerics_pathintegral_sparse_prune_threshold;
    // If new_nonzeros/nonzeros is below this value, the tensor is assumed to be fixed in size, increasing the densitychange_counter.
    double numerics_pathintegral_sparse_to_dense_threshold;
    // Propagator Cutoff; When iterated multiple times, M(t0->t1) may gain additional entries in between the RK iterations from t0 to t1. When set to true, the final non-zero matrix entries will be mapped onto the non-zero entries after
    // the first iteration, meaning any additional non-zero entries besides the ones created within the first iteration are lost.
    bool numerics_pathintegral_docutoff_propagator;
    // Force the use of a Dense Tensor evaluation
    bool numerics_pathintegral_force_dense;
    // Dynamic Cutoff; While true, the squared threshold will be increased or decreased until the number of ADM elements is approximately equal to the cutoff iterations set.
    long long numerics_pathintegral_dynamiccutoff_iterations_max; //=0
    // PI Density Change Coefficients
    int numerics_dynamic_densitychange_limit = 3; // If fillrate doesnt change for this number of iterations, switch tensor to dense

    bool numerics_pathintegral_use_qrt;
    bool numerics_pathintegral_set_couplings_zero;

    // Cache Switches
    // Enables caching of the PME Coefficients
    bool numerics_use_saved_coefficients;
    // Force Caching of the PME Coefficients even if use_saved_coefficients is false.
    bool numerics_enable_saving_coefficients;
    // Enables caching of the Hamilton Operators.
    bool numerics_use_saved_hamiltons;
    // Enables the caching of pulse- and chirp functions
    bool numerics_use_function_caching;

    bool detector_normalize_functions = true;

    // Unique Identifier for a given set of parameters for later evaluation
    std::vector<double> logfilecounter;

    // Interpolation Settings
    // Input String for interpolation Parameters
    std::string s_numerics_interpolate;
    // Interpolation Method for the Main time direction which will be output to files
    int numerics_interpolate_method_time;
    // Interpolation Method for the correlation grid
    int numerics_interpolate_method_tau;

    // Stepparameters for the RK45 Method
    Parameter numerics_rk_stepdelta, numerics_rk_stepmin, numerics_rk_stepmax, numerics_rk_order;
    std::vector<std::tuple<double, double>> numerics_rk_tol;
    // Enables the use of discrete timesteps using rk_stepdelta as an in- or decrease for dt
    bool numerics_rk_usediscrete_timesteps;

    // Time variables
    Parameter t_start, t_end, t_step, t_step_pathint;
    double numerics_subiterator_stepsize; // = 1E-12;

    // System Dimensions
    int maxStates;

    // Initial state
    // String state representation of the initial state
    std::string p_initial_state_s;
    // Matrix Index of the initial state
    int p_initial_state;

    // State Transition Delimiter
    std::string transition_delimiter = "=";

    // Environmental Coupling Constantes
    // QD-Cavity Coupling Rate
    Parameter p_omega_coupling;
    // Cavity-Environment Coupling / Cavity Decay Rate
    Parameter p_omega_cavity_loss;
    // Pure Dephasing Rate
    Parameter p_omega_pure_dephasing;
    // Radiative Decay Rate
    Parameter p_omega_decay;

    // Main Direction Trigger. Because some system and solver functions use e.g. cached variables, switching this trigger after the main t-direction is calculated enables those functions to adapt.
    // For example, the polaron functions will start to use cached and interpolated values as soon as this trigger is switched
    bool numerics_main_direction_done;

    // Custom Expectation Values
    std::vector<std::string> numerics_custom_expectation_values;

    // Constructor
    Parameters() = default;
    explicit Parameters( const std::vector<std::string> &arguments );

    /**
     * @brief A general input struct containing parameter maps. These maps either map to a string, a vector of strings, a parameter or a vector of parameters.
     * The struct also enables the general output of a given inputstruct for easier debugging.
     *
     */
    class universal_config {
       public:
        universal_config() = default;
        // String -> Parameter
        std::map<std::string, Parameter> property;
        // String -> String
        std::map<std::string, std::string> string;
        // String -> Parameter Vector
        std::map<std::string, std::vector<Parameter>> property_set;
        // String -> String Vector
        std::map<std::string, std::vector<std::string>> string_v;
        // Output Helper
        friend std::ostream &operator<<( std::ostream &os, const universal_config &is ) {
            os << "property Values:\n";
            for ( auto &p : is.property )
                os << p.first << " = " << p.second << "\n";
            os << std::endl;
            os << "String Values:\n";
            for ( auto &p : is.string )
                os << p.first << " = " << p.second << "\n";
            os << std::endl;
            if ( is.property_set.size() > 0 ) {
                os << "property Vector Values:\n";
                for ( auto &p : is.property_set ) {
                    os << p.first << " = ";
                    for ( auto &u : p.second )
                        os << u << " ";
                    os << std::endl;
                }
            }
            if ( is.string_v.size() > 0 ) {
                os << "String Vector Values:\n";
                for ( auto &p : is.string_v ) {
                    os << p.first << " = ";
                    for ( auto &u : p.second )
                        os << u << " ";
                    os << std::endl;
                }
            }
            return os;
        }

        const Parameter &get_value( const std::string &key ) const {
            return property.at( key );
        }
        const std::string &get_string( const std::string &key ) const {
            return string.at( key );
        }
        const std::vector<Parameter> &get_value_v( const std::string &key ) const {
            return property_set.at( key );
        }
        const std::vector<std::string> &get_string_v( const std::string &key ) const {
            return string_v.at( key );
        }
        Parameter &get_value( const std::string &key ) {
            return property.at( key );
        }
        std::string &get_string( const std::string &key ) {
            return string.at( key );
        }
        std::vector<Parameter> &get_value_v( const std::string &key ) {
            return property_set.at( key );
        }
        std::vector<std::string> &get_string_v( const std::string &key ) {
            return string_v.at( key );
        }
        /**
         * @brief Universal Get Method, when no index is provided, the property map is used, if index is provided, the property vector is used.
         *
         */
        static constexpr size_t no_index = 133713371337;
        Parameter &get( const std::string &key, const size_t &index = no_index ) {
            if ( index != no_index )
                return property_set.at( key ).at( index );
            return property.at( key );
        }
        const Parameter &get( const std::string &key, const size_t &index = no_index ) const {
            if ( index != no_index )
                return property_set.at( key ).at( index );
            return property.at( key );
        }
    };
    // Temporary Inputstrings
    std::string inputstring_electronic;
    std::string inputstring_photonic;
    std::string inputstring_pulse;
    std::string inputstring_chirp;
    std::string inputstring_spectrum, inputstring_indist, inputstring_conc, inputstring_gfunc, inputstring_wigner, inputstring_raman, inputstring_correlation_resolution, inputstring_SPconf;
    std::string inputstring_detector_time,inputstring_detector_spectral;
    std::string inputstring_densitymatrix_config;
    std::string inputstring_rk45_config;
    // Input Maps. The Temporary Inputstrings will get verwurstet into these maps, such that their parameters are available via their corresponding string key.
    std::map<std::string, universal_config> input_electronic;
    std::map<std::string, universal_config> input_photonic;
    std::map<std::string, universal_config> input_pulse;
    std::map<std::string, universal_config> input_chirp;
    std::map<std::string, std::vector<universal_config>> input_correlation; // Spectrum, Indist, Conc
    std::map<std::string, universal_config> input_correlation_resolution;   // G1/G2 correlation timesteps. length of g will be determined by gridres
    std::map<std::string, universal_config> input_conf;                     // Everything else

    /**
     * @brief Converts the inputstrings into input_maps. These maps will then be used to generate the operators and to output the inputsystem
     *
     */
    void parse_system();

    /**
     * @brief Log function; Uses log subclass (log.h). The initial state vector can be passed to be logged.
     *
     */
    void log( const Dense &initial_state_vector_ket );

    /**
     * @brief Parses the input from a string input vector. Uses the Parse_Parameters class from misc/commandlinearguments.h
     *
     */
    void parse_input( const std::vector<std::string> &arguments );

    /**
     * @brief Adjusts input before and after the main RK loop
     *
     */
    void pre_adjust_input();
    void post_adjust_input();
    void post_adjust_grids();
};
} // namespace QDLC