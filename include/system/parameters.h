#pragma once
// Dependencies
#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "misc/timer.h"
#include "system/parameter.h"

#define GLOBAL_PROGRAM_VERSION "3.3"
#define GLOBAL_PROGRAM_LASTCHANGE "Code Revamp"

class Parameters {
   public:
    // Numerical Parameters
    Parameter hbar, kb;
    std::string subfolder;
    bool numerics_use_interactionpicture, numerics_use_rwa, numerics_calculate_till_converged;
    // Also output electronic emission probabilities
    int numerics_order_timetrafo;
    int output_advanced_log, output_handlerstrings, output_operators, output_coefficients;
    int iterations_t_skip, iterations_tau_resolution, iterations_w_resolution;
    bool scale_parameters;
    double scale_value;
    // Runtime parameters and other stuff
    int iterations_t_max;
    std::vector<double> trace;
    bool output_full_dm, output_no_dm;
    int numerics_maximum_threads, numerics_phonon_approximation_markov1, numerics_phonon_approximation_order;
    double akf_deltaWmax;
    Parameter spectrum_frequency_center, spectrum_frequency_range;
    bool numerics_use_saved_coefficients, numerics_force_caching;
    bool numerics_output_raman_population;
    bool numerics_calculate_timeresolution_indistinguishability;
    int numerics_phonons_maximum_threads;
    bool numerics_use_saved_hamiltons;
    long unsigned int numerics_saved_coefficients_cutoff; // True: Only save last few coefficients (only viable for T-direction, not for G1/2)
    long unsigned int numerics_saved_coefficients_max_size;
    std::vector<int> logfilecounter;
    bool numerics_interpolate_outputs;
    std::string s_numerics_interpolate;
    int numerics_interpolate_method_time, numerics_interpolate_method_tau;
    Parameter numerics_rk_stepdelta, numerics_rk_stepmin, numerics_rk_stepmax, numerics_rk_order, numerics_rk_tol;
    bool numerics_rk_usediscrete_timesteps;
    bool numerics_pathint_partially_summed;
    bool numerics_phonon_nork45;
    bool numerics_use_function_caching;
    bool output_path;
    size_t numerics_groundstate;

    // Path Integral Numerics
    double numerics_pathintegral_stepsize_iterator;      // = 1E-12;
    double numerics_pathintegral_squared_threshold;      // = 1E-32;
    double numerics_pathintegral_sparse_prune_threshold; // = 1E-1;
    // Propagator Cutoff; When iterated multiple times, M(t0->t1) may gain additional entries in between the RK iterations from t0 to t1. When set to true, the final non-zero matrix entries will be mapped onto the non-zero entries after
    // the first iteration, meaning any additional non-zero entries besides the ones created within the first iteration are lost.
    bool numerics_pathintegral_docutoff_propagator;
    // Dynamic Cutoff; While true, the squared threshold will be increased or decreased until the number of ADM elements is approximately equal to the cutoff iterations set.
    long long numerics_pathintegral_dynamiccutoff_iterations_max; //=0

    // System Parameters
    // Time variables
    Parameter t_start, t_end, t_step, t_step_pathint;
    std::vector<double> grid_steps, grid_values; // Calculate T and Correlation Grid for these timesteps. maps Index -> t_t, t_step
    std::map<double, size_t> grid_value_indices; // maps t_t -> Index
    // System Dimensions
    int maxStates;
    // Starting state:
    std::string p_initial_state_s;
    int p_initial_state;

    Parameter p_omega_coupling;
    Parameter p_omega_cavity_loss;
    Parameter p_omega_pure_dephasing;
    Parameter p_omega_decay;
    Parameter p_phonon_b, p_phonon_alpha, p_phonon_wcutoff, p_phonon_T, p_phonon_tcutoff, p_phonon_pure_dephasing;
    bool p_phonon_adjust;
    int p_phonon_nc;

    // Constructor
    Parameters(){};
    Parameters( const std::vector<std::string> &arguments );

    // Variables for the inputstrings and input vectors
    struct input_s {
        std::map<std::string, Parameter> numerical;
        std::map<std::string, std::string> string;
        std::map<std::string, std::vector<Parameter>> numerical_v;
        std::map<std::string, std::vector<std::string>> string_v;
        friend std::ostream &operator<<( std::ostream &os, const input_s &is ) {
            os << "Numerical Values:\n";
            for ( auto &p : is.numerical )
                os << p.first << " = " << p.second << "\n";
            os << std::endl;
            os << "String Values:\n";
            for ( auto &p : is.string )
                os << p.first << " = " << p.second << "\n";
            os << std::endl;
            if ( is.numerical_v.size() > 0 ) {
                os << "Numerical Vector Values:\n";
                for ( auto &p : is.numerical_v ) {
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
        };
    };
    std::string inputstring_electronic;
    std::string inputstring_photonic;
    std::string inputstring_pulse;
    std::string inputstring_chirp;
    std::string inputstring_spectrum, inputstring_indist, inputstring_conc, inputstring_gfunc, inputstring_wigner, inputstring_raman, inputstring_correlation_resolution;
    std::map<std::string, input_s> input_electronic;
    std::map<std::string, input_s> input_photonic;
    std::map<std::string, input_s> input_pulse;
    std::map<std::string, input_s> input_chirp;
    std::map<std::string, input_s> input_correlation;            // spectrum, indist, conc
    std::map<std::string, input_s> input_correlation_resolution; // g1/g2 correlation timesteps. length of g will be determined by gridres
    // Converts the input strings into input vectors.
    // These Vectors will then be used to generate the operators and to output the inputsystem
    void parse_system();

    // Log function; Uses log subclass (log.h)
    // @param &info: Additional information this class does not have access too when created, e.g. basis
    void log( const Dense &initial_state_vector_ket );

    // Help function, output when --help is called
    static void help();

    // Parses the input from a string input vector. Uses the Parse_Parameters class from misc/commandlinearguments.h
    // @param &arguments: Vector of string arguments, e.g. from argv
    // @return Returns true if parsing was successful
    bool parseInput( const std::vector<std::string> &arguments );

    // Adjusting inputs
    // @return Returns true if successful
    bool adjustInput();

    // Scale inputs. Has to be called before anything else is calculated //TODO: finalize
    // @param scaling: Value to scale with, e.g. 1E12 for ps
    // @return Returns true if scaling was successful
    bool scaleInputs( const double scaling );

    // Scale single input
    // @param variable: Numerical value of variable to scale
    // @param scaling: Value to scale with, e.g. 1E12 for ps
    // @return Returns variable*scaling
    template <typename T>
    T scaleVariable( const T variable, const double scaling ) {
        if ( scale_parameters ) {
            return variable * scaling;
        }
        return variable;
    }
};