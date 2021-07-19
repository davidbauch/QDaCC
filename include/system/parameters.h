#pragma once
// Dependencies
#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "misc/timer.h"
#include "system/parameter.h"

#define GLOBAL_PROGRAM_VERSION "3.0"
#define GLOBAL_PROGRAM_LASTCHANGE "Procedural System"

class Parameters {
   public:
    // Numerical Parameters
    Parameter hbar, kb;
    std::string subfolder;
    bool numerics_calculate_spectrum, numerics_calculate_g2, numerics_use_simplified_g2, numerics_use_interactionpicture, numerics_use_rwa;
    // Also output electronic emission probabilities
    bool numerics_output_electronic_emission;
    int numerics_order_timetrafo;
    int numerics_order_t, numerics_order_tau, numerics_order_highest;
    int output_advanced_log, output_handlerstrings, output_operators, output_coefficients;
    int iterations_t_skip, iterations_tau_resolution, iterations_w_resolution;
    int iterations_wtau_skip;
    bool scale_parameters;
    double scale_value;
    // Runtime parameters and other stuff
    int iterations_t_max;
    std::vector<double> trace;
    bool output_full_dm, output_no_dm;
    int numerics_maximum_threads, numerics_phonon_approximation_markov1, numerics_phonon_approximation_order;
    double akf_deltaWmax;
    Parameter spectrum_frequency_center, spectrum_frequency_range;
    bool numerics_use_saved_coefficients;
    bool numerics_output_raman_population;
    bool numerics_calculate_timeresolution_indistinguishability;
    int numerics_phonons_maximum_threads;
    bool numerics_use_saved_hamiltons;
    long unsigned int numerics_saved_coefficients_cutoff; // True: Only save last few coefficients (only viable for T-direction, not for G1/2)
    long unsigned int numerics_saved_coefficients_max_size;
    std::vector<int> logfilecounter;
    bool numerics_stretch_correlation_grid, numerics_interpolate_outputs;
    Parameter numerics_rk_stepdelta, numerics_rk_stepmin, numerics_rk_stepmax, numerics_rk_order, numerics_rk_tol;
    bool numerics_rk_interpolate;

    // Path Integral Numerics
    double numerics_pathintegral_stepsize_iterator;      // = 1E-12;
    double numerics_pathintegral_squared_threshold;      // = 1E-32;
    double numerics_pathintegral_sparse_prune_threshold; // = 1E-1;
    // Propagator Cutoff; When iterated multiple times, M(t0->t1) may gain additional entries in between the RK iterations from t0 to t1. When set to true, the final non-zero matrix entries will be mapped onto the non-zero entries after
    // the first iteration, meaning any additional non-zero entries besides the ones created within the first iteration are lost.
    bool numerics_pathintegral_docutoff_propagator; //=false
    // Dynamic Cutoff; While true, the squared threshold will be increased or decreased until the number of ADM elements is approximately equal to the cutoff iterations set.
    long long numerics_pathintegral_dynamiccutoff_iterations_max; //=0

    // System Parameters
    // Time variables
    Parameter t_start, t_end, t_step, t_step_pathint;
    // System Dimensions
    int maxStates; //, p_max_photon_number; //REMOVE
    // Starting state:
    //int p_initial_state, p_initial_state_photon_h, p_initial_state_photon_v; //REMOVE
    //std::string p_initial_state_electronic;                                  //REMOVE
    std::string p_initial_state_s;
    int p_initial_state;
    bool startCoherent;

    // Non mandatory parameters, dependant on system chosen:
    // System Parameterss //TODO: die hier dann in einen vector<Parameter> mit den zugehörigen operator matrices. KEINE (!) neue subclass mit matrix und energie oder so, das is unnötig. einfach alles in vectoren schreiben.
    //Parameter p_omega_atomic_G_V; //TODO: die hier generieren in map p_omega_transitions oder so
    //Parameter p_omega_atomic_G_H; //TODO: die hier generieren in map p_omega_transitions oder so
    //Parameter p_omega_atomic_V_B; //TODO: die hier generieren in map p_omega_transitions oder so
    //Parameter p_omega_atomic_H_B; //TODO: die hier generieren in map p_omega_transitions oder so
    //Parameter p_omega_atomic_B;
    //Parameter p_omega_cavity_V;
    //Parameter p_omega_cavity_H;
    Parameter p_omega_coupling;
    Parameter p_omega_cavity_loss;
    Parameter p_omega_pure_dephasing;
    Parameter p_omega_decay;
    Parameter p_phonon_b, p_phonon_alpha, p_phonon_wcutoff, p_phonon_T, p_phonon_tcutoff, p_phonon_pure_dephasing;
    bool p_phonon_adjust;
    int p_phonon_nc;
    //Parameter p_deltaE;
    //Parameter p_biexciton_bindingenergy;

    // Calculated System properties:
    //Parameter init_detuning_G_H, init_detuning_G_V, init_detuning_H_B, init_detuning_V_B, max_detuning_G_H, max_detuning_G_V, max_detuning_H_B, max_detuning_V_B;
    //Parameter init_rabifrequenz_G_H, init_rabifrequenz_G_V, init_rabifrequenz_H_B, init_rabifrequenz_V_B, max_rabifrequenz_G_H, max_rabifrequenz_G_V, max_rabifrequenz_H_B, max_rabifrequenz_V_B;
    //Parameter init_detuning, max_detuning, init_rabifrequenz, max_rabifrequenz;

    // Chirp and Pulse properties:
    //std::vector<Parameter> pulse_center, pulse_amp, pulse_omega, pulse_sigma, pulse_omega_chirp;
    //std::vector<std::string> pulse_type;
    //std::vector<std::string> pulse_pol;
    //double chirp_total;
    //std::vector<Parameter> chirp_t, chirp_y, chirp_ddt;
    //std::string chirp_type;

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
            //os << "Numerical Values:\n";
            //for ( auto &p : is.numerical )
            //    os << p.first << " = " << p.second << "\n";
            //os << "Coupled_to Values:\n";
            //for ( auto &p : is.coupled_to )
            //    os << p << " ";
            //os << std::endl;
            //if ( is.numerical_v.size() > 0 ) {
            //    os << "Numerical Vector Values:\n";
            //    for ( auto &p : is.numerical_v ) {
            //        os << p.first << " = ";
            //        for ( auto &u : p.second )
            //            os << u << " ";
            //        os << std::endl;
            //    }
            //}
            return os;
        };
    };
    std::string inputstring_electronic; // = "G:0:H,V:0:0:0;H:1.365999eV:Z:1:1:1;V:1.366001eV:Z:1:1:1;Z:2.729eV:-:1:1:2";
    std::string inputstring_photonic;   // = "h:1.366eV:2:GH,HZ:1,1:1;v:1.366eV:2:GV,VZ:1,1:1";
    std::string inputstring_pulse;      //"p:GH,HZ:6:1.3645eV:4ps:30ps:0:gauss_pi"; //"h:GH,HZ:1,4:1.366eV,1.366eV:4ps,15ps:30ps,15ps:0,0:gauss,gauss;v:GV,VZ:1:1.366eV:4ps:15ps:0:gauss";
    std::string inputstring_chirp;      //"1:H,V,Z:1,1,2:0,1meV,0:0,100ps,200ps:0,0,0:monotone";
    //std::string inputstring_electronic = "G:0:X:0:0:0;X:1.366eV:-:1:1:1";
    //std::string inputstring_photonic = "c:1.366eV:1:GX:1:1";
    //std::string inputstring_pulse = "h:GX:1:1.366eV:4ps:30ps:0:gauss_pi";
    //std::string inputstring_chirp = ""; //1:X:1:0,1meV,0:0,100ps,200ps:0,0,0:monotone";
    std::string inputstring_spectrum, inputstring_indist, inputstring_conc, inputstring_gfunc, inputstring_wigner;
    std::map<std::string, input_s> input_electronic;
    std::map<std::string, input_s> input_photonic;
    std::map<std::string, input_s> input_pulse;
    std::map<std::string, input_s> input_chirp;
    std::map<std::string, input_s> input_correlation; //spectrum, indist, conc
    // Converts the input strings into input vectors.
    // These Vectors will then be used to generate the operators and to output the inputsystem
    void parse_system();

    // Log function; Uses log subclass (log.h)
    // @param &info: Additional information this class does not have access too when created, e.g. basis
    void log( const std::vector<std::string> &info );

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