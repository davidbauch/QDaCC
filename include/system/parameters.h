#pragma once
// Dependencies
#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "misc/timer.h"

template <typename T>
class Parameter {
    private:
        T val_si;
        T val_scaled;
        bool scaled;
        double scale_factor;
        void update() {
            if ( scale_factor != 0.0 ) {
                val_scaled = val_si*scale_factor;
                scaled = true;
            }
            else {
                val_scaled = val_si;
                scaled = false;
            }
        }
    public:
        Parameter( ) {};
        Parameter( T val ) : val_si(val), scale_factor(0.0), scaled(false) {};
        Parameter( T val, double scale_factor ) : val_si(val), scale_factor(scale_factor), scaled(true) {};
        T get() {
            return val_scaled;
        }
        T getSI() {
            return val_si;
        }
        T set( T new_val, double scale = 0.0 ) {
            if ( scale != 0.0 ) {
                scale_factor = scale;
            } 
            val_si = new_val;
            update();
            return val_scaled;
        }
        T setScale( double scale = 0.0 ) {
            return set(val_si, scale);
        }
        parameter =(T val) {
            set(val);
        }
        friend std::ostream& operator <<( std::ostream &stream, const Parameter &param ) {
            return stream << param.val_scaled;
        }
};

class Parameters {
   public:
        // Numerical Parameters
        Parameter<double> hbar, kb;
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
        double akf_deltaWmax, spectrum_frequency_center, spectrum_frequency_range;
        bool numerics_calculate_spectrum_H, numerics_calculate_spectrum_V;
        bool numerics_calculate_g2_H, numerics_calculate_g2_V, numerics_calculate_g2_C;
        bool numerics_use_saved_coefficients;
        bool numerics_output_raman_population;
        bool numerics_calculate_timeresolution_indistinguishability;
        int numerics_phonons_maximum_threads;
        bool numerics_use_saved_hamiltons;
        long unsigned int numerics_saved_coefficients_cutoff; // True: Only save last few coefficients (only viable for T-direction, not for G1/2)
        long unsigned int numerics_saved_coefficients_max_size;
        int logfilecounter;
        bool numerics_stretch_correlation_grid;

        // System Parameters
        // Time variables
        double t_start, t_end, t_step;
        // System Dimensions
        int maxStates, p_max_photon_number;
        // Starting state:
        int p_initial_state;
        bool startCoherent;

        // Non mandatory parameters, dependant on system chosen:
        // System Parameterss
        Parameter<double> p_omega_atomic_G_V;
        Parameter<double> p_omega_atomic_G_H;
        Parameter<double> p_omega_atomic_V_B;
        Parameter<double> p_omega_atomic_H_B;
        Parameter<double> p_omega_atomic_B;
        Parameter<double> p_omega_cavity_V;
        Parameter<double> p_omega_cavity_H;
        Parameter<double> p_omega_coupling;
        Parameter<double> p_omega_cavity_loss;
        Parameter<double> p_omega_pure_dephasing;
        Parameter<double> p_omega_decay;
        Parameter<double> p_phonon_b, p_phonon_alpha, p_phonon_wcutoff, p_phonon_T, p_phonon_tcutoff, p_phonon_pure_dephasing;
        bool p_phonon_adjust;
        Parameter<double> p_deltaE;
        Parameter<double> p_biexciton_bindingenergy;

        // Calculated System properties:
        double init_detuning_G_H, init_detuning_G_V, init_detuning_H_B, init_detuning_V_B, max_detuning_G_H, max_detuning_G_V, max_detuning_H_B, max_detuning_V_B;
        double init_rabifrequenz_G_H, init_rabifrequenz_G_V, init_rabifrequenz_H_B, init_rabifrequenz_V_B, max_rabifrequenz_G_H, max_rabifrequenz_G_V, max_rabifrequenz_H_B, max_rabifrequenz_V_B;
        double init_detuning, max_detuning, init_rabifrequenz, max_rabifrequenz;

        // Chirp and Pulse properties:
        std::vector<Parameter<double>> pulse_center, pulse_amp, pulse_omega, pulse_sigma, pulse_omega_chirp;
        std::vector<std::string> pulse_type;
        std::vector<std::string> pulse_pol;
        double chirp_total;
        std::vector<Parameter<double>> chirp_t, chirp_y, chirp_ddt;
        std::string chirp_type;

        // Constructor
        Parameters(){};
        Parameters( const std::vector<std::string> &arguments );

        // Log function; Uses log subclass (log.h)
        // @param &info: Additional information this class does not have access too when created, e.g. basis
        void log(const std::vector<std::string> &info);

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
        
        // Converts a named state vector, e.g. |g,0,0> into a matrix index
        // @param mode: Element of G,H,V,E
        // @param state: Photon number, 0,1,...
        // @return Returns matrix index corresponding to |X,n,m>
        int index_to_state( const char mode = 'E', const int state = 0 );

        // Approximates the ideal timestep for the initial system parameters chosen
        // @return Returns ideal timestep
        double getIdealTimestep();

};