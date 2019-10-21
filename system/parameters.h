#pragma once
// Dependencies
#include <stdlib.h>
#include "../misc/helperfunctions.h"
#include "../misc/log.h"
#include "../misc/timer.h"

class Parameters_Parent {
    // Mandatory Variables
   public:
    std::string subfolder;
    int numerics_calculate_spectrum, numerics_calculate_g2, numerics_use_interactionpicture, numerics_use_rwa, numerics_order_timetrafo;
    int numerics_order_t, numerics_order_tau, numerics_order_highest;
    int output_advanced_log, output_handlerstrings, output_operators, output_coefficients;
    int iterations_skips_t, iterations_tau_resolution, iterations_w_resolution;
    bool scale_parameters;
    double scale_value;
    // Functions
   public:
    // Constructor
    Parameters_Parent();
    Parameters_Parent( const std::vector<std::string> &arguments ) : Parameters_Parent(){};
    // Log function; Uses log subclass (log.h)
    virtual void log(){};
    // Help function, output when --help is called
    static void help();
    // @overwrite: Parsing inputs from string input vector. Uses parseinput functions
    virtual bool parseInput( const std::vector<std::string> &arguments ) { return false; };
    // @overwrite: Adjusting inputs
    virtual bool adjustInput() { return false; };
    // @overwrite: Scale inputs
    virtual bool scaleInputs( const double scaling ) { return false; };
    // @overrwrite: Scale single input
    virtual double scaleVariable( const double variable, const double scaling ) { return 0.0; };
    // To be called from child constructor
    void init( const std::vector<std::string> &arguments );
};