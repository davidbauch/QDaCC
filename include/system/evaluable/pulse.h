#pragma once
#include "global.h"
#include "system/parameters.h"

class Pulse {
   private:
    int size;
    int counter_evaluated, counter_returned;

   public:
    // Pulse values
    std::map<double, Scalar> pulsearray;
    std::map<double, Scalar> pulsearray_derivative;
    std::map<double, Scalar> pulsearray_integral;
    std::vector<Scalar> pulsearray_fourier;
    std::vector<double> fourier;
    std::vector<double> steps;
    double maximum;
    // Helper variables/functions
    Parameters::input_s inputs;

    Pulse(){};
    Pulse( Parameters::input_s& inputs, Parameters& p );
    void generate( double t_start, double t_end, double t_step, double omega_center = 0, double omega_range = 0, double dw = 0, double dt = 0 );
    void fileOutput( std::string filepath, double t_start, double t_end, double t_step );
    Scalar evaluate( double t );
    Scalar evaluate_derivative( double t, double dt );
    Scalar evaluate_integral( double t, double dt );
    Scalar get( double t, bool force_evaluate = false );
    Scalar derivative( double t, double t_step, bool force_evaluate = false );
    Scalar integral( double t, double t_step, bool force_evaluate = false );
    void log();
    static void fileOutput( std::string filepath, std::vector<Pulse> pulses, double t_start, double t_end, double t_step );
};