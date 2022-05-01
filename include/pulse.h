#pragma once
#include "global.h"
#include "system/parameters.h"

class Pulse {
   public:
    // class Inputs {
    //    public:
    //     std::vector<std::string> type;
    //     std::vector<double> center;
    //     std::vector<std::complex<double>> amp;
    //     std::vector<double> sigma;
    //     std::vector<double> omega;
    //     std::vector<double> omega_chirp;
    //     std::vector<double> super_amp;
    //     double t_start, t_end, t_step;
    //     int order;
    //     Inputs(){};
    //     Inputs( double t_start, double t_end, double t_step, int order ) : t_start( t_start ), t_end( t_end ), t_step( t_step ), order( order ) {}
    //     void add( double _center, double _amp, double _sigma, double _omega, double _omega_chirp, std::string _type );
    //     void add( std::vector<Parameter> &_center, std::vector<Parameter> &_amp, std::vector<Parameter> &_sigma, std::vector<Parameter> &_omega, std::vector<Parameter> &_omega_chirp, std::vector<Parameter> &_super_amp, std::vector<std::string> &_type, std::complex<double> amp_scaling = 1.0 );
    //     void add( std::vector<Parameter> &_center, std::vector<Parameter> &_amp, std::vector<Parameter> &_sigma, std::vector<Parameter> &_omega, std::vector<Parameter> &_omega_chirp, std::vector<Parameter> &_super_amp, std::vector<std::string> &_type, std::vector<std::string> &_filter, std::string to_match, std::complex<double> amp_scaling = 1.0 );
    // };

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