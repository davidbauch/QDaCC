#pragma once
#include "global.h"

class Pulse {
   public:
    class Inputs {
       public:
        std::vector<std::string> type;
        std::vector<double> center;
        std::vector<double> amp;
        std::vector<double> sigma;
        std::vector<double> omega;
        double t_start, t_end, t_step;
        int order;
        Inputs(){};
        Inputs( double t_start, double t_end, double t_step, int order ) : t_start( t_start ), t_end( t_end ), t_step( t_step ), order( order ) {}
        void add( double _center, double _amp, double _sigma, double _omega, std::string _type );
        void add( std::vector<double> &_center, std::vector<double> &_amp, std::vector<double> &_sigma, std::vector<double> &_omega, std::vector<std::string> &_type );
        void add( std::vector<double> &_center, std::vector<double> &_amp, std::vector<double> &_sigma, std::vector<double> &_omega, std::vector<std::string> &_type, std::vector<std::string> &_filter, std::string to_match );
    };

   private:
    int size;
    int counter_evaluated, counter_returned;
    void generate( Parameters &p );

   public:
    // Pulse values
    std::vector<dcomplex> pulsearray;
    // Saving the corresponding t doesn't hurt
    std::vector<double> timearray;
    std::vector<double> steps;
    // Helper variables/functions
    Inputs inputs;

    Pulse(){};
    Pulse( Inputs &inputs );
    void generate();
    void fileOutput( std::string filepath );
    dcomplex evaluate( double t );
    dcomplex get( double t, bool force_evaluate = false );
    void log();
    static void fileOutput( std::string filepath, std::vector<Pulse> pulses );
};