#pragma once
#include "global.h"
#include "misc/interpolant.cpp"

class Chirp {
   public:
    class Inputs {
       public:
        double t_start, t_end, t_step;
        std::vector<double> t, y, ddt;
        std::string type;
        int order;
        Inputs() {}
        Inputs( double t_start, double t_end, double t_step, std::string type, int order ) : t_start( t_start ), t_end( t_end ), t_step( t_step ), type( type ), order( order ) {}
        void add( double _t, double _y, double _ddt );
        void add( std::vector<double> &_t, std::vector<double> &_y, std::vector<double> &_ddt );
    };

   private:
    // Chirp values
    std::vector<double> chirparray;
    // Saving the corresponding t
    std::vector<double> timearray;
    // Helper variables/functions
    int size;
    std::vector<double> steps;
    Interpolant interpolant;
    Inputs inputs;

   public:
    Chirp(){};
    Chirp( Inputs &inputs );
    void generate();
    void fileOutput( std::string filepath );
    double get( double t ) const;
    double get( int i ) const;
};