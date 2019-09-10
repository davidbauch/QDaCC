#pragma once
#include "global.h"
#include "misc/interpolant.cpp"

class Chirp {
   private:
    class Inputs {
       public:
        Inputs() {}
        Inputs( Parameters &p ) {
        }
    };
    // Chirp values
    std::vector<double> chirparray;
    // Saving the corresponding t
    std::vector<double> timearray;
    // Helper variables/functions
    Inputs inputs;
    int size;
    double step;
    std::vector<double> steps;
    Interpolant interpolant;

   public:
    Chirp(){};
    Chirp( Parameters &p );
    void generate( Parameters &p );
    void fileOutput( std::string filepath, Parameters &p );
    double get( double t ) const;
    double get( int i ) const;
};