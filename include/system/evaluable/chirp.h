#pragma once
#include "global.h"
#include "misc/interpolant.h"
#include "system/parameters.h"
#include "system/evaluable/evaluable.h"

namespace QDLC {

class Chirp : public Evaluable {
   private:
    Interpolant interpolant;
    bool isSineChirp;

   public:
    Chirp() : Evaluable(){};
    Chirp( Parameters::universal_config& config, Parameters& p );
    Scalar evaluate( double t );
    Scalar evaluate_derivative( double t, double dt = 0 );
    Scalar evaluate_integral( double t, double dt = 0 );
    void calculate_fourier( Parameters& p );
    void log();
};
} // namespace QDLC