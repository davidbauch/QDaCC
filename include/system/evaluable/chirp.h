#pragma once
#include "global.h"
#include "misc/interpolant.h"
#include "system/parameters.h"
#include "system/evaluable/evaluable.h"

namespace QDACC {

class Chirp : public Evaluable {
   private:
    Interpolant interpolant;
    bool isSineChirp;

   public:
    Chirp() : Evaluable(){};
    Chirp( Parameters::universal_config& config, Parameters& p );
    Scalar evaluate( double t );
    void calculate_fourier( Parameters& p );
    void log();
};
} // namespace QDACC