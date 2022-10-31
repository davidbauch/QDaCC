#pragma once
#include "global.h"
#include "system/parameters.h"
#include "system/evaluable/evaluable.h"

namespace QDLC {

class Pulse : public Evaluable {
   public:
    Pulse() : Evaluable(){};
    Pulse( Parameters::input_s& config, Parameters& p );
    Scalar evaluate( double t );
    Scalar evaluate_derivative( double t, double dt = 0 );
    Scalar evaluate_integral( double t, double dt = 0 );
    void calculate_fourier( Parameters& p );
    void log();
};
} // namespace QDLC