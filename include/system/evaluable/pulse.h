#pragma once
#include "global.h"
#include "system/parameters.h"
#include "system/evaluable/evaluable.h"

namespace QDACC {

class Pulse : public Evaluable {
   public:
    Pulse() : Evaluable(){};
    Pulse( Parameters::universal_config& config, Parameters& p );
    Scalar evaluate( double t );
    void calculate_fourier( Parameters& p );
    void log();
};
} // namespace QDACC