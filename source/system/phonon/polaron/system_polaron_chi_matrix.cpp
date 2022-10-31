#include "system/system.h"

using namespace QDLC;

Sparse System::dgl_phonons_chi( const double t ) {
    Log::L3( "[System-PME] Generating Chi({})...\n", t );
    // QD-Cavity
    Sparse ret = operatorMatrices.polaron_factors[0];
    // QD-Pulse
    for ( int i = 0; i < pulse.size(); i++ ) {
        auto p = pulse[i].get( t );
        if ( parameters.numerics_use_rwa )
            ret += operatorMatrices.pulse_mat[2 * i + 1] * p;
        else
            ret += operatorMatrices.pulse_mat[2 * i + 1] * ( p + std::conj( p ) );
        Log::L3( "[System-PME]     Added Pulse with Omega(t) = {}\n", pulse[i].get( t ) );
    }
    Log::L3( "[System-PME]     Done, final untransformed Chi:\n{}\n", Dense( ret ).format( operatorMatrices.output_format ) );
    return dgl_timetrafo( ret, t );
}