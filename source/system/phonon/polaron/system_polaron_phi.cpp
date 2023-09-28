#include "system/system.h"

using namespace QDACC;

Scalar System::dgl_phonons_phi( const double tau ) {
    Scalar integral = 0;
    for ( double w = parameters.p_phonon_wcutoffdelta; w < 10.0 * parameters.p_phonon_wcutoff; w += parameters.p_phonon_wcutoffdelta ) {
        Scalar J = dgl_phonons_spectral_density( w ) / w / w;
        integral += parameters.p_phonon_wcutoffdelta * J * ( std::cos( w * tau ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) - 1.0i * std::sin( w * tau ) );
    }
    return integral;
}