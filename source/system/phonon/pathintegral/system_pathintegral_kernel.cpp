#include "system/system.h"

Scalar System::dgl_phonons_kernel( const double t, const double t_step ) {
    Scalar integral = 0;
    const double integral_upper_limit = 10.0 * parameters.p_phonon_wcutoff;
    for ( double w = 1E-10; w < integral_upper_limit; w += parameters.p_phonon_wcutoffdelta ) {
        double J = dgl_phonons_spectral_density(w);
        if ( t < t_step / 2.0 ) {
            const Scalar polaron_shift = -1.i * w * t_step;
            integral += parameters.p_phonon_wcutoffdelta * J * ( ( 1.0 - std::cos( w * t_step ) ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) + 1.0i * std::sin( w * t_step ) ); // + polaron_shift );
        } else {
            integral += parameters.p_phonon_wcutoffdelta * 2.0 * J * ( 1.0 - std::cos( w * t_step ) ) * ( std::cos( w * t ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) - 1.0i * std::sin( w * t ) );
        }
    }
    return integral;
}