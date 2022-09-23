#include "system/system.h"

double System::dgl_phonons_spectral_density( const double w ) {
    if ( parameters.p_phonon_qd_ae == 0.0 ) {
        const double gaussian_shape = std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff );
        return parameters.p_phonon_alpha * w * gaussian_shape;
    } else {
        const double gaussian_shape_electron = std::exp( -w * w * parameters.p_phonon_qd_ae * parameters.p_phonon_qd_ae / ( 4. * parameters.p_phonon_qd_cs * parameters.p_phonon_qd_cs ) );
        const double gaussian_shape_hole = std::exp( -w * w * parameters.p_phonon_qd_ae / parameters.p_phonon_qd_ratio * parameters.p_phonon_qd_ae / parameters.p_phonon_qd_ratio / ( 4. * parameters.p_phonon_qd_cs * parameters.p_phonon_qd_cs ) );
        const double normalized = ( 4. * 3.1415 * 3.1415 * parameters.p_phonon_qd_rho * std::pow( parameters.p_phonon_qd_cs, 5. ) );
        return w * parameters.hbar * std::pow( parameters.p_phonon_qd_de * gaussian_shape_electron - parameters.p_phonon_qd_dh * gaussian_shape_hole, 2. ) / normalized;
    }
}