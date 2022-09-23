#include "system/system.h"

Scalar System::dgl_phonons_greenf( double tau, const char mode ) {
    if ( not phi_vector.contains( tau ) )
        phi_vector[tau] = dgl_phonons_phi( tau );
    auto phi = phi_vector[tau];
    if ( mode == 'g' ) {
        return parameters.p_phonon_b * parameters.p_phonon_b * ( std::cosh( phi ) - 1.0 );
    }
    return parameters.p_phonon_b * parameters.p_phonon_b * std::sinh( phi );
}

Sparse &System::dgl_phonons_greenf_matrix( double tau, const char mode ) {
    if ( operatorMatrices.pme_greenfunction_matrix_cache_g.contains( tau ) ) {
        if ( mode == 'g' )
            return operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
        return operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
    } else {
        operatorMatrices.pme_greenfunction_matrix_cache_g[tau] = Sparse( operatorMatrices.polaron_phonon_coupling_matrix.cols(), operatorMatrices.polaron_phonon_coupling_matrix.rows() );
        operatorMatrices.pme_greenfunction_matrix_cache_u[tau] = Sparse( operatorMatrices.polaron_phonon_coupling_matrix.cols(), operatorMatrices.polaron_phonon_coupling_matrix.rows() );
    }
    if ( not phi_vector.contains( tau ) ) {
        phi_vector[tau] = dgl_phonons_phi( tau ); // std::cout << "Time " << t << " is not in phi vecotr\n";
    }
    auto phi = phi_vector[tau];
    auto &g = operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
    auto &u = operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
    for ( int k = 0; k < operatorMatrices.polaron_phonon_coupling_matrix.outerSize(); ++k ) {
        for ( Sparse::InnerIterator it( operatorMatrices.polaron_phonon_coupling_matrix, k ); it; ++it ) {
            auto row = it.row();
            auto col = it.col();
            auto val = it.value();
            g.coeffRef( row, col ) = ( std::cosh( phi * val ) - 1.0 );
            u.coeffRef( row, col ) = std::sinh( phi * val );
        }
    }
    Log::L3( "[System-PME] Calculated Phi Phonon Matrix cache for tau = {} -> Phi(tau) = {}\n", tau, phi );
    if ( mode == 'g' )
        return operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
    return operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
}