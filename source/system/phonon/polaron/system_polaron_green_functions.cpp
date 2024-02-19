#include "system/system.h"

using namespace QDACC;

Scalar System::dgl_phonons_greenf( double tau, const char mode ) {
    if ( not phi_vector.contains( tau ) )
        phi_vector[tau] = dgl_phonons_phi( tau );
    auto phi = phi_vector[tau];
    if ( mode == 'g' ) {
        return parameters.p_phonon_b * parameters.p_phonon_b * ( std::cosh( phi ) - 1.0 );
    }
    return parameters.p_phonon_b * parameters.p_phonon_b * std::sinh( phi );
}

MatrixMain &System::dgl_phonons_greenf_matrix( double tau, const char mode ) {
    if ( operatorMatrices.pme_greenfunction_matrix_cache_g.contains( tau ) ) {
        if ( mode == 'g' )
            return operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
        return operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
    } else {
        operatorMatrices.pme_greenfunction_matrix_cache_g[tau] = operatorMatrices.zero;
        operatorMatrices.pme_greenfunction_matrix_cache_u[tau] = operatorMatrices.zero;
    }
    if ( not phi_vector.contains( tau ) ) {
        phi_vector[tau] = dgl_phonons_phi( tau ); // std::cout << "Time " << t << " is not in phi vecotr\n";
    }
    auto phi = phi_vector[tau];
    auto &g = operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
    auto &u = operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
    //for ( int k = 0; k < operatorMatrices.polaron_phonon_coupling_matrix.outerSize(); ++k ) {
    //    for ( Sparse::InnerIterator it( operatorMatrices.polaron_phonon_coupling_matrix, k ); it; ++it ) {
    //        auto row = it.row();
    //        auto col = it.col();
    //        auto val = it.value();
    for( int row = 0; row < operatorMatrices.polaron_phonon_coupling_matrix.rows(); row++ ) {
        for( int col = 0; col < operatorMatrices.polaron_phonon_coupling_matrix.cols(); col++ ) {
            auto val = operatorMatrices.polaron_phonon_coupling_matrix.coeff( row, col );
            g.coeffRef( row, col ) = ( std::cosh( phi * val ) - 1.0 );
            u.coeffRef( row, col ) = std::sinh( phi * val );
        }
    }
    Log::L3( "[System-PME] Calculated Phi Phonon Matrix cache for tau = {} -> Phi(tau) = {}\n", tau, phi );
    if ( mode == 'g' )
        return operatorMatrices.pme_greenfunction_matrix_cache_g[tau];
    return operatorMatrices.pme_greenfunction_matrix_cache_u[tau];
}