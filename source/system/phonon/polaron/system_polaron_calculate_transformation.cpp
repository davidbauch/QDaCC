#include "system/system.h"

using namespace QDACC;

MatrixMain System::dgl_phonons_calculate_transformation( double t, double tau ) {
    // Backwards Integral
    if ( parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::BackwardsIntegral ) {
        // TODO
        // return QDACC::Numerics::calculate_definite_integral_vec( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, std::get<1>( parameters.numerics_rk_tol.front() ), parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() );
        auto chi = dgl_phonons_chi( t - tau );
        // auto func = [this](const MatrixMain &chi, const double t){return this->dgl_phonons_rungefunc(chi,t);};
        return QDACC::Numerics::calculate_definite_integral( chi, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, std::get<1>( parameters.numerics_rk_tol.front() ), parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() ).mat;
    }
    // Matrix Exponential
    else if ( parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::TransformationMatrix ) {
        auto chi_tau = dgl_phonons_chi( t - tau ); // TODO maybe cache these.
        #ifdef USE_SPARSE_MATRIX
        MatrixMain U = ( Dense( -1.0i * dgl_get_hamilton( t ) * tau ).exp() ).sparseView();
        #else
        MatrixMain U = Dense( -1.0i * dgl_get_hamilton( t ) * tau ).exp();
        #endif
        return ( U * chi_tau * U.adjoint() );
    }
    return dgl_phonons_chi( t - tau );
    // Backwards Integral if Threshold is met, else neglect transformation
    // else if ( parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::Mixed ) {
    //    double threshold = 0;
    //    for ( auto &p : pulse )
    //        threshold += std::abs( p.get( t ) / p.maximum );                                     // TODO: + photon zahl, wenn photon number > 0.1 oder so dann auch. für große kopplungen gibts sonst starke abweichungen. vil. number*g draufaddieren.
    //    if ( threshold > 1E-4 || ( not chirp.empty() and chirp.back().derivative( t ) != 0 ) ) { // TODO: threshold als parameter
    //        return QDACC::Numerics::calculate_definite_integral_vec( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, std::get<1>( parameters.numerics_rk_tol.front() ), parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() );
    //    }
    //}
}