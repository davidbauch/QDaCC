#include "system/system.h"

Sparse System::dgl_phonons_calculate_transformation( double t, double tau ) {
    // Backwards Integral
    if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_BACKWARDS_INTEGRAL ) {
        // TODO
        // return QDLC::Numerics::calculate_definite_integral_vec( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, std::get<1>( parameters.numerics_rk_tol.front() ), parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() );
        auto chi_tau = dgl_phonons_chi( t - tau );
        //auto func = [this](const Sparse &chi, const double t){return this->dgl_phonons_rungefunc(chi,t);};
        return QDLC::Numerics::calculate_definite_integral( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, std::get<1>( parameters.numerics_rk_tol.front() ), parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() ).mat;
    }
    // Matrix Exponential
    else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_TRANSFORMATION_MATRIX ) {
        auto chi_tau = dgl_phonons_chi( t - tau ); // TODO maybe cache these.
        Sparse U = ( Dense( -1.0i * dgl_get_hamilton( t ) * tau ).exp() ).sparseView();
        return ( U * chi_tau * U.adjoint() );
    }
    return dgl_phonons_chi( t - tau );
    // Backwards Integral if Threshold is met, else neglect transformation
    // else if ( parameters.numerics_phonon_approximation_order == PHONON_APPROXIMATION_MIXED ) {
    //    double threshold = 0;
    //    for ( auto &p : pulse )
    //        threshold += std::abs( p.get( t ) / p.maximum );                                     // TODO: + photon zahl, wenn photon number > 0.1 oder so dann auch. für große kopplungen gibts sonst starke abweichungen. vil. number*g draufaddieren.
    //    if ( threshold > 1E-4 || ( not chirp.empty() and chirp.back().derivative( t ) != 0 ) ) { // TODO: threshold als parameter
    //        return QDLC::Numerics::calculate_definite_integral_vec( chi_tau, std::bind( &System::dgl_phonons_rungefunc, this, std::placeholders::_1, std::placeholders::_2 ), t, std::max( t - tau, 0.0 ), parameters.numerics_subiterator_stepsize, std::get<1>( parameters.numerics_rk_tol.front() ), parameters.numerics_rk_stepmin, parameters.numerics_rk_stepmax, parameters.numerics_rk_usediscrete_timesteps ? parameters.numerics_rk_stepdelta.get() : 0.0, parameters.numerics_phonon_nork45 ? 4 : parameters.numerics_rk_order.get() );
    //    }
    //}
}