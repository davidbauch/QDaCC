#include "system/system.h"

using namespace QDLC;

/**
 * @brief Calculate the X_u and X_g matrices from the initial Chi(t) matrix.
 * using the formula X_g = Chi(t) + Chi(t)^dagger and X_u = i(Chi(t) - Chi(t)^dagger)
 * @param chi The initial Chi(t) matrix in the interaction frame.
 * @param chi_adjoint The adjoint of the initial Chi(t) matrix in the interaction frame.
 * @return A tuple containing the X_g and X_u matrices.
 */
static inline std::tuple<Sparse, Sparse> _chi_to_XG( const auto& chi, const auto& chi_adjoint ) {
    return std::make_tuple( chi + chi_adjoint /* g */, 1.i * ( chi - chi_adjoint ) /* u */ );
}

static inline int _get_tau_iterations( const double cutoff_time, const double iterator_stepsize, const double current_time ) {
    return std::min<int>( cutoff_time / iterator_stepsize, current_time / iterator_stepsize );
}

/**
 * @brief Get the cached coefficient matrix for the given time.
 * @param t The time for which the coefficient matrix should be returned.
 * @param saved_coefficients The map containing the saved coefficient matrices.
 * @return A tuple containing the X_g and X_u matrices.
 */
static inline std::tuple<Sparse, Sparse> _get_cached_coefficient( const double t, const auto& saved_coefficients ) {
    const auto& coeff = saved_coefficients.at( t ).at( 0.0 );
    return std::make_tuple( coeff.mat1 /* g */, coeff.mat2 /* u */ );
}

/**
 * @brief Interpolate between two cached coefficient matrices.
 * @param t The time for which the coefficient matrix should be returned.
 * @param saved_coefficients The map containing the saved coefficient matrices.
 * @return A tuple containing the X_g and X_u matrices.
 */
static inline std::tuple<Sparse, Sparse> _interpolate_cached_coefficient( const double t, const auto& saved_coefficients ) {
    // Hacky way to find [min,<t>,max], should work because this is an ordered map.
    const auto& greater_or_equal_than_t = saved_coefficients.lower_bound( t );
    const auto& smaller_than_t = std::prev( greater_or_equal_than_t );
    //  Interpolate
    const auto& XGTau = QDLC::Math::lerp( smaller_than_t->second.begin()->second.mat1, greater_or_equal_than_t->second.begin()->second.mat1, ( t - smaller_than_t->first ) / ( greater_or_equal_than_t->first - smaller_than_t->first ) );
    const auto& XUTau = QDLC::Math::lerp( smaller_than_t->second.begin()->second.mat2, greater_or_equal_than_t->second.begin()->second.mat2, ( t - smaller_than_t->first ) / ( greater_or_equal_than_t->first - smaller_than_t->first ) );
    Log::L3( "Returning interpolate coefficient for t = {} using t0 = {} and t1 = {}, where t1-t0 = {}\n", t, smaller_than_t->first * 1E12, greater_or_equal_than_t->first * 1E12, ( -smaller_than_t->first + greater_or_equal_than_t->first ) * 1E12 );
    return std::make_tuple( XGTau, XUTau );
}

static inline std::tuple<Sparse, Sparse> _get_thread_reduced_coefficients( const auto& threadmap_g, const auto& threadmap_u ) {
    return std::make_tuple(
        std::accumulate( threadmap_g.begin(), threadmap_g.end(), Sparse( threadmap_g.front().rows(), threadmap_g.front().cols() ) ),
        std::accumulate( threadmap_u.begin(), threadmap_u.end(), Sparse( threadmap_u.front().rows(), threadmap_u.front().cols() ) ) );
}

Sparse System::dgl_phonons_integrated_contribution( const double t, const Sparse& rho, const std::vector<QDLC::SaveState>& past_rhos ) {
    // Interger representation of the cutoff time
    int tau_max = _get_tau_iterations( parameters.p_phonon_tcutoff, parameters.numerics_subiterator_stepsize, t );
    // Phonon CPU Cores. Usually more than 8 threads are not beneficial.
    const auto phonon_iterator_threads = std::min<int>( parameters.numerics_maximum_secondary_threads, 16 );
    // Final Integrant
    Sparse integrant = Sparse( rho.rows(), rho.cols() );
    // Calculate the initial Chi(t) which is in the interaction frame and its adjoint.
    Sparse chi = dgl_phonons_chi( t );
    Sparse chi_adjoint = chi.adjoint();
    // Calculate generic X_u and X_g from the initial Chi(t)
    const auto& [XGT, XUT] = _chi_to_XG( chi, chi_adjoint );
    // Use markov approximation
    if ( parameters.numerics_phonon_approximation_markov1 ) {
        Sparse XGTau;
        Sparse XUTau;
        // Final Values for X_u and X_g
        // If we have saved coefficients, try to find the one for the current time, or interpolate between two saved ones.
        if ( parameters.numerics_use_saved_coefficients and savedCoefficients.contains( t ) ) {
            std::tie( XGTau, XUTau ) = _get_cached_coefficient( t, savedCoefficients );
            track_getcoefficient_read++;
        } else if ( parameters.numerics_use_saved_coefficients and savedCoefficients.size() > 2 and savedCoefficients.rbegin()->first > t ) {
            std::tie( XGTau, XUTau ) = _interpolate_cached_coefficient( t, savedCoefficients );
            track_getcoefficient_read_interpolated++;
        }
        // Not using saved_coefficients or index wasn't found or interpolation wasn't succesfull, recalculating.
        else {
            // Initialize temporary matrices to zero for threads to write to
            std::vector<Sparse> threadmap_g( phonon_iterator_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
            std::vector<Sparse> threadmap_u( phonon_iterator_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
            
            // Calculate backwards integral and sum it into threadmaps. Threadmaps will later be summed into one coefficient matrix.
#pragma omp parallel for schedule( dynamic ) num_threads( phonon_iterator_threads )
            for ( int tau_index = 0; tau_index < tau_max; tau_index++ ) {
                double tau = parameters.numerics_subiterator_stepsize * tau_index;
                Sparse chi_tau_back = dgl_phonons_calculate_transformation( t, tau );
                Sparse chi_tau_back_adjoint = chi_tau_back.adjoint();
                const auto& [x_tau_back_g, x_tau_back_u] = _chi_to_XG( chi_tau_back, chi_tau_back_adjoint );
                const auto thread = omp_get_thread_num();
                threadmap_g[thread] += x_tau_back_g.cwiseProduct( dgl_phonons_greenf_matrix( tau, 'g' ) );
                threadmap_u[thread] += x_tau_back_u.cwiseProduct( dgl_phonons_greenf_matrix( tau, 'u' ) );
            }
            // Sum all contributions from threadmaps into one coefficient
            std::tie( XGTau, XUTau ) = _get_thread_reduced_coefficients( threadmap_g, threadmap_u );
        }
        // Save coefficients
        if ( parameters.numerics_enable_saving_coefficients ) {
            savedCoefficients[t][0.0] = QDLC::SaveStateTau( XGTau, XUTau, t, 0 );
            track_getcoefficient_write++;
        }
        track_getcoefficient_calculate++;
        integrant = parameters.numerics_subiterator_stepsize * ( dgl_kommutator( XGT, XGTau * rho ) + dgl_kommutator( XUT, XUTau * rho ) );    
    }
    // Calculate phonon contributions from (saved/calculated) coefficients and rho(t)
    Sparse adjoint = integrant.adjoint();
    return -( integrant + adjoint );
}

//} else {
// TODO: FIXME
// std::vector<Sparse> threadmap_u( parameters.numerics_maximum_secondary_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
//// Dont use markov approximation; Full integral
// auto chis_transformed = dgl_phonons_calculate_transformation( chi, t, parameters.t_step * tau_max );
// #pragma omp parallel for ordered schedule( dynamic ) shared( savedCoefficients ) num_threads( parameters.numerics_maximum_secondary_threads )
// for ( int tau_index = 0; tau_index < tau_max; tau_index++ ) {
//     int rho_index = std::max( 0, (int)past_rhos.size() - 1 - tau_index ); // TODO: Das hier ist mit var timesteps der rhos auch grütze!
//     double tau = ( 1.0 * tau_index ) * parameters.t_step;
//
//    // Index not found, or no saving is used. Recalculate chi(t-tau).
//    auto &chi_tau_back = chis_transformed.at( std::min<int>( chis_transformed.size() - 1, tau_index ) ).mat;
//    Sparse chi_tau_back_adjoint = chi_tau_back.adjoint();
//    auto XGTau = chi_tau_back + chi_tau_back_adjoint;           // dgl_phonons_chiToX( chi_tau_back, 'u' );
//    auto XUTau = 1.i * ( chi_tau_back - chi_tau_back_adjoint ); // dgl_phonons_chiToX( chi_tau_back, 'g' );
//    // If saving is used, save current chi(t-tau).
//    if ( parameters.numerics_use_saved_coefficients or parameters.numerics_enable_saving_coefficients ) {
// #pragma omp critical
//        savedCoefficients[t][tau] = QDLC::SaveStateTau( XUTau, XGTau, t, tau );
//    }
//    Sparse integrant = dgl_phonons_greenf_matrix( tau, 'u' ).cwiseProduct( dgl_kommutator( XUT, XUTau * rho ) ); // FIXME: mit rk45 oder var gitter passt das rho hier nicht mehr! --> rho interpolieren!
//    integrant += dgl_phonons_greenf_matrix( tau, 'g' ).cwiseProduct( dgl_kommutator( XGT, XGTau * rho ) );       // FIXME: mit rk45 oder var gitter passt das rho hier nicht mehr! --> rho interpolieren!
//    Sparse adjoint = integrant.adjoint();
//    auto thread = omp_get_thread_num();
//    threadmap_u[thread] += ( integrant + adjoint ) * parameters.t_step;
//}
// ret -= std::accumulate( threadmap_u.begin(), threadmap_u.end(), Sparse( parameters.maxStates, parameters.maxStates ) );