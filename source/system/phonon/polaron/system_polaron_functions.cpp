#include "system/system.h"

using namespace QDLC;

void System::initialize_polaron_frame_functions() {
    if ( parameters.p_phonon_T >= 0 ) {
        Log::L2( "[System-PME] Initializing Polaron Frame Functions.\n" );
        // Initialize Phi(tau)
        double tau = 0.0;
        double last = 1.0;
        double first = std::abs( dgl_phonons_phi( 0.0 ) );
        while ( parameters.p_phonon_tcutoff < 0 ? true : ( tau < parameters.p_phonon_tcutoff ) ) {
            phi_vector[tau] = dgl_phonons_phi( tau );
            dgl_phonons_greenf_matrix( tau, 'g' );
            // FIXME: doesnt really work for stepsizes < 100fs
            if ( parameters.p_phonon_tcutoff < 0 and phi_vector.size() > 1 ) {
                double current = std::real( phi_vector[tau] );
                if ( std::abs( 1.0 - last / current ) < 1E-3 ) {
                    parameters.p_phonon_tcutoff = tau;
                    Log::L2( "[System-PME] Polaron t-cutoff was automatically determined to t_cutoff = {}\n", parameters.p_phonon_tcutoff );
                    break;
                }
                last = current;
            }
            tau += parameters.numerics_subiterator_stepsize; // Use iterator stepsize here to allow for dt = -1 (variable timestep)
        }

        // Output Phonon Functions.
        if ( parameters.output_dict.contains( "greenf" ) ) {
            auto &file = FileOutput::add_file( "phonon_greenfunctions" );
            file << fmt::format( "t\treal(phi(t))\timag(phi(t))\treal(g_u(t))\timag(g_u(t))\treal(g_g(t))\timag(g_g(t))\n" );
            for ( double t = parameters.t_start; t <= parameters.p_phonon_tcutoff; t += parameters.numerics_subiterator_stepsize ) {
                auto greenu = dgl_phonons_greenf( t, 'u' );
                auto greeng = dgl_phonons_greenf( t, 'g' );
                file << fmt::format( "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", t, std::real( phi_vector[t] ), std::imag( phi_vector[t] ), std::real( greenu ), std::imag( greenu ), std::real( greeng ), std::imag( greeng ) );
            }
        }
        if ( parameters.output_dict.contains( "phononJ" ) ) {
            auto &file = FileOutput::add_file( "phonon_spectral" );
            file << fmt::format( "omega\tJ(omega)\n" );
            for ( double w = parameters.p_phonon_wcutoffdelta; w < 10.0 * parameters.p_phonon_wcutoff; w += parameters.p_phonon_wcutoffdelta ) {
                file << fmt::format( "{}\t{}\n", w, std::real( dgl_phonons_spectral_density( w ) ) );
            }
        }
        if ( parameters.output_dict.contains( "phononcoefficients" ) ) {
            // fp_phonons = std::fopen( ( parameters.working_directory + "phonons_lb.txt" ).c_str(), "w" );
            // fmt::print( fp_phonons, "t\tL_a_+\tL_a_-\tL_c_+\tL_c_-\n" );
            // for ( double t = parameters.t_start; t < parameters.t_end; t += parameters.t_step ) {
            //     fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\n", t, dgl_phonons_lindblad_coefficients( t, 'L', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'L', -1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', 1.0 ), dgl_phonons_lindblad_coefficients( t, 'C', -1.0 ) );
            // }
            // std::fclose( fp_phonons );
        }
        Log::L2( "[System-PME] Done.\n" );
    }
}

// TODO: im RK45 fall chis einfach aus t-richtung cachen, für tau richtung dann interpolieren. d.h. wenn ich t einsetze und t zwischen t0,t1 liegt _(t0,t1 schon gespeichert) dann interpoliere einfach linear zwischen den beiden. mal testen.
//  t-richtung sollte eigentlich schon genug infos berechnen für chi. lineare interpolation sollte legit sein. NICHT neu integrieren für tau (auch im RK4 oder RK5 fall nicht, da hat man eh schon alle gecached)
//  im tau-fall eifnach direkt auslesen und ggf interpolieren. dazu interpolate_single funktion schreiben, in solver interpolator auch die nehmen. inlined.
// FIX: einfach immer mit tstep integrieren, auf schon für direction, dann immer interpolieren. lel. die rhos könnte der sich hier auch lokal linear interpolieren für ohne markov
// NEXT: das hier immer nur mit tStepIterator rechnen, dann interpolieren für returnvalue.
/**
 * @brief Calculates the Time Transformed Chi for the Polaron Master Equation. This function either evaluates the PME or returns an interpolated Chi if allowed.
 *
 * @param rho Current Rho Matrix
 * @param t Current Time
 * @param past_rhos All past Rhos. If the general Markov Approximation is disabled, this vector is used instead of the simple Rho
 * @return Sparse
 */

Sparse System::dgl_phonons_pmeq( const Sparse &rho, const double t, const std::vector<QDLC::SaveState> &past_rhos ) {
    track_getcoefficient_calcattempt++;
    // All Contributions will (finally) be reduced onto this return value matrix
    Sparse ret( rho.rows(), rho.cols() );
    // Most precise approximation used. Calculate polaron fram Chi by integrating backwards from t to t-tau.
    if ( parameters.numerics_phonon_approximation_order == QDLC::PhononApproximation::LindbladRates ) {
        return dgl_phonons_lindblad_contribution( t, rho );
    } else {
        // Calculate the initial Chi(t) which is in the interaction frame.
        Sparse chi = dgl_phonons_chi( t ); // dgl_timetrafo( dgl_phonons_chi( t ), t );
        Sparse chi_adjoint = chi.adjoint();
        // Calculate generic X_u and X_g from the initial Chi(t)
        Sparse XGT = chi + chi_adjoint;           // dgl_phonons_chiToX( chi, t, 'u' ); // X_G = Chi + H.c.
        Sparse XUT = 1.i * ( chi - chi_adjoint ); // dgl_phonons_chiToX( chi, t, 'g' ); // X_U = i(Chi - H.c.)
        Log::L3( "[System-PME] Calculating the PME contribution using the Polaron Transformation for t = {}. Initial Time Transformed Chi(t):\n{}\nXG(t):\n{}\nXU(t):\n{}\n", t, Dense( chi ).format( operatorMatrices.output_format ), Dense( XGT ).format( operatorMatrices.output_format ), Dense( XUT ).format( operatorMatrices.output_format ) );
        // Log::L2( "Tau max = {}, stepsize = {}\n", tau_max, parameters.numerics_subiterator_stepsize );
        // Use markov approximation
        if ( parameters.numerics_phonon_approximation_markov1 ) {
            // Temporary variables
            Sparse chi_tau_back_g;
            Sparse chi_tau_back_u;
            // Look if we already calculated coefficient sum for this specific t-value or any value for t > t-value, then interpolate to that specific matrix
            // std::pair<double, double> time_pair = std::make_pair( t, 0.0 );

            // Index was found, chi(t-tau) sum used from saved vector
            if ( parameters.numerics_use_saved_coefficients and savedCoefficients.contains( t ) ) {
                Log::L3( "[System-PME]     Thread #{} - Found {}\n", omp_get_thread_num(), t );
                const auto &coeff = savedCoefficients[t][0.0];
                chi_tau_back_u = coeff.mat1;
                chi_tau_back_g = coeff.mat2;
                track_getcoefficient_read++;
            }
            // Index wasn't found, try to interpolate
            else if ( parameters.numerics_use_saved_coefficients and savedCoefficients.size() > 2 and savedCoefficients.rbegin()->first > t ) {
                Log::L3( "[System-PME]     Element {} should be in here as max is {}, trying to find it...\n", t, savedCoefficients.rbegin()->first );
                // Hacky way to find [min,<t>,max], should work because this is an ordered map.
                const auto &greater_or_equal_than_t = savedCoefficients.lower_bound( t );
                const auto &smaller_than_t = std::prev( greater_or_equal_than_t );
                //  Interpolate
                Log::L3( "[System-PME]     Found Phonon index! Interpolating from {} - {} to {}\n", smaller_than_t->second.begin()->second.t, greater_or_equal_than_t->second.begin()->second.t, t );
                chi_tau_back_u = QDLC::Math::lerp( smaller_than_t->second.begin()->second.mat1, greater_or_equal_than_t->second.begin()->second.mat1, ( t - smaller_than_t->first ) / ( greater_or_equal_than_t->first - smaller_than_t->first ) );
                chi_tau_back_g = QDLC::Math::lerp( smaller_than_t->second.begin()->second.mat2, greater_or_equal_than_t->second.begin()->second.mat2, ( t - smaller_than_t->first ) / ( greater_or_equal_than_t->first - smaller_than_t->first ) );
                track_getcoefficient_read_interpolated++;
            }
            // Not using saved_coefficients or index wasn't found and interpolation wasn't succesfull, recalculating.
            else {
                // Index was not found, (re)calculate chi(t-tau) sum
                // Initialize temporary matrices to zero for threads to write to
                std::vector<Sparse> threadmap_u( parameters.numerics_maximum_secondary_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
                std::vector<Sparse> threadmap_g( parameters.numerics_maximum_secondary_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
                // Calculate backwards integral and sum it into threadmaps. Threadmaps will later be summed into one coefficient matrix.
                int tau_max = std::min<int>( parameters.p_phonon_tcutoff / parameters.numerics_subiterator_stepsize, t / parameters.numerics_subiterator_stepsize );
                Log::L3( "[System-PME]     Thread #{} - (Re)Calculating {} using {} steps in the PME integral\n", omp_get_thread_num(), t, tau_max );
#pragma omp parallel for schedule( dynamic ) num_threads( parameters.numerics_maximum_secondary_threads )
                for ( int tau_index = 0; tau_index < tau_max; tau_index++ ) {
                    double tau = parameters.numerics_subiterator_stepsize * tau_index;
                    Sparse chi_tau_back = dgl_phonons_calculate_transformation( t, tau );
                    Sparse chi_tau_back_adjoint = chi_tau_back.adjoint();
                    auto X_tau_back_g = chi_tau_back + chi_tau_back_adjoint;           // dgl_phonons_chiToX( chi_tau_back, 'g' );
                    auto X_tau_back_u = 1.i * ( chi_tau_back - chi_tau_back_adjoint ); // dgl_phonons_chiToX( chi_tau_back, 'u' );
                    auto thread = omp_get_thread_num();
                    threadmap_g[thread] += X_tau_back_g.cwiseProduct( dgl_phonons_greenf_matrix( tau, 'g' ) ); //*dgl_phonons_greenf( tau, 'g' )
                    threadmap_u[thread] += X_tau_back_u.cwiseProduct( dgl_phonons_greenf_matrix( tau, 'u' ) ); //*dgl_phonons_greenf( tau, 'u' )
                    Log::L3( "[System-PME]         Added Contributions for tau = {}. The sum of the contributions is {}, the phi value is {}\n", tau, threadmap_u[thread].sum() + threadmap_g[thread].sum(), phi_vector[tau] );
                }
// Sum all contributions from threadmaps into one coefficient
#pragma omp parallel sections
                {
#pragma omp section
                    {
                        chi_tau_back_u = std::accumulate( threadmap_u.begin(), threadmap_u.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
                    }
#pragma omp section
                    {
                        chi_tau_back_g = std::accumulate( threadmap_g.begin(), threadmap_g.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
                    }
                }
                // Save coefficients
                if ( parameters.numerics_enable_saving_coefficients ) {
                    savedCoefficients[t][0.0] = QDLC::SaveStateTau( chi_tau_back_u, chi_tau_back_g, t, 0 );
                    track_getcoefficient_write++;
                }
                track_getcoefficient_calculate++;
            }
            // Calculate phonon contributions from (saved/calculated) coefficients and rho(t)
            Sparse integrant = parameters.numerics_subiterator_stepsize * ( dgl_kommutator( XUT, chi_tau_back_u * rho ) + dgl_kommutator( XGT, chi_tau_back_g * rho ) );
            Sparse adjoint = integrant.adjoint();
            ret -= integrant + adjoint;
        } else {
            // TODO: FIXME
            // std::vector<Sparse> threadmap_u( parameters.numerics_maximum_secondary_threads, Sparse( parameters.maxStates, parameters.maxStates ) );
            //// Dont use markov approximation; Full integral
            // auto chis_transformed = dgl_phonons_calculate_transformation( chi, t, parameters.t_step * tau_max );
            //#pragma omp parallel for ordered schedule( dynamic ) shared( savedCoefficients ) num_threads( parameters.numerics_maximum_secondary_threads )
            // for ( int tau_index = 0; tau_index < tau_max; tau_index++ ) {
            //     int rho_index = std::max( 0, (int)past_rhos.size() - 1 - tau_index ); // TODO: Das hier ist mit var timesteps der rhos auch grütze!
            //     double tau = ( 1.0 * tau_index ) * parameters.t_step;
            //
            //    // Index not found, or no saving is used. Recalculate chi(t-tau).
            //    auto &chi_tau_back = chis_transformed.at( std::min<int>( chis_transformed.size() - 1, tau_index ) ).mat;
            //    Sparse chi_tau_back_adjoint = chi_tau_back.adjoint();
            //    auto chi_tau_back_g = chi_tau_back + chi_tau_back_adjoint;           // dgl_phonons_chiToX( chi_tau_back, 'u' );
            //    auto chi_tau_back_u = 1.i * ( chi_tau_back - chi_tau_back_adjoint ); // dgl_phonons_chiToX( chi_tau_back, 'g' );
            //    // If saving is used, save current chi(t-tau).
            //    if ( parameters.numerics_use_saved_coefficients or parameters.numerics_enable_saving_coefficients ) {
            //#pragma omp critical
            //        savedCoefficients[t][tau] = QDLC::SaveStateTau( chi_tau_back_u, chi_tau_back_g, t, tau );
            //    }
            //    Sparse integrant = dgl_phonons_greenf_matrix( tau, 'u' ).cwiseProduct( dgl_kommutator( XUT, chi_tau_back_u * rho ) ); // FIXME: mit rk45 oder var gitter passt das rho hier nicht mehr! --> rho interpolieren!
            //    integrant += dgl_phonons_greenf_matrix( tau, 'g' ).cwiseProduct( dgl_kommutator( XGT, chi_tau_back_g * rho ) );       // FIXME: mit rk45 oder var gitter passt das rho hier nicht mehr! --> rho interpolieren!
            //    Sparse adjoint = integrant.adjoint();
            //    auto thread = omp_get_thread_num();
            //    threadmap_u[thread] += ( integrant + adjoint ) * parameters.t_step;
            //}
            // ret -= std::accumulate( threadmap_u.begin(), threadmap_u.end(), Sparse( parameters.maxStates, parameters.maxStates ) );
        }
        Log::L3( "[System-PME]     Thread #{} - Adding ret value for {}. Matrix is: \n{}\n", omp_get_thread_num(), t, Dense( ret ).format( operatorMatrices.output_format ) );
    }
    Log::L3( "[System-PME] Thread #{} - Returning PME Contribution for {}:\n{}\n", omp_get_thread_num(), t, Dense( ret ).format( operatorMatrices.output_format ) );
    return ret;
}