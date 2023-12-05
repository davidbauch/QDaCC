#include "system/system.h"

using namespace QDACC;

void System::initialize_polaron_frame_functions() {
    if ( parameters.p_phonon_T >= 0 ) {
        Log::L2( "[System-PME] Initializing Polaron Frame Functions.\n" );
        // Initialize Phi(tau)
        double tau = 0.0;
        double last = 1.0;
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
            file << std::format( "t\treal(phi(t))\timag(phi(t))\treal(g_u(t))\timag(g_u(t))\treal(g_g(t))\timag(g_g(t))\n" );
            for ( double t = parameters.t_start; t <= parameters.p_phonon_tcutoff; t += parameters.numerics_subiterator_stepsize ) {
                auto greenu = dgl_phonons_greenf( t, 'u' );
                auto greeng = dgl_phonons_greenf( t, 'g' );
                file << std::format( "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", t, std::real( phi_vector[t] ), std::imag( phi_vector[t] ), std::real( greenu ), std::imag( greenu ), std::real( greeng ), std::imag( greeng ) );
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

Sparse System::dgl_phonons_pmeq( const Sparse &rho, const double t, const std::vector<QDACC::SaveState> &past_rhos ) {
    track_getcoefficient_calcattempt++;
    if ( parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::LindbladRates ) {
        return dgl_phonons_lindblad_contribution( t, rho );
    }
        return dgl_phonons_integrated_contribution( t, rho, past_rhos );
}