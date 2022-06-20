#include "solver/solver_ode.h"

bool QDLC::Numerics::ODESolver::calculate_wigner( System &s, const std::string &s_mode, const double x, const double y, const int resolution, const int skips ) {
    // Partially trace the mode to wigner:
    std::vector<Dense> reduced_rho;
    std::vector<Scalar> time;
    int base = s.operatorMatrices.el_states.count( s_mode ) != 0 ? s.operatorMatrices.el_states[s_mode].base : s.operatorMatrices.ph_states[s_mode].base;
    LOG2( "[Wigner] Calculating Wigner function for mode {}/{}\n", s_mode, base );
    // The Wigner function uses the interpolated savedStates, not the actually calculated states from the RK45 method.
    // TODO: interaction pic sollte man hier ausstellen können. am besten für normale rho ausgabe auch.
    bool int_trafo = s.parameters.input_conf["DMconfig"].string["interaction_picture"] == "int";
    LOG2( "[Wigner] Output Matrix format will be {}\n", int_trafo ? "in the interaction frame." : "in the Schrödinger frame." );
    for ( int i = 0; i < savedStates.size(); i += skips ) {
        reduced_rho.emplace_back( s.partial_trace( int_trafo ? get_rho_at( i ) : s.dgl_timetrafo( get_rho_at( i ), get_time_at( i ) ), base ) );
        time.emplace_back( get_time_at( i ) );
    }

    double g = std::sqrt( 2 );
    // Calculate Wigner functions over time:
    // Calculate Meshgrid
    auto [X_mat, Y_mat] = QDLC::Matrix::meshgrid( -x, -y, x, y, resolution );
    // Dense A = 0.5 * g * ( X_mat + 1.0i * Y_mat );
    Dense A = ( X_mat + 1.0i * Y_mat );

    std::vector<Dense> wigner( reduced_rho.size(), Dense::Zero( A.rows(), A.cols() ) );
    LOG2( "[Wigner] First Matrix for Wigner:\n{}\n", reduced_rho.front() );

    ProgressBar progressbar = ProgressBar();
    Timer &timer_w = Timers::create( "Wigner (" + s_mode + ")" );
    timer_w.start();
    // Iterate
    // Calculates the Wigner function at coordinate alpha
    auto Wigner = [&]( Dense &DM, Scalar alpha ) {
        double Wval = 0.0;
        for ( int m = 0; m < DM.rows(); m++ ) {
            if ( QDLC::Math::abs2( DM( m, m ) ) > 0.0 ) {
                Wval += std::real( DM( m, m ) * std::pow( -1., (double)m ) * std::assoc_laguerre( m, 0, 4. * QDLC::Math::abs2( alpha ) ) );
            }
        }

        for ( int m = 0; m < DM.rows(); m++ ) {
            for ( int n = m + 1; n < DM.rows(); n++ ) {
                if ( QDLC::Math::abs2( DM( m, n ) ) > 0.0 ) {
                    Wval += 2. * std::real( DM( m, n ) * std::pow( -1., (double)m ) * std::pow( 2. * alpha, (double)( n - m ) ) * std::sqrt( QDLC::Math::factorial( m ) / QDLC::Math::factorial( n ) ) * std::assoc_laguerre( m, n - m, 4. * QDLC::Math::abs2( alpha ) ) );
                }
            }
        }
        return ( 1 / ( 2. * 3.1415 ) * std::exp( -2. * QDLC::Math::abs2( alpha ) ) * Wval );
    };
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < reduced_rho.size(); i++ ) {
        Dense W = Dense::Zero( A.rows(), A.cols() );
        for ( int k = 0; k < A.rows(); k++ ) {
            for ( int l = 0; l < A.cols(); l++ ) {
                W( k, l ) = Wigner( reduced_rho[i], A( k, l ) );
            }
        }
        wigner.at( i ) = W;
        timer_w.iterate();
        Timers::outputProgress( timer_w, progressbar, timer_w.getTotalIterationNumber(), reduced_rho.size(), "Wigner (" + s_mode + "): " );
    }
    timer_w.end();
    Timers::outputProgress( timer_w, progressbar, timer_w.getTotalIterationNumber(), reduced_rho.size(), "Wigner (" + s_mode + "): ", Timers::PROGRESS_FORCE_OUTPUT );
    // Add to Fileoutput:
    if ( to_output["Wigner"].size() == 0 )
        to_output["Wigner"]["Time"] = time;
    to_output_m["Wigner"][s_mode] = wigner;
    to_output_m["Wigner"]["rho_" + s_mode] = reduced_rho;
    return true;
}