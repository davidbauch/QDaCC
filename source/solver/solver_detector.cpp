#include "solver/solver_ode.h"

void QDLC::Numerics::ODESolver::apply_detector_function( System &s, Dense &mat, const Dense &timemat, const std::string &purpose ) {
    if ( s.parameters.input_conf["Detector"].numerical_v["time_center"].size() > 0 ) {
        if ( detector_temporal_mask.rows() == 0 ) {
            detector_temporal_mask = Dense( mat.rows(), mat.cols() );
            for ( int c = 0; c < s.parameters.input_conf["Detector"].numerical_v["time_center"].size(); c++ ) {
                double detector_t_center = s.parameters.input_conf["Detector"].numerical_v["time_center"][c];
                double detector_t_range = s.parameters.input_conf["Detector"].numerical_v["time_range"][c];
                double detector_power = s.parameters.input_conf["Detector"].numerical_v["power_amplitude"][c];
                for ( int i = 0; i < mat.rows(); i++ ) {
                    double time = std::real( timemat( i, 0 ) );
                    for ( int j = 0; j < mat.cols() - i; j++ ) {
                        double tau = std::imag( timemat( i, j ) );
                        detector_temporal_mask( i, j ) += std::exp( -std::pow( ( time - detector_t_center ) / detector_t_range, detector_power ) ) * std::exp( -std::pow( ( tau - detector_t_center ) / detector_t_range, detector_power ) );
                    }
                }
            }
        }
        mat = mat.cwiseProduct( detector_temporal_mask );
    }
    if ( s.parameters.input_conf["Detector"].numerical_v["spectral_range"].size() > 0 ) {
        Timer &timer = Timers::create( "Detector (" + purpose + ")" );
        ProgressBar progressbar = ProgressBar();
        timer.start();
        // Disgusting piece of code
        if ( detector_frequency_mask.size() == 0 ) {
            Log::L2( "[PhotonStatistics] Calculating Detector Spectral Mask... \n" );
            for ( int c = 0; c < s.parameters.input_conf["Detector"].numerical_v["spectral_range"].size(); c++ ) {
                double center = s.parameters.input_conf["Detector"].numerical_v["spectral_center"][c];
                double range = s.parameters.input_conf["Detector"].numerical_v["spectral_range"][c];
                double points = s.parameters.input_conf["Detector"].numerical_v["spectral_number_points"][c]; // Points-per-sigma
                double power = s.parameters.input_conf["Detector"].numerical_v["spectral_power_amplitude"][c];
                double delta_omega = range / points;
                double w = 0;
                while ( true ) {
                    double amplitude = std::exp( -0.5 * std::pow( w / points, power ) );
                    detector_frequency_mask.emplace_back( std::make_tuple( center - w / points * range, amplitude, delta_omega ) );
                    detector_frequency_mask.emplace_back( std::make_tuple( center + w / points * range, amplitude, delta_omega ) );
                    w += 1;
                    if ( amplitude < 1E-8 )
                        break;
                }
            }
            std::sort( detector_frequency_mask.begin(), detector_frequency_mask.end(), []( std::tuple<double, double, double> &first, std::tuple<double, double, double> &second ) { return std::get<0>( first ) > std::get<0>( second ); } );
            Log::L2( "[PhotonStatistics] Calculating Detector Spectral Mask done, Mask size is {} \n", detector_frequency_mask.size() );
        }
        Log::L2( "[PhotonStatistics] Applying Detector Mask for {}\n", purpose );
        // Calculate main fourier transform integral with spectral amplitude in tau direction.
        Dense mat_transformed = Dense( mat.rows(), detector_frequency_mask.size() );
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
        for ( long unsigned int i = 0; i < mat.rows(); i++ ) {
            for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
                auto &[frequency_w, frequency_amp, frequency_delta] = detector_frequency_mask[w];
                for ( int j = 0; j < mat.cols() - i; j++ ) { // -i because of triangular grid
                    double tau = std::imag( timemat( i, j ) ) - std::real( timemat( i, j ) );
                    double dtau = Numerics::get_taudelta( timemat, i, j );
                    mat_transformed( i, w ) += frequency_amp * std::exp( -1.0i * frequency_w * tau ) * mat( i, j ) * dtau / 2.0 / 3.1415;
                }
            }
            Timers::outputProgress( timer, progressbar, i, mat.rows(), "Detector FT (" + purpose + "): " );
            timer.iterate();
        }
        // Output. TODO: only if wanted
        /*
        FILE *f_gfunc = std::fopen( ( s.parameters.working_directory + purpose + "_mFT.txt" ).c_str(), "w" );
        fmt::print( f_gfunc, "Time\tTau\tAbs\tReal\tImag\n" );
        for ( int k = 0; k < mat_transformed.rows(); k++ ) {
            for ( int l = 0; l < mat_transformed.cols(); l++ ) {
                fmt::print( f_gfunc, "{:d}\t{:d}\t{:.8e}\t{:.8e}\t{:.8e}\n", k, l, std::abs( mat_transformed( k, l ) ), std::real( mat_transformed( k, l ) ), std::imag( mat_transformed( k, l ) ) );
            }
            fmt::print( f_gfunc, "\n" );
        }
        std::fclose( f_gfunc );
        */
        // Transform Back in tau direction.
        mat = Dense::Zero( mat.rows(), mat.cols() );
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
        for ( long unsigned int i = 0; i < mat.rows(); i++ ) {
            for ( int j = 0; j < mat.cols() - i; j++ ) { // -i because of triangular grid
                double t = std::real( timemat( i, j ) );
                double tau = std::imag( timemat( i, j ) ) - t;
                for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
                    auto &[frequency_w, frequency_amp, frequency_delta] = detector_frequency_mask[w];
                    mat( i, j ) += std::exp( 1.0i * tau * frequency_w ) * mat_transformed( i, w ) * frequency_delta;
                }
            }
            Timers::outputProgress( timer, progressbar, i, mat.rows(), "Detector iFT (" + purpose + "): " );
            timer.iterate();
        }

        // Do everythogn again in t-direction
        /*
        Log::L2( "[PhotonStatistics] Applying Detector Mask for t direction {}\n", purpose );
        // Calculate main fourier transform integral with spectral amplitude in t direction.
        mat_transformed = Dense( mat.rows(), detector_frequency_mask.size() );
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
        for ( long unsigned int j = 0; j < mat.rows(); j++ ) {
            for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
                auto &[frequency_w, frequency_amp, frequency_delta] = detector_frequency_mask[w];
                for ( int i = 0; i < mat.cols() - j; i++ ) { // -i because of triangular grid
                    double t = std::real( timemat( i, j ) );
                    double dt = Numerics::get_tdelta( timemat, i, j );
                    mat_transformed( j, w ) += frequency_amp * std::exp( -1.0i * frequency_w * t ) * mat( i, j ) * dt / 2.0 / 3.1415;
                }
            }
            Timers::outputProgress( timer, progressbar, j, mat.rows(), "Detector FT (" + purpose + "): " );
            timer.iterate();
        }
        // Output. TODO: only if wanted

        // FILE *f_gfunc = std::fopen( ( s.parameters.working_directory + purpose + "_mFT.txt" ).c_str(), "w" );
        // fmt::print( f_gfunc, "Time\tTau\tAbs\tReal\tImag\n" );
        // for ( int k = 0; k < mat_transformed.rows(); k++ ) {
        //     for ( int l = 0; l < mat_transformed.cols(); l++ ) {
        //         fmt::print( f_gfunc, "{:d}\t{:d}\t{:.8e}\t{:.8e}\t{:.8e}\n", k, l, std::abs( mat_transformed( k, l ) ), std::real( mat_transformed( k, l ) ), std::imag( mat_transformed( k, l ) ) );
        //     }
        //     fmt::print( f_gfunc, "\n" );
        // }
        // std::fclose( f_gfunc );

        // Transform Back in t direction.
        mat = Dense::Zero( mat.rows(), mat.cols() );
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
        for ( long unsigned int j = 0; j < mat.rows(); j++ ) {
            for ( int i = 0; i < mat.cols() - j; i++ ) { // -i because of triangular grid
                double t = std::real( timemat( i, j ) );
                for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
                    auto &[frequency_w, frequency_amp, frequency_delta] = detector_frequency_mask[w];
                    mat( i, j ) += std::exp( 1.0i * t * frequency_w ) * mat_transformed( j, w ) * frequency_delta;
                }
            }
            Timers::outputProgress( timer, progressbar, j, mat.rows(), "Detector iFT (" + purpose + "): " );
            timer.iterate();
        }
*/
        timer.end();
        Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), timer.getTotalIterationNumber(), "Detector (" + purpose + ")", Timers::PROGRESS_FORCE_OUTPUT );
    }
}
