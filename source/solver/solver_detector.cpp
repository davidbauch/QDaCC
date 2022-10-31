#include "solver/solver_ode.h"

void QDLC::Numerics::ODESolver::apply_detector_function( System &s, CacheMatrix &mat ) {
    const auto &purpose = mat.get_name();
    const auto dim = mat.dim();
    if ( s.parameters.input_conf["Detector"].numerical_v["time_center"].size() > 0 ) {
        if ( detector_temporal_mask.rows() == 0 ) {
            Log::L2( "[PhotonStatistics] Calculating Detector Temporal Mask... \n" );
            detector_temporal_mask = mat.empty();
            for ( int c = 0; c < s.parameters.input_conf["Detector"].numerical_v["time_center"].size(); c++ ) {
                double detector_t_center = s.parameters.input_conf["Detector"].numerical_v["time_center"][c];
                double detector_t_range = s.parameters.input_conf["Detector"].numerical_v["time_range"][c];
                double detector_power = s.parameters.input_conf["Detector"].numerical_v["time_power_amplitude"][c];
                Log::L2( "[PhotonStatistics] Adding Temporal Mask for t_0 = {}, delta = {}, amp = {}\n", detector_t_center, detector_t_range, detector_power );
                for ( int i = 0; i < dim; i++ ) {
                    double time = mat.t( i );
                    for ( int j = 0; j < dim - i; j++ ) {
                        double tau = mat.tau( j, i );
                        detector_temporal_mask( i, j ) += std::exp( -std::pow( ( time - detector_t_center ) / detector_t_range, detector_power ) ) * std::exp( -std::pow( ( tau - detector_t_center ) / detector_t_range, detector_power ) );
                    }
                }
            }
            Log::L2( "[PhotonStatistics] Saving Detector Matrix to detector_temporal_mask.txt...\n" );
            auto &f_detector = FileOutput::add_file( "detector_temporal_mask" );
            f_detector << "Time_index\ttau_index\tD(t)*D(t+tau)\n";
            for ( int k = 0; k < detector_temporal_mask.rows(); k++ ) {
                const auto t = mat.t( k );
                for ( int l = 0; l < detector_temporal_mask.cols(); l++ ) {
                    const auto tau = mat.tau( l );
                    const auto val = std::real( detector_temporal_mask( k, l ) );
                    f_detector << fmt::format( "{}\t{}\t{:.8e}\n", t, tau, val );
                }
                f_detector << "\n";
            }
        }
        Log::L2( "[PhotonStatistics] Applying Detector Temporal Mask to {}... \n", purpose );
        mat.get() = mat.get().cwiseProduct( detector_temporal_mask );
    }
    if ( s.parameters.input_conf["Detector"].numerical_v["spectral_range"].size() > 0 ) {
        Timer &timer = Timers::create( "Detector (" + purpose + ")" ).start();
        auto progressbar = ProgressBar();
        // Disgusting piece of code
        if ( detector_frequency_mask.empty() ) {
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
                    if ( w != 0.0 )
                        detector_frequency_mask.emplace_back( center - w / points * range, amplitude, delta_omega );
                    detector_frequency_mask.emplace_back( center + w / points * range, amplitude, delta_omega );
                    w += 1;
                    if ( amplitude < 1E-4 )
                        break;
                }
            }
            std::ranges::sort( detector_frequency_mask.begin(), detector_frequency_mask.end(), []( std::tuple<double, double, double> &first, std::tuple<double, double, double> &second ) { return std::get<0>( first ) < std::get<0>( second ); } );
            for ( int i = 0; i < detector_frequency_mask.size() - 1; i++ ) {
                auto &[w1, a, dw1] = detector_frequency_mask[i];
                auto &[w2, c, dw2] = detector_frequency_mask[i + 1];
                double dw = w2 - w1;
                dw1 = dw;
                dw2 = dw;
            }
            Log::L2( "[PhotonStatistics] Calculating Detector Spectral Mask done, Mask size is {} \n", detector_frequency_mask.size() );
        }
        Log::L2( "[PhotonStatistics] Applying Detector Mask for {}\n", purpose );
        Log::L2( "[PhotonStatistics] Calculating FT({})...\n", purpose );
        // Calculate main fourier transform integral with spectral amplitude in tau direction.
        // Transform Tau Direction for every t_i
        Dense mat_transformed = Dense::Zero( dim, detector_frequency_mask.size() );
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
        for ( int i = 0; i < dim; i++ ) {
            for ( int v = 0; v < detector_frequency_mask.size(); v++ ) {
                auto &[frequency_v, frequency_amp, frequency_delta] = detector_frequency_mask[v];
                for ( int j = 0; j < dim; j++ ) {   // -i because of triangular grid
                    double tau = mat.tau( j, i );   // std::imag( timemat( i, j ) ) - std::real( timemat( i, j ) );
                    double dtau = mat.dtau( j, i ); // Numerics::get_taudelta( timemat, i, j );
                    mat_transformed( i, v ) += frequency_amp * std::exp( -1.0i * frequency_v * tau ) * mat( i, j ) * dtau / 2.0 / 3.1415;
                }
            }
            Timers::outputProgress( timer, progressbar, i, dim, "Detector FT tau (" + purpose + "): " );
            timer.iterate();
        }
        // Output.
        if ( s.parameters.output_dict.contains( "detectortrafo" ) ) {
            Log::L2( "[PhotonStatistics] Outputting FT({})...\n", purpose );
            FILE *f_gfunc = std::fopen( ( s.parameters.working_directory + purpose + "_mFT.txt" ).c_str(), "w" );
            fmt::print( f_gfunc, "Time\tOmega_tau\tAbs\tReal\tImag\n" );
            for ( int k = 0; k < dim; k++ ) {
                for ( int l = 0; l < detector_frequency_mask.size(); l++ ) {
                    const auto &[frequency_w_l, frequency_amp_l, frequency_delta_l] = detector_frequency_mask[l];
                    fmt::print( f_gfunc, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", std::real( mat.t( k ) ), frequency_w_l, std::abs( mat_transformed( k, l ) ), std::real( mat_transformed( k, l ) ), std::imag( mat_transformed( k, l ) ) );
                }
                fmt::print( f_gfunc, "\n" );
            }
            std::fclose( f_gfunc );
        }
        // Transform Back
        Log::L2( "[PhotonStatistics] Calculating iFT({})...\n", purpose );
        mat.set_empty();
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
        for ( int i = 0; i < dim; i++ ) {
            for ( int j = 0; j < dim; j++ ) {
                double tau = mat.tau( j, i );
                // double dtau = Numerics::get_taudelta( timemat, i, j );
                for ( int v = 0; v < detector_frequency_mask.size(); v++ ) {
                    const auto &[frequency_v, frequency_amp, frequency_delta] = detector_frequency_mask[v];
                    mat( i, j ) += std::exp( 1.0i * frequency_v * tau ) * mat_transformed( i, v ) * frequency_delta;
                }
            }
            Timers::outputProgress( timer, progressbar, i, mat.dim(), "Detector iFT t (" + purpose + "): " );
            timer.iterate();
        }

        // Tau summation
        /*  Dense mat_transformed_tau = Dense::Zero( mat.rows(), detector_frequency_mask.size() );
 #pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
         for ( int i = 0; i < mat.rows(); i++ ) {
             for ( int v = 0; v < detector_frequency_mask.size(); v++ ) {
                 auto &[frequency_v, frequency_amp, frequency_delta] = detector_frequency_mask[v];
                 for ( int j = 0; j < mat.cols() - i; j++ ) { // -i because of triangular grid
                     double tau = std::imag( timemat( i, j ) ) - std::real( timemat( i, j ) );
                     double dtau = Numerics::get_taudelta( timemat, i, j );
                     mat_transformed_tau( i, v ) += std::exp( -1.0i * frequency_v * tau ) * mat( i, j ) * dtau / 2.0 / 3.1415;
                 }
             }
             Timers::outputProgress( timer, progressbar, i, mat.rows(), "Detector FT tau (" + purpose + "): " );
             timer.iterate();
         }
         // t summation
         Dense mat_transformed = Dense::Zero( detector_frequency_mask.size(), detector_frequency_mask.size() );
 #pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
         for ( int v = 0; v < detector_frequency_mask.size(); v++ ) {
             for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
                 auto &[frequency_w, frequency_amp, frequency_delta] = detector_frequency_mask[w];
                 for ( int i = 0; i < mat.cols(); i++ ) { // -i because of triangular grid
                     double t = std::real( timemat( i, 0 ) );
                     double dt = Numerics::get_tdelta( timemat, i, 0 );
                     mat_transformed( w, v ) += std::exp( -1.0i * frequency_w * t ) * mat_transformed_tau( i, v ) * dt / 2.0 / 3.1415;
                 }
             }
             Timers::outputProgress( timer, progressbar, v, detector_frequency_mask.size(), "Detector FT t (" + purpose + "): " );
             timer.iterate();
         }
         // Output. TODO: only if wanted
         FILE *f_gfunc = std::fopen( ( s.parameters.working_directory + purpose + "_mFT.txt" ).c_str(), "w" );
         fmt::print( f_gfunc, "Time\tTau\tAbs\tReal\tImag\n" );
         for ( int k = 0; k < mat_transformed.rows(); k++ ) {
             auto &[frequency_w_k, frequency_amp_k, frequency_delta_k] = detector_frequency_mask[k];
             for ( int l = 0; l < mat_transformed.cols(); l++ ) {
                 auto &[frequency_w_l, frequency_amp_l, frequency_delta_l] = detector_frequency_mask[l];
                 fmt::print( f_gfunc, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", frequency_w_k, frequency_w_l, std::abs( mat_transformed( k, l ) ), std::real( mat_transformed( k, l ) ), std::imag( mat_transformed( k, l ) ) );
             }
             fmt::print( f_gfunc, "\n" );
         }
         std::fclose( f_gfunc );
         // Transform Back in t direction.
         Dense mat_v = Dense::Zero( mat.cols(), detector_frequency_mask.size() );
 #pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
         for ( int v = 0; v < detector_frequency_mask.size(); v++ ) {
             for ( int i = 0; i < timemat.rows(); i++ ) {
                 double t = std::real( timemat( i, 0 ) );
                 double dt = Numerics::get_tdelta( timemat, i, 0 );
                 for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
                     auto &[frequency_w, frequency_amp, frequency_delta] = detector_frequency_mask[w];
                     mat_v( i, v ) += std::exp( 1.0i * frequency_w * t ) * mat_transformed( w, v ) * dt;
                 }
             }
             Timers::outputProgress( timer, progressbar, v, detector_frequency_mask.size(), "Detector iFT t (" + purpose + "): " );
             timer.iterate();
         }
         // tau summation
         mat = Dense::Zero( mat.cols(), mat.rows() );
 #pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
         for ( int i = 0; i < mat.cols(); i++ ) {
             for ( int j = 0; j < mat.rows(); j++ ) {
                 double tau = std::imag( timemat( i, j ) ) - std::real( timemat( i, j ) );
                 double dtau = Numerics::get_taudelta( timemat, i, j );
                 for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
                     auto &[frequency_w, frequency_amp, frequency_delta] = detector_frequency_mask[w];
                     mat( i, j ) += std::exp( 1.0i * frequency_w * tau ) * mat_v( i, w ) * dtau;
                 }
             }
             Timers::outputProgress( timer, progressbar, i, mat.rows(), "Detector iFT tau (" + purpose + "): " );
             timer.iterate();
         } */
        // Transform back in t direction

        /*
        // Calculate main 2d fourier transform integral with spectral amplitude in tau direction.
        Dense mat_transformed = Dense::Zero( detector_frequency_mask.size(), detector_frequency_mask.size() );
#pragma omp parallel for collapse( 2 ) schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
        for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
            for ( int v = 0; v < detector_frequency_mask.size(); v++ ) {
                auto &[frequency_w, frequency_amp_w, frequency_delta_w] = detector_frequency_mask[w];
                auto &[frequency_v, frequency_amp_v, frequency_delta_v] = detector_frequency_mask[v];
                for ( long unsigned int i = 0; i < mat.rows(); i++ ) {
                    double t = std::real( timemat( i, 0 ) );
                    double dt = Numerics::get_tdelta( timemat, i, 0 );
                    for ( int j = 0; j < mat.cols() - i; j++ ) { // -i because of triangular grid
                        double tau = std::imag( timemat( i, j ) ) - t;
                        double dtau = Numerics::get_taudelta( timemat, i, j );
                        mat_transformed( w, v ) += std::exp( -1.0i * ( frequency_w * t + frequency_v * tau ) ) * mat( i, j ) * dtau * dt / ( 2.0 * 3.1415 ) / ( 2.0 * 3.1415 );
                    }
                }
                Timers::outputProgress( timer, progressbar, w * detector_frequency_mask.size() + v, detector_frequency_mask.size() * detector_frequency_mask.size(), "Detector FT (" + purpose + "): " );
                timer.iterate();
            }
        }
        // Output. TODO: only if wanted
        FILE *f_gfunc = std::fopen( ( s.parameters.working_directory + purpose + "_mFT.txt" ).c_str(), "w" );
        fmt::print( f_gfunc, "Time\tTau\tAbs\tReal\tImag\n" );
        for ( int k = 0; k < mat_transformed.rows(); k++ ) {
            for ( int l = 0; l < mat_transformed.cols(); l++ ) {
                fmt::print( f_gfunc, "{:d}\t{:d}\t{:.8e}\t{:.8e}\t{:.8e}\n", k, l, std::abs( mat_transformed( k, l ) ), std::real( mat_transformed( k, l ) ), std::imag( mat_transformed( k, l ) ) );
            }
            fmt::print( f_gfunc, "\n" );
        }
        std::fclose( f_gfunc );

        // Transform Back in tau direction.
        mat = Dense::Zero( mat.rows(), mat.cols() );

#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
        for ( long unsigned int i = 0; i < mat.rows(); i++ ) {
            double t = std::real( timemat( i, 0 ) );
            double dt = Numerics::get_tdelta( timemat, i, 0 );
            for ( int j = 0; j < mat.cols() - i; j++ ) { // -i because of triangular grid
                double tau = std::imag( timemat( i, j ) ) - t;
                double dtau = Numerics::get_taudelta( timemat, i, j );
                for ( int w = 0; w < detector_frequency_mask.size(); w++ ) {
                    auto &[frequency_w, frequency_amp_w, frequency_delta_w] = detector_frequency_mask[w];
                    for ( int v = 0; v < detector_frequency_mask.size(); v++ ) {
                        auto &[frequency_v, frequency_amp_v, frequency_delta_v] = detector_frequency_mask[v];
                        mat( i, j ) += std::exp( 1.0i * ( frequency_w * t + frequency_v * tau ) ) * mat_transformed( w, v ) * frequency_delta_w * frequency_delta_v;
                    }
                }
                Timers::outputProgress( timer, progressbar, i * mat.rows() + j, mat.rows() * mat.rows(), "Detector iFT (" + purpose + "): " );
                timer.iterate();
            }
        }
*/
        timer.end();
        Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), timer.getTotalIterationNumber(), "Detector (" + purpose + ")", Timers::PROGRESS_FORCE_OUTPUT );
    }
}
