#include "solver/solver_ode.h"

/*
Apply Mask if mode is contained in mat_mode
Example: Let mat_mode = G1-hbd-bd
Let mode be G1: Apply Mask
Let mode be G2: Do not apply Mask
Let Mode be G: Apply Mask
Let Mode be G1-hbd: Apply Mask
Let Mode be G1-hbd-bd: Apply Mask
Let Mode be G1-hbd-vb: Do not apply Mask
...
*/

void QDLC::Numerics::ODESolver::initialize_detector_functions( System &s, CacheMatrix &mat ) {
    const auto dim = mat.dim();
    auto &config_time = s.parameters.input_conf["DetectorTime"];
    // Initialize Temporal Mask
    for ( const auto &mode : config_time.string_v["time_mode"] ) {
        if ( not detector_temporal_mask.contains( mode ) ) {
            Log::L2( "[PhotonStatistics] Initializing Detector Temporal Mask for mode {}...\n", mode );
            detector_temporal_mask[mode] = dDense::Zero( dim, dim );
        }
    }
    // Calculate Temporal Mask
    for ( int d = 0; d < config_time.string_v["time_mode"].size(); d++ ) {
        const auto mode = config_time.string_v["time_mode"][d];
        const double t_center = config_time.property_set["time_center"][d];
        const double t_range = config_time.property_set["time_range"][d];
        const double t_amp = config_time.property_set["time_amplitude"][d];
        const double power = config_time.property_set["time_power_amplitude"][d];
        Log::L2( "[PhotonStatistics] Adding Temporal Mask for mode {}: t_0 = {}, delta = {}, amp = {}\n", mode, t_center, t_range, power );
        auto &current_detector_temporal_mask = detector_temporal_mask[mode];
        for ( int i = 0; i < dim; i++ ) {
            double time = mat.t( i );
            for ( int j = 0; j < dim - i; j++ ) {
                double tau = mat.tau( j, i );
                current_detector_temporal_mask( i, j ) += t_amp * std::exp( -std::pow( ( time - t_center ) / t_range, power ) ) * std::exp( -std::pow( ( tau - t_center ) / t_range, power ) );
            }
        }
    }
    // Normalize each mask to 1
    if ( s.parameters.detector_normalize_functions )
        for ( auto &[mode, mask] : detector_temporal_mask ) {
            const auto max_coeff = mask.maxCoeff();
            const auto min_coeff = mask.minCoeff();
            // mask = ( mask.array() - min_coeff ) / ( max_coeff - min_coeff );
            mask /= max_coeff;
        }

    auto &config_spectral = s.parameters.input_conf["DetectorSpectral"];
    // We cache the minimum delta omega to clean up the detector mask later.
    double domega_min = std::numeric_limits<double>::max();
    // Calculate Spectral Mask
    std::map<std::string, std::vector<double>> omega_array;
    // First, Calculate omega array only
    for ( int d = 0; d < config_spectral.property_set["spectral_range"].size(); d++ ) {
        const auto mode = config_spectral.string_v["spectral_mode"][0];
        const double center = config_spectral.property_set["spectral_center"][d];
        const double range = config_spectral.property_set["spectral_range"][d];
        const double points = config_spectral.property_set["spectral_number_points"][d]; // Points-per-sigma
        const double power = config_spectral.property_set["spectral_power_amplitude"][d];
        const double amp = config_spectral.property_set["spectral_amplitude"][d];
        // This loop ensures the spectral envelope is always high enough resolution
        double w = 0;
        while ( true ) {
            double amplitude = amp * std::exp( -0.5 * std::pow( w / points, power ) );
            if ( w != 0.0 )
                omega_array[mode].emplace_back( center - w / points * range );
            omega_array[mode].emplace_back( center + w / points * range );
            w += 1;
            if ( std::abs( amplitude ) < 1E-4 )
                break;
        }
    }
    for ( auto &[mode, omega] : omega_array ) {
        // Now that we have a combined omega array, we can calculate the amplitude
        // We sort the array beforehand to make sure we get the correct delta omega.
        std::ranges::sort( omega );
        // Remove dublicates from all omega_arrays
        auto dublicates = std::unique( omega.begin(), omega.end() );
        omega.erase( dublicates, omega.end() );
        // Calculate amps. Very inefficient, but whatever.
        auto &current_detector_frequency_mask = detector_frequency_mask[mode];
        for ( int i = 0; i < omega.size(); i++ ) {
            double amplitude = 0;
            for ( int d = 0; d < config_spectral.property_set["spectral_range"].size(); d++ ) {
                const auto current_mode = config_spectral.string_v["spectral_mode"][0];
                if ( current_mode.compare( mode ) != 0 )
                    continue;
                const double center = config_spectral.property_set["spectral_center"][d];
                const double range = config_spectral.property_set["spectral_range"][d];
                const double points = config_spectral.property_set["spectral_number_points"][d]; // Points-per-sigma
                const double power = config_spectral.property_set["spectral_power_amplitude"][d];
                const double amp = config_spectral.property_set["spectral_amplitude"][d];
                auto next = amp * std::exp( -0.5 * std::pow( ( omega[i] - center ) / range, power ) );
                if ( not std::isnan( next ) ) {
                    amplitude += next;
                }
            }
            const double delta_omega = i == 0 ? omega[i] - omega[i - 1] : omega[i + 1] - omega[i];
            current_detector_frequency_mask.emplace_back( omega[i], amplitude, delta_omega );
        }
    }
    // Normalize the amplitude to 1.
    for ( auto &[mode, current_detector_frequency_mask] : detector_frequency_mask ) {
        if ( s.parameters.detector_normalize_functions ) {
            auto [minimum, maximum] = std::ranges::minmax_element( current_detector_frequency_mask, []( std::tuple<double, double, double> &first, std::tuple<double, double, double> &second ) { return std::get<1>( first ) < std::get<1>( second ); } );
            auto min = std::get<1>( *minimum );
            auto max = std::get<1>( *maximum );
            // std::ranges::for_each( current_detector_frequency_mask, [min, max]( std::tuple<double, double, double> &tuple ) { std::get<1>( tuple ) = ( std::get<1>( tuple ) - min ) / ( max - min ); } );
            std::ranges::for_each( current_detector_frequency_mask, [min, max]( std::tuple<double, double, double> &tuple ) { std::get<1>( tuple ) = std::get<1>( tuple ) /= max; } );
        }
        Log::L2( "[PhotonStatistics] Calculating Detector Spectral Mask for mode {} done, Mask size is {}.\n", mode, current_detector_frequency_mask.size() );
    }
}

void QDLC::Numerics::ODESolver::apply_detector_function( System &s, CacheMatrix &mat, const std::string &mat_mode ) {
    // If both masks are empty, calculate them
    if ( detector_temporal_mask.empty() and detector_frequency_mask.empty() )
        initialize_detector_functions( s, mat );

    const auto &purpose = mat.get_name();
    const auto dim = mat.dim();
    auto &config = s.parameters.input_conf["DetectorTime"];
    if ( config.property_set["time_center"].size() > 0 ) {
        // Current mode
        const std::string mode = config.string_v["time_mode"][0];
        auto &current_detector_temporal_mask = detector_temporal_mask[mode];
        // Apply Mask if mode is contained in mat_mode
        Log::L2( "[PhotonStatistics] Applying Detector Temporal Mask to {}... \n", purpose );
        if ( mat_mode.find( mode ) != std::string::npos )
            mat.get() = mat.get().cwiseProduct( detector_temporal_mask[mode] );
    }

    auto &spectralconfig = s.parameters.input_conf["DetectorSpectral"];
    if ( spectralconfig.property_set["spectral_range"].size() > 0 ) {
        Timer &timer = Timers::create( "Detector (" + purpose + ")" ).start();
        auto progressbar = ProgressBar();
        // Current mode
        for ( const auto &[mode,current_detector_frequency_mask] : detector_frequency_mask ) {
            // If mode is not contained in mat_mode, leave
            if ( mat_mode.find( mode ) == std::string::npos )
                return;
            // Else, apply mask
            Log::L2( "[PhotonStatistics] Applying Detector Mask '{}' for correlation function '{}'\n", mode, purpose );
            Log::L2( "[PhotonStatistics] Calculating FT({})...\n", purpose );
            // Calculate main fourier transform integral with spectral amplitude in tau direction.
            // Transform Tau Direction for every t_i
            Dense mat_transformed = Dense::Zero( dim, current_detector_frequency_mask.size() );
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_primary_threads )
            for ( int i = 0; i < dim; i++ ) {
                for ( int v = 0; v < current_detector_frequency_mask.size(); v++ ) {
                    auto &[frequency_v, frequency_amp, frequency_delta] = current_detector_frequency_mask[v];
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
                FILE *f_gfunc = std::fopen( ( s.parameters.working_directory + purpose + "_FTm.txt" ).c_str(), "w" );
                fmt::print( f_gfunc, "Time\tOmega_tau\tAbs\tReal\tImag\n" );
                for ( int k = 0; k < dim; k++ ) {
                    for ( int l = 0; l < current_detector_frequency_mask.size(); l++ ) {
                        const auto &[frequency_w_l, frequency_amp_l, frequency_delta_l] = current_detector_frequency_mask[l];
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
                    for ( int v = 0; v < current_detector_frequency_mask.size(); v++ ) {
                        const auto &[frequency_v, frequency_amp, frequency_delta] = current_detector_frequency_mask[v];
                        mat( i, j ) += std::exp( 1.0i * frequency_v * tau ) * mat_transformed( i, v ) * frequency_delta;
                    }
                }
                Timers::outputProgress( timer, progressbar, i, mat.dim(), "Detector iFT t (" + purpose + "): " );
                timer.iterate();
            }
        }
        timer.end();
        Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), timer.getTotalIterationNumber(), "Detector (" + purpose + ")", Timers::PROGRESS_FORCE_OUTPUT );
    }
}
