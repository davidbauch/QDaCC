#include "system/evaluable/pulse.h"

using namespace QDLC;

Pulse::Pulse( Parameters::input_s &config, Parameters &p ) : Evaluable( config ) {
    Log::L2( "[System-Pulse] Creating total pulse with {} individual pulses...\n", config.numerical_v.at( "Amplitude" ).size() );
    generate( p );
    Log::L2( "[System-Pulse] Done!\n" );
}

Scalar Pulse::evaluate( double t ) {
    Scalar ret = 0;
    const auto &config = get_inputs();
    for ( int i = 0; i < (int)config.numerical_v.at( "Amplitude" ).size(); i++ ) { // TODO: zip this using cpp23
        // SUPER modulation
        const double centered_time = t - config.numerical_v.at( "Center" )[i];
        double main_freq = config.numerical_v.at( "Frequency" )[i] * centered_time;
        if ( config.numerical_v.at( "SUPERFreq" )[i] != 0.0 )
            main_freq += config.numerical_v.at( "SUPERDelta" )[i] / config.numerical_v.at( "SUPERFreq" )[i] * std::cos( config.numerical_v.at( "SUPERFreq" )[i] * centered_time );
        if ( config.string_v.at( "Type" )[i].compare( "cw" ) == 0 && t >= config.numerical_v.at( "Center" )[i] ) {
            ret += config.numerical_v.at( "Amplitude" )[i] * std::exp( -1.0i * ( main_freq + config.numerical_v.at( "ChirpRate" )[i] * std::pow( centered_time, 2.0 ) ) );
        } else if ( config.string_v.at( "Type" )[i].compare( "gauss" ) == 0 ) {
            const double amplitude = std::sqrt( std::pow( config.numerical_v.at( "ChirpRate" )[i] / config.numerical_v.at( "Width" )[i], 2.0 ) + std::pow( config.numerical_v.at( "Width" )[i], 2.0 ) );
            const double gaussian_shape = std::exp( -0.5 * std::pow( centered_time / amplitude, config.numerical_v.at( "GaussAmp" )[i] ) );
            // Chirped frequency
            double freq = config.numerical_v.at( "ChirpRate" )[i] / ( std::pow( config.numerical_v.at( "ChirpRate" )[i], 2.0 ) + std::pow( config.numerical_v.at( "Width" )[i], 4.0 ) );
            if ( config.numerical_v.at( "CutoffDelta" )[i] != 0 ) {
                ret += config.numerical_v.at( "Amplitude" )[i] * gaussian_shape * std::exp( -1.0i * ( main_freq - config.numerical_v.at( "CutoffDelta" )[i] * centered_time + 0.5 * freq * std::pow( centered_time, 2.0 ) ) );
                ret += config.numerical_v.at( "Amplitude" )[i] * gaussian_shape * std::exp( -1.0i * ( main_freq + config.numerical_v.at( "CutoffDelta" )[i] * centered_time + 0.5 * freq * std::pow( centered_time, 2.0 ) ) );
            } else {
                ret += config.numerical_v.at( "Amplitude" )[i] * gaussian_shape * std::exp( -1.0i * ( main_freq + 0.5 * freq * std::pow( centered_time, 2.0 ) ) );
            }
        }
    }
    return ret;
}

// TODO: implement analytically
Scalar Pulse::evaluate_derivative( double t, double dt ) {
    return 0;
}
// TODO: implement analytically
Scalar Pulse::evaluate_integral( double t, double dt ) {
    return 0;
}

// Generate array of energy-values corresponding to the Pulse
void Pulse::calculate_fourier( Parameters &p ) {
    // Log::L3( "generating type " + inputs.pulse_type + "... " );
    double omega_center = p.input_conf["PulseConf"].numerical["Center"];
    double omega_range = p.input_conf["PulseConf"].numerical["Range"];
    double dw = p.input_conf["PulseConf"].numerical["Res"];
    double dt = p.input_conf["PulseConf"].numerical["dt"];

    Log::L2( "[System-Pulse] Pulsearray.size() = {}... \n", size() );
    const auto &config = get_inputs();
    if ( config.numerical_v.at( "Frequency" ).size() > 0 ) {
        Log::L2( "[System-Pulse] Calculating pulse fourier transformation...\n" );
        Log::L2( "[System-Pulse] Initial w_center = {}, w_range = {}, dw = {} (dt = {})\n", omega_center, omega_range, dw, dt );
        // Finding Spectral integration window
        if ( omega_center == 0 ) {
            for ( auto &input : config.numerical_v.at( "Frequency" ) ) {
                omega_center += input;
            }
            omega_center /= (double)( config.numerical_v.at( "Frequency" ).size() );
        }
        if ( omega_range == 0 ) {
            for ( const auto &input : config.numerical_v.at( "Width" ) ) {
                if ( input != 0.0 )
                    omega_range += 25.0 / input;
            }
            for ( const auto &input : config.numerical_v.at( "SUPERFreq" ) ) {
                omega_range += 2.0 * input;
            }
            omega_range /= (double)( config.numerical_v.at( "Frequency" ).size() );
        }
        if ( dw == 0 )
            dw = 2.0 * omega_range / double( size() );
        // Finding Temporal integration window
        double pulse_begin = 0.0;
        double pulse_end = 0.0;
        for ( int i = 0; i < config.numerical_v.at( "Width" ).size(); i++ ) {
            pulse_begin = std::min<double>( pulse_begin, config.numerical_v.at( "Center" )[i] - 8.0 * config.numerical_v.at( "Width" )[i] );
            pulse_end = std::max<double>( pulse_begin, config.numerical_v.at( "Center" )[i] + 8.0 * config.numerical_v.at( "Width" )[i] );
        }
        pulse_begin = std::max<double>( pulse_begin, p.t_start );
        pulse_end = std::min<double>( pulse_end, p.t_end );
        Log::L2( "[System-Pulse] Using w_center = {}, w_range = {}, dw = {} (dt = {}) in time window {} to {}\n", omega_center, omega_range, dw, dt, pulse_begin, pulse_end );
        for ( double w = omega_center - omega_range; w < omega_center + omega_range; w += dw ) {
            Scalar spectral_amp = 0;
            for ( double t = pulse_begin; t < pulse_end; t += dt ) {
                spectral_amp += get( t ) * std::exp( 1.i * w * t ) * dt;
            }
            set_fourier_value( w, spectral_amp );
        }
        Log::L2( "[System-Pulse] Fourier pulsearray.size() = {}...\n", get_fourier_size() );
    }
}

void Pulse::log() {
    Evaluable::log( "Pulse" );
};