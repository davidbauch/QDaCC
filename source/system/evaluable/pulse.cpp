#include "system/evaluable/pulse.h"

using namespace QDLC;

Pulse::Pulse( Parameters::universal_config &config, Parameters &p ) : Evaluable( config ) {
    Log::L2( "[System-Pulse] Creating total pulse with {} individual pulses...\n", config.property_set.at( "Amplitude" ).size() );
    generate( p );
    Log::L2( "[System-Pulse] Done!\n" );
}

Scalar Pulse::evaluate( double t ) { 
    Scalar ret = 0;
    const auto &config = get_inputs();
    for ( int i = 0; i < (int)config.property_set.at( "Amplitude" ).size(); i++ ) { // TODO: zip this using cpp23
        // Phase
        double phase = config.property_set.at( "Phase" )[i];
        // SUPER modulation
        const double centered_time = t - config.property_set.at( "Center" )[i];
        double main_freq = config.property_set.at( "Frequency" )[i] * centered_time;
        if ( config.property_set.at( "SUPERFreq" )[i] != 0.0 )
            main_freq += config.property_set.at( "SUPERDelta" )[i] / config.property_set.at( "SUPERFreq" )[i] * std::cos( config.property_set.at( "SUPERFreq" )[i] * centered_time + phase );
        if ( config.string_v.at( "Type" )[i] == "cw" && t >= config.property_set.at( "Center" )[i] ) {
            ret += config.property_set.at( "Amplitude" )[i] * std::exp( -1.0i * ( main_freq + phase + config.property_set.at( "ChirpRate" )[i] * std::pow( centered_time, 2.0 ) ) );
        } else {
            double shape = 0;
            if ( config.string_v.at( "Type" )[i] == "gauss" ) {
                const double amplitude = std::sqrt( std::pow( config.property_set.at( "ChirpRate" )[i] / config.property_set.at( "Width" )[i], 2.0 ) + std::pow( config.property_set.at( "Width" )[i], 2.0 ) );
                shape = std::exp( -0.5 * std::pow( centered_time / amplitude, config.property_set.at( "GaussAmp" )[i] ) );
            } else if ( config.string_v.at( "Type" )[i] == "sech" ) {
                const double amplitude = std::sqrt( std::pow( config.property_set.at( "ChirpRate" )[i] / config.property_set.at( "Width" )[i], 2.0 ) + std::pow( config.property_set.at( "Width" )[i], 2.0 ) );
                shape = 1./std::cosh(  centered_time / amplitude );
            } else if (config.string_v.at( "Type" )[i] == "lorentz") {
                const auto x = centered_time / config.property_set.at( "Width" )[i];
                shape = 1./(1. + std::pow( x, config.property_set.at( "GaussAmp" )[i] ));
            } 
            // Chirped frequency
            double freq = config.property_set.at( "ChirpRate" )[i] / ( std::pow( config.property_set.at( "ChirpRate" )[i], 2.0 ) + std::pow( config.property_set.at( "Width" )[i], 4.0 ) );
            if ( config.property_set.at( "CutoffDelta" )[i] != 0 ) {
                ret += config.property_set.at( "Amplitude" )[i] * shape * std::exp( -1.0i * ( main_freq + phase - config.property_set.at( "CutoffDelta" )[i] * centered_time + 0.5 * freq * std::pow( centered_time, 2.0 ) ) );
                ret += config.property_set.at( "Amplitude" )[i] * shape * std::exp( -1.0i * ( main_freq + phase + config.property_set.at( "CutoffDelta" )[i] * centered_time + 0.5 * freq * std::pow( centered_time, 2.0 ) ) );
            } else {
                ret += config.property_set.at( "Amplitude" )[i] * shape * std::exp( -1.0i * ( main_freq + phase + 0.5 * freq * std::pow( centered_time, 2.0 ) ) );
            }
        }
    }
    return ret; 
}  

// Generate array of energy-values corresponding to the Pulse
void Pulse::calculate_fourier( Parameters &p ) {
    // Log::L3( "generating type " + inputs.pulse_type + "... " );
    double omega_center = p.input_conf["PulseConf"].property["Center"];
    double omega_range = p.input_conf["PulseConf"].property["Range"];
    double dw = p.input_conf["PulseConf"].property["Res"];
    double dt = p.input_conf["PulseConf"].property["dt"]; //TODO: dt can be negative... 
    if (dt <= 0) {
        dt = get_approximated_dt();
    }
    Log::L2( "[System-Pulse] Pulsearray.size() = {}... \n", size() );
    const auto &config = get_inputs();
    if ( config.property_set.at( "Frequency" ).size() > 0 ) {
        Log::L2( "[System-Pulse] Calculating pulse fourier transformation...\n" );
        Log::L2( "[System-Pulse] Initial w_center = {}, w_range = {}, dw = {} (dt = {})\n", omega_center, omega_range, dw, dt );
        // Finding Spectral integration window
        if ( omega_center == 0 ) {
            for ( auto &input : config.property_set.at( "Frequency" ) ) {
                omega_center += input;
            }
            omega_center /= (double)( config.property_set.at( "Frequency" ).size() );
        }
        if ( omega_range == 0 ) {
            for ( const auto &input : config.property_set.at( "Width" ) ) {
                if ( input != 0.0 )
                    omega_range += 25.0 / input;
            }
            for ( const auto &input : config.property_set.at( "SUPERFreq" ) ) {
                omega_range += 2.0 * input;
            }
            omega_range /= (double)( config.property_set.at( "Frequency" ).size() );
        }
        if ( dw == 0 )
            dw = 2.0 * omega_range / double( size() );
        // Finding Temporal integration window
        double pulse_begin = 0.0;
        double pulse_end = 0.0;
        for ( int i = 0; i < config.property_set.at( "Width" ).size(); i++ ) {
            pulse_begin = std::min<double>( pulse_begin, config.property_set.at( "Center" )[i] - 8.0 * config.property_set.at( "Width" )[i] );
            pulse_end = std::max<double>( pulse_begin, config.property_set.at( "Center" )[i] + 8.0 * config.property_set.at( "Width" )[i] );
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