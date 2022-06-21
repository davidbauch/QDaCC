#include "pulse.h"
#include "solver/solver.h"
#include "system/fileoutput.h"

Pulse::Pulse( Parameters::input_s &inputs, Parameters &p ) : inputs( inputs ) {
    Log::L2( "[System-Pulse] Creating total pulse with {} individual pulses...\n", inputs.numerical_v["Amplitude"].size() );
    counter_evaluated = 0;
    counter_returned = 0;
    maximum = 0;
    int n = (int)( ( p.t_end - p.t_start ) / p.t_step * 6.0 + 5 );
    if ( p.numerics_rk_order > 5 )
        steps = { QDLC::Numerics::RKCoefficients::a1 * p.t_step, QDLC::Numerics::RKCoefficients::a2 * p.t_step, QDLC::Numerics::RKCoefficients::a3 * p.t_step, QDLC::Numerics::RKCoefficients::a4 * p.t_step, QDLC::Numerics::RKCoefficients::a5 * p.t_step };
    else
        steps = { 0, 0.5 * p.t_step };
    Log::L2( "[System-Pulse] Done initializing class, creating precalculated pulse...\n" );
    generate( p.t_start, p.t_end, p.t_step, p.input_conf["PulseConf"].numerical["Center"], p.input_conf["PulseConf"].numerical["Range"], p.input_conf["PulseConf"].numerical["Res"], p.input_conf["PulseConf"].numerical["dt"] );
    Log::L2( "[System-Pulse] Done!\n" );
}

Scalar Pulse::evaluate( double t ) {
    Scalar ret = 0;
    for ( int i = 0; i < (int)inputs.numerical_v["Amplitude"].size(); i++ ) {
        if ( inputs.string_v["Type"][i].compare( "cw" ) == 0 && t >= inputs.numerical_v["Center"][i] ) {
            ret += inputs.numerical_v["Amplitude"][i] * std::exp( -1.0i * ( inputs.numerical_v["Frequency"][i] * ( t - inputs.numerical_v["Center"][i] ) + inputs.numerical_v["ChirpRate"][i] * std::pow( ( t - inputs.numerical_v["Center"][i] ), 2.0 ) ) );
        } else if ( inputs.string_v["Type"][i].compare( "gauss" ) == 0 ) {
            double amp = std::sqrt( std::pow( inputs.numerical_v["ChirpRate"][i] / inputs.numerical_v["Width"][i], 2.0 ) + std::pow( inputs.numerical_v["Width"][i], 2.0 ) );
            // Chirped frequency
            // FIXME: ?? chirp is doch so falsch??
            double freq = inputs.numerical_v["ChirpRate"][i] / ( std::pow( inputs.numerical_v["ChirpRate"][i], 2.0 ) + std::pow( inputs.numerical_v["Width"][i], 4.0 ) );
            // SUPER modulation
            double main_freq = inputs.numerical_v["Frequency"][i] * ( t - inputs.numerical_v["Center"][i] );
            if ( inputs.numerical_v["SUPERFreq"][i] != 0.0 )
                main_freq += inputs.numerical_v["SUPERDelta"][i] / inputs.numerical_v["SUPERFreq"][i] * std::cos( inputs.numerical_v["SUPERFreq"][i] * ( t - inputs.numerical_v["Center"][i] ) );
            if ( inputs.numerical_v["CutoffDelta"][i] != 0 ) {
                ret += inputs.numerical_v["Amplitude"][i] * std::exp( -0.5 * std::pow( ( t - inputs.numerical_v["Center"][i] ) / amp, inputs.numerical_v["GaussAmp"][i] ) - 1.0i * ( main_freq - inputs.numerical_v["CutoffDelta"][i] * ( t - inputs.numerical_v["Center"][i] ) + 0.5 * freq * std::pow( ( t - inputs.numerical_v["Center"][i] ), 2.0 ) ) );
                ret += inputs.numerical_v["Amplitude"][i] * std::exp( -0.5 * std::pow( ( t - inputs.numerical_v["Center"][i] ) / amp, inputs.numerical_v["GaussAmp"][i] ) - 1.0i * ( main_freq + inputs.numerical_v["CutoffDelta"][i] * ( t - inputs.numerical_v["Center"][i] ) + 0.5 * freq * std::pow( ( t - inputs.numerical_v["Center"][i] ), 2.0 ) ) );
            } else {
                ret += inputs.numerical_v["Amplitude"][i] * std::exp( -0.5 * std::pow( ( t - inputs.numerical_v["Center"][i] ) / amp, inputs.numerical_v["GaussAmp"][i] ) - 1.0i * ( main_freq + 0.5 * freq * std::pow( ( t - inputs.numerical_v["Center"][i] ), 2.0 ) ) );
            }
        }
    }
    return ret;
}

// DEPRECATED
Scalar Pulse::evaluate_derivative( double t, double dt ) {
    return ( evaluate( t ) - evaluate( t - dt ) ) / dt;
}
// DEPRECATED
Scalar Pulse::evaluate_integral( double t, double dt ) {
    if ( pulsearray_integral.size() == 0 || t == 0.0 || pulsearray_integral.count( t ) == 0 )
        return evaluate( t ) * dt;
    return integral( t, dt ) + evaluate( t ) * dt;
}

// Generate array of energy-values corresponding to the Pulse
void Pulse::generate( double t_start, double t_end, double t_step, double omega_center, double omega_range, double dw, double dt ) {
    // Log::L3( "generating type " + inputs.pulse_type + "... " );
    double t;
    for ( double t1 = t_start; t1 < t_end + t_step * steps.size(); t1 += t_step ) {
        for ( int i = 0; i < (int)steps.size(); i++ ) {
            t = t1 + steps[i];
            Scalar val = get( t );
            // pulsearray[t] = val;
            if ( std::abs( val ) > maximum )
                maximum = std::abs( val );
        }
        pulsearray_derivative[t1] = derivative( t1, t_step );
        pulsearray_integral[t1] = 0.0; // integral( t1 ); // FIXME: segmentation fault, just integrate properly.
    }

    size = pulsearray.size();
    Log::L2( "[System-Pulse] Pulsearray.size() = {}... \n", size );
    if ( inputs.numerical_v["Frequency"].size() > 0 ) {
        Log::L2( "[System-Pulse] Calculating pulse fourier transformation...\n" );
        Log::L2( "[System-Pulse] Initial w_center = {}, w_range = {}, dw = {} (dt = {})\n", omega_center, omega_range, dw, dt );
        // Finding Spectral integration window
        if ( omega_center == 0 ) {
            for ( auto &input : inputs.numerical_v["Frequency"] ) {
                omega_center += input;
            }
            omega_center /= (double)( inputs.numerical_v["Frequency"].size() );
        }
        if ( omega_range == 0 ) {
            for ( auto &input : inputs.numerical_v["Width"] ) {
                if ( input != 0.0 )
                    omega_range += 25.0 / input;
            }
            for ( auto &input : inputs.numerical_v["SUPERFreq"] ) {
                omega_range += 2.0 * input;
            }
            omega_range /= (double)( inputs.numerical_v["Frequency"].size() );
        }
        if ( dw == 0 )
            dw = 2.0 * double( steps.size() ) * omega_range / double( size );
        // Finding Temporal integration window
        double pulse_begin = 0.0;
        double pulse_end = 0.0;
        for ( int i = 0; i < inputs.numerical_v["Width"].size(); i++ ) {
            pulse_begin = std::min<double>( pulse_begin, inputs.numerical_v["Center"][i] - 8.0 * inputs.numerical_v["Width"][i] );
            pulse_end = std::max<double>( pulse_begin, inputs.numerical_v["Center"][i] + 8.0 * inputs.numerical_v["Width"][i] );
        }
        pulse_begin = std::max<double>( pulse_begin, t_start );
        pulse_end = std::min<double>( pulse_end, t_end );
        Log::L2( "[System-Pulse] Using w_center = {}, w_range = {}, dw = {} (dt = {}) in time window {} to {}\n", omega_center, omega_range, dw, dt, pulse_begin, pulse_end );
        for ( double w = omega_center - omega_range; w < omega_center + omega_range; w += dw ) {
            fourier.emplace_back( w );
            Scalar spectral_amp = 0;
            for ( double t = pulse_begin; t < pulse_end; t += dt ) {
                spectral_amp += get( t ) * std::exp( 1.i * w * t ) * dt;
            }
            pulsearray_fourier.emplace_back( spectral_amp );
        }
        Log::L2( "[System-Pulse] Fourier pulsearray.size() = {}...\n", pulsearray_fourier.size() );
    }
}

void Pulse::fileOutput( std::string filepath, double t_start, double t_end, double t_step ) {
    Log::L2( "Outputting Pulse Array to file...\n" );
    auto &pulsefile = FileOutput::add_file( "pulse" );
    pulsefile << "t\tabs(Omega(t))\treal(Omega(t)))\tw\treal(FT(Omega(t)))\tabs(FT(Omega(t)))\n";
    int i = 0;
    for ( double t = t_start; t < t_end + t_step; t += t_step ) {
        if ( i < pulsearray_fourier.size() )
            pulsefile << fmt::format( "{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}\n", t, std::abs( pulsearray[t] ), std::real( pulsearray[t] ), fourier[i], std::abs( pulsearray_fourier[i] ), std::real( pulsearray_fourier[i] ) );
        else
            pulsefile << fmt::format( "{:.10e}\t{:.10e}\t{:.10e}\n", t, std::abs( pulsearray[t] ), std::real( pulsearray[t] ) );
        i++;
    }
    while ( i < pulsearray_fourier.size() ) {
        pulsefile << fmt::format( "NaN\tNaN\tNaN\t{:.10e}\t{:.10e}\t{:.10e}\n", fourier[i], std::abs( pulsearray_fourier[i] ), std::real( pulsearray_fourier[i] ) );
        i++;
    }
}

Scalar Pulse::get( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate( t );
    }
    if ( pulsearray.count( t ) == 0 ) {
        counter_evaluated++;
        Scalar val = evaluate( t );
#pragma omp critical
        pulsearray[t] = val;
        return val;
    }
    counter_returned++;
    return pulsearray[t];
}

Scalar Pulse::derivative( double t, double t_step, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_derivative( t, t_step );
    }

    if ( pulsearray_derivative.count( t ) == 0 ) {
        counter_evaluated++;
        Scalar val = evaluate_derivative( t, t_step );
#pragma omp critical
        pulsearray_derivative[t] = val;
        return val;
    }
    counter_returned++;
    return pulsearray_derivative[t];
}

Scalar Pulse::integral( double t, double t_step, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_integral( t, t_step );
    }

    if ( pulsearray_integral.count( t ) == 0 ) {
        counter_evaluated++;
        Scalar val = evaluate_integral( t, t_step );
#pragma omp critical
        pulsearray_integral[t] = val;
        return val;
    }
    counter_returned++;
    return pulsearray_integral[t];
}

void Pulse::log() {
    Log::L2( "Pulse evaluations/returns: {}/{}\n", counter_evaluated, counter_returned );
}

void Pulse::fileOutput( std::string filepath, std::vector<Pulse> pulses, double t_start, double t_end, double t_step ) {
    auto &pulsefile = FileOutput::add_file( "pulse" );
    pulsefile << "t";
    for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
        pulsefile << fmt::format( "\tRe(Omega_{0}(t))\tIm(Omega_{0}(t)))", p );
    }
    for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
        pulsefile << fmt::format( "\tw\tRe(FT(Omega_{0}(t)))\tIm(FT(Omega_{0}(t)))\t", p );
    }
    pulsefile << "\n";
    int i = 0;
    for ( double t = t_start; t < t_end + t_step; t += t_step ) {
        pulsefile << fmt::format( "{:.8e}\t", t );
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            pulsefile << fmt::format( "{:.8e}\t{:.8e}\t", std::real( pulses.at( p ).pulsearray[t] ), std::imag( pulses.at( p ).pulsearray[t] ) );
        }
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            if ( i < pulses.at( p ).pulsearray_fourier.size() )
                pulsefile << fmt::format( "{:.10e}\t{:.10e}\t{:.10e}\t", pulses.at( p ).fourier[i], std::real( pulses.at( p ).pulsearray_fourier[i] ), std::imag( pulses.at( p ).pulsearray_fourier[i] ) );
            else
                pulsefile << "\t\t";
        }
        i++;
        pulsefile << "\n";
    }
    while ( i < pulses.at( 0 ).pulsearray_fourier.size() ) {
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            pulsefile << "NaN\tNaN\tNaN\t";
        }
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            pulsefile << fmt::format( "{:.10e}\t{:.10e}\t{:.10e}\t", pulses.at( p ).fourier[i], std::real( pulses.at( p ).pulsearray_fourier[i] ), std::imag( pulses.at( p ).pulsearray_fourier[i] ) );
        }
        pulsefile << "\n";
        i++;
    }
}