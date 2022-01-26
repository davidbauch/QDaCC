#include "pulse.h"
#include "solver/solver.h"

Pulse::Pulse( Parameters::input_s &inputs, Parameters &p ) : inputs( inputs ) {
    Log::L2( "Creating total pulse with {} individual pulses...\n", inputs.numerical_v["Amplitude"].size() );
    counter_evaluated = 0;
    counter_returned = 0;
    maximum = 0;
    int n = (int)( ( p.t_end - p.t_start ) / p.t_step * 6.0 + 5 );
    if ( p.numerics_rk_order > 5 )
        steps = { QDLC::Numerics::RKCoefficients::a1 * p.t_step, QDLC::Numerics::RKCoefficients::a2 * p.t_step, QDLC::Numerics::RKCoefficients::a3 * p.t_step, QDLC::Numerics::RKCoefficients::a4 * p.t_step, QDLC::Numerics::RKCoefficients::a5 * p.t_step };
    else
        steps = { 0, 0.5 * p.t_step };
    Log::L2( "Done initializing class, creating precalculated pulse...\n" );
    generate( p.t_start, p.t_end, p.t_step );
    Log::L2( "Done!\n" );
}

Scalar Pulse::evaluate( double t ) {
    Scalar ret = 0;
    for ( int i = 0; i < (int)inputs.numerical_v["Amplitude"].size(); i++ ) {
        if ( inputs.string_v["Type"][i].compare( "cw" ) == 0 && t >= inputs.numerical_v["Center"][i] ) {
            ret += inputs.numerical_v["Amplitude"][i] * std::exp( -1.0i * ( inputs.numerical_v["Frequency"][i] * ( t - inputs.numerical_v["Center"][i] ) + inputs.numerical_v["ChirpRate"][i] * std::pow( ( t - inputs.numerical_v["Center"][i] ), 2.0 ) ) );
        } else if ( inputs.string_v["Type"][i].compare( "gauss" ) == 0 ) {
            double amp = std::sqrt( std::pow( inputs.numerical_v["ChirpRate"][i] / inputs.numerical_v["Width"][i], 2.0 ) + std::pow( inputs.numerical_v["Width"][i], 2.0 ) );
            // Chirped frequency
            double freq = inputs.numerical_v["ChirpRate"][i] / ( std::pow( inputs.numerical_v["ChirpRate"][i], 2.0 ) + std::pow( inputs.numerical_v["Width"][i], 4.0 ) );
            // SUPER modulation
            double main_freq = inputs.numerical_v["Frequency"][i] + inputs.numerical_v["SUPERDelta"][i] * std::sin( inputs.numerical_v["SUPERFreq"][i] * t );
            if ( inputs.string_v["Type"][i].compare( "cutoff" ) == 0 ) {
                ret += inputs.numerical_v["Amplitude"][i] * std::exp( -0.5 * std::pow( ( t - inputs.numerical_v["Center"][i] ) / amp, inputs.numerical_v["GaussAmp"][i] ) - 1.0i * ( ( main_freq - inputs.numerical_v["CutoffDelta"][i] ) * ( t - inputs.numerical_v["Center"][i] ) + 0.5 * freq * std::pow( ( t - inputs.numerical_v["Center"][i] ), 2.0 ) ) );
                ret += inputs.numerical_v["Amplitude"][i] * std::exp( -0.5 * std::pow( ( t - inputs.numerical_v["Center"][i] ) / amp, inputs.numerical_v["GaussAmp"][i] ) - 1.0i * ( ( main_freq + inputs.numerical_v["CutoffDelta"][i] ) * ( t - inputs.numerical_v["Center"][i] ) + 0.5 * freq * std::pow( ( t - inputs.numerical_v["Center"][i] ), 2.0 ) ) );
            } else {
                ret += inputs.numerical_v["Amplitude"][i] * std::exp( -0.5 * std::pow( ( t - inputs.numerical_v["Center"][i] ) / amp, inputs.numerical_v["GaussAmp"][i] ) - 1.0i * ( main_freq * ( t - inputs.numerical_v["Center"][i] ) + 0.5 * freq * std::pow( ( t - inputs.numerical_v["Center"][i] ), 2.0 ) ) );
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
void Pulse::generate( double t_start, double t_end, double t_step ) {
    // Log::L3( "generating type " + inputs.pulse_type + "... " );
    double t;
    for ( double t1 = t_start; t1 < t_end + t_step * steps.size(); t1 += t_step ) {
        for ( int i = 0; i < (int)steps.size(); i++ ) {
            t = t1 + steps[i];
            Scalar val = get( t );
            pulsearray[t] = val;
            if ( std::abs( val ) > maximum )
                maximum = std::abs( val );
        }
        pulsearray_derivative[t1] = derivative( t1, t_step );
        pulsearray_integral[t1] = 0.0; // integral( t1 ); // FIXME: segmentation fault, just integrate properly.
    }

    size = pulsearray.size();
    Log::L2( "Pulsearray.size() = {}... \n", size );
    if ( inputs.numerical_v["Frequency"].size() > 0 ) {
        Log::L2( "Calculating pulse fourier transformation...\n" );
        double omega_center = 0;
        double omega_range = 0;
        for ( auto &input : inputs.numerical_v["Frequency"] ) {
            omega_center += input;
        }
        for ( auto &input : inputs.numerical_v["Width"] ) {
            omega_range += 25.0 / input;
        }
        for ( auto &input : inputs.numerical_v["SUPERFreq"] ) {
            omega_range += 2.0 * input;
        }
        omega_center /= (double)( inputs.numerical_v["Frequency"].size() );
        omega_range /= (double)( inputs.numerical_v["Frequency"].size() );
        double dw = omega_range / 500.0;
        for ( double w = omega_center - omega_range; w < omega_center + omega_range; w += dw ) {
            fourier.emplace_back( w );
            Scalar spectral_amp = 0;
            for ( double t = t_start; t < t_end + t_step * steps.size(); t += t_step ) { // bad for bigger timesteps
                spectral_amp += get( t ) * std::exp( 1.i * w * t ) * t_step;
            }
            pulsearray_fourier.emplace_back( spectral_amp );
        }
        Log::L2( "Fourier pulsearray.size() = {}...\n", pulsearray_fourier.size() );
    }
}

void Pulse::fileOutput( std::string filepath, double t_start, double t_end, double t_step ) {
    Log::L2( "Outputting Pulse Array to file...\n" );
    FILE *pulsefile = std::fopen( filepath.c_str(), "w" );
    if ( !pulsefile ) {
        Log::L2( "Failed to open outputfile for Pulse!\n" );
        return;
    }
    fmt::print( pulsefile, "t\tabs(Omega(t))\treal(Omega(t)))\tw\tabs(FT(Omega(t)))\n" );
    int i = 0;
    for ( double t = t_start; t < t_end + t_step; t += t_step ) {
        if ( i < pulsearray_fourier.size() )
            fmt::print( pulsefile, "{:.10e}\t{:.10e}\t{:.10e}\t{:.1e}\t{:.10e}\n", t, std::abs( pulsearray[t] ), std::real( pulsearray[t] ), fourier[i], std::abs( pulsearray_fourier[i] ) );
        else
            fmt::print( pulsefile, "{:.10e}\t{:.10e}\t{:.10e}\n", t, std::abs( pulsearray[t] ), std::real( pulsearray[t] ) );
        i++;
    }
    std::fclose( pulsefile );
    Log::L2( "Done!\n" );
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
    FILE *pulsefile = std::fopen( filepath.c_str(), "w" );
    if ( !pulsefile ) {
        Log::L2( "Failed to open outputfile for Pulse!\n" );
        return;
    }
    fmt::print( pulsefile, "t" );
    for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
        fmt::print( pulsefile, "\tabs(Omega_{0}(t))\treal(Omega_{0}(t)))", p );
    }
    for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
        fmt::print( pulsefile, "\tw\tabs(FT(Omega_{0}(t)))\t", p );
    }
    fmt::print( pulsefile, "\n" );
    int i = 0;
    for ( double t = t_start; t < t_end + t_step; t += t_step ) {
        fmt::print( pulsefile, "{:.8e}\t", t );
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            fmt::print( pulsefile, "{:.8e}\t{:.8e}\t", std::abs( pulses.at( p ).pulsearray[t] ), std::real( pulses.at( p ).pulsearray[t] ) );
        }
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            if ( i < pulses.at( p ).pulsearray_fourier.size() )
                fmt::print( pulsefile, "{:.10e}\t{:.10e}\t", pulses.at( p ).fourier[i], std::abs( pulses.at( p ).pulsearray_fourier[i] ) );
            else
                fmt::print( pulsefile, "\t\t" );
        }
        i++;
        fmt::print( pulsefile, "\n" );
    }
    std::fclose( pulsefile );
}