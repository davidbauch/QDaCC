#include "pulse.h"

Pulse::Pulse( Pulse::Inputs &inputs ) : inputs( inputs ) {
    Log::L2( "Creating total pulse with {} individual pulses...\n", inputs.amp.size() );
    counter_evaluated = 0;
    counter_returned = 0;
    maximum = 0;
    int n = (int)( ( inputs.t_end - inputs.t_start ) / inputs.t_step * 6.0 + 5 );
    if ( inputs.order == 5 )
        steps = {0, 1. / 5. * inputs.t_step, 3. / 10. * inputs.t_step, 1. / 2. * inputs.t_step, 4. / 5. * inputs.t_step, 8. / 9. * inputs.t_step};
    else
        steps = {0, 0.5 * inputs.t_step};
    Log::L2( "Done initializing class, creating precalculated pulse...\n" );
    generate();
    Log::L2( "Done!\n" );
}

Scalar Pulse::evaluate( double t ) {
    Scalar ret = 0;
    for ( int i = 0; i < (int)inputs.amp.size(); i++ ) {
        if ( inputs.type.at( i ).compare( "cw" ) == 0 && t >= inputs.center.at( i ) ) {
            ret += inputs.amp.at( i ) * std::exp( -1i * ( inputs.omega.at( i ) * ( t - inputs.center.at( i ) ) + inputs.omega_chirp.at( i ) * std::pow( ( t - inputs.center.at( i ) ), 2.0 ) ) );
        } else if ( inputs.type.at( i ).compare( "gauss" ) == 0 ) {
            double amp = std::sqrt( std::pow( inputs.omega_chirp.at( i ) / inputs.sigma.at( i ), 2.0 ) + std::pow( inputs.sigma.at( i ), 2.0 ) );
            double freq = inputs.omega_chirp.at( i ) / ( std::pow( inputs.omega_chirp.at( i ), 2.0 ) + std::pow( inputs.sigma.at( i ), 4.0 ) );
            ret += inputs.amp.at( i ) * std::exp( -0.5 * std::pow( ( t - inputs.center.at( i ) ) / amp, 2. ) - 1i * ( inputs.omega.at( i ) * ( t - inputs.center.at( i ) ) + 0.5 * freq * std::pow( ( t - inputs.center.at( i ) ), 2.0 ) ) );
        } else if ( inputs.type.at( i ).compare( "cutoff" ) == 0 ) {
            ret += inputs.amp.at( i ) * std::exp( -0.5 * std::pow( ( t - inputs.center.at( i ) ) / inputs.sigma.at( i ), 2. ) ) * ( std::exp( -1i * ( ( inputs.omega.at( i ) - inputs.omega_chirp.at( i ) ) * ( t - inputs.center.at( i ) ) ) ) + std::exp( -1i * ( ( inputs.omega.at( i ) + inputs.omega_chirp.at( i ) ) * ( t - inputs.center.at( i ) ) ) ) );
        }
    }
    return ret;
}

Scalar Pulse::evaluate_derivative( double t ) {
    return ( evaluate( t ) - evaluate( t - inputs.t_step ) ) / inputs.t_step;
}

Scalar Pulse::evaluate_integral( double t ) {
    if ( pulsearray_integral.size() == 0 || t == 0.0 || pulsearray_integral.count( t ) == 0 )
        return evaluate( t ) * inputs.t_step;
    return integral( t ) + evaluate( t ) * inputs.t_step;
}

// Generate array of energy-values corresponding to the Pulse
void Pulse::generate() {
    //Log::L3( "generating type " + inputs.pulse_type + "... " );
    double t;
    for ( double t1 = inputs.t_start; t1 < inputs.t_end + inputs.t_step * steps.size(); t1 += inputs.t_step ) {
        for ( int i = 0; i < (int)steps.size(); i++ ) {
            t = t1 + steps[i];
            Scalar val = evaluate( t );
            pulsearray[t] = val;
            pulsearray_derivative[t] = evaluate_derivative( t );
            pulsearray_integral[t] = evaluate_integral( t );
            if ( std::abs( val ) > maximum )
                maximum = std::abs( val );
        }
    }

    size = pulsearray.size();
    Log::L2( "Pulsearray.size() = {}... \n", size );
    if ( inputs.omega.size() > 0 ) {
        Log::L2( "Calculating pulse fourier transformation...\n" );
        double omega_center = 0;
        double omega_range = 0;
        for ( auto &input : inputs.omega ) {
            omega_center += input;
        }
        for ( auto &input : inputs.sigma ) {
            omega_range += 25.0 / input;
        }
        omega_center /= (double)( inputs.omega.size() );
        omega_range /= (double)( inputs.omega.size() );
        double dw = omega_range / 500.0;
        for ( double w = omega_center - omega_range; w < omega_center + omega_range; w += dw ) {
            fourier.emplace_back( w );
            Scalar spectral_amp = 0;
            for ( double t = inputs.t_start; t < inputs.t_end + inputs.t_step * steps.size(); t += inputs.t_step ) { // bad for bigger timesteps
                spectral_amp += get( t ) * std::exp( 1.i * w * t ) * inputs.t_step;
            }
            pulsearray_fourier.emplace_back( spectral_amp );
        }
        Log::L2( "Fourier pulsearray.size() = {}...\n", pulsearray_fourier.size() );
    }
}

void Pulse::fileOutput( std::string filepath ) {
    Log::L2( "Outputting Pulse Array to file...\n" );
    FILE *pulsefile = std::fopen( filepath.c_str(), "w" );
    if ( !pulsefile ) {
        Log::L2( "Failed to open outputfile for Pulse!\n" );
        return;
    }
    fmt::print( pulsefile, "t\tabs(Omega(t))\treal(Omega(t)))\tw\tabs(FT(Omega(t)))\n" );
    int i = 0;
    for ( double t = inputs.t_start; t < inputs.t_end + inputs.t_step; t += inputs.t_step ) {
        if ( i < pulsearray_fourier.size() )
            fmt::print( pulsefile, "{:.10e}\t{:.10e}\t{:.10e}\t{:.1e}\t{:.10e}\n", t, std::abs( pulsearray[t] ), std::real( pulsearray[t] ), fourier.at( i ), std::abs( pulsearray_fourier.at( i ) ) );
        else
            fmt::print( pulsefile, "{:.10e}\t{:.10e}\t{:.10e}\n", t, std::abs( pulsearray[t] ), std::real( pulsearray[t] ) );
        i++;
    }
    std::fclose( pulsefile );
    Log::L2( "Done!\n" );
}

void Pulse::Inputs::add( double _center, double _amp, double _sigma, double _omega, double _omega_chirp, std::string _type ) { //,double chirp = 0) {
    center.emplace_back( _center );
    amp.emplace_back( _amp );
    sigma.emplace_back( _sigma );
    omega.emplace_back( _omega );
    omega_chirp.emplace_back( _omega_chirp );
    type.emplace_back( _type );
    Log::L2( "Added Pulse with parameters: center = {}, amp = {}, sigma = {}, omega = {}, chirp = {}, type = {}. No filter was used.\n", _center, _amp, _sigma, _omega, _omega_chirp, _type );
}

void Pulse::Inputs::add( std::vector<Parameter> &_center, std::vector<Parameter> &_amp, std::vector<Parameter> &_sigma, std::vector<Parameter> &_omega, std::vector<Parameter> &_omega_chirp, std::vector<std::string> &_type, std::complex<double> amp_scaling ) {
    if ( !( _center.size() == _amp.size() && _sigma.size() == _omega.size() && _amp.size() == _sigma.size() && _sigma.size() == _type.size() ) ) {
        Log::L2( "Input arrays don't have the same length! No Vectors are created, initializing pulse will fail!\n" );
        return;
    }
    for ( int i = 0; i < (int)_amp.size(); i++ ) {
        center.emplace_back( _center.at( i ) );
        amp.emplace_back( _amp.at( i ) );
        sigma.emplace_back( _sigma.at( i ) );
        omega.emplace_back( _omega.at( i ) );
        omega_chirp.emplace_back( _omega_chirp.at( i ) );
        type.emplace_back( _type.at( i ) );
        Log::L2( "Added Pulse with parameters: center = {}, amp = {}, sigma = {}, omega = {}, chirp = {}, type = {}. No filter was used.\n", _center.at( i ), _amp.at( i ), _sigma.at( i ), _omega.at( i ), _omega_chirp.at( i ), _type.at( i ) );
    }
}
void Pulse::Inputs::add( std::vector<Parameter> &_center, std::vector<Parameter> &_amp, std::vector<Parameter> &_sigma, std::vector<Parameter> &_omega, std::vector<Parameter> &_omega_chirp, std::vector<std::string> &_type, std::vector<std::string> &_filter, std::string to_match, std::complex<double> amp_scaling ) {
    if ( !( _center.size() == _amp.size() && _sigma.size() == _omega.size() && _amp.size() == _sigma.size() && _sigma.size() == _type.size() && _type.size() == _filter.size() && _filter.size() == _omega_chirp.size() ) ) {
        Log::L2( "Input arrays don't have the same length! No Vectors are created, initializing pulse will fail!\n" );
        return;
    }
    for ( int i = 0; i < (int)_amp.size(); i++ ) {
        if ( to_match.compare( _filter.at( i ) ) == 0 ) {
            center.emplace_back( _center.at( i ) );
            amp.emplace_back( _amp.at( i ) * amp_scaling );
            sigma.emplace_back( _sigma.at( i ) );
            omega.emplace_back( _omega.at( i ) );
            omega_chirp.emplace_back( _omega_chirp.at( i ) );
            type.emplace_back( _type.at( i ) );
            Log::L2( "Added Pulse with parameters: center = {}, amp = {}, sigma = {}, omega = {}, chirp = {}, type = {}. Filter {} was used.\n", _center.at( i ), _amp.at( i ), _sigma.at( i ), _omega.at( i ), _omega_chirp.at( i ), _type.at( i ), to_match );
        } else {
            //Log::L3( "Failed to add Pulse with parameters: center = {}, amp = {}, sigma = {}, omega = {}, type = {}. Mismatched Filter {} was used.\n", _center.at( i ), _amp.at( i ), _sigma.at( i ), _omega.at( i ), _type.at( i ), to_match );
        }
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
        pulsearray[t] = val;
        return val;
    }
    counter_returned++;
    return pulsearray[t];
}

Scalar Pulse::derivative( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_derivative( t );
    }

    if ( pulsearray_derivative.count( t ) == 0 ) {
        counter_evaluated++;
        Scalar val = evaluate_derivative( t );
        pulsearray_derivative[t] = val;
        return val;
    }
    counter_returned++;
    return pulsearray_derivative[t];
}

Scalar Pulse::integral( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_integral( t );
    }

    if ( pulsearray_integral.count( t ) == 0 ) {
        counter_evaluated++;
        Scalar val = evaluate_integral( t );
        pulsearray_integral[t] = val;
        return val;
    }
    counter_returned++;
    return pulsearray_integral[t];
}

void Pulse::log() {
    Log::L2( "Pulse evaluations/returns: {}/{}\n", counter_evaluated, counter_returned );
}

void Pulse::fileOutput( std::string filepath, std::vector<Pulse> pulses ) {
    FILE *pulsefile = std::fopen( filepath.c_str(), "w" );
    if ( !pulsefile ) {
        Log::L2( "Failed to open outputfile for Pulse!\n" );
        return;
    }
    int i = 0;
    for ( double t = pulses.front().inputs.t_start; t < pulses.front().inputs.t_end + pulses.front().inputs.t_step; t += pulses.front().inputs.t_step ) {
        fmt::print( pulsefile, "{:.8e}\t", t );
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            fmt::print( pulsefile, "{:.8e}\t{:.8e}\t", std::abs( pulses.at( p ).pulsearray[t] ), std::real( pulses.at( p ).pulsearray[t] ) );
        }
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            if ( i < pulses.at( p ).pulsearray_fourier.size() )
                fmt::print( pulsefile, "{:.10e}\t{:.10e}", pulses.at( p ).fourier.at( i ), std::abs( pulses.at( p ).pulsearray_fourier.at( i ) ) );
            else
                fmt::print( pulsefile, "\t\t" );
        }
        i++;
        fmt::print( pulsefile, "\n" );
    }
    std::fclose( pulsefile );
}