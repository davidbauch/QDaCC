#pragma once
#include "pulse.h"

Pulse::Pulse( Pulse::Inputs &inputs ) : inputs( inputs ) {
    logs.level2( "Creating total pulse with {} individual pulses... ", inputs.amp.size() );
    counter_evaluated = 0;
    counter_returned = 0;
    int n = (int)( ( inputs.t_end - inputs.t_start ) / inputs.t_step * 6.0 + 5 );
    pulsearray.reserve( n );
    timearray.reserve( n );
    if ( inputs.order == 5 )
        steps = {0, 1. / 5. * inputs.t_step, 3. / 10. * inputs.t_step, 1. / 2. * inputs.t_step, 4. / 5. * inputs.t_step, 8. / 9. * inputs.t_step};
    else
        steps = {0, 0.5 * inputs.t_step};
    logs.level2( "Done initializing class, creating precalculated chirp... " );
    generate();
    logs.level2( "Done!\n" );
}

dcomplex Pulse::evaluate( double t ) {
    dcomplex ret = 0;
    for ( int i = 0; i < (int)inputs.amp.size(); i++ ) {
        double amp = std::sqrt( std::pow( inputs.omega_chirp.at( i ) / inputs.sigma.at( i ), 2.0 ) + std::pow( inputs.sigma.at( i ), 2.0 ) );
        double freq = inputs.omega_chirp.at( i ) / ( std::pow( inputs.omega_chirp.at( i ), 2.0 ) + std::pow( inputs.sigma.at( i ), 4.0 ) );
        if ( inputs.type.at( i ).compare( "cw" ) == 0 && t >= inputs.center.at( i ) )
            ret += inputs.amp.at( i ) * std::exp( -1i * ( inputs.omega.at( i ) * ( t - inputs.center.at( i ) ) + inputs.omega_chirp.at( i ) * std::pow( ( t - inputs.center.at( i ) ), 2.0 ) ) );
        else if ( inputs.type.at( i ).compare( "gauss" ) == 0 )
            ret += inputs.amp.at( i ) * std::exp( -0.5 * std::pow( ( t - inputs.center.at( i ) ) / amp, 2. ) - 1i * ( inputs.omega.at( i ) * ( t - inputs.center.at( i ) ) + 0.5 * freq * std::pow( ( t - inputs.center.at( i ) ), 2.0 ) ) );
    }
    return ret;
}

dcomplex Pulse::evaluate_derivative( double t ) {
    return ( evaluate( t ) - evaluate( t - inputs.t_step ) ) / inputs.t_step;
}

dcomplex Pulse::evaluate_integral( double t ) {
    if ( pulsearray_integral.size() == 0 || t == 0.0 )
        return evaluate( t ) * inputs.t_step;
    return pulsearray_integral.at( std::max( (int)std::floor( t / inputs.t_step ), (int)pulsearray_integral.size() - 1 ) ) + evaluate( t ) * inputs.t_step;
}

// Generate array of energy-values corresponding to the Pulse
void Pulse::generate() {
    //logs.level2( "generating type " + inputs.pulse_type + "... " );
    double t;
    for ( double t1 = inputs.t_start; t1 < inputs.t_end + inputs.t_step * steps.size(); t1 += inputs.t_step ) {
        for ( int i = 0; i < (int)steps.size(); i++ ) {
            t = t1 + steps[i];
            dcomplex val = evaluate( t );
            pulsearray.push_back( val );
            pulsearray_derivative.push_back( evaluate_derivative( t ) );
            pulsearray_integral.push_back( evaluate_integral( t ) );
            timearray.push_back( t );
        }
    }

    pulsearray.shrink_to_fit();
    timearray.shrink_to_fit();
    size = pulsearray.size();
    logs.level2( "pulsearray.size() = {}... ", size );
}

void Pulse::fileOutput( std::string filepath ) {
    FILE *pulsefile = std::fopen( filepath.c_str(), "w" );
    if ( !pulsefile ) {
        logs.level2( "Failed to open outputfile for Pulse!\n" );
        return;
    }
    for ( long unsigned int i = 0; i < timearray.size() - steps.size(); i += steps.size() ) {
        fmt::print( pulsefile, "{:.10e}\t{:.10e}\n", timearray.at( i ), std::abs( pulsearray.at( i ) ) );
    }
    std::fclose( pulsefile );
}

void Pulse::Inputs::add( double _center, double _amp, double _sigma, double _omega, double _omega_chirp, std::string _type ) { //,double chirp = 0) {
    center.emplace_back( _center );
    amp.emplace_back( _amp );
    sigma.emplace_back( _sigma );
    omega.emplace_back( _omega );
    omega_chirp.emplace_back( _omega_chirp );
    type.emplace_back( _type );
    logs.level2( "Added Pulse with parameters: center = {}, amp = {}, sigma = {}, omega = {}, chirp = {}, type = {}. No filter was used.\n", _center, _amp, _sigma, _omega, _omega_chirp, _type );
}

void Pulse::Inputs::add( std::vector<double> &_center, std::vector<double> &_amp, std::vector<double> &_sigma, std::vector<double> &_omega, std::vector<double> &_omega_chirp, std::vector<std::string> &_type, std::complex<double> amp_scaling ) {
    if ( !( _center.size() == _amp.size() && _sigma.size() == _omega.size() && _amp.size() == _sigma.size() && _sigma.size() == _type.size() ) ) {
        logs.level2( "Input arrays don't have the same length! No Vectors are created, initializing pulse will fail!\n" );
        return;
    }
    for ( int i = 0; i < (int)_amp.size(); i++ ) {
        center.emplace_back( _center.at( i ) );
        amp.emplace_back( _amp.at( i ) );
        sigma.emplace_back( _sigma.at( i ) );
        omega.emplace_back( _omega.at( i ) );
        omega_chirp.emplace_back( _omega_chirp.at( i ) );
        type.emplace_back( _type.at( i ) );
        logs.level2( "Added Pulse with parameters: center = {}, amp = {}, sigma = {}, omega = {}, chirp = {}, type = {}. No filter was used.\n", _center.at( i ), _amp.at( i ), _sigma.at( i ), _omega.at( i ), _omega_chirp.at( i ), _type.at( i ) );
    }
}
void Pulse::Inputs::add( std::vector<double> &_center, std::vector<double> &_amp, std::vector<double> &_sigma, std::vector<double> &_omega, std::vector<double> &_omega_chirp, std::vector<std::string> &_type, std::vector<std::string> &_filter, std::string to_match, std::complex<double> amp_scaling ) {
    if ( !( _center.size() == _amp.size() && _sigma.size() == _omega.size() && _amp.size() == _sigma.size() && _sigma.size() == _type.size() && _type.size() == _filter.size() && _filter.size() == _omega_chirp.size() ) ) {
        logs.level2( "Input arrays don't have the same length! No Vectors are created, initializing pulse will fail!\n" );
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
            logs.level2( "Added Pulse with parameters: center = {}, amp = {}, sigma = {}, omega = {}, chirp = {}, type = {}. Filter {} was used.\n", _center.at( i ), _amp.at( i ), _sigma.at( i ), _omega.at( i ), _omega_chirp.at( i ), _type.at( i ), to_match );
        } else {
            //logs.level2( "Failed to add Pulse with parameters: center = {}, amp = {}, sigma = {}, omega = {}, type = {}. Mismatched Filter {} was used.\n", _center.at( i ), _amp.at( i ), _sigma.at( i ), _omega.at( i ), _type.at( i ), to_match );
        }
    }
}

dcomplex Pulse::get( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate( t );
    }

    int i = std::max( 0.0, std::floor( t / inputs.t_step - 1 ) ) * steps.size();
    while ( timearray.at( i ) < t ) {
        i++;
    }
    if ( std::abs( timearray.at( i ) - t ) > 1E-14 ) {
        counter_evaluated++;
        return evaluate( t );
    }
    if ( i < 0 || i >= size ) {
        logs.level2( "!! Warning: requested pulsevalue at index {} is out of range! pulsearray.size() = {}\n", i, pulsearray.size() );
        counter_evaluated++;
        return evaluate( t );
    }
    counter_returned++;
    return pulsearray.at( i );
}

dcomplex Pulse::derivative( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_derivative( t );
    }

    int i = std::max( 0.0, std::floor( t / inputs.t_step - 1 ) ) * steps.size();
    while ( timearray.at( i ) < t ) {
        i++;
    }
    if ( std::abs( timearray.at( i ) - t ) > 1E-14 ) {
        counter_evaluated++;
        return evaluate_derivative( t );
    }
    if ( i < 0 || i >= size ) {
        logs.level2( "!! Warning: requested pulsevalue at index {} is out of range! pulsearray.size() = {}\n", i, pulsearray.size() );
        counter_evaluated++;
        return evaluate_derivative( t );
    }
    counter_returned++;
    return pulsearray_derivative.at( i );
}

dcomplex Pulse::integral( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_integral( t );
    }

    int i = std::max( 0.0, std::floor( t / inputs.t_step - 1 ) ) * steps.size();
    while ( timearray.at( i ) < t ) {
        i++;
    }
    if ( std::abs( timearray.at( i ) - t ) > 1E-14 ) {
        counter_evaluated++;
        return evaluate_integral( t );
    }
    if ( i < 0 || i >= size ) {
        logs.level2( "!! Warning: requested pulsevalue at index {} is out of range! pulsearray.size() = {}\n", i, pulsearray.size() );
        counter_evaluated++;
        return evaluate_integral( t );
    }
    counter_returned++;
    return pulsearray_integral.at( i );
}

void Pulse::log() {
    logs.level2( "Pulse evaluations/returns: {}/{}\n", counter_evaluated, counter_returned );
}

void Pulse::fileOutput( std::string filepath, std::vector<Pulse> pulses ) {
    FILE *pulsefile = std::fopen( filepath.c_str(), "w" );
    if ( !pulsefile ) {
        logs.level2( "Failed to open outputfile for Pulse!\n" );
        return;
    }
    for ( long unsigned int i = 0; i < pulses.at( 0 ).timearray.size() - pulses.at( 0 ).steps.size(); i += pulses.at( 0 ).steps.size() ) {
        fmt::print( pulsefile, "{:.8e}\t", pulses.at( 0 ).timearray.at( i ) );
        for ( long unsigned int p = 0; p < pulses.size(); p++ ) {
            if ( i < pulses.at( p ).timearray.size() ) {
                fmt::print( pulsefile, "{:.8e}\t", std::real( pulses.at( p ).pulsearray.at( i ) ) );
            }
        }
        fmt::print( pulsefile, "\n" );
    }
    std::fclose( pulsefile );
}