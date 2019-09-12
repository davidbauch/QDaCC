#pragma once
#include "pulse.h"

Pulse::Pulse( Pulse::Inputs &inputs ) : inputs( inputs ) {
    logs.level2( "Creating total pulse with {} individual pulses... ", inputs.amp.size() );
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

// Generate array of energy-values corresponding to the Pulse
void Pulse::generate() {
    //logs.level2( "generating type " + inputs.pulse_type + "... " );
    for ( int i = 0; i < (int)inputs.amp.size(); i++ ) {
        if ( inputs.type.at( i ).compare( "gauss_pi" ) == 0 ) {
            inputs.type.at( i ) = "gauss";
            // Adjust amplitute in case this wasn't done in parameters adjustInput() function before:
            inputs.amp.at( i ) = M_PI / ( std::sqrt( 2.0 * M_PI ) * inputs.sigma.at( i ) );
        }
    }
    double t;
    for ( double t1 = inputs.t_start; t1 < inputs.t_end + inputs.t_step * steps.size(); t1 += inputs.t_step ) {
        for ( int i = 0; i < (int)steps.size(); i++ ) {
            std::complex<double> val = 0;
            t = t1 + steps[i];
            for ( int i = 0; i < (int)inputs.amp.size(); i++ ) {
                if ( inputs.type.at( i ).compare( "cw" ) == 0 && t >= inputs.center.at( i ) )
                    val += inputs.amp.at( i ) * std::exp( -1i * inputs.omega.at( i ) * ( t - inputs.center.at( i ) ) );
                else if ( inputs.type.at( i ).compare( "gauss" ) == 0 )
                    val += inputs.amp.at( i ) * std::exp( -0.5 * std::pow( ( t - inputs.center.at( i ) ) / inputs.sigma.at( i ), 2. ) - 1i * inputs.omega.at( i ) * ( t - inputs.center.at( i ) ) );
            }
            if ( std::abs( val ) > 1E-10 )
                pulsearray.push_back( val );
            else
                pulsearray.push_back( 0 );
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
        std::fprintf( pulsefile, "%.10e\t%.10e\n", timearray.at( i ), std::abs( pulsearray.at( i ) ) );
    }
    std::fclose( pulsefile );
}

void Pulse::Inputs::add( double _center, double _amp, double _sigma, double _omega, std::string _type ) { //,double chirp = 0) {
    center.emplace_back( _center );
    amp.emplace_back( _amp );
    sigma.emplace_back( _sigma );
    omega.emplace_back( _omega );
    type.emplace_back( _type );
}

void Pulse::Inputs::add( std::vector<double> &_center, std::vector<double> &_amp, std::vector<double> &_sigma, std::vector<double> &_omega, std::vector<std::string> &_type ) {
    if ( !( _center.size() == _amp.size() && _sigma.size() == _omega.size() && _amp.size() == _sigma.size() && _sigma.size() == _type.size() ) ) {
        logs.level2( "Input arrays don't have the same length! No Vectors are created, initializing pulse will fail!\n" );
        return;
    }
    for ( int i = 0; i < (int)_amp.size(); i++ ) {
        center.emplace_back( _center.at( i ) );
        amp.emplace_back( _amp.at( i ) );
        sigma.emplace_back( _sigma.at( i ) );
        omega.emplace_back( _omega.at( i ) );
        type.emplace_back( _type.at( i ) );
    }
}

std::complex<double> Pulse::get( double t ) const {
    int i = std::floor( t / inputs.t_step - 1 ) * steps.size();
    while ( timearray.at( i ) < t ) {
        i++;
    }
    if ( i < 0 || i >= size ) {
        logs.level2( "!! Warning: requested pulsevalue at index {} is out of range! pulsearray.size() = {}\n", i, pulsearray.size() );
        i = 0;
    }
    //logs.level2("Requested t = {:.10e}, returned t = {:.10e}\n",t,timearray.at(i));
    return pulsearray.at( i );
    //return lerp<std::complex<double>>(pulsearray.at(j),pulsearray.at(j+1),delta);
}

std::complex<double> Pulse::get( int i ) const {
    if ( i < 0 || i >= size ) {
        logs.level2( "!! Warning: requested pulsevalue at index {} is out of range!\n", i );
        i = 0;
    }
    return pulsearray.at( i );
}