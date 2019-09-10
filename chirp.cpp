#include "chirp.h"
//#include "misc/spline.h"

Chirp::Chirp( Parameters &p ) {
    int n = (int)( ( p.t_end - p.t_start ) / p.t_step * 2.0 + 5 );
    chirparray.reserve( n );
    timearray.reserve( n );
    step = p.t_step;
    //if (p.numerics_order_t == 5 || p.numerics_order_tau == 5)
    //    steps = { 0, 1./5.*p.t_step, 3./10.*p.t_step, 1./2.*p.t_step, 4./5.*p.t_step, 8./9.*p.t_step };
    //else
    steps = {0, 0.5 * p.t_step};
    interpolant = Interpolant( p.chirp_t, p.chirp_y, p.chirp_ddt, p.chirp_type );
    generate( p );
}

/* Generate array of energy-values corresponding to the chirp */
void Chirp::generate( Parameters &p ) {
    logs.level2( "generating... " );
    for ( double t1 = p.t_start; t1 < p.t_end + 10 * p.t_step; t1 += p.t_step / 2.0 ) {
        chirparray.push_back( interpolant.evaluate( t1 ) );
        timearray.push_back( t1 );
    }
    chirparray.shrink_to_fit();
    timearray.shrink_to_fit();
    size = chirparray.size();
    logs.level2( "chirparray.size() = {}... ", size );
}

void Chirp::fileOutput( std::string filepath, Parameters &p ) {
    FILE *chirpfile = std::fopen( filepath.c_str(), "w" );
    if ( !chirpfile ) {
        logs.level2( "Failed to open outputfile for chirp!\n" );
        return;
    }
    for ( long unsigned int i = 0; i < chirparray.size(); i++ ) {
        std::fprintf( chirpfile, "%.10e\t%.10e\n", timearray.at( i ), chirparray.at( i ) );
    }
    std::fclose( chirpfile );
}

// Per index and rest is lerp delta
// TODO: option for new evaluation instead of precalculating
double Chirp::get( double t ) const {
    int i = std::floor( t / step - 1 ) * steps.size();
    while ( timearray.at( i ) < t ) {
        i++;
    }
    if ( i < 0 || i >= size ) {
        logs.level2( "!! Warning: requested chirpvalue at index {} is out of range! chirparray.size() = {}\n", i, chirparray.size() );
        i = 0;
    }
    //logs.level2("Requested t = {:.10e}, returned t = {:.10e}\n",t,timearray.at(i));
    return chirparray.at( i );
    //return lerp(chirparray.at(j),chirparray.at(j+1),delta);
}

// Per index
double Chirp::get( int i ) const {
    if ( i < 0 || i >= size ) {
        logs.level2( "!! Warning: requested chirpvalue at index {} is out of range!\n", i );
        i = 0;
    }
    return chirparray.at( i );
}