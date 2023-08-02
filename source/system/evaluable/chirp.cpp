#include "system/evaluable/chirp.h"
#include "solver/solver.h"
#include "system/fileoutput.h"

using namespace QDLC;

Chirp::Chirp( Parameters::universal_config &config, Parameters &p ) : Evaluable( config ) {
    const std::string type = config.string.at( "Type" );
    Log::L2( "[System-Chirp] Creating Chirp with {} points of type {}...\n", config.property_set.at( "Times" ).size(), type );
    if ( type != "sine" ) {
        Log::L2( "[System-Chirp] Done initializing class, creating interpolant...\n" );
        interpolant = Interpolant( config.property_set.at( "Times" ), config.property_set.at( "Amplitude" ), config.property_set.at( "ddt" ), config.string.at( "Type" ) );
        isSineChirp = false;
        Log::L2( "[System-Chirp] Done creating interpolant, generating chirp...\n" );
    } else {
        Log::L2( "[System-Chirp] No interpolant class used, using sine chirp instead.\n" );
        isSineChirp = true;
    }
    generate( p );
    Log::L2( "[System-Chirp] Done!\n" );
}

Scalar Chirp::evaluate( double t ) {  
    if ( isSineChirp ) { 
        Scalar ret = 0; 
        const auto &config = get_inputs();
        for ( long unsigned int i = 0; i < config.property_set.at( "Times" ).size(); i++ ) {
            ret += config.property_set.at( "Amplitude" ).at( i ) * std::sin( 6.283 * t * config.property_set.at( "ddt" ).at( i ) + 6.283 * config.property_set.at( "Times" ).at( i ) );
        }
        return ret;
    }
    return interpolant.evaluate( t );
}

void Chirp::log() {
    Evaluable::log( "Chirp" );
}

void Chirp::calculate_fourier( Parameters &p ){};