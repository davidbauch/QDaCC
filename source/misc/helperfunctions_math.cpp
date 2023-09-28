#include "misc/helperfunctions_math.h"

int QDACC::Math::delta( int i, int j ) {
    return ( i == j ) ? 1 : 0;
}

double QDACC::Math::Hz_to_eV( double hz, double scaling ) {
    return hz * QDACC::Math::ev_conversion * scaling;
}

double QDACC::Math::eV_to_Hz( double eV ) {
    return eV / QDACC::Math::ev_conversion;
}

double QDACC::Math::Hz_to_wavelength( double hz, double scaling ) {
    return 299792458.0 * 2.0 * QDACC::Math::PI / hz * scaling;
}

double QDACC::Math::wavelength_to_Hz( double wl ) {
    return 299792458.0 * 2.0 * QDACC::Math::PI / wl;
}

double QDACC::Math::rabiFrequency( double deltaE, double g, double n ) {
    return std::sqrt( std::pow( deltaE, 2 ) + 4. * std::pow( g, 2 ) * n );
}

double QDACC::Math::factorial( double n ) {
    return std::tgamma( n + 1.0 ); //( n == 1.0 || n == 0.0 ) ? 1.0 : QDACC::Math::factorial( n - 1.0 ) * n;
}

double QDACC::Math::factorial_range( double upper, double lower ) {
    return ( upper == 1.0 || upper == 0.0 || upper <= lower ) ? 1.0 : QDACC::Math::factorial_range( upper - 1.0, lower ) * upper;
}

QDACC::Type::Scalar QDACC::Math::getCoherent( QDACC::Type::Scalar alpha, double N ) {
    return std::exp( -std::pow( std::abs(alpha), 2.0 ) / 2.0 ) * std::pow( alpha, N ) / std::sqrt(QDACC::Math::factorial( N ));
}


QDACC::Type::Scalar QDACC::Math::getThermal( 
    QDACC::Type::Scalar alpha, double N ) {
    return std::pow( alpha, N ) / std::pow( 1.0 + alpha, N + 1 );
}

QDACC::Type::Scalar QDACC::Math::getSqueezed( double r, double phi, double N ) {
    return 1.0 / std::sqrt( std::cosh( r ) ) * std::pow( -std::exp( 1.0i * phi ) * std::tanh( r ), N ) * std::sqrt( factorial_range( 2 * N, N ) * factorial( N ) ) / factorial( N ) / std::pow( 2.0, N );
}

bool QDACC::Math::is_number( const std::string &s ) {
    return !s.empty() && std::find_if( s.begin(), s.end(), []( char c ) { return ( !( std::isdigit( c ) || c == 'e' || c == 'E' || c == '-' || c == '+' || c == '.' ) ); } ) == s.end();
}