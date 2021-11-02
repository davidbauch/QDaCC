#include "misc/helperfunctions_math.h"

int QDLC::Math::delta( int i, int j ) {
    return ( i == j ) ? 1 : 0;
}

double QDLC::Math::Hz_to_eV( double hz, double scaling ) {
    return hz * 6.582119516885722624E-16 * scaling;
}

double QDLC::Math::eV_to_Hz( double eV ) {
    return eV / 6.582119516885722624E-16;
}

double QDLC::Math::Hz_to_wavelength( double hz, double scaling ) {
    return 299792458.0 * 2.0 * QDLC::Math::PI / hz * scaling;
}

double QDLC::Math::rabiFrequency( double deltaE, double g, double n ) {
    return std::sqrt( std::pow( deltaE, 2 ) + 4. * std::pow( g, 2 ) * n );
}

double QDLC::Math::factorial( double n ) {
    return std::tgamma( n + 1.0 );//( n == 1.0 || n == 0.0 ) ? 1.0 : QDLC::Math::factorial( n - 1.0 ) * n;
}

double QDLC::Math::factorial_range( double upper, double lower ) {
    return ( upper == 1.0 || upper == 0.0 || upper <= lower ) ? 1.0 : QDLC::Math::factorial_range( upper - 1.0, lower ) * upper;
}

double QDLC::Math::getCoherent( double alpha, double N ) {
    return std::exp( -std::pow( alpha, 2.0 ) ) * std::pow( std::pow( alpha, 2.0 ), N ) / QDLC::Math::factorial( N );
}
QDLC::Type::Scalar QDLC::Math::getSqueezed( double r, double phi, double N ) {
    return 1.0 / std::sqrt( std::cosh( r ) ) * std::pow(  -std::exp( 1.0i * phi ) * std::tanh( r ) , N ) * std::sqrt( factorial_range( 2 * N, N ) * factorial( N ) ) / factorial( N ) / std::pow( 2.0, N );
}

bool QDLC::Math::is_number( const std::string &s ) {
    return !s.empty() && std::find_if( s.begin(), s.end(), []( char c ) { return ( !( std::isdigit( c ) || c == 'e' || c == 'E' || c == '-' || c == '+' || c == '.' ) ); } ) == s.end();
}