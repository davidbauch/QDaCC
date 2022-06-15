#pragma once

#include <string>
#include <complex>
#include <cmath>
#include "misc/helperfunctions_math.h"
#include "typedef.h"

namespace QDLC {

namespace Math {

const double PI = 3.14159265358979323846264338327950288;
const double vLight = 299792458.0;
const double ev_conversion = 6.582119516885722624E-16;

inline int delta( int i, int j );

double Hz_to_eV( double hz, double scaling = 1.0 );

double eV_to_Hz( double eV );

double Hz_to_wavelength( double hz, double scaling = 1E9 );
double wavelength_to_Hz( double wl );

double rabiFrequency( double deltaE, double g, double n );

template <typename T>
T lerp( T a, T b, double c ) {
    assert( c < 0 || c > 1 );
    return ( 1.0 - c ) * a + c * b;
}

bool is_number( const std::string &s );

double factorial( double n );

double factorial_range( double upper, double lower );

double getCoherent( double alpha, double N );

double getThermal( double alpha, double N );

QDLC::Type::Scalar getSqueezed( double r, double phi, double N );

bool is_number( const std::string &s );
bool is_number( const char &s );

template <typename T>
double abs2( const T &v1 ) {
    return std::real( v1 ) * std::real( v1 ) + std::imag( v1 ) * std::imag( v1 );
}

} // namespace Math

} // namespace QDLC