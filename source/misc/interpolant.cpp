#include "misc/interpolant.h"
#include <format>
// 1D

// Takes X,Y[,Z], Values to return evaluations at
Interpolant::Interpolant( const std::vector<double> &interpolationPointsX, const std::vector<double> &interpolationPointsY, const std::string t ) {
    X = "[";
    Y = "[";
    for ( unsigned long i = 0; i < interpolationPointsX.size(); i++ ) {
        X += std::format( "{:.15e}", interpolationPointsX.at( i ) ) + ( i < interpolationPointsX.size() - 1 ? "," : "" );
        Y += std::format( "{:.15e}", interpolationPointsY.at( i ) ) + ( i < interpolationPointsY.size() - 1 ? "," : "" );
    }
    X += "]";
    Y += "]";
    Z = "[0]";
    generate( t );
}

Interpolant::Interpolant( const std::vector<double> &interpolationPointsX, const std::vector<double> &interpolationPointsY, const std::vector<double> &interpolationPointsZ, const std::string t ) {
    // TODO: verify inputs (length, etc)
    // Interpolant
    // Generate input strings for interpolant
    X = "[";
    Y = "[";
    Z = "[";
    for ( unsigned long i = 0; i < interpolationPointsX.size(); i++ ) {
        X += std::format( "{:.15e}", interpolationPointsX.at( i ) ) + ( i < interpolationPointsX.size() - 1 ? "," : "" );
        Y += std::format( "{:.15e}", interpolationPointsY.at( i ) ) + ( i < interpolationPointsX.size() - 1 ? "," : "" );
        Z += std::format( "{:.15e}", interpolationPointsZ.at( i ) ) + ( i < interpolationPointsX.size() - 1 ? "," : "" );
    }
    X += "]";
    Y += "]";
    Z += "]";
    // std::cout << "X = " << X << "\nY = " << Y << "\nZ = " << Z << "\nType = " << t <<"\n";
    generate( t );
}

void Interpolant::generate( std::string t = "linear" ) {
    type = t;
    alglib::real_1d_array x = X.c_str();
    alglib::real_1d_array y = Y.c_str();
    alglib::real_1d_array z = Z.c_str();
    if ( !type.compare( "cubic" ) )
        alglib::spline1dbuildcubic( x, y, p );
    else if ( !type.compare( "hermite" ) )
        alglib::spline1dbuildhermite( x, y, z, p );
    else if ( !type.compare( "akima" ) )
        alglib::spline1dbuildakima( x, y, p );
    else if ( !type.compare( "catmullrom" ) )
        alglib::spline1dbuildcatmullrom( x, y, p );
    else if ( !type.compare( "monotone" ) )
        alglib::spline1dbuildmonotone( x, y, p );
    else
        alglib::spline1dbuildlinear( x, y, p );
}

double Interpolant::evaluate( double x ) const {
    return (double)alglib::spline1dcalc( p, x );
}

std::vector<double> Interpolant::evaluate( std::vector<double> &xar ) {
    std::vector<double> ret;
    ret.reserve( xar.size() );
    for ( double x : xar ) {
        ret.emplace_back( (double)alglib::spline1dcalc( p, x ) );
    }
    return ret;
}

// 2D
template <class T>
Interpolant2d<T>::Interpolant2d( std::vector<double> &interpolationPointsX, std::vector<double> &interpolationPointsY, T &interpolationPointsZ ) {
    X = "[";
    Y = "[";
    Z = "[";
    for ( unsigned long i = 0; i < interpolationPointsX.size(); i++ ) {
        X += std::format( "{:.15e}", interpolationPointsX.at( i ) ) + ( i < interpolationPointsX.size() - 1 ? "," : "" );
        Y += std::format( "{:.15e}", interpolationPointsY.at( i ) ) + ( i < interpolationPointsY.size() - 1 ? "," : "" );
    }
    for ( unsigned long i = 0; i < interpolationPointsY.size(); i++ ) {
        for ( unsigned long j = 0; j < interpolationPointsX.size(); j++ ) {
            Z += std::format( "{:.15e},", interpolationPointsZ( j, i ) );
        }
    }
    X += "]";
    Y += "]";
    Z.back() = "]";
    sx = interpolationPointsX.size();
    sy = interpolationPointsY.size();
    generate();
}

template <class T>
void Interpolant2d<T>::generate() {
    alglib::real_1d_array x = X.c_str();
    alglib::real_1d_array y = Y.c_str();
    alglib::real_1d_array z = Z.c_str();
    alglib::spline2dbuildbilinearv( x, sx, y, sy, z, 1, p );
}

template <class T>
double Interpolant2d<T>::evaluate( double x, double y ) const {
    return (double)alglib::spline2dcalc( p, x, y );
}

template <class T>
T Interpolant2d<T>::evaluate( std::vector<double> &xar, std::vector<double> &yar ) {
    T ret;
    T::Zero( xar.size(), yar.size() );
    int i = 0;
    for ( double x : xar ) {
        int j = 0;
        for ( double y : yar ) {
            ret( i, j ) = (double)alglib::spline2dcalc( p, x, y );
        }
    }
    return ret;
}

// g++ .\interpolant.cpp .\*.o
// set yrange[-1.1:1.1]; plot "out.txt" u 1:2 w l t 'linear',"out.txt" u 1:3 w l t 'cubic',"out.txt" u 1:4 w l t 'hermite',"out.txt" u 1:5 w l t 'dakima',"out.txt" u 1:6 w l t 'catmullrom',"out.txt" u 1:7 w l t 'monotone'