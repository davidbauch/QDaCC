#include "misc/helperfunctions.h"

int delta( int i, int j ) {
    return ( i == j ) ? 1 : 0;
}

double Hz_to_eV( double hz, double scaling ) {
    return hz * 6.582119516885722624E-16 * scaling;
}

double eV_to_Hz( double eV ) {
    return eV / 6.582119516885722624E-16;
}

double Hz_to_wavelength( double hz, double scaling ) {
    return 299792458 * 2 * M_PI / hz * scaling;
}

double rabiFrequency( double deltaE, double g, double n ) {
    return std::sqrt( std::pow( deltaE, 2 ) + 4. * std::pow( g, 2 ) * n );
}

//int in( char **arr, int len, char *target ) {
//    int i;
//    for ( i = 0; i < len; i++ ) {
//        if ( std::strncmp( arr[i], target, strlen( target ) ) == 0 ) {
//            return i;
//        }
//    }
//    return -1;
//}

bool is_number( const std::string &s ) {
    //for (int i = 0; i < s.size(); i++) {
    //    if (!(std::isdigit(c) || c == 'e' || c == 'E' || c == '-' || c == '+' || c == '.'))
    //}
    return !s.empty() && std::find_if( s.begin(), s.end(), []( char c ) { return ( !( std::isdigit( c ) || c == 'e' || c == 'E' || c == '-' || c == '+' || c == '.' ) ); } ) == s.end();
}

std::vector<std::string> getNextInputVectorString( const std::vector<std::string> &arguments, const std::string name, int &index ) {
    return QDLC::Misc::String::str_to_vec( arguments.at( index++ ) );
}

int vec_find_str( std::string toFind, const std::vector<std::string> &input ) {
    for ( int i = 0; i < (int)input.size(); i++ ) {
        if ( input.at( i ).compare( toFind ) == 0 )
            return i;
    }
    return -1;
}

std::vector<std::string> argv_to_vec( int argc, char **argv ) {
    std::vector<std::string> ret;
    ret.reserve( argc );
    for ( int i = 0; i < argc; i++ )
        ret.push_back( std::string( argv[i] ) );
    return ret;
}

inline double factorial( double n ) {
    return ( n == 1.0 || n == 0.0 ) ? 1.0 : factorial( n - 1.0 ) * n;
}
inline double factorial_range( double upper, double lower ) {
    return ( upper == 1.0 || upper == 0.0 || upper <= lower ) ? 1.0 : factorial_range( upper - 1.0, lower ) * upper;
}

double getCoherent( double alpha, double N ) {
    return std::exp( -std::pow( alpha, 2.0 ) ) * std::pow( std::pow( alpha, 2.0 ), N ) / factorial( N );
}
std::complex<double> getSqueezed( double r, double phi, double N ) {
    return 1.0 / std::sqrt( std::cosh( r ) ) * std::pow(  -std::exp( 1i * phi ) * std::tanh( r ) , N ) * std::sqrt( factorial_range( 2 * N, N ) * factorial( N ) ) / factorial( N ) / std::pow( 2.0, N );
}

//std::vector<std::string> get_parameter_vector( const std::string &key, const std::string &subkey ) {
//    std::string arg = QDLC::CommandlineArguments::get( key, subkey );
//    // Check if input is not a vector, then output will be a new vector with only one element
//    if ( arg.at( 0 ) != '[' ) {
//        return std::vector<std::string>( { arg } );
//    }
//    return str_to_vec( arg );
//}

Eigen::MatrixXcd project_matrix( const Eigen::MatrixXcd &input ) {
    Eigen::MatrixXcd ret = Eigen::MatrixXcd::Zero( input.rows(), input.cols() );
    for ( int i = 0; i < ret.rows(); i++ ) {
        for ( int j = 0; j < ret.cols(); j++ ) {
            if ( real( input( i, j ) ) != 0 || imag( input( i, j ) ) != 0 )
                ret( i, j ) = 1.0;
        }
    }
    return ret;
}

std::vector<std::string> splitline( const std::string &input, const char splitter ) {
    std::string token, subtoken;
    std::vector<std::string> set;
    std::istringstream iss( input );
    while ( std::getline( iss, subtoken, '\t' ) ) {
        std::istringstream subiss( subtoken );
        while ( std::getline( subiss, token, splitter ) )
            if ( token.size() > 0 ) {
                set.push_back( token );
            }
    }
    return set;
}

std::string getNextInputString( const std::vector<std::string> &arguments, const std::string name, int &index ) {
    Log::L2( "Trying to convert input named '{}' at index '{}' to '{}'... done\n", name, index, arguments.at( index ) );
    return arguments.at( index++ );
}