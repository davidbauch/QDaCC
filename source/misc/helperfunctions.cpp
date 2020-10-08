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

int in( char **arr, int len, char *target ) {
    int i;
    for ( i = 0; i < len; i++ ) {
        if ( std::strncmp( arr[i], target, strlen( target ) ) == 0 ) {
            return i;
        }
    }
    return -1;
}

int instr( const std::string &arr, const std::string tofind, int start ) {
    bool found = true;
    for ( int i = start; i < (int)arr.size() + 1 - (int)tofind.size(); i++ ) {
        found = true;
        for ( int j = 0; j < (int)tofind.size(); j++ ) {
            //fmt::print("comparing {} and {}... ",arr.at(i+j),tofind.at(j));
            if ( tofind.at( j ) != arr.at( i + j ) ) {
                found = false;
                j = (int)tofind.size();
            }
        }
        if ( found ) {
            //fmt::print("found at index {}\n",i);
            return i;
        }
    }
    return -1;
}

std::vector<std::string> str_to_vec( std::string input ) {
    std::vector<std::string> ret;
    if ( (int)input.size() < 3 ) {
        ret = {};
        return ret;
    }
    int s = 1; // Starting index
    int e = 1; // End
    while ( ( e = instr( input, ",", s ) ) != -1 ) {
        ret.emplace_back( input.substr( s, e - s ) );
        s = e + 1;
    }
    ret.emplace_back( input.substr( s, input.size() - s - 1 ) );
    //for (std::string el : ret)
    //    std::cout << el << "\n";
    return ret;
}

bool is_number( const std::string &s ) {
    //for (int i = 0; i < s.size(); i++) {
    //    if (!(std::isdigit(c) || c == 'e' || c == 'E' || c == '-' || c == '+' || c == '.'))
    //}
    return !s.empty() && std::find_if( s.begin(), s.end(), []( char c ) { return ( !( std::isdigit( c ) || c == 'e' || c == 'E' || c == '-' || c == '+' || c == '.' ) ); } ) == s.end();
}

std::vector<std::string> getNextInputVectorString( const std::vector<std::string> &arguments, const std::string name, int &index ) {
    return str_to_vec( arguments.at( index++ ) );
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

double getCoherent( double alpha, double N ) {
    return std::exp( -std::pow( alpha, 2.0 ) ) * std::pow( std::pow( alpha, 2.0 ), N ) / factorial( N );
}

std::string get_parameter( const std::string& key, const std::string& subkey ) {
    std::string arg = CommandlineArguments::cla.get(key,subkey);
    return arg;
}

std::vector<std::string> get_parameter_vector( const std::string& key, const std::string& subkey ) {
    std::string arg = CommandlineArguments::cla.get(key,subkey);
    // Check if input is not a vector, then output will be a new vector with only one element
    if (arg.at(0) != '[') {
        return std::vector<std::string>({arg});
    }
    return str_to_vec(arg);
}

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

std::vector<std::string> splitline( const std::string &input ) {
    std::string token, subtoken;
    std::vector<std::string> set;
    std::istringstream iss( input );
    while ( std::getline( iss, subtoken, '\t' ) ) {
        std::istringstream subiss( subtoken );
        while ( std::getline( subiss, token, ' ' ) )
            if ( token.size() > 0 ) {
                set.push_back( token );
            }
    }
    return set;
}

bool get_parameter_passed(const std::string& key, const std::string& subkey) {
    return CommandlineArguments::cla.get(key,subkey).toBool();
}

std::string getNextInputString( const std::vector<std::string> &arguments, const std::string name, int &index ) {
    logs.level2( "Trying to convert input named '{}' at index '{}' to '{}'... done\n", name, index, arguments.at( index ) );
    return arguments.at( index++ );
}