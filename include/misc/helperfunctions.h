#pragma once
#include "global.h"

#ifndef M_PI
#    define M_PI 3.1415926535
#endif

int delta( int i, int j );

double Hz_to_eV( double hz, double scaling = 1.0 );

double eV_to_Hz( double eV );

double Hz_to_wavelength( double hz, double scaling = 1E9 );

double rabiFrequency( double deltaE, double g, double n );

// is target in arr of length len?
// returns index where target was found, -1 if target is not in arr
int in( char **arr, int len, char *target );

int instr( const std::string &arr, const std::string tofind, int start = 0 );

// Converts an input string ( like e.g. [1,2,3,4] ) into a vector<string>
std::vector<std::string> str_to_vec( std::string input = "[]" );

#define toStr( x ) std::to_string( x )

template <typename T>
T lerp( T a, T b, double c ) {
    if ( c < 0 || c > 1 ) {
        //logs("Warning: Invalid lerp delta! delta = {:e}\n",c);
        return a;
    }
    return ( 1.0 - c ) * a + c * b;
}

bool is_number( const std::string &s );

// Valid conversions from ns,ps,Hz,eV,meV,mueV to corresponding SI unit (not scaled)
template <typename T>
T convertParam( const std::string input ) {
    double value = 0;
    double conversion = 1;
    int index;
    logs.level2( "Attempting to convert '{}'... ", input );
    if ( -1 != ( index = instr( input, "eV" ) ) ) {
        // Found 'eV' as unit (energy), now check for scaling
        if ( input.at( index - 1 ) == 'm' ) {
            // meV
            logs.level2( "from meV to Hz..." );
            value = eV_to_Hz( std::stod( input.substr( 0, index - 1 ) ) );
            conversion = 1E-3;
        } else if ( input.at( index - 2 ) == 'm' && input.at( index - 1 ) == 'u' ) {
            // mueV
            logs.level2( "from mueV to Hz..." );
            value = eV_to_Hz( std::stod( input.substr( 0, index - 2 ) ) );
            conversion = 1E-6;
        } else if ( is_number( input.substr( index - 1, 1 ) ) ) {
            // eV
            logs.level2( "from eV to Hz..." );
            value = eV_to_Hz( std::stod( input.substr( 0, index ) ) );
            conversion = 1.0;
        } else {
            logs.level2( "Conversion of input '{}' from eV failed!\n", input );
            return (T)0.0;
        }
    } else if ( -1 != ( index = instr( input, "s" ) ) ) {
        // Found 's' as unit (time)
        //fmt::print("\n {} {} {} {}\n",index, input.at(index-1)=='n',input.compare(index-1,1,"n") ,input.at(index-1));
        if ( input.at( index - 1 ) == 'n' ) {
            // ns
            logs.level2( "from ns to s..." );
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1E-9;
        } else if ( input.at( index - 1 ) == 'p' ) {
            // ps
            logs.level2( "from ps to s..." );
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1E-12; //fmt::print("{} {} ... ", value, conversion);
        } else if ( is_number( input.substr( index - 1, 1 ) ) ) {
            // s
            logs.level2( "from s to s..." );
            value = std::stod( input.substr( 0, index ) );
            conversion = 1.0;
        } else {
            logs.level2( "Conversion from input '{}' from time failed!\n", input );
            return (T)0.0;
        }
    } else if ( -1 != ( index = instr( input, "Hz" ) ) ) {
        // Found 'Hz' as unit (Frequency)
        logs.level2( "from Hz to Hz..." );
        if ( is_number( input.substr( index - 1, 1 ) ) ) {
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1.0;
        } else {
            logs.level2( "Conversion from input '{}' from frequency failed!\n", input );
            return (T)0.0;
        }
    } else if ( -1 != ( index = instr( input, "pi" ) ) ) {
        // Found 'Hz' as unit (Frequency)
        logs.level2( "from Xpi to rad..." );
        if ( is_number( input.substr( index - 1, 1 ) ) ) {
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1.0;
        } else {
            logs.level2( "Conversion from input '{}' from frequency failed!\n", input );
            return (T)0.0;
        }
    } else if ( is_number( input ) ) {
        // Assuming Frequency input
        logs.level2( "no conversion..." );
        value = std::stod( input );
    } else {
        // Input type is unknown
        logs.level2( "Input Type of input '{}' is unkown!\n", input );
        return (T)0.0;
    }
    logs.level2( "done, final value = {}!\n", value * conversion );
    return ( T )( value * conversion );
}

template <typename T>
std::vector<T> convertParam( const std::vector<std::string> input ) {
    std::vector<T> ret;
    for ( std::string in : input ) {
        ret.emplace_back( convertParam<T>( in ) );
    }
    return ret;
}

template <typename T>
T getNextInput( const std::vector<std::string> &arguments, const std::string name, int &index ) {
    logs.level2( "Trying to convert input named '{}' at index '{}' to '{}'...", name, index, arguments.at( index ) );
    return convertParam<T>( arguments.at( index++ ) );
}

std::string getNextInputString( const std::vector<std::string> &arguments, const std::string name, int &index );

template <typename T>
std::vector<T> getNextInputVector( const std::vector<std::string> &arguments, const std::string name, int &index ) {
    logs.level2( "Trying to convert input named '{}' at index '{}' to '{}'...", name, index, arguments.at( index ) );
    return convertParam<T>( str_to_vec( arguments.at( index++ ) ) );
}

std::vector<std::string> getNextInputVectorString( const std::vector<std::string> &arguments, const std::string name, int &index );

template <typename T>
T vec_max( std::vector<T> input, bool norm = true ) {
    if ( input.size() == 0 )
        return ( T )( 0.0 );
    T ret = input.at( 0 );
    for ( T element : input )
        if ( ( norm && ( std::abs( element ) > std::abs( ret ) ) ) || ( !norm && ( element > ret ) ) )
            ret = element;
    return ret;
}
template <typename T>
T vec_min( std::vector<T> input, bool norm = true ) {
    if ( input.size() == 0 )
        return ( T )( 0.0 );
    T ret = input.at( 0 );
    for ( T element : input )
        if ( ( norm && ( std::abs( element ) < std::abs( ret ) ) ) || ( !norm && ( element < ret ) ) )
            ret = element;
    return ret;
}

// finds a string in
int vec_find_str( std::string toFind, const std::vector<std::string> &input );

std::vector<std::string> argv_to_vec( int argc, char **argv );

double factorial( double n );

double getCoherent( double alpha, double N );

std::string get_parameter(const std::string& key, const std::string& subkey = "");

template <class T>
T get_parameter(const std::string& key, const std::string& subkey = "") {
    std::string arg = CommandlineArguments::cla.get(key,subkey);
    return convertParam<T>(arg);
}

bool get_parameter_passed(const std::string& key, const std::string& subkey = "");

template <class T>
std::vector<T> get_parameter_vector(const std::string& key, const std::string& subkey = "") {
    std::string arg = CommandlineArguments::cla.get(key,subkey);
    // Check if input is not a vector, then output will be a new vector with only one element
    if (arg.at(0) != '[') {
        T elem = convertParam<T>(arg);
        return std::vector<T>({elem});
    }
    return convertParam<T>( str_to_vec(arg) );
}

std::vector<std::string> get_parameter_vector(const std::string& key, const std::string& subkey = "");

// Map a vector onto comp, meaning that if in is shorter than comp, in will be filled with standard value
template <typename T1, typename T2, typename T3> 
void map_vector_to_standard(const std::vector<T1>& comp, std::vector<T2>& in, const T3 standard) {
    while (comp.size() > in.size()) {
        in.emplace_back((T2)standard);
    }
}

class Parse_Parameters {
   private:
    std::vector<std::string> parameters;

   public:
    Parse_Parameters(){};
    Parse_Parameters( const std::vector<std::string> &arguments, const std::vector<std::string> &flags, const std::vector<int> &number_of_parameters_per_flag, std::string name = "" ) {
        if ( flags.size() != number_of_parameters_per_flag.size() ) {
            logs.level2( "Error: Flag input vector and number vector have different sizes!\n" );
        }
        if ( name.size() > 0 )
            logs.level2( "Parsing input block '{}'\n", name );
        int i = 0;
        for ( auto flag : flags ) {
            int index = vec_find_str( flag, arguments );
            for ( int j = 0; j < number_of_parameters_per_flag.at( i ); j++ ) {
                if ( index == -1 )
                    parameters.push_back( "_standard" );
                else if ( arguments.at( index ).at( 1 ) != '-' )
                    parameters.push_back( "_found" );
                else
                    parameters.push_back( arguments.at( index + j + 1 ) );
            }
            i++;
        }
        //logs.level2( "Parameters to read: {}\nRead Parameters: {}\nReturned Parameters: {}\n", fmt::join( flags, ", " ), fmt::join( arguments, ", " ), fmt::join( parameters, ", " ) );
    }
    // Returns a parsed parameter at position [at]. If this parameter could not be parsed, the passed standard parameter is returned instead.
    std::string get( int at, std::string standard ) {
        if ( parameters.at( at ).compare( "_standard" ) )
            return parameters.at( at );
        return standard;
    }
    // If a vector of ats is passed, the last non-standard parameter is returned.
    std::string get( std::vector<int> ats, std::string standard ) {
        std::string ret = standard;
        for ( int at : ats ) {
            if ( parameters.at( at ).compare( "_standard" ) )
                ret = parameters.at( at );
        }
        return ret;
    }
    template <typename T>
    T get( std::vector<int> ats, std::string standard ) {
        std::string ret = standard;
        for ( int at : ats ) {
            if ( parameters.at( at ).compare( "_standard" ) )
                ret = parameters.at( at );
        }
        return convertParam<T>( ret );
    }
    template <typename T>
    T get( int at, std::string standard = "" ) {
        if ( parameters.at( at ).compare( "_standard" ) )
            return convertParam<T>( parameters.at( at ) );
        return convertParam<T>( standard );
    }
    // Special syntax: If no standard value is passed, the get function returns a boolean wether the string was found or not, instead of a value
    bool get( int at ) {
        if ( parameters.at( at ).compare( "_standard" ) )
            return true;
        return false;
    }
};

Eigen::MatrixXcd project_matrix( const Eigen::MatrixXcd &input );

template <typename T>
Eigen::SparseMatrix<T> project_matrix_sparse( const Eigen::SparseMatrix<T> &input ) {
    Eigen::SparseMatrix<T> ret = Eigen::SparseMatrix<T>( input.rows(), input.cols() );
    std::vector<Eigen::Triplet<T>> ret_v;
    for ( int k = 0; k < input.outerSize(); ++k ) {
        for ( Eigen::SparseMatrix<std::complex<double>>::InnerIterator it( input, k ); it; ++it ) {
            ret_v.emplace_back( it.row(), it.col(), 1.0 );
        }
    }
    // Generate new Matrix from triplet list
    ret.setFromTriplets( ret_v.begin(), ret_v.end() );
    return ret;
}

std::vector<std::string> splitline( const std::string &input = "" );

template <typename T>
void init_sparsevector( std::vector<Eigen::SparseMatrix<T>> &mat, int dim, int count ) {
    mat.clear();
    Eigen::SparseMatrix<T> fill = Eigen::SparseMatrix<T>( dim, dim );
    for ( int i = 0; i < count; i++ ) {
        mat.emplace_back( fill );
    }
}