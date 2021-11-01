#pragma once
#include <cmath>
#include <complex>
#include <string>
#include <vector>
// #define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std::complex_literals;
#include "misc/log.h"
#include "misc/helperfunctions_string.h"
#include "misc/helperfunctions_math.h"

namespace QDLC {

namespace Misc {

// Valid conversions from ns,ps,Hz,eV,meV,mueV to corresponding SI unit (not scaled)
template <typename T>
T convertParam( const std::string input ) {
    double value = 0;
    double conversion = 1;
    int index;
    Log::L2( "Attempting to convert '{}'... ", input );
    if ( -1 != ( index = QDLC::String::instr( input, "eV" ) ) ) {
        // Found 'eV' as unit (energy), now check for scaling
        if ( input.at( index - 1 ) == 'm' ) {
            // meV
            Log::L2( "from meV to Hz..." );
            value = QDLC::Math::eV_to_Hz( std::stod( input.substr( 0, index - 1 ) ) );
            conversion = 1E-3;
        } else if ( input.at( index - 2 ) == 'm' && input.at( index - 1 ) == 'u' ) {
            // mueV
            Log::L2( "from mueV to Hz..." );
            value = QDLC::Math::eV_to_Hz( std::stod( input.substr( 0, index - 2 ) ) );
            conversion = 1E-6;
        } else if ( QDLC::Math::is_number( input.substr( index - 1, 1 ) ) ) {
            // eV
            Log::L2( "from eV to Hz..." );
            value = QDLC::Math::eV_to_Hz( std::stod( input.substr( 0, index ) ) );
            conversion = 1.0;
        } else {
            Log::L2( "Conversion of input '{}' from eV failed!\n", input );
            return (T)0.0;
        }
    } else if ( -1 != ( index = QDLC::String::instr( input, "s" ) ) ) {
        // Found 's' as unit (time)
        //fmt::print("\n {} {} {} {}\n",index, input.at(index-1)=='n',input.compare(index-1,1,"n") ,input.at(index-1));
        if ( input.at( index - 1 ) == 'n' ) {
            // ns
            Log::L2( "from ns to s..." );
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1E-9;
        } else if ( input.at( index - 1 ) == 'p' ) {
            // ps
            Log::L2( "from ps to s..." );
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1E-12; //fmt::print("{} {} ... ", value, conversion);
        } else if ( QDLC::Math::is_number( input.substr( index - 1, 1 ) ) ) {
            // s
            Log::L2( "from s to s..." );
            value = std::stod( input.substr( 0, index ) );
            conversion = 1.0;
        } else {
            Log::L2( "Conversion from input '{}' from time failed!\n", input );
            return (T)0.0;
        }
    } else if ( -1 != ( index = QDLC::String::instr( input, "Hz" ) ) ) {
        // Found 'Hz' as unit (Frequency)
        Log::L2( "from Hz to Hz..." );
        if ( QDLC::Math::is_number( input.substr( index - 1, 1 ) ) ) {
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1.0;
        } else {
            Log::L2( "Conversion from input '{}' from frequency failed!\n", input );
            return (T)0.0;
        }
    } else if ( -1 != ( index = QDLC::String::instr( input, "pi" ) ) ) {
        // Found 'Hz' as unit (Frequency)
        Log::L2( "from Xpi to rad..." );
        if ( QDLC::Math::is_number( input.substr( index - 1, 1 ) ) ) {
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1.0;
        } else {
            Log::L2( "Conversion from input '{}' from frequency failed!\n", input );
            return (T)0.0;
        }
    } else if ( QDLC::Math::is_number( input ) ) {
        // Assuming Frequency input
        Log::L2( "no conversion..." );
        value = std::stod( input );
    } else {
        // Input type is unknown
        Log::L2( "Input Type of input '{}' is unkown!\n", input );
        return (T)0.0;
    }
    Log::L2( "done, final value = {}!\n", value * conversion );
    return (T)( value * conversion );
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
T vec_max( std::vector<T> input, bool norm = true ) {
    if ( input.size() == 0 )
        return (T)( 0.0 );
    T ret = input.at( 0 );
    for ( T element : input )
        if ( ( norm && ( std::abs( element ) > std::abs( ret ) ) ) || ( !norm && ( element > ret ) ) )
            ret = element;
    return ret;
}
template <typename T>
T vec_min( std::vector<T> input, bool norm = true ) {
    if ( input.size() == 0 )
        return (T)( 0.0 );
    T ret = input.at( 0 );
    for ( T element : input )
        if ( ( norm && ( std::abs( element ) < std::abs( ret ) ) ) || ( !norm && ( element < ret ) ) )
            ret = element;
    return ret;
}

// Map a vector onto comp, meaning that if in is shorter than comp, in will be filled with standard value
template <typename T1, typename T2, typename T3>
void map_vector_to_fixed_size( const std::vector<T1> &comp, std::vector<T2> &in, const T3 standard ) {
    while ( comp.size() > in.size() ) {
        in.emplace_back( (T2)standard );
    }
}

template <typename T>
std::string vec_to_str( const std::vector<T> &input ) {
    std::string ret = "[";
    for ( auto &element : input ) {
        ret += std::to_string( element ) + ",";
    }
    ret.back() = ']';
    return ret;
}

} // namespace Misc

} // namespace QDLC