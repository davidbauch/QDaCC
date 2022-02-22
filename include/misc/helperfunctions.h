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
T convertParam( std::string input ) {
    T value = (T)0;
    T conversion = (T)1;
    bool complex = false;
    if ( input.front() == 'i' and QDLC::Math::is_number( input.substr( 1, 1 ) ) ) {
        complex = true;
        input = input.substr( 1 );
    }
    int index;
    Log::L3( "Attempting to convert '{}'...\n", input );
    // "Meter" Scale
    if ( input.back() == 'm' ) {
        index = input.size() - 1;
        if ( input.at( index - 1 ) == 'm' ) {
            Log::L3( "\tfrom mm to Hz...\n" );
            value = QDLC::Math::wavelength_to_Hz( 1E-3 * std::stod( input.substr( 0, index - 1 ) ) );
        }
        if ( input.at( index - 2 ) == 'm' and input.at( index - 1 ) == 'u' ) {
            Log::L3( "\tfrom mum to Hz...\n" );
            value = QDLC::Math::wavelength_to_Hz( 1E-6 * std::stod( input.substr( 0, index - 2 ) ) );
        }
        if ( input.at( index - 1 ) == 'n' ) {
            Log::L3( "\tfrom nm to Hz...\n" );
            value = QDLC::Math::wavelength_to_Hz( 1E-9 * std::stod( input.substr( 0, index - 1 ) ) );
        }
    }
    // Electron Volt Scale
    else if ( -1 != ( index = QDLC::String::instr( input, "eV" ) ) ) {
        // Found 'eV' as unit (energy), now check for scaling
        if ( input.at( index - 1 ) == 'm' ) {
            // meV
            Log::L3( "\tfrom meV to Hz...\n" );
            value = QDLC::Math::eV_to_Hz( std::stod( input.substr( 0, index - 1 ) ) );
            conversion = 1E-3;
        } else if ( int( index ) > 1 and input.at( index - 2 ) == 'm' and input.at( index - 1 ) == 'u' ) {
            // mueV
            Log::L3( "\tfrom mueV to Hz...\n" );
            value = QDLC::Math::eV_to_Hz( std::stod( input.substr( 0, index - 2 ) ) );
            conversion = 1E-6;
        } else if ( QDLC::Math::is_number( input.substr( index - 1, 1 ) ) ) {
            // eV
            Log::L3( "\tfrom eV to Hz...\n" );
            value = QDLC::Math::eV_to_Hz( std::stod( input.substr( 0, index ) ) );
            conversion = 1.0;
        } else {
            Log::L3( "Conversion of input '{}' from eV failed!\n", input );
            return (T)0.0;
        }
    }
    // Second Scale
    else if ( -1 != ( index = QDLC::String::instr( input, "s" ) ) ) {
        // Found 's' as unit (time)
        // fmt::print("\n {} {} {} {}\n",index, input.at(index-1)=='n',input.compare(index-1,1,"n") ,input.at(index-1));
        if ( input.at( index - 1 ) == 'n' ) {
            // ns
            Log::L3( "\tfrom ns to s...\n" );
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1E-9;
        } else if ( input.at( index - 1 ) == 'p' ) {
            // ps
            Log::L3( "\tfrom ps to s...\n" );
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1E-12; // fmt::print("{} {} ... ", value, conversion);
        } else if ( input.at( index - 1 ) == 'f' ) {
            // fs
            Log::L3( "\tfrom fs to s...\n" );
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1E-15; // fmt::print("{} {} ... ", value, conversion);
        } else if ( QDLC::Math::is_number( input.substr( index - 1, 1 ) ) ) {
            // s
            Log::L3( "\tfrom s to s...\n" );
            value = std::stod( input.substr( 0, index ) );
            conversion = 1.0;
        } else {
            Log::L3( "Conversion from input '{}' from time failed!\n", input );
            return (T)0.0;
        }
    } else if ( -1 != ( index = QDLC::String::instr( input, "Hz" ) ) ) {
        // Found 'Hz' as unit (Frequency)
        Log::L3( "\tfrom Hz to Hz...\n" );
        if ( QDLC::Math::is_number( input.substr( index - 1, 1 ) ) ) {
            value = std::stod( input.substr( 0, index - 1 ) );
            conversion = 1.0;
        } else {
            Log::L3( "Conversion from input '{}' from frequency failed!\n", input );
            return (T)0.0;
        }
    } else if ( -1 != ( index = QDLC::String::instr( input, "pi" ) ) ) {
        // Found 'Hz' as unit (Frequency)
        Log::L3( "\tfrom Xpi to rad...\n" );
        if ( QDLC::Math::is_number( input.substr( index - 1, 1 ) ) ) {
            value = std::stod( input.substr( 0, index ) );
            conversion = 1.0;
        } else {
            Log::L3( "Conversion from input '{}' from frequency failed!\n", input );
            return (T)0.0;
        }
    } else if ( QDLC::Math::is_number( input ) ) {
        // Assuming Frequency input
        Log::L3( "\tno conversion...\n" );
        value = std::stod( input );
    } else {
        // Input type is unknown
        Log::L3( "Input Type of input '{}' is unkown!\n", input );
        return (T)0.0;
    }
    Log::L3( "Done! Final value = {}{}\n", value * conversion, complex ? ", is imag!" : "!" );
    // if ( complex )
    //     return (T)( 1.i * value * conversion );
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

template <typename T, typename L>
T vec_filter( const std::vector<T> &input, const L &lambda ) {
    if ( input.size() == 0 )
        return (T)( 0.0 );
    T ret = input.at( 0 );
    for ( T element : input )
        if ( lambda( element, ret ) )
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