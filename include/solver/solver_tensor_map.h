#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ostream>
#include <set>
#include <unordered_map>
#include <map>
#include <functional>
#include <array>

#include "typedef.h"

namespace QDLC {

namespace Numerics {

template <class T>
class Tensor {
   public:
    // https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
    struct vector_hash : std::unary_function<QDLC::Type::iVector, size_t> {
        static std::hash<QDLC::Type::iVector::Scalar> hasher;
        std::size_t operator()( QDLC::Type::iVector const &vec ) const {
            size_t seed = 0;
            for ( size_t i = 0; i < vec.size(); ++i ) {
                auto elem = *( vec.data() + i );
                seed ^= hasher( elem ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
            }
            return seed;
        }
    };
    struct vector_compare {
        bool operator()( const QDLC::Type::iVector &A, const QDLC::Type::iVector &B ) const {
            return vector_hash()( A ) < vector_hash()( B );
        }
    };

    typedef std::unordered_map<QDLC::Type::iVector, std::unordered_map<QDLC::Type::iVector, T, vector_hash>, vector_hash> TensorMap;

    // Ordered map including SparseIndex;Value pairs.
    std::vector<TensorMap> values; // index is [ValueVector]
    std::vector<QDLC::Type::iVector> indices;

    // Counter for which indices/values vector is active
    int current_value_vector = 0;
    int next_value_vector = 1;

    void permute( QDLC::Type::iVector current, const std::vector<int> &dimensions, int index, T value ) {
        if ( index >= dimensions.size() ) {
            indices.emplace_back( current );
            values[0][QDLC::Type::iVector::Zero( dimensions.size() )][current] = value;
        } else {
            for ( int c = 0; c < dimensions[index]; c++ ) {
                current( index ) = c;
                permute( current, dimensions, index + 1, value );
            }
        }
    }

   public:
    Tensor(){};
    Tensor( const std::vector<int> &dimensions, T init_value = (T)0 ) {
        values = std::vector<TensorMap>( 2 );
        permute( QDLC::Type::iVector::Zero( dimensions.size() ), dimensions, 0, init_value );
        // for ( auto &[key, val] : values[0][QDLC::Type::iVector::Zero( dimensions.size() )] ) {
        //     values[0][key] = values[0][QDLC::Type::iVector::Zero( dimensions.size() )];
        // }
        for ( const auto &indx : indices )
            for ( const auto &indy : indices ) {
                values[0][indx][indy] = init_value;
                values[1][indx][indy] = init_value;
            }
        // values[1] = values[0];
    }
    Tensor( const Tensor &other ) : values( other.values ),
                                    indices( other.indices ),
                                    next_value_vector( other.next_value_vector ),
                                    current_value_vector( other.current_value_vector ) {}

    TensorMap &getCurrentValues() {
        return values[current_value_vector];
    }
    TensorMap &getNextValues() {
        return values[next_value_vector];
    }
    std::vector<QDLC::Type::iVector> &getIndices() {
        return indices;
    }

    void setTriplet( const QDLC::Type::iVector &indicesX, const QDLC::Type::iVector &indicesY, const T &value ) {
        values[next_value_vector][indicesX][indicesY] = value;
    }

    void addToTriplet( const QDLC::Type::iVector &indicesX, const QDLC::Type::iVector &indicesY, const T &value ) {
        values[next_value_vector][indicesX][indicesY] += value;
    }

    T &getTriplet( const QDLC::Type::iVector &indicesX, const QDLC::Type::iVector &indicesY ) {
        return values[current_value_vector][indicesX][indicesY];
    }

    T &getNextTriplet( const QDLC::Type::iVector &indicesX, const QDLC::Type::iVector &indicesY ) {
        return values[next_value_vector][indicesX][indicesY];
    }

    // void addTriplet( QDLC::Type::iVector indicesX, QDLC::Type::iVector indicesY, const QDLC::Type::Scalar &value, int thread, const int i_n = -1, const int j_n = -1, const int gi_n = -1, const int gj_n = -1 ) {
    //     if ( i_n != -1 && j_n != -1 ) {
    //         for ( int i = indicesX.size() - 1; i > 0; i-- ) {
    //             indicesX( i ) = indicesX( i - 1 );
    //             indicesY( i ) = indicesY( i - 1 );
    //         }
    //         indicesX( 0 ) = i_n;
    //         indicesY( 0 ) = j_n;
    //     }
    //     if ( gi_n != -1 && gj_n != -1 ) {
    //         indicesX( 1 ) = gi_n;
    //         indicesY( 1 ) = gj_n;
    //     }
    //     values[current_value_vector][indicesX][indicesY] += value;
    // }

    long long int nonZeros() {
        long long s = 0;
        for ( auto &[xIndices, map] : values[current_value_vector] )
            s += map.size();
        return s;
    }

    void swap() {
        next_value_vector = current_value_vector;
        current_value_vector = ( current_value_vector + 1 ) % 2;
    }

    int size() {
        return nonZeros() * sizeof( T );
    }
};

} // namespace Numerics

} // namespace QDLC