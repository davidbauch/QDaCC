#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ostream>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>
#include <functional>
#include <array>

#include "typedef.h"

namespace QDLC {

namespace Numerics {

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
struct tuple_vector_hash : std::unary_function<std::tuple<QDLC::Type::iVector, QDLC::Type::iVector>, size_t> {
    static std::hash<QDLC::Type::iVector::Scalar> hasher;
    std::size_t operator()( std::tuple<QDLC::Type::iVector, QDLC::Type::iVector> const &tup ) const {
        size_t seed = 0;
        auto &[vec1, vec2] = tup;
        for ( size_t i = 0; i < vec1.size(); ++i ) {
            auto elem = *( vec1.data() + i );
            seed ^= hasher( elem ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
        }
        for ( size_t i = 0; i < vec2.size(); ++i ) {
            auto elem = *( vec2.data() + i );
            seed ^= hasher( elem ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
        }
        return seed;
    }
};
struct tuple_vector_compare {
    bool operator()( const std::tuple<QDLC::Type::iVector, QDLC::Type::iVector> &A, const std::tuple<QDLC::Type::iVector, QDLC::Type::iVector> &B ) const {
        return tuple_vector_hash()( A ) < tuple_vector_hash()( B );
    }
};

template <typename T>
class Tensor {
   public:
    const static int TYPE_DENSE = 0;
    const static int TYPE_SPARSE = 1;

   public:
    typedef std::unordered_map<QDLC::Type::iVector, std::unordered_map<QDLC::Type::iVector, T, vector_hash>, vector_hash> TensorMap;

   private:
    T zero = (T)0.0;
    std::vector<int> dimensions;
    // Ordered map including SparseIndex;Value pairs.
    std::vector<TensorMap> values;
    std::vector<QDLC::Type::iVector> indices_single;
    std::vector<std::tuple<QDLC::Type::iVector, QDLC::Type::iVector>> indices;
    std::unordered_set<std::tuple<QDLC::Type::iVector, QDLC::Type::iVector>, tuple_vector_hash> unique_indices;

    // Counter for which indices/values vector is active
    int current_value_vector = 0;
    int next_value_vector = 1;
    int tensor_type;

    void permute( QDLC::Type::iVector current, const std::vector<int> &dimensions, int index, T value ) {
        if ( index >= dimensions.size() ) {
            indices_single.emplace_back( current );
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
    Tensor( const std::vector<int> &dimensions, int type = TYPE_DENSE, T init_value = (T)0 ) : dimensions( dimensions ), tensor_type( type ) {
        zero = init_value;
        values = std::vector<TensorMap>( 2 );
        if ( type == TYPE_DENSE ) {
            permute( QDLC::Type::iVector::Zero( dimensions.size() ), dimensions, 0, init_value );
            for ( const auto &indx : indices_single )
                for ( const auto &indy : indices_single ) {
                    values[0][indx][indy] = init_value;
                    values[1][indx][indy] = init_value;
                    // indices.emplace_back( std::make_tuple( indx, indy ) );
                    addIndex( indx, indy );
                }
        }
    }
    Tensor( const TensorMap &other ) : dimensions( other.dimensions ),
                                       values( other.values ),
                                       indices( other.indices ),
                                       next_value_vector( other.next_value_vector ),
                                       current_value_vector( other.current_value_vector ),
                                       tensor_type( other.tensor_type ) {}

    void addIndex( const QDLC::Type::iVector &x, const QDLC::Type::iVector &y ) {
        indices.emplace_back( std::make_tuple( x, y ) );
    }
    void addUniqueIndex( const QDLC::Type::iVector &x, const QDLC::Type::iVector &y ) {
        unique_indices.insert( std::make_tuple( x, y ) );
    }

    void convertToSparse() {
        unique_indices.clear();
        indices.clear();
        tensor_type = TYPE_SPARSE;
    }
    void convertToDense() {
        // indices.clear();
        // for ( auto &[indx, map] : getCurrentValues() )
        //     for ( auto &[indy, v] : map )
        //         indices.emplace_back( std::make_tuple( indx, indy ) );
        indices = std::vector<std::tuple<QDLC::Type::iVector, QDLC::Type::iVector>>( unique_indices.begin(), unique_indices.end() );
        for ( auto &[x, y] : indices ) {
            if ( not getCurrentValues().contains( x ) or not getCurrentValues()[x].contains( y ) )
                getCurrentValues()[x][y] = 0.0;
            if ( not getNextValues().contains( x ) or not getNextValues()[x].contains( y ) )
                getNextValues()[x][y] = 0.0;
        }

        tensor_type = TYPE_DENSE;
    }

    void switchType() {
        if ( tensor_type == TYPE_DENSE ) {
            convertToSparse();
        } else {
            convertToDense();
        }
    }

    inline bool isDenseTensor() {
        return tensor_type == TYPE_DENSE;
    }
    inline bool isSparseTensor() {
        return tensor_type == TYPE_SPARSE;
    }

    inline TensorMap &getCurrentValues() {
        return values[current_value_vector];
    }
    inline TensorMap &getNextValues() {
        return values[next_value_vector];
    }
    inline std::vector<std::tuple<QDLC::Type::iVector, QDLC::Type::iVector>> &getIndices() {
        return indices;
    }

    inline void setTriplet( const QDLC::Type::iVector &indicesX, const QDLC::Type::iVector &indicesY, const T &value ) {
        getNextValues()[indicesX][indicesY] = value;
    }

    void addToTriplet( QDLC::Type::iVector indicesX, QDLC::Type::iVector indicesY, const T &value, const int i_n = -1, const int j_n = -1, const int gi_n = -1, const int gj_n = -1 ) {
        if ( i_n != -1 && j_n != -1 ) {
            for ( int i = indicesX.size() - 1; i > 0; i-- ) {
                indicesX( i ) = indicesX( i - 1 );
                indicesY( i ) = indicesY( i - 1 );
            }
            indicesX( 0 ) = i_n;
            indicesY( 0 ) = j_n;
        }
        if ( gi_n != -1 && gj_n != -1 ) {
            indicesX( 1 ) = gi_n;
            indicesY( 1 ) = gj_n;
        }
        auto &X = values[next_value_vector][indicesX];
        if ( X.contains( indicesY ) ) {
            X[indicesY] += value;
        } else {
            X[indicesY] = value;
        }
        if ( indicesX.size() == dimensions.size() )
            addUniqueIndex( indicesX, indicesY );
    }

    T &getTriplet( const QDLC::Type::iVector &indicesX, const QDLC::Type::iVector &indicesY ) {
        if ( not getCurrentValues().contains( indicesX ) or not getCurrentValues()[indicesX].contains( indicesY ) )
            return zero;
        return getCurrentValues()[indicesX][indicesY];
    }

    T &getNextTriplet( const QDLC::Type::iVector &indicesX, const QDLC::Type::iVector &indicesY ) {
        return values[next_value_vector][indicesX][indicesY];
    }

    u_int64_t nonZeros() {
        u_int64_t s = 0;
        for ( auto &[xIndices, map] : getCurrentValues() ) {
            s += map.size();
        }
        return s;
    }

    void swap() {
        // next_value_vector = current_value_vector;
        // current_value_vector = ( current_value_vector + 1 ) % 2;
        std::swap( next_value_vector, current_value_vector );
        if ( tensor_type == TYPE_SPARSE ) {
            for ( auto &[indx, map] : values[next_value_vector] ) {
                map.clear(); // Could also just set to zero i guess
            }
        }
    }

    int size() {
        return nonZeros() * sizeof( T );
    }
};

} // namespace Numerics

} // namespace QDLC