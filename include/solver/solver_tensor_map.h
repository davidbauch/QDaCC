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
    template <typename T>
    struct vector_compare {
        bool operator()( const T &A, const T &B ) const {
            return vector_hash()( A ) < vector_hash()( B );
        }
    };

    typedef std::unordered_map<QDLC::Type::iVector, std::unordered_map<QDLC::Type::iVector, QDLC::Type::Scalar, vector_hash>, vector_hash> ValueMap; // Save Map of (ivec,jvec,value) triplets
    typedef std::unordered_map<QDLC::Type::Index, std::unordered_map<QDLC::Type::Index, ValueMap>> TensorMap;                                        // inside of map sorted by their first two indices (i0,j0)

    // Ordered map including SparseIndex;Value pairs.
    std::vector<TensorMap> values; // index is [CacheVector]
    std::vector<TensorMap> cache;  // index is [Thread]

    // Counter for which indices/values QDLC::Type::iVector is active
    int current_value_vector = 0;

   public:
    Tensor(){};
    Tensor( int max_threads = 1 ) {
        auto filler = TensorMap();
        values.emplace_back( filler );
        values.emplace_back( filler );
        for ( auto i = 0; i < max_threads; i++ )
            cache.emplace_back( filler );
    }
    Tensor( const Tensor &other ) : values( other.values ),
                                    current_value_vector( other.current_value_vector ),
                                    cache( other.cache ) {}

    TensorMap &get() {
        return values[current_value_vector];
    }

    void addTriplet( QDLC::Type::iVector indicesX, QDLC::Type::iVector indicesY, const QDLC::Type::Scalar &value, int thread, const int i_n = -1, const int j_n = -1, const int gi_n = -1, const int gj_n = -1 ) {
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
        //TODO: create NxN matrix of omp_lock_t and lock this before every &X=..., then no pragma omp critical is needed. or execute by master only, idk
        auto &X = cache[thread][indicesX( 0 )][indicesY( 0 )][indicesX];
        if ( X.count( indicesY ) > 0 ) {
            X[indicesY] += value;
        } else {
            X[indicesY] = value;
        }
    }

    long long int nonZeros() {
        long long s = 0;
        for ( auto &[i0, a0] : values[current_value_vector] )
            for ( auto &[j0, a1] : a0 )
                for ( auto &[indx, b0] : a1 )
                    s += b0.size();
        return s;
    }

    void swapAndClearCache( int max_index, int threads ) {
        current_value_vector = ( current_value_vector + 1 ) % 2;
#pragma omp parallel for collapse( 2 ) num_threads( threads )
        for ( int i = 0; i < max_index; i++ ) {
            for ( int j = 0; j < max_index; j++ ) {
                values[current_value_vector][i][j].clear();
                for ( auto t = 0; t < cache.size(); t++ ) {
                    values[current_value_vector][i][j].merge( cache[t][i][j] );
                    // Manually insert preaxisting elements
                    if ( cache[t][i][j].size() > 0 ) {
                        for ( auto &[sparse_index_x, inner] : cache[t][i][j] ) {
                            for ( auto &[sparse_index_y, value] : inner ) {
#pragma omp critical
                                { values[current_value_vector][i][j][sparse_index_x][sparse_index_y] += value; }
                            }
                        }
                        //std::cout << "Manually moved values in cache[" << t << "][" << i << "][" << j << "]: " << cache[t][i][j].size() << " to values\n";
                        cache[t][i][j].clear();
                    }
                }
            }
        }
    }

    int size() {
        return nonZeros() * sizeof( QDLC::Type::Scalar );
    }
};

} // namespace Numerics

} // namespace QDLC