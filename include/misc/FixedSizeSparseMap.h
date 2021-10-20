#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ostream>
#include <set>
#include <unordered_map>
#include <map>
#include <functional>
#include <array>

template <typename Scalar>
class FixedSizeSparseMap {
   public:
    // Index Type definition
    typedef uint64_t Index;
    typedef Eigen::VectorXi Vector;

    // https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
    template <typename T>
    struct vector_hash : std::unary_function<T, size_t> {
        static std::hash<typename T::Scalar> hasher;
        std::size_t operator()( T const &matrix ) const {
            size_t seed = 0;
            for ( size_t i = 0; i < matrix.size(); ++i ) {
                auto elem = *( matrix.data() + i );
                seed ^= hasher( elem ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
            }
            return seed;
        }
    };
    template <typename T>
    struct vector_compare {
        bool operator()( const T &A, const T &B ) const {
            return vector_hash<Vector>()( A ) < vector_hash<Vector>()( B );
        }
    };

    typedef std::unordered_map<Vector, std::unordered_map<Vector, Scalar, vector_hash<Vector>>, vector_hash<Vector>> ValueMap; // Save Map of (ivec,jvec,value) triplets
    typedef std::unordered_map<Index, std::unordered_map<Index, ValueMap>> TensorMap;                                          // inside of map sorted by their first two indices (i0,j0)

    // Ordered map including SparseIndex;Value pairs.
    std::vector<TensorMap> values; // Index is [CacheVector]
    std::vector<TensorMap> cache;  // Index is [Thread]

    // Counter for which indices/values Vector is active
    int current_value_vector = 0;

   public:
    FixedSizeSparseMap(){};
    FixedSizeSparseMap( int max_threads = 1 ) {
        auto filler = TensorMap();
        values.emplace_back( filler );
        values.emplace_back( filler );
        for ( auto i = 0; i < max_threads; i++ )
            cache.emplace_back( filler );
    }
    FixedSizeSparseMap( const FixedSizeSparseMap &other ) : values( other.values ),
                                                            current_value_vector( other.current_value_vector ),
                                                            cache( other.cache ) {}

    TensorMap &get() {
        return values[current_value_vector];
    }

    void addTriplet( Vector indicesX, Vector indicesY, const Scalar &value, int thread, const int i_n = -1, const int j_n = -1, const int gi_n = -1, const int gj_n = -1 ) {
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
        return nonZeros() * sizeof( Scalar );
    }
};