#pragma once

#include <vector>
#include <functional>

namespace QDLC::Numerics {

/**
 * @brief Hash function for std::vector
 * Required to store keys of type std::vector in std::unordered_map or std::map
 * After https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
*/
template <typename T>
struct vector_hash {
    static std::hash<T> hasher;
    std::size_t operator()( std::vector<T> const &vec ) const {
        size_t seed = 0;
        for ( const auto &el : vec ) {
            seed ^= hasher( el ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
        }
        return seed;
    }
};

/**
 * @brief Compare function for std::vector
 * Required to store keys of type std::vector in std::map
*/
template <typename T>
struct vector_compare {
    bool operator()( const std::vector<T> &A, const std::vector<T> &B ) const {
        return vector_hash<T>()( A ) < vector_hash<T>()( B );
    }
};

} // namespace QDLC::Numerics