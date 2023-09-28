#pragma once

#include <vector>
#include <functional>
#include "typedef.h"

namespace QDACC::Numerics {

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


template <typename T>
struct iVector_hash {
    static std::hash<T> hasher;
    std::size_t operator()( Eigen::Matrix<T, -1, 1> const &vec ) const {
        size_t seed = 0;
        for ( const auto &el : vec ) {
            seed ^= hasher( el ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
        }
        return seed;
    }
};


/**
 * @brief Compare function for Eigen Vectors
 * Required to store keys of type Vector in std::map
*/
template <typename T>
struct iVector_compare {
    bool operator()( const Eigen::Matrix<T, -1, 1> &A, const Eigen::Matrix<T, -1, 1> &B ) const {
        return iVector_hash<T>()( A ) < iVector_hash<T>()( B );
    }
};



} // namespace QDACC::Numerics