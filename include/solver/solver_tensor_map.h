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
#include "misc/log.h"

// Windows Workaround for unsigned int struct
#ifndef u_int64_t
#    define u_int64_t uint64_t
#endif

namespace QDLC {

namespace Numerics {

// TODO:
//  Switch iVector index to plain c++ vector
//  Implement uint_64t index instead of vector if partially summed is used -> much less memory
//  And feasable, because e.g. 36^2*2^(2*7) (NC = 8) = 21233664 which easyly fits into an integer.
// Maybe even limit to uint64_t index, because all other cases are useless.
// AND: have 1 vector of indices that map the integer onto the index vector!

// That way indices can even be restored, because the main index vector is just 1,2,3....,max (length of index map)

// Other TODO: correlation func for PI reimplementieren, dann einen punkt für 09 nachrechnen.
// Rechnungen auf cluster meddeln

// https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
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
template <typename T>
struct vector_compare {
    bool operator()( const std::vector<T> &A, const std::vector<T> &B ) const {
        return vector_hash<T>()( A ) < vector_hash<T>()( B );
    }
};
// struct tuple_vector_hash {
//     static std::hash<QDLC::Type::iVector::QDLC::Type::Scalar> hasher;
//     std::size_t operator()( std::tuple<QDLC::Type::iVector, QDLC::Type::iVector> const &tup ) const {
//         size_t seed = 0;
//         auto &[vec1, vec2] = tup;
//         for ( size_t i = 0; i < vec1.size(); ++i ) {
//             auto elem = *( vec1.data() + i );
//             seed ^= hasher( elem ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
//         }
//         for ( size_t i = 0; i < vec2.size(); ++i ) {
//             auto elem = *( vec2.data() + i );
//             seed ^= hasher( elem ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
//         }
//         return seed;
//     }
// };
// struct tuple_vector_compare {
//     bool operator()( const std::tuple<QDLC::Type::iVector, QDLC::Type::iVector> &A, const std::tuple<QDLC::Type::iVector, QDLC::Type::iVector> &B ) const {
//         return tuple_vector_hash()( A ) < tuple_vector_hash()( B );
//     }
// };

// TODO: Implementation in .cpp file

class Tensor {
   public:
    using IndexFlat = uint64_t;
    using Index = uint16_t;
    // Tensor Rank indices (N0,N1,...,NN,M0,M1,...,MM)
    using IndexVector = std::vector<Index>;
    // QDLC::Type::Scalar QDLC::Type::Scalar type
    using ValueVector = std::vector<QDLC::Type::Scalar>;
    using IndexMap = std::map<IndexVector, IndexFlat, vector_compare<Index>>;
    using IndexRefVector = std::vector<IndexVector>;

   private:
    // Current Tensor Values. Creating a new Tensor will allocate this as new
    ValueVector values;
    // Index Maps can be static because we only need them once in memory.
    // Save IndexVector - to - IndexFlat map
    static IndexMap index_vector_to_index_flat_struct;
    // Save References of the keys in vector to allow for fast acces / conversion of IndexFlat - to - IndexVector
    // using IndexRefVector = std::vector<std::reference_wrapper<IndexVector>>;
    static IndexRefVector index_flat_to_index_vector_struct;

    void permute( IndexVector current, const IndexVector &dimensions, int index = 0 ) {
        if ( index == dimensions.size() ) {
            const auto index = index_vector_to_index_flat_struct.size();
            index_vector_to_index_flat_struct[current] = index;
            if ( current.size() > 0 )
                index_flat_to_index_vector_struct[index] = current;
        } else {
            for ( Index c = 0; c < dimensions[index]; c++ ) {
                current[index] = c;
                permute( current, dimensions, index + 1 );
            }
        }
    }

   public:
    Tensor() = default;
    Tensor( const IndexFlat num, const QDLC::Type::Scalar &init_value = 0 ) {
        values = ValueVector( num, init_value );
    }
    Tensor( const IndexVector &dimensions, QDLC::Type::Scalar init_value = 0 ) {
        // Reserve Memory and initialize value vector
        auto max_elements = 1;
        index_vector_to_index_flat_struct.clear();
        index_flat_to_index_vector_struct.clear();
        std::ranges::for_each( dimensions, [&]( const auto &el ) { max_elements *= el; } );
        index_flat_to_index_vector_struct = IndexRefVector( max_elements ); // FIXME: max_elements ist für Biex nicht gleich der hinzugefügten elemente?????????????
        values = ValueVector( max_elements, init_value );
        // Calculate and Save all possible index permutations
        permute( IndexVector( dimensions.size(), 0 ), dimensions );
        // Save References to map keys in vector for fast lookup
        // std::ranges::for_each( index_vector_to_index_flat_struct, []( const auto &el ) { index_flat_to_index_vector_struct[el.second] = el.first; } );
        // Lazy fix for upper fixme:
        for ( int i = index_flat_to_index_vector_struct.size() - 1; i > 0; i-- )
            if ( index_flat_to_index_vector_struct[i].size() == 0 )
                index_flat_to_index_vector_struct.erase( index_flat_to_index_vector_struct.begin() + i );
        Log::L2( "[PathIntegral] Added {} elements to the dimension vector ({} elements to the inverse map).\n", index_vector_to_index_flat_struct.size(), index_flat_to_index_vector_struct.size() );
    }
    Tensor( const Tensor &other ) = default; //: values( other.values ){};

    // Access Operator
    inline QDLC::Type::Scalar &operator()( const IndexVector &index ) {
        return values[to_flat( index )];
    }
    inline QDLC::Type::Scalar &operator()( const IndexFlat &index ) {
        return values[index];
    }

    // Converts the indexvector into a flat index using the precalculated map.
    inline IndexFlat &to_flat( const IndexVector &index ) {
        return index_vector_to_index_flat_struct[index];
    }
    // Converts the flat intex into an indexvector using the precalculated map
    inline IndexVector &to_index( const IndexFlat &index ) {
        return index_flat_to_index_vector_struct[index];
    }

    inline IndexRefVector &get_indices() {
        return index_flat_to_index_vector_struct;
    }

    inline size_t nonZeros() const {
        return values.size();
    }

    inline size_t size() const {
        return values.size() * ( sizeof( QDLC::Type::Scalar ) + sizeof( Index ) * index_flat_to_index_vector_struct.front().size() ); // / 1024. / 1024;
    }
};

} // namespace Numerics

} // namespace QDLC