#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ostream>
#include <set>
#include <unordered_map>
#include <functional>
//using int128 = __uint128_t;
//std::ostream &operator<<( std::ostream &os, const int128 i ) noexcept {
//    std::ostream::sentry s( os );
//    if ( s ) {
//        unsigned __int128 tmp = i < 0 ? -i : i;
//        char buffer[128];
//        char *d = std::end( buffer );
//        do {
//            --d;
//            *d = "0123456789"[tmp % 10];
//            tmp /= 10;
//        } while ( tmp != 0 );
//        if ( i < 0 ) {
//            --d;
//            *d = '-';
//        }
//        int len = std::end( buffer ) - d;
//        if ( os.rdbuf()->sputn( d, len ) != len ) {
//            os.setstate( std::ios_base::badbit );
//        }
//    }
//    return os;
//}

template <typename Scalar>
class FixedSizeSparseMap {
   public:
    // Index Type definition
    typedef long long Index;
    // Overwrite << stream operator
    friend std::ostream &operator<<( std::ostream &os, const Index i ) noexcept {
        os << (long long)i;
        return os;
    }

    class DoubleIndex {
       public:
        Index x, y;
        DoubleIndex() : x( 0 ), y( 0 ){};
        DoubleIndex( const DoubleIndex &other ) : x( other.x ), y( other.y ){};
        DoubleIndex( const Index &_x, const Index &_y ) : x( _x ), y( _y ){};
        bool operator==( const DoubleIndex &other ) const {
            return x == other.x && y == other.y;
        }
    };

    struct hash {
        static std::hash<Index> hasher;
        inline size_t operator()( const DoubleIndex &val ) const {
            size_t seed = 0;
            seed ^= hasher( val.x ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
            seed ^= hasher( val.y ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
            return seed;
        }
    };

    class SparseMapTriplet {
       public:
        std::vector<int> indicesX, indicesY;
        DoubleIndex sparse_index;
        Scalar value;
        bool created_zero = false;
        static DoubleIndex getSparseIndex( const std::vector<int> &indicesX, const std::vector<int> &indicesY, const std::vector<Index> &dimensions_scaled ) {
            DoubleIndex sparse_index( indicesX[0], indicesY[0] );
            for ( int i = 1; i < indicesX.size(); i++ ) {
                sparse_index.x += indicesX[i] * dimensions_scaled[i - 1];
                sparse_index.y += indicesY[i] * dimensions_scaled[i - 1];
            }
            return sparse_index;
        }
        SparseMapTriplet() : value( 0 ), created_zero( true ){};
        SparseMapTriplet( const SparseMapTriplet &other ) : indicesX( other.indicesX ), indicesY( other.indicesY ), sparse_index( other.sparse_index ), value( other.value ){};
        SparseMapTriplet( const Scalar &value, const std::vector<int> &indicesX, const std::vector<int> &indicesY, const std::vector<Index> &dimensions_scaled ) : value( value ), indicesX( indicesX ), indicesY( indicesY ) {
            sparse_index = getSparseIndex( indicesX, indicesY, dimensions_scaled );
        }
        SparseMapTriplet( const Scalar &value, const std::vector<int> &indicesX, const std::vector<int> &indicesY, DoubleIndex sparse_index, const std::vector<Index> &dimensions_scaled ) : value( value ), indicesX( indicesX ), indicesY( indicesY ), sparse_index( sparse_index ){};
        void addWithCheckForExistence( const Scalar &_value, const std::vector<int> &_indicesX, const std::vector<int> &_indicesY, DoubleIndex _sparse_index ) {
            value += _value;
            if ( created_zero ) {
                indicesX = _indicesX;
                indicesY = _indicesY;
                sparse_index = _sparse_index;
                created_zero = false;
            }
        }
        void addWithCheckForExistence( const SparseMapTriplet &other ) {
            value += other.value;
            if ( created_zero ) {
                indicesX = other.indicesX;
                indicesY = other.indicesY;
                sparse_index = other.sparse_index;
                created_zero = false;
            }
        }
        bool operator==( const SparseMapTriplet &other ) const {
            return sparse_index == other.sparse_index;
        }
    };

   private:
    std::vector<int> dimensions;
    std::vector<Index> dimensions_scaled;
    Index sparse_matrix_dimension;

    // Unordered map including SparseIndex;Value pairs.
    std::vector<std::vector<std::unordered_map<DoubleIndex, SparseMapTriplet, hash>>> values;
    // Counter for which indices/values vector is active. Use at least 2 and shift current working/getting vector
    int current_value_vector = 0;
    int current_cache_vector = 1;

   public:
    FixedSizeSparseMap(){};
    FixedSizeSparseMap( const std::vector<int> &init_dimensions, int max_threads = 1 ) {
        sparse_matrix_dimension = 1;
        // X
        for ( int dim : init_dimensions ) {
            sparse_matrix_dimension *= dim;
            dimensions_scaled.emplace_back( sparse_matrix_dimension );
        }
        dimensions = init_dimensions;
        std::cout << "Sparse Elements: " << (long long)( sparse_matrix_dimension ) << std::endl;
        for ( int i = 0; i < 2; i++ ) {
            auto filler = std::vector<std::unordered_map<DoubleIndex, SparseMapTriplet, hash>>( max_threads );
            values.emplace_back( filler );
        }
    }

    std::vector<std::unordered_map<DoubleIndex, SparseMapTriplet, hash>> &get() {
        return values[current_value_vector];
    }

    void addTriplet( std::vector<int> indicesX, std::vector<int> indicesY, const Scalar &value, int thread, const long long int i_n = -1, const long long int j_n = -1 ) {
        if ( i_n != -1 && j_n != -1 ) {
            for ( int i = indicesX.size() - 1; i > 0; i-- ) {
                indicesX[i] = indicesX[i - 1];
                indicesY[i] = indicesY[i - 1];
            }
            indicesX[0] = i_n;
            indicesY[0] = j_n;
        }
        auto sparse_index = SparseMapTriplet::getSparseIndex( indicesX, indicesY, dimensions_scaled );
        values[current_cache_vector][thread][sparse_index].addWithCheckForExistence( value, indicesX, indicesY, sparse_index );
    }

    long long int nonZeros() {
        long long s = values[current_value_vector][0].size();
        for ( int i = 1; i < values[current_value_vector].size(); i++ )
            s += values[current_value_vector][i].size();
        return s;
    }

    void reduceDublicates() {
        std::swap( current_value_vector, current_cache_vector );
        for ( auto &vec : values[current_cache_vector] )
            vec.clear();
        //if ( values[current_value_vector].size() > 1 ) {
        //    int current_thread_vector = 0;
        //    // Merge all vectors into vector 0
        //    for ( int cv = 0; cv < values[current_value_vector].size(); cv++ ) {
        //        for ( auto &keyv : values[current_value_vector][cv] ) {
        //            values[current_cache_vector][0][keyv.first].addWithCheckForExistence( keyv.second );
        //        }
        //        values[current_value_vector][cv].clear();
        //    }
        //    // Redistribute by moving each element to current_thread_vector, then incrementing current_thread_vector
        //    for ( auto &keyv : values[current_cache_vector][0] ) {
        //        values[current_value_vector][current_thread_vector][keyv.first].addWithCheckForExistence( keyv.second );
        //        current_thread_vector = ( current_thread_vector + 1 ) % values[current_value_vector].size();
        //    }
        //    values[current_cache_vector][0].clear();
        //}
    }

    int getSizeOfTensor() {
        return nonZeros() * ( sizeof( Scalar ) + sizeof( Index ) + sizeof( int ) * dimensions_scaled.size() );
    }
    std::string getSizesOfCache() {
        std::string cachevs = "";
        return cachevs;
    }
};