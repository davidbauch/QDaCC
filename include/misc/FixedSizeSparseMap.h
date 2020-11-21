#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

using int128 = __uint128_t;

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
    class SparseMapTriplet {
       public:
        std::vector<int> indicesX, indicesY;
        int128 sparse_index;
        Scalar value;
        SparseMapTriplet(){};
        SparseMapTriplet( const SparseMapTriplet &other ) : indicesX( other.indicesX ), indicesY( other.indicesY ), sparse_index( other.sparse_index ), value( other.value ){};
        SparseMapTriplet( const Scalar &value, const std::vector<int> &indicesX, const std::vector<int> &indicesY, const std::vector<int128> &dimensions_scaled ) : value( value ), indicesX( indicesX ), indicesY( indicesY ) {
            //if ( std::real( value ) != 0 || std::imag( value ) != 0 ) {
            sparse_index = indicesX[0];
            for ( int i = 1; i < indicesX.size(); i++ ) {
                sparse_index += indicesX[i] * dimensions_scaled[i - 1];
            }
            for ( int i = 0; i < indicesY.size(); i++ ) {
                sparse_index += indicesY[i] * dimensions_scaled[indicesX.size() + i];
            }
            //}
        }
        bool operator>( const SparseMapTriplet &other ) const {
            return sparse_index > other.sparse_index;
        }
        bool operator>=( const SparseMapTriplet &other ) const {
            return sparse_index >= other.sparse_index;
        }
        bool operator<( const SparseMapTriplet &other ) const {
            return sparse_index < other.sparse_index;
        }
        bool operator<=( const SparseMapTriplet &other ) const {
            return sparse_index <= other.sparse_index;
        }
        bool operator==( const SparseMapTriplet &other ) const {
            return sparse_index == other.sparse_index;
            //bool first = true;
            //bool second = true;
            //if ( sparse_index != other.sparse_index ) first = false;
            //for ( int i = 0; i < indicesX.size(); i++ ) {
            //    if ( indicesX[i] != other.indicesX[i] ) {
            //        second = false;
            //        if ( first != second ) std::cout << "FUCK\n";
            //        return false;
            //    }
            //    if ( indicesY[i] != other.indicesY[i] ) {
            //        second = false;
            //        if ( first != second ) std::cout << "FUCK\n";
            //        return false;
            //    }
            //}
            //if ( first != second ) std::cout << "FUCK\n";
            //return true;
        }
    };

   private:
    std::vector<SparseMapTriplet> cache;
    std::vector<SparseMapTriplet> triplets;
    std::vector<int> dimensions;
    std::vector<int128> dimensions_scaled;
    int128 sparse_matrix_dimension;

   public:
    FixedSizeSparseMap(){};
    FixedSizeSparseMap( const std::vector<int> &init_dimensions ) {
        sparse_matrix_dimension = 1;
        // X
        for ( int dim : init_dimensions ) {
            sparse_matrix_dimension *= dim;
            dimensions_scaled.emplace_back( sparse_matrix_dimension );
        }
        // Y
        for ( int dim : init_dimensions ) {
            sparse_matrix_dimension *= dim;
            dimensions_scaled.emplace_back( sparse_matrix_dimension );
        }
        std::cout << "Sparse Elements: " << (long long)( sparse_matrix_dimension ) << std::endl;
        dimensions = init_dimensions;
    }

    std::vector<SparseMapTriplet> &get() {
        return cache;
    }

    void addTriplet( const std::vector<int> &indicesX, const std::vector<int> &indicesY, const Scalar &value, const int i_n = -1, const int j_n = -1 ) {
        // Convert index array to total matrix index
        if ( i_n != -1 && j_n != -1 ) {
            auto newIndexX = indicesX;
            auto newIndexY = indicesY;
            newIndexX[0] = i_n;
            newIndexY[0] = j_n;
            for ( int i = 1; i < newIndexX.size(); i++ ) {
                newIndexX[i] = indicesX[i - 1];
                newIndexY[i] = indicesY[i - 1];
            }
            triplets.emplace_back( value, newIndexX, newIndexY, dimensions_scaled );
        } else {
            triplets.emplace_back( value, indicesX, indicesY, dimensions_scaled );
        }
    }

    void addTripletTo( const std::vector<int> &indicesX, const std::vector<int> &indicesY, const Scalar &value, std::vector<SparseMapTriplet> &vec, const long long int i_n = -1, const long long int j_n = -1 ) {
        if ( i_n != -1 && j_n != -1 ) {
            auto newIndexX = indicesX;
            auto newIndexY = indicesY;
            newIndexX[0] = i_n;
            newIndexY[0] = j_n;
            for ( int i = 1; i < newIndexX.size(); i++ ) {
                newIndexX[i] = indicesX[i - 1];
                newIndexY[i] = indicesY[i - 1];
            }
            vec.emplace_back( value, newIndexX, newIndexY, dimensions_scaled );
        } else {
            vec.emplace_back( value, indicesX, indicesY, dimensions_scaled );
        }
    }

    std::vector<Eigen::Triplet<Scalar>> &getTriplets() {
        return triplets;
    }

    long long int nonZeros() {
        return cache.size();
    }

    void set() {
        //cache.clear();
        //cache = triplets;
        reduceDublicates( triplets );
        auto size = triplets.size();
        triplets.clear();
        triplets.reserve( size );
    }

    void reduceDublicates( std::vector<SparseMapTriplet> &triplets ) {
        int128 limit = triplets.size();
        cache.clear();
        std::sort( triplets.begin(), triplets.end() );
        cache.emplace_back( triplets.front() );

        for ( int128 i = 1; i < limit; i++ ) {
            if ( triplets[i] == cache.back() ) {
                cache.back().value += triplets[i].value;
            } else {
                cache.emplace_back( triplets[i] );
            }
        }
    }

    int getSizeOfCache() {
        return nonZeros() * sizeof( Scalar );
    }
};