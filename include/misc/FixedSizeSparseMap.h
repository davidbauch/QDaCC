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
        static int128 getSparseIndex( const std::vector<int> &indicesX, const std::vector<int> &indicesY, const std::vector<int128> &dimensions_scaled ) {
            int128 sparse_index = indicesX[0];
            for ( int i = 1; i < indicesX.size(); i++ ) {
                sparse_index += indicesX[i] * dimensions_scaled[i - 1];
            }
            for ( int i = 0; i < indicesY.size(); i++ ) {
                sparse_index += indicesY[i] * dimensions_scaled[indicesX.size() + i];
            }
            return sparse_index;
        }
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
        // Sparse Index Comparisons
        bool operator>( const int128 &other ) const {
            return sparse_index > other;
        }
        bool operator>=( const int128 &other ) const {
            return sparse_index >= other;
        }
        bool operator<( const int128 &other ) const {
            return sparse_index < other;
        }
        bool operator<=( const int128 &other ) const {
            return sparse_index <= other;
        }
        // Scalar equal operator
        bool operator==( const Scalar &other ) const {
            return value == other;
        }
    };

   private:
    std::vector<std::vector<SparseMapTriplet>> value_triplets;
    std::vector<std::vector<SparseMapTriplet>> cache_triplets;
    std::vector<int> dimensions;
    std::vector<int128> dimensions_scaled;
    int128 sparse_matrix_dimension;

   public:
    FixedSizeSparseMap(){};
    FixedSizeSparseMap( const std::vector<int> &init_dimensions, int max_threads = 1 ) {
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
        for ( int thread = 0; thread < max_threads; thread++ ) {
            value_triplets.emplace_back( std::vector<SparseMapTriplet>() );
            cache_triplets.emplace_back( std::vector<SparseMapTriplet>() );
        }
    }

    std::vector<std::vector<SparseMapTriplet>> &get() {
        return value_triplets;
    }

    void addTriplet( std::vector<int> indicesX, std::vector<int> indicesY, const Scalar &value, int thread, const long long int i_n = -1, const long long int j_n = -1, bool sortedInsertion = false ) {
        if ( i_n != -1 && j_n != -1 ) {
            for ( int i = indicesX.size() - 1; i > 0; i-- ) {
                indicesX[i] = indicesX[i - 1];
                indicesY[i] = indicesY[i - 1];
            }
            indicesX[0] = i_n;
            indicesY[0] = j_n;
        }
        if ( sortedInsertion ) {
            auto sparse_index = SparseMapTriplet::getSparseIndex( indicesX, indicesY, dimensions_scaled );
            // Find sparse_index
            auto lower_bound = std::lower_bound( cache_triplets[thread].begin(), cache_triplets[thread].end(), sparse_index );
            // Not found, insert
            if ( lower_bound == cache_triplets[thread].end() ) {
                cache_triplets[thread].emplace_back( value, indicesX, indicesY, dimensions_scaled );
                return;
            }
            // Check if found index is equal to the to-inserted index, if yes add them, if no insert new element here.
            auto index = lower_bound - cache_triplets[thread].begin();
            // Equal indices
            if ( cache_triplets[thread][index] == sparse_index ) {
                cache_triplets[thread][index].value += value;
                return;
            } else {
                cache_triplets[thread].insert( lower_bound, {value, indicesX, indicesY, dimensions_scaled} );
            }
        }
        // Insert element
        cache_triplets[thread].emplace_back( value, indicesX, indicesY, dimensions_scaled );
    }

    long long int nonZeros() {
        long long int nonz = 0;
        for ( int i = 0; i < value_triplets.size(); i++ ) {
            nonz += value_triplets[i].size();
        }
        return nonz;
    }

    void reduceDublicates( bool clear_equals_setzero = false ) {
        int threads = cache_triplets.size();
#pragma omp parallel for num_threads( threads )
        for ( int thread = 0; thread < threads; thread++ ) {
            value_triplets[thread].clear();
            std::sort( cache_triplets[thread].begin(), cache_triplets[thread].end() );
            value_triplets[thread].emplace_back( cache_triplets[thread].front() );
            // If a value is still zero (after adm advancing), then just remove it. It can be re-added next iteration
            if ( clear_equals_setzero ) {
                cache_triplets[thread].erase( std::remove( cache_triplets[thread].begin(), cache_triplets[thread].end(), 0.0 ), cache_triplets[thread].end() );
            }
            // Add all dublicates together, or add new element if element is unique
            for ( int128 i = 1; i < cache_triplets[thread].size(); i++ ) {
                if ( cache_triplets[thread][i] == value_triplets[thread].back() ) {
                    value_triplets[thread].back().value += cache_triplets[thread][i].value;
                } else {
                    value_triplets[thread].emplace_back( cache_triplets[thread][i] );
                }
                if ( clear_equals_setzero ) {
                    cache_triplets[thread][i].value = 0.0;
                }
            }
            if ( !clear_equals_setzero ) {
                auto size = cache_triplets[thread].size();
                cache_triplets[thread].clear();
                cache_triplets[thread].reserve( size );
            }
        }
    }

    int getSizeOfTensor() {
        return nonZeros() * ( sizeof( Scalar ) + sizeof( int128 ) + sizeof( int ) * dimensions_scaled.size() );
    }
};