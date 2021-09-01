#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ostream>
#include <set>
#include <unordered_map>
#include <map>
#include <functional>
#include <array>
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

//template <>
//struct std::less<Eigen::VectorXi> {
//    bool operator()( const Eigen::VectorXi &a, const Eigen::VectorXi &b ) const {
//        for ( size_t i = 0; i < a.size(); ++i ) {
//            if ( a[i] > b[i] ) return false;
//        }
//        return true;
//    }
//};

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

    static Eigen::VectorXi dimensions_scaled;
    //struct hash {
    //    static std::hash<Index> hasher;
    //    inline size_t operator()( const DoubleIndex &val ) const {
    //        size_t seed = 0;
    //        seed ^= hasher( val.x ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
    //        seed ^= hasher( val.y ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
    //        return seed;
    //    }
    //};

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
            return A.dot( dimensions_scaled ) < B.dot( dimensions_scaled );
        }
    };

    typedef std::unordered_map<Eigen::VectorXi, Scalar, vector_hash<Eigen::VectorXi>> SubTensorMapEigen;
    typedef std::unordered_map<Eigen::VectorXi, SubTensorMapEigen, vector_hash<Eigen::VectorXi>> TensorMapEigen;
    //typedef std::map<Eigen::VectorXi, Scalar, vector_compare<Eigen::VectorXi>> SubTensorMapEigen;
    //typedef std::map<Eigen::VectorXi, SubTensorMapEigen, vector_compare<Eigen::VectorXi>> TensorMapEigen;

   private:
    std::vector<int> dimensions;
    Index sparse_matrix_dimension;

    // Ordered map including SparseIndex;Value pairs.
    std::vector<TensorMapEigen> values;

    // Counter for which indices/values vector is active. Use at least 2 and shift current working/getting vector
    int current_value_vector = 0;
    int current_cache_vector = 1;

   public:
    FixedSizeSparseMap(){};
    FixedSizeSparseMap( const std::vector<int> &init_dimensions, int max_threads = 1 ) {
        sparse_matrix_dimension = 1;
        dimensions_scaled = Eigen::VectorXi::Zero( init_dimensions.size() );
        // X
        dimensions_scaled( 0 ) = 1;
        for ( int i = 0; i < init_dimensions.size(); i++ ) {
            sparse_matrix_dimension *= init_dimensions.at( i );
            if ( i + 1 < init_dimensions.size() )
                dimensions_scaled( i + 1 ) = sparse_matrix_dimension;
        }
        //vector_compare<Eigen::VectorXi>::ds = dimensions_scaled;
        dimensions = init_dimensions;
        std::cout << "Tensor Size: ( ";
        std::copy( std::begin( init_dimensions ), std::end( init_dimensions ), std::ostream_iterator<int>( std::cout, ", " ) );
        std::cout << ") - " << (long long)( sparse_matrix_dimension * sparse_matrix_dimension ) << " total." << std::endl;
        for ( int i = 0; i < 2; i++ ) {
            //auto filler = std::vector<std::map<DoubleIndex, SparseMapTriplet, hash>>( max_threads );
            auto filler = TensorMapEigen();
            values.emplace_back( filler );
        }
    }

    //std::vector<std::map<DoubleIndex, SparseMapTriplet, hash>> &get() {
    TensorMapEigen &get() {
        return values[current_value_vector];
    }

    size_t getDim() {
        return dimensions.front();
    }

    void addTriplet( Eigen::VectorXi indicesX, Eigen::VectorXi indicesY, const Scalar &value, int thread, const int i_n = -1, const int j_n = -1, const int gi_n = -1, const int gj_n = -1 ) {
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
#pragma omp critical
        {
            auto &curX = values[current_cache_vector][indicesX];
            if ( curX.count( indicesY ) > 0 ) {
                curX[indicesY] += value;
            } else {
                curX[indicesY] = value;
            }
        }
        //values[current_cache_vector][sparse_index_x][sparse_index_y].addWithCheckForExistence( value, indicesX, indicesY, sparse_index_x, sparse_index_y );
    }

    long long int nonZeros() {
        long long s = 0;
        for ( auto &[indxX, outer] : values[current_value_vector] )
            s += outer.size();
        return s;
    }

    std::vector<std::pair<Eigen::VectorXi, SubTensorMapEigen>> getIndices() {
        std::vector<std::pair<Eigen::VectorXi, SubTensorMapEigen>> ret( values[current_value_vector].size() );
        std::move( values[current_value_vector].begin(), values[current_value_vector].end(), ret.begin() );
        return ret;
    }

    //TODO: damit mehrere threads worken können müsste man heir halt alles in eine map schreiben... sonst arbeiten verschiedene threads an gleichen summanden
    void reduceDublicates() {
        std::swap( current_value_vector, current_cache_vector );
        values[current_cache_vector].clear();
        //FIXME
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
    int size() {
        return nonZeros() * sizeof( Scalar );
    }
};