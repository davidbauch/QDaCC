#pragma once

#include "typedef.h"
#include <ranges>

namespace QDACC::Numerics {

/**
 * This function is supposed to work like the CacheMatrix class, but for multidimensional matrices.
 * Instead of saving a 2D matrix by default, we save a nD matrix in a std::vector. get({i,j,...}) and set({i,j,...})
 * will then calculate a reduced index from the multidimensional index and use that to access the std::vector.
 */
class MultidimensionalCacheMatrix {
   public:
    typedef size_t Index;

   private:
    std::string m_name;
    std::vector<Scalar> m_matrix;
    Index max_number_of_values;
    using IndexVector = std::vector<size_t>;
    IndexVector m_reduction_map;
    bool fft = false;
    int order;

    /**
     * Given the reduction_map and an index vector, this function calculates the reduced index.
     * The reduction map is a vector of the form { dim[0], dim[1]*dim[0], dim[2]*dim[1]*dim[0], ...}
     * The index vector is a vector of the form { i[0], i[1], i[2], ...}
     * The reduced index is then calculated by summing { i[0]*reduction_map[0], ..., i[N]*reduction_map[N] }
     */
    inline Index reduceIndex( const IndexVector& vector_index ) const {
        // Calculate the reduced index by summing { vector_index[0]*reduction_map[0], ..., vector_index[N]*reduction_map[N] }
        const auto index = std::inner_product( vector_index.begin(), vector_index.end(), m_reduction_map.begin(), 0 );
        return index;
    }
    /**
     * Reverse function.
     */
    inline IndexVector expandIndex( Index index ) const {
        IndexVector vector_index( m_reduction_map.size(), 0 );
        for ( size_t i = m_reduction_map.size(); i > 0; --i ) {
            vector_index[i - 1] = index / m_reduction_map[i - 1];
            index -= vector_index[i - 1] * m_reduction_map[i - 1];
        }
        return vector_index;
    }

   public:
    MultidimensionalCacheMatrix() = default;
    MultidimensionalCacheMatrix( const std::vector<size_t> dimensions, const std::string& name, bool initialize = true, const Scalar initial_value = { 0, 0 } ) : m_name( name ) {
        order = dimensions.size();

        m_reduction_map.push_back( 1 ); // { dim[0], ...}
        for ( const size_t& el : dimensions | std::views::drop( 1 ) ) {
            m_reduction_map.push_back( m_reduction_map.back() * el ); // { dim[0], dim[1]*dim[0], dim[2]*dim[1]*dim[0], ...}
        }
        max_number_of_values = dimensions.front() * m_reduction_map.back();

        // We test the maximum size of this vector. If it is too large, we promt the user with a warning.
        // We do not throw an error, because large matrices may as well be intended, especially on HPC systems.
        // The threshold we choose is 10GB, which is feasable even for a consumer PC, and well in the range
        // of what a HPC system can handle and is equivalent to a 3D 880x880x880 Grid.
        const auto memory_size_mb = max_number_of_values * sizeof( Scalar ) / 1024 / 1024;
        if ( memory_size_mb > 1000 ) {
            Log::Warning( "[MultidimensionalCacheMatrix] Warning: You are trying to allocate {} MB of memory for the matrix {}. This may be intended, but is usually a mistake.\n", memory_size_mb, m_name );
        }

        m_matrix.reserve( max_number_of_values );
        // Initialize the vector if desired. Only initialize the "upper tringular matrix"
        if ( initialize ) {
            m_matrix = std::vector<Scalar>( max_number_of_values, initial_value );
        }

        Log::L2( "[MultidimensionalCacheMatrix] Created matrix {} with {} elements ({} MB).\n", m_name, max_number_of_values, memory_size_mb );
    }

    Scalar& get( const std::vector<size_t>& current_vector_index ) {
        const auto index = reduceIndex( current_vector_index );
        return m_matrix[index];
    }
    Scalar get( const std::vector<size_t>& current_vector_index ) const {
        const auto index = reduceIndex( current_vector_index );
        return m_matrix.at(index);
    }

    Scalar at( const std::vector<size_t>& current_vector_index ) {
        const auto index = reduceIndex( current_vector_index );
        return m_matrix.at( index );
    }

    void set( const std::vector<size_t>& current_vector_index, const Scalar& value ) {
        const auto index = reduceIndex( current_vector_index );
        m_matrix[index] = value;
    }

    // Special Get and Set methods for 2D case, which is G1 and G2:
    Scalar& get( const size_t row_index, const size_t col_index ) {
        const auto index = reduceIndex( { row_index, col_index } );
        return m_matrix[index];
    }
    Scalar get( const size_t row_index, const size_t col_index ) const {
        const auto index = reduceIndex( { row_index, col_index } );
        return m_matrix.at( index );
    }
    Scalar at( const size_t row_index, const size_t col_index ) const {
        const auto index = reduceIndex( { row_index, col_index } );
        return m_matrix.at( index );
    }
    void set( const size_t row_index, const size_t col_index, const Scalar& value ) {
        const auto index = reduceIndex( { row_index, col_index } );
        m_matrix[index] = value;
    }

    // Single size_t Get and Set methods without calling reduceIndex
    Scalar& get( const size_t index ) { return m_matrix[index]; }
    Scalar get( const size_t index ) const { return m_matrix[index]; }
    IndexVector get_index( const size_t index ) const { return expandIndex( index ); }

    // Some helper functions
    inline size_t dim() const { return m_reduction_map.at( 1 ); }
    inline Index size() const { return m_matrix.size(); }

    inline std::string get_name() const { return m_name; }
    inline Scalar maxCoeff() const {
        return *std::ranges::max_element( m_matrix, []( const Scalar& a, const Scalar& b ) { return std::real( a ) < std::real( b ); } );
    }
    inline Scalar minCoeff() const {
        return *std::ranges::min_element( m_matrix, []( const Scalar& a, const Scalar& b ) { return std::real( a ) < std::real( b ); } );
    }
    inline Scalar maxAbsCoeff() const {
        return *std::ranges::max_element( m_matrix, []( const Scalar& a, const Scalar& b ) { return std::abs( a ) < std::abs( b ); } );
    }
    inline Scalar minAbsCoeff() const {
        return *std::ranges::min_element( m_matrix, []( const Scalar& a, const Scalar& b ) { return std::abs( a ) < std::abs( b ); } );
    }
    inline Scalar sum() const { return std::accumulate( m_matrix.begin(), m_matrix.end(), Scalar{ 0, 0 } ); }

    // If sum(index) > dim, we are outside the triangular matrix
    inline bool insideTriangular( const IndexVector index ) const { return std::accumulate( index.begin(), index.end(), 0 ) < dim(); }

    // Math Helper functions
    template <class callable>
    void applyFunction( callable&& func ) {
        std::transform( m_matrix.begin(), m_matrix.end(), m_matrix.begin(), func );
    }
    template <class callable>
    void transformFunction( callable&& func, const MultidimensionalCacheMatrix& other ) {
        if ( other.size() != size() ) {
            Log::Error( "[MultidimensionalCacheMatrix] Error: You are trying to combine two matrices of different size. This is not allowed and will crash the program.\n" );
            throw std::runtime_error( "MultidimensionalCacheMatrix: Matrix size mismatch" );
        }
#pragma omp parallel for
        for ( size_t i = 0; i < size(); ++i ) {
            m_matrix[i] = func( m_matrix[i], other.get( i ) );
        }
    }
    void addScalar( const Scalar to_add ) {
        applyFunction( [to_add]( const Scalar& a ) { return a + to_add; } );
    }
    void multiplyScalar( const Scalar to_multiply ) {
        applyFunction( [to_multiply]( const Scalar& a ) { return a * to_multiply; } );
    }
    void addMatrix( const MultidimensionalCacheMatrix& other_matrix ) {
        transformFunction( []( const Scalar& a, const Scalar& b ) { return a + b; }, other_matrix );
    }
    void multiplyMatrix( const MultidimensionalCacheMatrix& other_matrix ) {
        transformFunction( []( const Scalar& a, const Scalar& b ) { return a * b; }, other_matrix );
    }

    // FFT tracker
    bool hasBeenFourierTransformed() {
        if ( not fft ) {
            fft = true;
            return false;
        }
        return true;
    }

    // Order
    int getOrder() const { return order; }
};

} // namespace QDACC::Numerics