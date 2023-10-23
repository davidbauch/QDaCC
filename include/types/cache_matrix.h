#pragma once

#include <Eigen/Dense>
#include "typedef.h"
#include <ranges>

namespace QDACC::Numerics {

/**
 * @brief Matrix Wrapper Class to save time-value pairs
 * get() and set() functions for two-dimensional time and timedeltas are provided
 * DEPRECATED: Use MultidimensionalCacheMatrix instead
 */
class CacheMatrix {
   private:
    std::string m_name;
    Dense m_matrix;
    Dense m_time;
    bool fft = false;

   public:
    // No Copying. This object is ususally only copied by mistake
    CacheMatrix( const CacheMatrix& other ) = delete;
    CacheMatrix() = default;
    CacheMatrix( const Dense& matrix, const Dense& time, const std::string& name ) : m_matrix( matrix ), m_time( time ), m_name( name ){};
    CacheMatrix( const Dense::Index& dim, const std::string& name ) : m_name( name ) {
        Log::L2( "[CacheMatrix] Created Empty Cache Matrix {0} with dimensions {1}x{1}\n", name, dim );
        m_matrix = Dense::Zero( dim, dim );
        m_time = Dense::Zero( dim, dim );
    }

    /**
     * @brief Getter function for the matrix element. Also available through the () operator
     *
     * @param i Row
     * @param j Column
     * @return Dense&
     */
    Scalar& get( const Dense::Index i, const Dense::Index j ) {
        // No range check because Eigen does that for us
        return m_matrix( i, j );
    }
    const Scalar& get( const Dense::Index i, const Dense::Index j ) const {
        // No range check because Eigen does that for us
        return m_matrix( i, j );
    }
    Scalar& operator()( const Dense::Index i, const Dense::Index j ) {
        return get( i, j );
    }
    const Scalar& operator()( const Dense::Index i, const Dense::Index j ) const {
        return get( i, j );
    }
    Scalar& get_time( const Dense::Index i, const Dense::Index j ) {
        // No range check because Eigen does that for us
        return m_time( i, j );
    }
    const Scalar& get_time( const Dense::Index i, const Dense::Index j ) const {
        // No range check because Eigen does that for us
        return m_time( i, j );
    }
    /**
     * @brief Getter function for the whole matrix.
     *
     * @return Dense&
     */
    Dense& get() {
        return m_matrix;
    }
    const Dense& get() const {
        return m_matrix;
    }
    Dense& get_time() {
        return m_time;
    }
    const Dense& get_time() const {
        return m_time;
    }

    bool hasBeenFourierTransformed() {
        if (not fft) {
            fft = true;
            return false;
        }
        return true;
    }

    /**
     * @brief Getter for the t-time at a specific index
     *
     * @param index Row == Time
     * @param fixed Column
     * @return Scalar
     */
    double t( const Dense::Index& index, const Dense::Index& fixed = 0 ) const {
        const auto& el = m_time( index, fixed );
        return std::real( el );
    }
    double dt( const Dense::Index& index, const Dense::Index& fixed = 0 ) const {
        if ( index == 0 )
            return std::real( m_time( index + 1, fixed ) - m_time( index, fixed ) );
        return std::real( m_time( index, fixed ) - m_time( index - 1, fixed ) );
    }
    Dense t() {
        return m_time.real();
    }
    Dense t() const {
        return m_time.real();
    }
    /**
     * @brief Getter for the tau-time at a specific index
     *
     * @param index Col == Tau
     * @param fixed Row
     * @return Scalar
     */
    double tau( const Dense::Index& index, const Dense::Index& fixed = 0 ) const {
        const auto& el = m_time( fixed, index );
        return std::imag( el ) - std::real( el );
    }
    double dtau( const Dense::Index& index, const Dense::Index& fixed = 0 ) const {
        if ( index == 0 )
            return std::imag( m_time( fixed, index + 1 ) - m_time( fixed, index ) );
        return std::imag( m_time( fixed, index ) - m_time( fixed, index - 1 ) );
    }
    Dense tau() {
        return m_time.imag();
    }
    Dense tau() const {
        return m_time.imag();
    }

    /**
     * @brief Getter for name
     *
     * @return std::string&
     */
    std::string& get_name() {
        return m_name;
    }
    const std::string& get_name() const {
        return m_name;
    }
    /**
     * @brief Return Inner Matrix Dimensions. Matrices in this case are always square.
     *
     * @return Dense::Index
     */
    Dense::Index dim() const {
        return m_matrix.rows();
    }
    Dense empty() const {
        const auto rc = dim();
        return Dense::Zero( rc, rc );
    }
    void set_empty() {
        m_matrix = empty();
    }

    /**
     * @brief Returns this object with a conjugated matrix. This will leave the time matrix unchanged.
     *
     * @return CacheMatrix
     */
    CacheMatrix conjugate( const std::string& purpose = "" ) {
        std::string new_purpose = purpose.size() > 1 ? purpose : m_name;
        return { m_matrix.conjugate(), m_time, new_purpose };
    }
};

/**
* This function is supposed to work like the CacheMatrix class, but for multidimensional matrices.
* Instead of saving a 2D matrix by default, we save a nD matrix in a std::vector. get({i,j,...}) and set({i,j,...}) 
* will then calculate a reduced index from the multidimensional index and use that to access the std::vector. 
*/
class MultidimensionalCacheMatrix {
    public:
     typedef u_int64_t Index;
    private:
    std::string m_name;
    std::vector<Scalar> m_matrix;
    std::vector<Scalar> m_time;
    Index max_number_of_values;
    using IndexVector = std::vector<int>;
    IndexVector m_reduction_map;
    bool fft = false;
    int order;

    /**
     * Given the reduction_map and an index vector, this function calculates the reduced index.
     * The reduction map is a vector of the form { dim[0], dim[1]*dim[0], dim[2]*dim[1]*dim[0], ...}
     * The index vector is a vector of the form { i[0], i[1], i[2], ...}
     * The reduced index is then calculated by summing { i[0]*reduction_map[0], ..., i[N]*reduction_map[N] }
    */
    inline Index reduceIndex(const IndexVector& vector_index) const {
        // Calculate the reduced index by summing { vector_index[0]*reduction_map[0], ..., vector_index[N]*reduction_map[N] }
        const auto index = std::inner_product(vector_index.begin(), vector_index.end(), m_reduction_map.begin(), 0);
        return index;
    }
    /**
     * Reverse function.
    */
    inline IndexVector expandIndex(Index index) const {
        IndexVector vector_index(m_reduction_map.size(), 0);
        for (int i = m_reduction_map.size() - 1; i >= 0; --i) {
            vector_index[i] = index / m_reduction_map[i];
            index -= vector_index[i] * m_reduction_map[i];
        }
        return vector_index;
    }

    public:
    MultidimensionalCacheMatrix() = default;
    MultidimensionalCacheMatrix( const std::vector<int> dimensions, const std::string& name, bool initialize = true, const Scalar initial_value = {0,0} ) : m_name( name ) {
        order = dimensions.size();
        
        m_reduction_map.push_back( 1 ); // { dim[0], ...}
        for (const int& el : dimensions | std::views::drop(1)) {
            m_reduction_map.push_back( m_reduction_map.back() * el ); // { dim[0], dim[1]*dim[0], dim[2]*dim[1]*dim[0], ...}
        }
        max_number_of_values = dimensions.front() * m_reduction_map.back();
        
        // We test the maximum size of this vector. If it is too large, we promt the user with a warning.
        // We do not throw an error, because large matrices may as well be intended, especially on HPC systems.
        // The threshold we choose is 10GB, which is feasable even for a consumer PC, and well in the range
        // of what a HPC system can handle and is equivalent to a 3D 880x880x880 Grid.
        const auto memory_size_mb = max_number_of_values * sizeof( Scalar ) / 1024 / 1024;
        if ( memory_size_mb > 10000 ) {
            Log::Warning( "[MultidimensionalCacheMatrix] Warning: You are trying to allocate {0} MB of memory for the matrix {1}. This may be intended, but is usually a mistake.\n", memory_size_mb, m_name );
        }

        m_matrix.reserve( max_number_of_values );
        m_time.reserve( max_number_of_values );
        // Initialize the vector if desired. Only initialize the "upper tringular matrix"
        if (initialize) {
            m_matrix = std::vector<Scalar>( max_number_of_values, initial_value );
            m_time = std::vector<Scalar>( max_number_of_values, initial_value );
        }
    }

    Scalar& get( const std::vector<int>& current_vector_index ) {
        const auto index = reduceIndex(current_vector_index);
        return m_matrix[index];
    }
    Scalar& get_time( const std::vector<int>& current_vector_index ) {
        const auto index = reduceIndex(current_vector_index);
        return m_time[index];
    }
    void set( const std::vector<int>& current_vector_index, const Scalar& value ) {
        const auto index = reduceIndex(current_vector_index);
        m_matrix[index] = value;
    }
    void set_time( const std::vector<int>& current_vector_index, const Scalar& value ) {
        const auto index = reduceIndex(current_vector_index);
        m_time[index] = value;
    }

    // Special Get and Set methods for 2D case, which is G1 and G2:
    Scalar& get( const int row_index, const int col_index ) {
        const auto index = reduceIndex({row_index, col_index});
        return m_matrix[index];
    }
    Scalar& get_time( const int row_index, const int col_index ) {
        const auto index = reduceIndex({row_index, col_index});
        return m_time[index];
    }
    Scalar get( const int row_index, const int col_index ) const {
        const auto index = reduceIndex({row_index, col_index});
        return m_matrix[index];
    }
    Scalar get_time( const int row_index, const int col_index ) const {
        const auto index = reduceIndex({row_index, col_index});
        return m_time[index];
    }
    void set( const int row_index, const int col_index, const Scalar& value ) {
        const auto index = reduceIndex({row_index, col_index});
        m_matrix[index] = value;
    }
    void set_time( const int row_index, const int col_index, const Scalar& value ) {
        const auto index = reduceIndex({row_index, col_index});
        m_time[index] = value;
    }

    // Single int Get and Set methods without calling reduceIndex
    Scalar& get( const int index ) {
        return m_matrix[index];
    }
    Scalar& get_time( const int index ) {
        return m_time[index];
    }
    Scalar get( const int index ) const {
        return m_matrix[index];
    }
    Scalar get_time( const int index ) const {
        return m_time[index];
    }
    IndexVector get_index( const int index ) const {
        return expandIndex(index);
    }

    // Some helper functions
    inline int dim() const {
        return m_reduction_map.at(1);
    }
    inline Index size() const {
        return m_matrix.size();
    }
    inline Index size_time() const {
        return m_time.size();
    }
    inline std::string get_name() const {
        return m_name;
    }
    inline Scalar maxCoeff() const {
        return *std::ranges::max_element(m_matrix, [](const Scalar& a, const Scalar& b) { return std::real(a) < std::real(b); });
    }
    inline Scalar minCoeff() const {
        return *std::ranges::min_element(m_matrix, [](const Scalar& a, const Scalar& b) { return std::real(a) < std::real(b); });
    }
    // If sum(index) > dim, we are outside the triangular matrix
    inline bool insideTriangular(const IndexVector index) const {
        return std::accumulate(index.begin(), index.end(), 0) < dim();
    }

    // TODO: Move these to the system class
    // TODO: these dont work; fix them and remove m_time from the program.
    // maybe write a getTime and getDeltaTime function in the system class that just regenerates the correct time
    inline double getTimeOf(int index, const IndexVector& index_vector, const std::vector<double>& time_vector) {
        const int total_index = std::accumulate(index_vector.begin(), index_vector.begin() + index, 0);
        // Edge case, we are at the end of the time vector. here, we just extrapolate the time
        const auto final_index = total_index+index_vector[index];
        if (final_index >= time_vector.size())
            return NAN; // ... for now. after moving this function to the system class, it can just extrapolate the time
        return time_vector[total_index+index_vector[index]] - time_vector[total_index];
    }
    inline double getDeltaTimeOf(int index, const IndexVector& index_vector, const std::vector<double>& delta_time_vector) {
        int total_index = std::accumulate(index_vector.begin(), index_vector.end()-1, 0);
        total_index = std::max<int>(total_index , delta_time_vector.size()-1);
        return delta_time_vector[total_index];
    }

    // Math Helper functions
    template<class callable>
    void applyFunction(callable&& func) {
        std::transform(m_matrix.begin(), m_matrix.end(), m_matrix.begin(), func);
    }
    template <class callable>
    void transformFunction(callable&& func, const MultidimensionalCacheMatrix& other) {
        if (other.size() != size()) {
            Log::Error("[MultidimensionalCacheMatrix] Error: You are trying to combine two matrices of different size. This is not allowed and will crash the program.\n");
            throw std::runtime_error("MultidimensionalCacheMatrix: Matrix size mismatch");
        }
        #pragma omp parallel for
        for (int i = 0; i < size(); ++i) {
            m_matrix[i] = func(m_matrix[i], other.get(i));
        }
    }
    void addScalar(const Scalar to_add) {
        applyFunction([to_add](const Scalar& a) { return a + to_add; });
    }
    void multiplyScalar(const Scalar to_multiply) {
        applyFunction([to_multiply](const Scalar& a) { return a * to_multiply; });
    }
    void addMatrix(const MultidimensionalCacheMatrix& other_matrix) {
        transformFunction([](const Scalar& a, const Scalar& b) { return a + b; }, other_matrix);
    }
    void multiplyMatrix(const MultidimensionalCacheMatrix& other_matrix) {
        transformFunction([](const Scalar& a, const Scalar& b) { return a * b; }, other_matrix);
    }

    /**
     * @brief Getter for the t-time at a specific index
     *
     * @param index Row == Time
     * @param fixed Column
     * @return Scalar
     */
    double t( const int index, const int fixed = 0 ) const {
        const auto& el = get_time( index, fixed );
        return std::real( el );
    }
    double dt( const int index, const int fixed = 0 ) const {
        if ( index == 0 )
            return std::real( get_time( index + 1, fixed ) - get_time( index, fixed ) );
        return std::real( get_time( index, fixed ) - get_time( index - 1, fixed ) );
    }
    /**
     * @brief Getter for the tau-time at a specific index
     *
     * @param index Col == Tau
     * @param fixed Row
     * @return Scalar
     */
    double tau( const int index, const int fixed = 0 ) const {
        const auto& el = get_time( fixed, index );
        return std::imag( el ) - std::real( el );
    }
    double dtau( const int index, const int fixed = 0 ) const {
        if ( index == 0 )
            return std::imag( get_time( fixed, index + 1 ) - get_time( fixed, index ) );
        return std::imag( get_time( fixed, index ) - get_time( fixed, index - 1 ) );
    }

    // FFT tracker
    bool hasBeenFourierTransformed() {
        if (not fft) {
            fft = true;
            return false;
        }
        return true;
    }

    // Order
    int getOrder() const {
        return order;
    }
};

} // namespace QDACC::Numerics