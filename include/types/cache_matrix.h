#pragma once

#include <Eigen/Dense>
#include "typedef.h"

namespace QDLC::Numerics {

/**
 * @brief Matrix Wrapper Class to save time-value pairs
 * get() and set() functions for two-dimensional time and timedeltas are provided
 *
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

} // namespace QDLC::Numerics