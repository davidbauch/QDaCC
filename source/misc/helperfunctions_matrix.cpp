#include "misc/helperfunctions_matrix.h"

QDLC::Type::Dense QDLC::Matrix::dense_projector( const QDLC::Type::Dense &input ) {
    QDLC::Type::Dense ret = QDLC::Type::Dense::Zero( input.rows(), input.cols() );
    for ( int i = 0; i < ret.rows(); i++ ) {
        for ( int j = 0; j < ret.cols(); j++ ) {
            ret( i, j ) = 1.0;
        }
    }
    // Generate new Matrix from triplet list
    return ret;
}

QDLC::Type::Sparse QDLC::Matrix::sparse_projector( const QDLC::Type::Sparse &input ) {
    QDLC::Type::Sparse ret = QDLC::Type::Sparse( input.rows(), input.cols() );
    std::vector<QDLC::Type::Triplet> ret_v;
    for ( int k = 0; k < input.outerSize(); ++k ) {
        for ( QDLC::Type::Sparse::InnerIterator it( input, k ); it; ++it ) {
            ret_v.emplace_back( it.row(), it.col(), 1.0 );
        }
    }
    // Generate new Matrix from triplet list
    ret.setFromTriplets( ret_v.begin(), ret_v.end() );
    return ret;
}

// Calculates the tensor product of matrices a and b
// @param &a,&b: Input matrices
// @return Returns a x b where x is the tensor product
QDLC::Type::Dense QDLC::Matrix::tensor( const QDLC::Type::Dense &a, const QDLC::Type::Dense &b ) {
    assert( a.rows() == a.cols() && b.rows() == b.cols() && "Only Square matrices accepted" );
    QDLC::Type::Dense ret = QDLC::Type::Dense::Zero( a.cols() * b.cols(), a.rows() * b.rows() );
    for ( int i = 0; i < a.rows(); i++ )
        for ( int j = 0; j < a.cols(); j++ )
            for ( int k = 0; k < b.rows(); k++ )
                for ( int l = 0; l < b.cols(); l++ ) {
                    ret( i * b.rows() + k, j * b.cols() + l ) = a( i, j ) * b( k, l );
                }
    return ret;
}

// Chains the tensor products
QDLC::Type::Dense QDLC::Matrix::tensor( const std::vector<QDLC::Type::Dense> &m ) {
    if ( m.size() < 2 )
        return m.front();
    else if ( m.size() == 2 )
        return tensor( m[0], m[1] );
    auto md = m;
    md.pop_back();
    return tensor( tensor( md ), m.back() );
}

// Determines the resulting base of the tensor product of multiple matrices
// @param &a,&b input string vectors containing the named base
// @return Returns a string vector containing the combined basis vector names
std::vector<std::string> QDLC::Matrix::tensor( const std::vector<std::string> &a, const std::vector<std::string> &b ) {
    std::vector<std::string> ret;
    for ( int i = 0; i < (int)a.size(); i++ )
        for ( int k = 0; k < (int)b.size(); k++ ) {
            ret.emplace_back( a.at( i ) + ":" + b.at( k ) );
        }
    return ret;
}

// Chains the string tensor products
std::vector<std::string> QDLC::Matrix::tensor( const std::vector<std::vector<std::string>> &m ) {
    if ( m.size() < 2 )
        return m.front();
    else if ( m.size() == 2 )
        return tensor( m[0], m[1] );
    auto md = m;
    md.pop_back();
    return tensor( tensor( md ), m.back() );
}

// Meshgrid
// Creates a XY-tuple containing x- and y values of a meshgrid containing [N] values between [xmin],[ymin] and [xmax],[ymax].
// If the endpoints (xmax,ymax) are not to be included
std::pair<QDLC::Type::Dense, QDLC::Type::Dense> QDLC::Matrix::meshgrid( double x_min, double y_min, double x_max, double y_max, int N, bool include_endpoint ) {
    assert( x_min < x_max && y_min < y_max );
    QDLC::Type::Dense X = QDLC::Type::Dense::Zero( N, N );
    QDLC::Type::Dense Y = QDLC::Type::Dense::Zero( N, N );
    for ( int i = 0; i < N; i++ )
        for ( int j = 0; j < N; j++ ) {
            X( i, j ) = x_min + (double)i / (double)( N - ( include_endpoint ? 1 : 0 ) ) * ( x_max - x_min );
            Y( i, j ) = y_min + (double)j / (double)( N - ( include_endpoint ? 1 : 0 ) ) * ( y_max - y_min );
        }
    return std::make_pair( X, Y );
}