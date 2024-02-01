#include "misc/helperfunctions_matrix.h"

QDACC::Type::MatrixMain QDACC::Matrix::projector( const QDACC::Type::MatrixMain &input ) {
#ifdef USE_SPARSE_MATRIX
    return sparse_projector( input );
#else
    return dense_projector( input );
#endif
}


QDACC::Type::Dense QDACC::Matrix::dense_projector( const QDACC::Type::Dense &input ) {
    QDACC::Type::Dense ret = QDACC::Type::Dense::Zero( input.rows(), input.cols() );
    for ( int i = 0; i < ret.rows(); i++ ) {
        for ( int j = 0; j < ret.cols(); j++ ) {
            ret( i, j ) = 1.0;
        }
    }
    // Generate new Matrix from triplet list
    return ret;
}

QDACC::Type::Sparse QDACC::Matrix::sparse_projector( const QDACC::Type::Sparse &input ) {
    QDACC::Type::Sparse ret = QDACC::Type::Sparse( input.rows(), input.cols() );
    std::vector<QDACC::Type::Triplet> ret_v;
    for ( int k = 0; k < input.outerSize(); ++k ) {
        for ( QDACC::Type::Sparse::InnerIterator it( input, k ); it; ++it ) {
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
QDACC::Type::MatrixMain QDACC::Matrix::tensor( const QDACC::Type::MatrixMain &a, const QDACC::Type::MatrixMain &b ) {
    assert( a.rows() == a.cols() && b.rows() == b.cols() && "Only Square matrices accepted" );
    QDACC::Type::MatrixMain ret( a.cols() * b.cols(), a.rows() * b.rows() );
    ret.setZero();
    for ( int i = 0; i < a.rows(); i++ )
        for ( int j = 0; j < a.cols(); j++ )
            for ( int k = 0; k < b.rows(); k++ )
                for ( int l = 0; l < b.cols(); l++ ) {
                    ret.coeffRef( i * b.rows() + k, j * b.cols() + l ) = a.coeff( i, j ) * b.coeff( k, l );
                }
    return ret;
}

// Chains the tensor products
QDACC::Type::MatrixMain QDACC::Matrix::tensor( const std::vector<QDACC::Type::MatrixMain> &m ) {
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
std::vector<std::string> QDACC::Matrix::tensor( const std::vector<std::string> &a, const std::vector<std::string> &b ) {
    std::vector<std::string> ret;
    for ( int i = 0; i < (int)a.size(); i++ )
        for ( int k = 0; k < (int)b.size(); k++ ) {
            ret.emplace_back( a.at( i ) + ":" + b.at( k ) );
        }
    return ret;
}

// Chains the string tensor products
std::vector<std::string> QDACC::Matrix::tensor( const std::vector<std::vector<std::string>> &m ) {
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
std::pair<QDACC::Type::Dense, QDACC::Type::Dense> QDACC::Matrix::meshgrid( double x_min, double y_min, double x_max, double y_max, int N, bool include_endpoint ) {
    assert( x_min < x_max && y_min < y_max );
    QDACC::Type::Dense X = QDACC::Type::Dense::Zero( N, N );
    QDACC::Type::Dense Y = QDACC::Type::Dense::Zero( N, N );
    for ( int i = 0; i < N; i++ )
        for ( int j = 0; j < N; j++ ) {
            X( i, j ) = x_min + (double)i / (double)( N - ( include_endpoint ? 1 : 0 ) ) * ( x_max - x_min );
            Y( i, j ) = y_min + (double)j / (double)( N - ( include_endpoint ? 1 : 0 ) ) * ( y_max - y_min );
        }
    return std::make_pair( X, Y );
}