#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

template <typename Scalar>
class FixedSizeSparseMap {
   private:
    Eigen::SparseMatrix<Scalar> cache;
    std::vector<int> dimensions;
    std::vector<int> dimensions_scaled;
    int sparse_matrix_dimension;
    std::vector<Eigen::Triplet<Scalar>> triplets;

   public:
    FixedSizeSparseMap(){};
    FixedSizeSparseMap( const std::vector<int> &init_dimensions ) {
        sparse_matrix_dimension = 1;
        for ( int dim : init_dimensions ) {
            sparse_matrix_dimension *= dim;
            dimensions_scaled.emplace_back( sparse_matrix_dimension );
        }
        dimensions = init_dimensions;
        cache = Eigen::SparseMatrix<Scalar>( sparse_matrix_dimension, sparse_matrix_dimension );
    }

    Eigen::SparseMatrix<Scalar> &mat() {
        return cache;
    }

    std::vector<int> indexToIndices( int i ) {
        //std::vector<int> ret;
        //ret.reserve( dimensions.size() );
        //std::cout << "Input index = " << i << " -> ";
        //for ( int k = dimensions.size() - 1; k > 0; k-- ) {
        //    int current = std::floor( 1.0 * i / dimensions_scaled[k - 1] );
        //    ret.emplace_back( current );
        //    i -= current * dimensions_scaled[k];
        //    //std::cout << current << " ";
        //}
        //ret.emplace_back( i );
        ////std::cout << std::endl;
        //std::reverse( ret.begin(), ret.end() );
        auto ii = i;
        int i3 = std::floor( i / dimensions_scaled[2] );
        i -= i3 * dimensions_scaled[2];
        int i2 = std::floor( i / dimensions_scaled[1] );
        i -= i2 * dimensions_scaled[1];
        int i1 = std::floor( i / dimensions_scaled[0] );
        i -= i1 * dimensions_scaled[0];
        return {i, i1, i2, i3};
    }

    // For now, hardcode get(numvar)
    Scalar get( int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3 ) {
        // Convert index array to total matrix index
        int indexX = i0 + i1 * dimensions_scaled[0] + i2 * dimensions_scaled[1] + i3 * dimensions_scaled[2];
        int indexY = j0 + j1 * dimensions_scaled[0] + j2 * dimensions_scaled[1] + j3 * dimensions_scaled[2];
        assert( indexX < sparse_matrix_dimension && indexY < sparse_matrix_dimension && "Index mismatch" );
        // Return copy of cache entry
        return cache.coeff( indexX, indexY );
    }

    Scalar &getRef( int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3 ) { // Convert index array to total matrix index
        int indexX = i0 + i1 * dimensions_scaled[0] + i2 * dimensions_scaled[1] + i3 * dimensions_scaled[2];
        int indexY = j0 + j1 * dimensions_scaled[0] + j2 * dimensions_scaled[1] + j3 * dimensions_scaled[2];
        // Return reference of cache entry
        return cache.coeffRef( indexX, indexY );
    }

    void makeCompressed() {
        cache.makeCompressed();
    }

    void setFromTripletList( bool clear = true ) {
        cache.setFromTriplets( triplets.begin(), triplets.end() );
        if ( clear ) {
            clearTripletList();
        }
        makeCompressed();
    }
    void addTriplet( int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, const Scalar &value ) {
        int indexX = i0 + i1 * dimensions_scaled[0] + i2 * dimensions_scaled[1] + i3 * dimensions_scaled[2];
        int indexY = j0 + j1 * dimensions_scaled[0] + j2 * dimensions_scaled[1] + j3 * dimensions_scaled[2];
        triplets.emplace_back( indexX, indexY, value );
    }

    void addTripletTo( int i0, int i1, int i2, int i3, int j0, int j1, int j2, int j3, const Scalar &value, std::vector<Eigen::Triplet<Scalar>> &vec ) {
        int indexX = i0 + i1 * dimensions_scaled[0] + i2 * dimensions_scaled[1] + i3 * dimensions_scaled[2];
        int indexY = j0 + j1 * dimensions_scaled[0] + j2 * dimensions_scaled[1] + j3 * dimensions_scaled[2];
        vec.emplace_back( indexX, indexY, value );
    }

    std::vector<Eigen::Triplet<Scalar>> &getTriplets() {
        return triplets;
    }

    long unsigned int nonZeros() {
        return cache.nonZeros();
    }

    void clearTripletList() {
        triplets.clear();
    }

    double getFillRatio() {
        return nonZeros() / ( 1.0 * sparse_matrix_dimension * sparse_matrix_dimension );
    }
    int getSizeOfCache() {
        return nonZeros() * sizeof( Scalar );
    }
    long unsigned int getCacheDimensions() {
        return sparse_matrix_dimension;
    }

    //Scalar &operator() {}
};