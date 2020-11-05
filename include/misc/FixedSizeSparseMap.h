#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

template <typename Scalar>
class FixedSizeSparseMap {
   private:
    Eigen::SparseMatrix<Scalar> cache;
    std::vector<int> dimensions;
    std::vector<int> dimensions_scaled;
    long long int sparse_matrix_dimension;
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

    std::vector<long long int> indexToIndices( long long int i ) {
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

        std::vector<long long int> ret;
        ret.reserve( dimensions_scaled.size() );
        for ( int k = dimensions_scaled.size() - 2; k >= 0; k-- ) {
            auto val = std::floor( i / dimensions_scaled[k] );
            ret.emplace_back( val );
            i -= val * dimensions_scaled[k];
        }
        ret.emplace_back( i );
        std::reverse( ret.begin(), ret.end() );
        return ret;

        //auto ii = i;
        //int i3 = std::floor( i / dimensions_scaled[2] );
        //i -= i3 * dimensions_scaled[2];
        //int i2 = std::floor( i / dimensions_scaled[1] );
        //i -= i2 * dimensions_scaled[1];
        //int i1 = std::floor( i / dimensions_scaled[0] );
        //i -= i1 * dimensions_scaled[0];
        //return {i, i1, i2, i3};
    }

    std::vector<long long int> shiftVector( const Scalar &i, std::vector<long long int> &input ) {
        std::vector<long long int> ret{i};
        ret.insert( ret.begin(), input.begin(), input.end() - 1 );
        std::cout << "input = ";
        for ( auto &a : input )
            std::cout << a << " ";
        std::cout << " -> Output = ";
        for ( auto &a : ret )
            std::cout << a << " ";
        std::cout << std::endl;
        exit( 0 );
        return ret;
    }

    // For now, hardcode get(numvar)template<typename... IndexTypes>
    //inline const Scalar &operator()( Index firstIndex, Index secondIndex, IndexTypes... otherIndices ) const {
    //    // The number of indices used to access a tensor coefficient must be equal to the rank of the tensor.
    //    EIGEN_STATIC_ASSERT( sizeof...( otherIndices ) + 2 == NumIndices, YOU_MADE_A_PROGRAMMING_MISTAKE )
    //    return this->operator()( array<Index, NumIndices>{{firstIndex, secondIndex, otherIndices...}} );
    //}
    Scalar get( const std::vector<long long int> &indicesX, const std::vector<long long int> &indicesY ) {
        // Convert index array to total matrix index const Index index = i4 + m_dimensions[4] * (i3 + m_dimensions[3] * (i2 + m_dimensions[2] * (i1 + m_dimensions[1] * i0)));
        int indexX = indicesX[0]; //i0 + i1 * dimensions_scaled[0] + i2 * dimensions_scaled[1] + i3 * dimensions_scaled[2];
        int indexY = indicesY[0]; //j0 + j1 * dimensions_scaled[0] + j2 * dimensions_scaled[1] + j3 * dimensions_scaled[2];
        for ( int i = 1; i < indicesX.size(); i++ ) {
            indexX += indicesX[i] * dimensions_scaled[i - 1];
            indexY += indicesY[i] * dimensions_scaled[i - 1];
        }
        assert( indexX < sparse_matrix_dimension && indexY < sparse_matrix_dimension && "Index mismatch" );
        // Return copy of cache entry
        return cache.coeff( indexX, indexY );
    }

    Scalar &getRef( const std::vector<long long int> &indicesX, const std::vector<long long int> &indicesY ) {
        // Convert index array to total matrix index
        int indexX = indicesX[0];
        int indexY = indicesY[0];
        for ( int i = 1; i < indicesX.size(); i++ ) {
            indexX += indicesX[i] * dimensions_scaled[i - 1];
            indexY += indicesY[i] * dimensions_scaled[i - 1];
        }
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
    void addTriplet( const std::vector<long long int> &indicesX, const std::vector<long long int> &indicesY, const Scalar &value ) {
        // Convert index array to total matrix index
        int indexX = indicesX[0];
        int indexY = indicesY[0];
        for ( int i = 1; i < indicesX.size(); i++ ) {
            indexX += indicesX[i] * dimensions_scaled[i - 1];
            indexY += indicesY[i] * dimensions_scaled[i - 1];
        }
        triplets.emplace_back( indexX, indexY, value );
    }

    void addTripletTo( const std::vector<long long int> &indicesX, const std::vector<long long int> &indicesY, const Scalar &value, std::vector<Eigen::Triplet<Scalar>> &vec ) {
        // Convert index array to total matrix index
        int indexX = indicesX[0];
        int indexY = indicesY[0];
        for ( int i = 1; i < indicesX.size(); i++ ) {
            indexX += indicesX[i] * dimensions_scaled[i - 1];
            indexY += indicesY[i] * dimensions_scaled[i - 1];
        }
        vec.emplace_back( indexX, indexY, value );
    }

    std::vector<Eigen::Triplet<Scalar>> &getTriplets() {
        return triplets;
    }

    long long int nonZeros() {
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
    long long int getCacheDimensions() {
        return sparse_matrix_dimension;
    }

    //Scalar &operator() {}
};