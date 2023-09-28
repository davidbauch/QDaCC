#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ostream>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>
#include <functional>
#include <array>

#include "typedef.h"
#include "misc/log.h"
#include "solver/solver_hash.h"

// Das hier ist rekursiv. TODO: solver_tensor_map als ode solver klasse. Oder die funktion unten als template. aber eig kann die implementation hiervon einfach in ne .cpp.
// Der inlcude von fileoutput ist dann ok
//#include "system/fileoutput.h"

// Windows Workaround for unsigned int struct
#ifndef u_int64_t
#    define u_int64_t uint64_t
#endif

namespace QDACC::Numerics {

// TODO:
//  Switch iVector index to plain c++ vector
//  Implement uint_64t index instead of vector if partially summed is used -> much less memory
//  And feasable, because e.g. 36^2*2^(2*7) (NC = 8) = 21233664 which easyly fits into an integer.
// Maybe even limit to uint64_t index, because all other cases are useless.
// AND: have 1 vector of indices that map the integer onto the index vector!

// That way indices can even be restored, because the main index vector is just 1,2,3....,max (length of index map)

// Other TODO: correlation func for PI reimplementieren, dann einen punkt für 09 nachrechnen.
// Rechnungen auf cluster meddeln



// TODO: Implementation in .cpp file

class Tensor {
   public:
    using IndexFlat = uint64_t;
    using Index = uint16_t;
    // Tensor Rank indices (N0,N1,...,NN,M0,M1,...,MM)
    using IndexVector = std::vector<Index>;
    // QDACC::Type::Scalar QDACC::Type::Scalar type
    using Value = QDACC::Type::Scalar;
    using ValueVector = std::vector<Value>;
    using IndexMap = std::map<IndexVector, IndexFlat, vector_compare<Index>>;
    using IndexRefVector = std::vector<IndexVector>;

   private:
    // Current Tensor Values. Creating a new Tensor will allocate this as new
    ValueVector values;
    // When loading the tensor, save the non-zero indices here. Used to approximate sparse behaviour.
    static IndexRefVector index_flat_to_index_vector_struct_pruned;
    // Index Maps can be static because we only need them once in memory.
    // Save IndexVector - to - IndexFlat map
    static IndexMap index_vector_to_index_flat_struct;
    // Save References of the keys in vector to allow for fast acces / conversion of IndexFlat - to - IndexVector
    // using IndexRefVector = std::vector<std::reference_wrapper<IndexVector>>;
    static IndexRefVector index_flat_to_index_vector_struct;
    static int different_dimensions;
    static int tensor_dimensions;
    static int tensor_size;
    /**
     * @brief Recursive function to calculate all possible index permutations
     * @param current Current index permutation
     * @param dimensions Dimensions of the tensor
     * @param index Current index
     */
    void permute( IndexVector current, const IndexVector &dimensions, int index = 0 );

   public:
    Tensor() = default;
    Tensor( const IndexFlat num, const Value &init_value = 0 );
    Tensor( const IndexVector &dimensions, Value init_value = 0 );
    Tensor( const Tensor &other ) = default; //: values( other.values ){};

    // Access Operator
    inline Value &operator()( const IndexVector &index ) {
        return values[to_flat( index )];
    }
    inline Value &operator()( const IndexFlat &index ) {
        return values[index];
    }

    // Converts the indexvector into a flat index using the precalculated map.
    inline IndexFlat &to_flat( const IndexVector &index ) {
        return index_vector_to_index_flat_struct[index];
    }
    // Converts the flat intex into an indexvector using the precalculated map
    inline IndexVector &to_index( const IndexFlat &index ) {
        return index_flat_to_index_vector_struct[index];
    }

    inline IndexRefVector &get_indices(bool sparse = false) {
        sparse = not index_flat_to_index_vector_struct_pruned.empty();
        if (sparse)
            return index_flat_to_index_vector_struct_pruned;
        return index_flat_to_index_vector_struct;
    }

    inline size_t nonZeros() const {
        return values.size();
    }

    inline size_t size() const {
        return values.size() * ( sizeof( Value ) + sizeof( Index ) * index_flat_to_index_vector_struct.front().size() ); // / 1024. / 1024;
    }

    inline int primary_dimensions() const {
        return tensor_dimensions;
    }
    inline int secondary_dimensions() const {
        return different_dimensions;
    }

    inline int index_size() const {
        return index_flat_to_index_vector_struct.front().size();
    }

    void save_to_file(const int index);
    void load_from_file(const int index);
    void make_indices_sparse(const double eps = 1E-15);

    bool is_pruned() {
        return not index_flat_to_index_vector_struct_pruned.empty();
    }

    static inline int max_size() {
        return tensor_size;
    }
};

} // namespace QDACC::Numerics