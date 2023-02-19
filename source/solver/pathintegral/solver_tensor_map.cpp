#include "solver/solver_tensor_map.h"
#include "system/fileoutput.h"

using namespace QDLC::Numerics;

QDLC::Numerics::Tensor::IndexMap QDLC::Numerics::Tensor::index_vector_to_index_flat_struct{};
QDLC::Numerics::Tensor::IndexRefVector QDLC::Numerics::Tensor::index_flat_to_index_vector_struct{};
int QDLC::Numerics::Tensor::tensor_dimensions{};
int QDLC::Numerics::Tensor::different_dimensions{};
int QDLC::Numerics::Tensor::tensor_size{};

Tensor::Tensor( const IndexFlat num, const QDLC::Type::Scalar &init_value ) {
    values = Tensor::ValueVector( num, init_value );
}

Tensor::Tensor( const IndexVector &dimensions, QDLC::Type::Scalar init_value ) {
    // Set dimensions
    tensor_dimensions = dimensions[0];
    different_dimensions = dimensions[2];
    Log::L2("[PathIntegral] Adding Tensor with primary dimension: {} and secondary dimension: {}\n", tensor_dimensions, different_dimensions);
    // Reserve Memory and initialize value vector
    auto max_elements = 1;
    index_vector_to_index_flat_struct.clear();
    index_flat_to_index_vector_struct.clear();
    std::ranges::for_each( dimensions, [&]( const auto &el ) { max_elements *= el; } );
    index_flat_to_index_vector_struct = IndexRefVector( max_elements ); // FIXME: max_elements ist für Biex nicht gleich der hinzugefügten elemente?????????????
    values = ValueVector( max_elements, init_value );
    tensor_size = max_elements;
    // Calculate and Save all possible index permutations
    permute( IndexVector( dimensions.size(), 0 ), dimensions );
    // Save References to map keys in vector for fast lookup
    // std::ranges::for_each( index_vector_to_index_flat_struct, []( const auto &el ) { index_flat_to_index_vector_struct[el.second] = el.first; } );
    // Lazy fix for upper fixme:
    for ( int i = index_flat_to_index_vector_struct.size() - 1; i > 0; i-- )
        if ( index_flat_to_index_vector_struct[i].size() == 0 )
            index_flat_to_index_vector_struct.erase( index_flat_to_index_vector_struct.begin() + i );
    Log::L2( "[PathIntegral] Added {} elements to the dimension vector ({} elements to the inverse map).\n", index_vector_to_index_flat_struct.size(), index_flat_to_index_vector_struct.size() );
}

static inline std::string _get_tensor_name( const int index ) {
    return "tensor_" + std::to_string( index );
}

void Tensor::save_to_file( const int index ) {
    const std::string file_name = _get_tensor_name( index );
    auto &f_adm = QDLC::FileOutput::add_file( file_name, "adm" );
    for ( auto i = 0; i < values.size(); i++ ) {
        const auto real = values[i].real();
        const auto imag = values[i].real();
        if ( real != 0 || imag != 0)
            f_adm << i << " " << real << " " << imag << "\n";
    }
    f_adm.close();
}

void Tensor::load_from_file( const int index ) {
    const std::string file_name = _get_tensor_name( index );
    Log::L2( "[PathIntegral] Loading tensor from file {}.adm.\n", file_name );
    auto f_adm = FileOutput::load_file( file_name, "adm" ); 
    // Read all lines in file:
    int non_zeros = 0;
    int line;
    double real, imag;
    //for ( auto i = 0; i < values.size(); i++ ) {
    while ( f_adm >> line >> real >> imag ) {
        f_adm >> line >> real >> imag;
        values[line] = Tensor::Value( real, imag );
        non_zeros++;
    }
    f_adm.close();
    Log::L2( "[PathIntegral] Loaded {} non-zero elements from file {}.\n", non_zeros, file_name );
}

void Tensor::permute( IndexVector current, const IndexVector &dimensions, int index ) {
    if ( index == dimensions.size() ) {
        const auto index = index_vector_to_index_flat_struct.size();
        index_vector_to_index_flat_struct[current] = index;
        if ( current.size() > 0 )
            index_flat_to_index_vector_struct[index] = current;
    } else {
        for ( Index c = 0; c < dimensions[index]; c++ ) {
            current[index] = c;
            permute( current, dimensions, index + 1 );
        }
    }
}