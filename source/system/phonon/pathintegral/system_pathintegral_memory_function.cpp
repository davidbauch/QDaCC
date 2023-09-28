#include "system/system.h"

using namespace QDACC;

Scalar System::dgl_phonon_memory_function( const int t_delta, const int i_n, const int j_n, const int i_nd, const int j_nd ) {
    Scalar result = 0.0;
    const int coupling_index_upper = i_n;
    const int coupling_index_lower = j_n;
    const Scalar coupling_value_upper = operatorMatrices.phonon_group_index_to_coupling_value[coupling_index_upper];
    const Scalar coupling_value_lower = operatorMatrices.phonon_group_index_to_coupling_value[coupling_index_lower];
    const Scalar coupling_value = phi_vector_int[t_delta];
    const Scalar coupling_value_conj = std::conj( coupling_value );
    if ( i_n == i_nd )
        result -= coupling_value * coupling_value_upper;
    if ( j_n == j_nd )
        result -= coupling_value_conj * coupling_value_lower;
    if ( i_n == j_nd )
        result += coupling_value_conj * coupling_value_upper;
    if ( j_n == i_nd )
        result += coupling_value * coupling_value_lower;
    return result;
}