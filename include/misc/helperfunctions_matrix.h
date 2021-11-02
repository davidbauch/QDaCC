#pragma once

#include "typedef.h"

namespace QDLC {

namespace Matrix {

QDLC::Type::Dense dense_projector( const QDLC::Type::Dense &input );

QDLC::Type::Sparse sparse_projector( const QDLC::Type::Sparse &input );

void init_sparsevector( std::vector<QDLC::Type::Sparse> &vec, int dim, int count );

QDLC::Type::Dense tensor( const QDLC::Type::Dense &a, const QDLC::Type::Dense &b );

QDLC::Type::Dense tensor( const std::vector<QDLC::Type::Dense> &m );

std::vector<std::string> tensor( const std::vector<std::string> &a, const std::vector<std::string> &b );

std::vector<std::string> tensor( const std::vector<std::vector<std::string>> &m );

std::pair<QDLC::Type::Dense, QDLC::Type::Dense> meshgrid( double x_min, double y_min, double x_max, double y_max, int N, bool include_endpoint = true );

} // namespace Matrix

} // namespace QDLC