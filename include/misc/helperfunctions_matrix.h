#pragma once

#include "typedef.h"

namespace QDACC {

namespace Matrix {

QDACC::Type::MatrixMain projector( const QDACC::Type::MatrixMain &input );

QDACC::Type::Dense dense_projector( const QDACC::Type::Dense &input );

QDACC::Type::Sparse sparse_projector( const QDACC::Type::Sparse &input );

QDACC::Type::MatrixMain tensor( const QDACC::Type::MatrixMain &a, const QDACC::Type::MatrixMain &b );

QDACC::Type::MatrixMain tensor( const std::vector<QDACC::Type::MatrixMain> &m );

std::vector<std::string> tensor( const std::vector<std::string> &a, const std::vector<std::string> &b );

std::vector<std::string> tensor( const std::vector<std::vector<std::string>> &m );

std::pair<QDACC::Type::Dense, QDACC::Type::Dense> meshgrid( double x_min, double y_min, double x_max, double y_max, int N, bool include_endpoint = true );

} // namespace Matrix

} // namespace QDACC