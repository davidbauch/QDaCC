#pragma once
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std::complex_literals;

namespace QDACC::Type {

using Scalar = std::complex<double>;
using Dense = Eigen::MatrixXcd;
using dDense = Eigen::MatrixXd;
using Sparse = Eigen::SparseMatrix<Scalar>;
using dSparse = Eigen::SparseMatrix<double>;
using iVector = Eigen::VectorXi;
using dVector = Eigen::VectorXd;
using Vector = Eigen::VectorXcd;
using dTriplet = Eigen::Triplet<double>;
using Triplet = Eigen::Triplet<Scalar>;
using Index = uint64_t;

//#define USE_SPARSE_MATRIX

#ifdef USE_SPARSE_MATRIX
using MatrixMain = Sparse;
using dMatrixMain = dSparse;
#else
using MatrixMain = Dense;
using dMatrixMain = dDense;
#endif

template <typename T>
using NestedVector = std::vector<std::vector<T>>;
} // namespace QDACC::Type

namespace QDACC::Message::Prefix {

// TODO: remove handler output, deprecated!
const std::string PERCENT = "@#PERCENT#@";
const std::string PERCENT_TIME = "@#PERCENTTIME#@";
const std::string PERCENT_TIME_FINAL = "@#PERCENTTIMEFINAL#@";
const std::string WARNING = "@#WARNING#@";
const std::string PERROR = "@#ERROR#@";
const std::string DEBUG = "@#DEBUG#@";
const std::string OUTPUT = "@#OUTPUT#@";
const std::string SUFFIX = "@#SUFFIX#@";

} // namespace QDACC::Message::Prefix

// namespace QDACC