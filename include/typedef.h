#pragma once
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std::complex_literals;

namespace QDLC {

namespace Type {

typedef std::complex<double> Scalar;
typedef Eigen::MatrixXcd Dense;
typedef Eigen::MatrixXd dDense;
typedef Eigen::SparseMatrix<Scalar> Sparse;
typedef Eigen::SparseMatrix<double> dSparse;
typedef Eigen::VectorXi iVector;
typedef Eigen::VectorXd dVector;
typedef Eigen::VectorXcd Vector;
typedef Eigen::Triplet<double> dTriplet;
typedef Eigen::Triplet<Scalar> Triplet;

} // namespace Type

} // namespace QDLC