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
typedef uint64_t Index;

} // namespace Type

namespace Message {

namespace Prefix {

const std::string PERCENT = "@#PERCENT#@";
const std::string PERCENT_TIME = "@#PERCENTTIME#@";
const std::string PERCENT_TIME_FINAL = "@#PERCENTTIMEFINAL#@";
const std::string WARNING = "@#WARNING#@";
const std::string ERROR = "@#ERROR#@";
const std::string DEBUG = "@#DEBUG#@";
const std::string OUTPUT = "@#OUTPUT#@";
const std::string SUFFIX = "@#SUFFIX#@";

} // namespace Prefix

const std::string global_normaltermination = "[Programm Terminated normally]";
const std::string global_error_divergent = "[FATAL ERROR: PROCESS DIVERGED, PROGRAM TERMINATED]";
const std::string global_error_wrong_number_input = "[FATAL ERROR: WRONG NUMBER OF INPUT PARAMETERS, PROGRAM TERMINATED]";

} // namespace Message

} // namespace QDLC