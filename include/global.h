#pragma once

#include <cmath>
#include <complex>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
//#include <stdarg.h>
#include <string>
//#define EIGEN_SPARSEMATRIX_PLUGIN "matrixextension.h"
#include <fmt/core.h> // -DFMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <vector>
// #define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>
#include <omp.h> // -fopenmp
#include <iostream>
#include <sstream>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <map>
#include <filesystem>
//#include <experimental/filesystem>
//namespace std {
//namespace filesystem = std::experimental::filesystem::v1;
//}

extern std::string PREFIX_PERCENT;
extern std::string PREFIX_PERCENT_TIME;
extern std::string PREFIX_PERCENT_TIME_FINAL;
extern std::string PREFIX_WARNING;
extern std::string PREFIX_ERROR;
extern std::string PREFIX_DEBUG;
extern std::string PREFIX_OUTPUT;
extern std::string PREFIX_SUFFIX;

#include "misc/commandlinearguments.h"
#include "misc/log.h"
#include "misc/ProgressBar.h"
#include "misc/helperfunctions.h"
#include "misc/timer.h"
#include "misc/sysinfo.h"
#include "misc/FixedSizeSparseMap.h"

//using Eigen::MatrixXcd;
//using Eigen::MatrixXd;

using namespace std::complex_literals;

extern std::string global_message_normaltermination;
extern std::string global_message_error_divergent;
extern std::string global_message_error_wrong_number_input;

#define DIR_T 1
#define DIR_TAU 2
#define DIR_W 3

#define TIMETRANSFORMATION_ANALYTICAL 0
#define TIMETRANSFORMATION_MATRIXEXPONENTIAL 1

#define OPERATOR_PHOTONIC_CREATE 0
#define OPERATOR_PHOTONIC_ANNIHILATE 1
#define OPERATOR_PHOTONIC_STATE 2

#define PHONON_APPROXIMATION_BACKWARDS_INTEGRAL 0
#define PHONON_APPROXIMATION_TRANSFORMATION_MATRIX 1
#define PHONON_APPROXIMATION_TIMETRANSFORMATION 2
#define PHONON_APPROXIMATION_LINDBLAD_FULL 3
#define PHONON_APPROXIMATION_MIXED 4
#define PHONON_PATH_INTEGRAL 5

typedef std::complex<double> Scalar;
typedef Eigen::MatrixXcd Dense;
typedef Eigen::SparseMatrix<Scalar> Sparse;

// Vector for mat/time, tuple
class SaveState {
   public:
    Sparse mat;
    double t;
    SaveState( const Sparse &mat, const double time ) : mat( mat ), t( time ){};
};
class SaveStateTau {
   public:
    Sparse mat1, mat2;
    double t, tau;
    SaveStateTau( const Sparse &mat1, const Sparse &mat2, const double t, const double tau ) : mat1( mat1 ), mat2( mat2 ), t( t ), tau( tau ){};
    SaveStateTau( const Sparse &mat, const double t ) : mat1( mat ), t( t ), tau( 0 ){};
    SaveStateTau(){};
};
class SaveScalar {
   public:
    Scalar scalar;
    double t, tau;
    //SaveScalarTuple( const Sparse &mat1, const Sparse &mat2, const double t, const double tau ) : mat1( mat1 ), mat2( mat2 ), t( t ), tau( tau ){};
    SaveScalar( const Scalar &scalar, const double t, const double tau = 0.0 ) : scalar( scalar ), t( t ), tau( tau ){};
};

bool Save_State_sort_t( const SaveStateTau &ss1, const SaveStateTau &ss2 );
bool Save_State_sort_tau( const SaveStateTau &ss1, const SaveStateTau &ss2 );