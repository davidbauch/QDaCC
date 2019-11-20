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
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h> // -fopenmp
#include <iostream>
#include <functional>
#include <limits>

std::string PREFIX_PERCENT = "@#PERCENT#@";
std::string PREFIX_PERCENT_TIME = "@#PERCENTTIME#@";
std::string PREFIX_PERCENT_TIME_FINAL = "@#PERCENTTIMEFINAL#@";
std::string PREFIX_WARNING = "@#WARNING#@";
std::string PREFIX_ERROR = "@#ERROR#@";
std::string PREFIX_DEBUG = "@#DEBUG#@";
std::string PREFIX_OUTPUT = "@#OUTPUT#@";
std::string PREFIX_SUFFIX = "@#SUFFIX#@";

#include "misc/log.h"
#include "misc/ProgressBar.h"
#include "misc/helperfunctions.h"
#include "misc/timer.h"

//using Eigen::MatrixXcd;
//using Eigen::MatrixXd;

using namespace std::complex_literals;

std::string global_message_normaltermination = "[Programm Terminated normally]";
std::string global_message_error_divergent = "[FATAL ERROR: PROCESS DIVERGED, PROGRAM TERMINATED]";
std::string global_message_error_wrong_number_input = "[FATAL ERROR: WRONG NUMBER OF INPUT PARAMETERS, PROGRAM TERMINATED]";

#define DIR_T 1
#define DIR_TAU 2
#define DIR_W 3

#define TIMETRANSFORMATION_ANALYTICAL 0
#define TIMETRANSFORMATION_MATRIXEXPONENTIAL 1

#define OPERATOR_PHOTONIC_CREATE 0
#define OPERATOR_PHOTONIC_ANNIHILATE 1

#define PHONON_APPROXIMATION_BACKWARDS_INTEGRAL 0
#define PHONON_APPROXIMATION_TRANSFORMATION_MATRIX 1
#define PHONON_APPROXIMATION_TIMETRANSFORMATION 2
#define PHONON_APPROXIMATION_LINDBLAD_FULL 3
#define PHONON_APPROXIMATION_MULITPLY_HAMILTON 4

typedef std::complex<double> dcomplex;
typedef Eigen::SparseMatrix<dcomplex> SparseMat;
typedef Eigen::MatrixXcd DenseMat;
typedef Eigen::SparseMatrix<dcomplex> MatType;

// Vector for mat/time, tuple
class SaveState {
   public:
    MatType mat;
    double t;
    SaveState( const MatType &mat, const double time ) : mat( mat ), t( time ){};
};
class SaveStateTau {
   public:
    MatType mat1, mat2;
    double t, tau;
    SaveStateTau( const MatType &mat1, const MatType &mat2, const double t, const double tau ) : mat1( mat1 ), mat2( mat2 ), t( t ), tau( tau ){};
    SaveStateTau( const MatType &mat, const double t ) : mat1( mat ), t( t ), tau( 0 ){};
};

bool Save_State_sort_t(const SaveStateTau &ss1, const SaveStateTau &ss2) {
    return (ss1.t < ss2.t);
}
bool Save_State_sort_tau(const SaveStateTau &ss1, const SaveStateTau &ss2) {
    return (ss1.tau < ss2.tau);
}