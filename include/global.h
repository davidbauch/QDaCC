#pragma once

#include <cmath>
#include <complex>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fmt/core.h> // -DFMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <vector>
template <class T>
std::ostream &operator<<( std::ostream &os, const std::vector<T> &v ) {
    os << "[";
    // std::for_each(v.begin(), v.end(), [os](T& el){os << el;});
    for ( typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii ) {
        os << *ii;
        if ( ii != v.end() )
            os << " ";
    }
    os << "]";
    return os;
}
// #define EIGEN_DEFAULT_DENSE_INDEX_TYPE int
//#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h> // -fopenmp
#include <iostream>
#include <sstream>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <map>
#include <set>
#include <unordered_map>
#include <filesystem>
//#include <experimental/filesystem>
// namespace std {
// namespace filesystem = std::experimental::filesystem::v1;
//}

#include "misc/commandlinearguments.h"
#include "misc/log.h"
#include "misc/ProgressBar.h"
#include "misc/helperfunctions_string.h"
#include "misc/helperfunctions_matrix.h"
#include "misc/helperfunctions_math.h"
#include "misc/helperfunctions.h"
#include "misc/timer.h"
#include "misc/sysinfo.h"
#include "solver/solver_tensor_map.h"

#include "typedef.h"
using namespace QDLC::Type;

#include "system/savestate.h"

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