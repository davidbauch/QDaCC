#pragma once

#include <cmath>
#include <complex>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
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
//#define EIGEN_USE_BLAS
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
#include <ranges>
#include <format>
template<typename T>
struct std::formatter<std::complex<T>> : std::formatter<std::string> {
    auto format(const std::complex<T>& c, std::format_context& ctx) {
        return std::formatter<std::string>::format(
            std::format("({},{})", c.real(), c.imag()), ctx);
    }
};
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
using namespace QDACC::Type;

#include "system/savestate.h"

// Global Enums
namespace QDACC {

// Time Transformation order including analytical and matrix exponential
enum class TransformationOrder {
    // The time transformation will be carried out analytically
    Analytical,
    // The time transformation will be carried out by using Eigen's matrix exponential exp(-i*H_0)
    MatrixExponential
};

enum class PhotonicOperator {
    Create,
    Annihilate,
    State
};

enum class PhononApproximation {
    BackwardsIntegral,
    TransformationMatrix,
    Timetransformation,
    LindbladRates,
    Mixed,
    PathIntegral
};

} // namespace QDACC