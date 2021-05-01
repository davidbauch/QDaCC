#pragma once
#include "global.h"
#include "misc/interpolant.h"

namespace Solver {
// Command variables
extern const unsigned int CHANGE_TO_SINGLETHREADED_MAINPROGRAM;

template <typename T>
Sparse iterate_definite_integral( const Sparse &rho, T rungefunction, const double t, const double step ) {
    Sparse rk1 = rungefunction( rho, t );
    Sparse rk2 = rungefunction( rho + step * 0.5 * rk1, t + step * 0.5 );
    Sparse rk3 = rungefunction( rho + step * 0.5 * rk2, t + step * 0.5 );
    Sparse rk4 = rungefunction( rho + step * rk3, t + step );
    // Dichtematrix
    return rho + step / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
}

// Description: Integrates rho from t0 to t1 via Rungefunc.
// Type: Solver public function
// @return: [vector<SaveState>] Vector of save state tuples (matrix, time)
std::vector<SaveState> calculate_definite_integral_vec( Sparse rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step );

// Description: Integrates rho from t0 to t1 via Rungefunc.
// Type: Solver public function
// @return: [SaveState] Save state tuple (matrix, time)
SaveState calculate_definite_integral( Sparse rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step );

// Description: Uses the interpolation class to monotone-cubic spline interpolate a given vector of saved states, resulting in much smoother output. Should probably not be used
// Type: Solver public function
// @return Returns a vector of interpolated SaveStates
std::vector<SaveState> calculate_smooth_curve( const std::vector<SaveState> &input, double t_start, double t_end, int num_of_points, bool output_handler = true );
}; // namespace Solver