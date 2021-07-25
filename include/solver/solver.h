#pragma once
#include "global.h"
#include "misc/interpolant.h"

namespace Solver {
// Command variables
extern const unsigned int CHANGE_TO_SINGLETHREADED_MAINPROGRAM;

// RK 4&5 coefficients (Dormand–Prince method)
const double a1 = 0;
const double a2 = 1. / 5.;
const double a3 = 3. / 10.;
const double a4 = 4. / 5.;
const double a5 = 8. / 9.;
const double a6 = 1.0;

const double b11 = 1. / 5.;
const double b21 = 3. / 40.;
const double b31 = 44. / 45.;
const double b41 = 19372. / 6561.;
const double b51 = 9017. / 3168.;
const double b61 = 35. / 384.;
const double b22 = 9. / 40.;
const double b32 = -56. / 15.;
const double b42 = -25360. / 2187.;
const double b52 = -355. / 33.;
const double b33 = 32. / 9.;
const double b43 = 64448. / 6561.;
const double b53 = 46732. / 5247.;
const double b63 = 500. / 1113.;
const double b44 = -212. / 729.;
const double b54 = 49. / 176.;
const double b64 = 125. / 192.;
const double b55 = -5103 / 18656.;
const double b65 = -2187. / 6784.;
const double b66 = 11. / 84.;

const double e1 = 5179. / 57600.; //71. / 57600.;
const double e2 = 0;
const double e3 = 7571. / 16695.;    //-71. / 16695.;
const double e4 = 393. / 640.;       //71. / 1920.;
const double e5 = -92097. / 339200.; //-17253. / 339200.;
const double e6 = 187. / 2100.;      //22. / 525.;
const double e7 = 1. / 40.;          //-1. / 40.;

template <typename T>
std::pair<Sparse, double> iterate_definite_integral( const Sparse &rho, T rungefunction, const double t, const double t_step ) {
    //Sparse rk1 = rungefunction( rho, t );
    //Sparse rk2 = rungefunction( rho + step * 0.5 * rk1, t + step * 0.5 );
    //Sparse rk3 = rungefunction( rho + step * 0.5 * rk2, t + step * 0.5 );
    //Sparse rk4 = rungefunction( rho + step * rk3, t + step );
    //// Dichtematrix
    //return rho + step / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
    // k1-6 ausrechnen
    Sparse k1 = rungefunction( rho, t );
    Sparse k2 = rungefunction( rho + t_step * Solver::b11 * k1, t + Solver::a2 * t_step );
    Sparse k3 = rungefunction( rho + t_step * ( Solver::b21 * k1 + Solver::b22 * k2 ), t + Solver::a3 * t_step );
    Sparse k4 = rungefunction( rho + t_step * ( Solver::b31 * k1 + Solver::b32 * k2 + Solver::b33 * k3 ), t + Solver::a4 * t_step );
    Sparse k5 = rungefunction( rho + t_step * ( Solver::b41 * k1 + Solver::b42 * k2 + Solver::b43 * k3 + Solver::b44 * k4 ), t + Solver::a5 * t_step );
    Sparse k6 = rungefunction( rho + t_step * ( Solver::b51 * k1 + Solver::b52 * k2 + Solver::b53 * k3 + Solver::b54 * k4 + Solver::b55 * k5 ), t + Solver::a6 * t_step );

    Sparse drho = b61 * k1 + b63 * k3 + b64 * k4 + b65 * k5 + b66 * k6;
    Sparse ret = rho + t_step * drho;
    // Error
    Sparse k7 = rungefunction( ret, t + t_step );
    dSparse errmat = ( drho - ( k1 * Solver::e1 + k3 * Solver::e3 + k4 * Solver::e4 + k5 * Solver::e5 + k6 * Solver::e6 + k7 * Solver::e7 ) ).cwiseAbs2();
    double err = errmat.sum() / drho.cwiseAbs2().sum();

    // Dichtematrix
    return std::make_pair( ret, err );
}

// Description: Integrates rho from t0 to t1 via Rungefunc.
// Type: Solver public function
// @return: [vector<SaveState>] Vector of save state tuples (matrix, time)
std::vector<SaveState> calculate_definite_integral_vec( const Sparse &rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step, const double tolerance );

// Description: Integrates rho from t0 to t1 via Rungefunc.
// Type: Solver public function
// @return: [SaveState] Save state tuple (matrix, time)
SaveState calculate_definite_integral( Sparse rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step, const double tolerance );

// Description: Uses the interpolation class to monotone-cubic spline interpolate a given vector of saved states, resulting in much smoother output. Should probably not be used
// Type: Solver public function
// @return Returns a vector of interpolated SaveStates
std::vector<SaveState> calculate_smooth_curve( const std::vector<SaveState> &input, double t_start, double t_end, int num_of_points, bool output_handler = true );
}; // namespace Solver