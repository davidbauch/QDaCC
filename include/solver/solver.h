#pragma once
#include "global.h"

namespace QDACC {

namespace Numerics {

namespace RKCoefficients {

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

const double e1 = 5179. / 57600.; // 71. / 57600.;
const double e2 = 0;
const double e3 = 7571. / 16695.;    //-71. / 16695.;
const double e4 = 393. / 640.;       // 71. / 1920.;
const double e5 = -92097. / 339200.; //-17253. / 339200.;
const double e6 = 187. / 2100.;      // 22. / 525.;
const double e7 = 1. / 40.;          //-1. / 40.;

} // namespace RKCoefficients

std::pair<Sparse, double> iterate_definite_integral( const Sparse &rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t, const double t_step, const int order );

// Description: Integrates rho from t0 to t1 via Rungefunc.
// Type: Solver public function
// @return: [vector<QDACC::SaveState>] Vector of save state tuples (matrix, time)
std::vector<QDACC::SaveState> calculate_definite_integral_vec( const Sparse &rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step, const double tolerance, const double stepmin, const double stepmax, const double stepdelta, const int order );

// Description: Integrates rho from t0 to t1 via Rungefunc.
// Type: Solver public function
// @return: [QDACC::SaveState] Save state tuple (matrix, time)
QDACC::SaveState calculate_definite_integral( Sparse rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step, const double tolerance, const double stepmin, const double stepmax, const double stepdelta, const int order );

// Description: Uses the interpolation class to monotone-cubic spline interpolate a given vector of saved states, resulting in much smoother output. Should probably not be used
// Type: Solver public function
// @return Returns a vector of interpolated QDACC::SaveStates
std::vector<QDACC::SaveState> interpolate_curve( const std::vector<QDACC::SaveState> &input, double t_start, double t_end, double t_step, int threads, int order = 0, bool output_handler = false );
std::vector<QDACC::SaveState> interpolate_curve( const std::vector<QDACC::SaveState> &input, double t_start, double t_end, const std::vector<double> &t_values, const std::vector<double> &t_steps, const std::map<double, size_t> &t_index, int order = 0, bool output_handler = false );

double get_tdelta( const Dense &gmat_time, size_t fixed_index, size_t var_index );
double get_taudelta( const Dense &gmat_time, size_t fixed_index, size_t var_index );
double get_tdelta( const std::vector<SaveState> &savedStates, size_t var_index );

} // namespace Numerics

} // namespace QDACC