#pragma once
#include "global.h"

namespace Solver {
    // Command variables
    extern const unsigned int CHANGE_TO_SINGLETHREADED_MAINPROGRAM;

    template <typename T>
    MatType iterate_definite_integral( const MatType &rho, T rungefunction, const double t, const double step ) {
        MatType rk1 = rungefunction( rho, t );
        MatType rk2 = rungefunction( rho + step * 0.5 * rk1, t + step * 0.5 );
        MatType rk3 = rungefunction( rho + step * 0.5 * rk2, t + step * 0.5 );
        MatType rk4 = rungefunction( rho + step * rk3, t + step );
        // Dichtematrix
        return rho + step / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
    }

    // Description: Integrates rho from t0 to t1 via Rungefunc.
    // Type: Solver public function
    // @return: [vector<SaveState>] Vector of save state tuples (matrix, time)
    std::vector<SaveState> calculate_definite_integral_vec( MatType rho, std::function<MatType( const MatType &, const double )> const &rungefunction, const double t0, const double t1, const double step );

    // Description: Integrates rho from t0 to t1 via Rungefunc.
    // Type: Solver public function
    // @return: [SaveState] Save state tuple (matrix, time)
    SaveState calculate_definite_integral( MatType rho, std::function<MatType( const MatType &, const double )> const &rungefunction, const double t0, const double t1, const double step );
};