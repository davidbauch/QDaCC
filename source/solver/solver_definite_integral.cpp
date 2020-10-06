#include "solver/solver.h"

std::vector<SaveState> Solver::calculate_definite_integral_vec( MatType rho, std::function<MatType( const MatType &, const double )> const &rungefunction, const double t0, const double t1, const double step ) { //std::function<MatrixXcd( const MatrixXcd &, const double )>
    std::vector<SaveState> ret;
    //ret.reserve( std::ceil( std::abs( t1 - t0 ) / std::abs( step ) ) );
    ret.emplace_back( SaveState( rho, t0 ) );
    if ( step > 0 )
        for ( double t = t0; t < t1; t += step ) {
            rho = iterate_definite_integral( rho, rungefunction, t, step );
            ret.emplace_back( SaveState( rho, t ) );
        }
    else if ( step < 0 )
        for ( double t = t0; t > t1; t += step ) {
            rho = iterate_definite_integral( rho, rungefunction, t, step );
            ret.emplace_back( SaveState( rho, t ) );
        }
    return ret;
}

SaveState Solver::calculate_definite_integral( MatType rho, std::function<MatType( const MatType &, const double )> const &rungefunction, const double t0, const double t1, const double step ) { //std::function<MatrixXcd( const MatrixXcd &, const double )>
    //ret.reserve( std::ceil( std::abs( t1 - t0 ) / std::abs( step ) ) );
    if ( step > 0 )
        for ( double t = t0; t < t1; t += step ) {
            rho = iterate_definite_integral( rho, rungefunction, t, step );
        }
    else if ( step < 0 )
        for ( double t = t0; t > t1; t += step ) {
            rho = iterate_definite_integral( rho, rungefunction, t, step );
        }
    return SaveState( rho, t1 );
}

const unsigned int Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM = 1111111111;