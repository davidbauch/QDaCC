#include "solver/solver.h"

std::vector<SaveState> Solver::calculate_definite_integral_vec( const Sparse &rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step, const double tolerance ) { //std::function<MatrixXcd( const MatrixXcd &, const double )>
    double t_step = step;
    std::vector<SaveState> ret;
    ret.emplace_back( SaveState( rho, t0 ) );

    double t_t = t0;
    while ( t0 < t1 ? t_t < t1 : t_t > t1 ) {
        bool accept = false;
        // Hit endpoint exactly
        if ( t0 < t1 ? t_t + t_step > t1 : t_t + t_step < t1 ) {
            t_step = t1 - t_t;
            accept = true;
        }
        // Runge-Kutta iteration
        auto rkret = iterate_definite_integral( rho, rungefunction, t_t, t_step );
        double error = rkret.second;
        double dh = std::pow( tolerance / 2. / error, 0.25 );
        if ( error < tolerance ) {
            accept = true;
        }
        if ( accept ) {
            t_t += t_step;
            ret.emplace_back( SaveState( rkret.first, t_t ) );
        }
        t_step = t_step * dh;
    }
    return ret;
}

SaveState Solver::calculate_definite_integral( Sparse rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step, const double tolerance ) { //std::function<MatrixXcd( const MatrixXcd &, const double )>
    double t_step = step;
    double t_t = t0;
    int iterations = 0;
    double maxerror = 0;
    Log::L3( "Calculating Definite Integral from t0 = {} to t1 = {} at initial step = {}\n", t0, t1, step );
    Log::L3( "Beginning Rho:\n{}\n", Dense( rho ) );
    while ( t0 < t1 ? t_t < t1 : t_t > t1 ) {
        bool accept = false;
        // Hit endpoint exactly
        if ( t0 < t1 ? t_t + t_step > t1 : t_t + t_step < t1 ) {
            Log::L3( " - - - - Timestep would overshoot to t1 = {} (dt = {}), adjusting to {}\n", t1 + t_step, t_step, t_t - t1 );
            t_step = t1 - t_t;
            accept = true;
        }
        // Runge-Kutta iteration
        auto rkret = iterate_definite_integral( rho, rungefunction, t_t, t_step );
        double error = rkret.second;
        double dh = std::pow( tolerance / 2. / error, 0.25 );
        Log::L3( " - [t = {}] - - - Local error: {} - dh = {}, current timestep is: {}, new timestep will be: {}\n", t_t, error, dh, t_step, t_step * dh );
        if ( error < tolerance ) {
            accept = true;
        }
        if ( accept ) {
            Log::L3( " - [t = {}] - Accepdet step - Local error: {} - current timestep: {}, dh = {}\n", t_t, error, t_step, dh );
            t_t += t_step;
            rho = rkret.first;
            maxerror = std::max( maxerror, error );
        }
        t_step = t_step * dh;
        iterations++;
    }
    Log::L3( "Calculating Definite Integral from t0 = {} to t1 = {} at initial step = {} -- Done.\n", t0, t1, step );
    Log::L3( "Done Chi integrating, (t0 = {}) t1 was supposted to be {} and ended up {}, maxerror was {}, did {}/{} iterations (var/const dt) \n", t0, t1, t_t, maxerror, iterations, std::floor( t1 - t0 / step ) );
    Log::L3( "End Rho:\n{}\n", Dense( rho ) );
    return SaveState( rho, t1 );
}

const unsigned int Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM = 1111111111;