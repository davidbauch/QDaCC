#include "solver/solver.h"

std::pair<Sparse, double> Solver::iterate_definite_integral( const Sparse &rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t, const double t_step, const int order ) {
    if ( order == 4 ) {
        Sparse rk1 = rungefunction( rho, t );
        Sparse rk2 = rungefunction( rho + t_step * 0.5 * rk1, t + t_step * 0.5 );
        Sparse rk3 = rungefunction( rho + t_step * 0.5 * rk2, t + t_step * 0.5 );
        Sparse rk4 = rungefunction( rho + t_step * rk3, t + t_step );
        // Dichtematrix
        Sparse ret = rho + t_step / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
        return std::make_pair( ret, 0 );
    }
    // k1-6 ausrechnen
    Sparse k1 = rungefunction( rho, t );
    Sparse k2 = rungefunction( rho + t_step * Solver::b11 * k1, t + Solver::a2 * t_step );
    Sparse k3 = rungefunction( rho + t_step * ( Solver::b21 * k1 + Solver::b22 * k2 ), t + Solver::a3 * t_step );
    Sparse k4 = rungefunction( rho + t_step * ( Solver::b31 * k1 + Solver::b32 * k2 + Solver::b33 * k3 ), t + Solver::a4 * t_step );
    Sparse k5 = rungefunction( rho + t_step * ( Solver::b41 * k1 + Solver::b42 * k2 + Solver::b43 * k3 + Solver::b44 * k4 ), t + Solver::a5 * t_step );
    Sparse k6 = rungefunction( rho + t_step * ( Solver::b51 * k1 + Solver::b52 * k2 + Solver::b53 * k3 + Solver::b54 * k4 + Solver::b55 * k5 ), t + Solver::a6 * t_step );

    Sparse drho = b61 * k1 + b63 * k3 + b64 * k4 + b65 * k5 + b66 * k6;
    Sparse ret = rho + t_step * drho;
    if ( order == 5 )
        return std::make_pair( ret, 0 );
    // Error
    Sparse k7 = rungefunction( ret, t + t_step );
    dSparse errmat = ( drho - ( k1 * Solver::e1 + k3 * Solver::e3 + k4 * Solver::e4 + k5 * Solver::e5 + k6 * Solver::e6 + k7 * Solver::e7 ) ).cwiseAbs2();
    double err = errmat.sum() / drho.cwiseAbs2().sum();

    // Dichtematrix
    return std::make_pair( ret, err );
}

std::vector<SaveState> Solver::calculate_definite_integral_vec( const Sparse &rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step, const double tolerance, const int order ) { //std::function<MatrixXcd( const MatrixXcd &, const double )>
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
        auto rkret = iterate_definite_integral( rho, rungefunction, t_t, t_step, order );
        double error = rkret.second;
        double dh = std::pow( tolerance / 2. / error, 0.25 );
        if ( error < tolerance or order != 45 ) {
            accept = true;
        }
        if ( accept ) {
            t_t += t_step;
            ret.emplace_back( SaveState( rkret.first, t_t ) );
        }
        if ( order == 45 )
            t_step = t_step * dh;
    }
    return ret;
}

SaveState Solver::calculate_definite_integral( Sparse rho, std::function<Sparse( const Sparse &, const double )> const &rungefunction, const double t0, const double t1, const double step, const double tolerance, const int order ) { //std::function<MatrixXcd( const MatrixXcd &, const double )>
    double t_step = step;
    double t_t = t0;
    //int iterations = 0;
    //double maxerror = 0;
    //Log::L3( "Calculating Definite Integral from t0 = {} to t1 = {} at initial step = {}\n", t0, t1, step );
    //Log::L3( "Beginning Rho:\n{}\n", Dense( rho ) );
    while ( t0 < t1 ? t_t < t1 : t_t > t1 ) {
        bool accept = false;
        // Hit endpoint exactly
        if ( t0 < t1 ? t_t + t_step > t1 : t_t + t_step < t1 ) {
            //Log::L3( " - - - - Timestep would overshoot to t1 = {} (dt = {}), adjusting to {}\n", t1 + t_step, t_step, t_t - t1 );
            t_step = t1 - t_t;
            accept = true;
        }
        // Runge-Kutta iteration
        auto rkret = iterate_definite_integral( rho, rungefunction, t_t, t_step, order );
        double error = rkret.second;
        double dh = std::pow( tolerance / 2. / error, 0.25 );
        //Log::L3( " - [t = {}] - - - Local error: {} - dh = {}, current timestep is: {}, new timestep will be: {}\n", t_t, error, dh, t_step, t_step * dh );
        if ( error < tolerance or order != 45 ) {
            accept = true;
        }
        if ( accept ) {
            //Log::L3( " - [t = {}] - Accepdet step - Local error: {} - current timestep: {}, dh = {}\n", t_t, error, t_step, dh );
            t_t += t_step;
            rho = rkret.first;
            //maxerror = std::max( maxerror, error );
        }
        if ( order == 45 )
            t_step = t_step * dh;
        //iterations++;
    }
    //Log::L3( "Calculating Definite Integral from t0 = {} to t1 = {} at initial step = {} -- Done.\n", t0, t1, step );
    //Log::L3( "Done Chi integrating, (t0 = {}) t1 was supposted to be {} and ended up {}, maxerror was {}, did {}/{} iterations (var/const dt) \n", t0, t1, t_t, maxerror, iterations, std::floor( t1 - t0 / step ) );
    //Log::L3( "End Rho:\n{}\n", Dense( rho ) );
    return SaveState( rho, t1 );
}

//const unsigned int Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM = 1111111111;