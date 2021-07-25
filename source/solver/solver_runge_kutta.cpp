#include "solver/solver_ode.h"
#include "solver/solver.h"

Sparse ODESolver::iterateRungeKutta4( const Sparse &rho, System &s, const double t, const double t_step, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-4 ausrechnen
    Sparse H_calc_k1 = getHamilton( s, t );
    Sparse H_calc_k23 = getHamilton( s, t + t_step * 0.5 );
    Sparse H_calc_k4 = getHamilton( s, t + t_step );
    // k1-4 ausrechnen
    Sparse rk1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    Sparse rk2 = s.dgl_rungeFunction( rho + t_step * 0.5 * rk1, H_calc_k23, t + t_step * 0.5, savedStates );
    Sparse rk3 = s.dgl_rungeFunction( rho + t_step * 0.5 * rk2, H_calc_k23, t + t_step * 0.5, savedStates );
    Sparse rk4 = s.dgl_rungeFunction( rho + t_step * rk3, H_calc_k4, t + t_step, savedStates );
    // Dichtematrix
    return rho + t_step / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
}

Sparse ODESolver::iterateRungeKutta5( const Sparse &rho, System &s, const double t, const double t_step, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-6 ausrechnen
    Sparse H_calc_k1 = getHamilton( s, t );
    Sparse H_calc_k2 = getHamilton( s, t + Solver::a2 * t_step );
    Sparse H_calc_k3 = getHamilton( s, t + Solver::a3 * t_step );
    Sparse H_calc_k4 = getHamilton( s, t + Solver::a4 * t_step );
    Sparse H_calc_k5 = getHamilton( s, t + Solver::a5 * t_step );
    Sparse H_calc_k6 = getHamilton( s, t + Solver::a6 * t_step );
    // k1-6 ausrechnen
    Sparse k1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    Sparse k2 = s.dgl_rungeFunction( rho + t_step * Solver::b11 * k1, H_calc_k2, t + Solver::a2 * t_step, savedStates );
    Sparse k3 = s.dgl_rungeFunction( rho + t_step * ( Solver::b21 * k1 + Solver::b22 * k2 ), H_calc_k3, t + Solver::a3 * t_step, savedStates );
    Sparse k4 = s.dgl_rungeFunction( rho + t_step * ( Solver::b31 * k1 + Solver::b32 * k2 + Solver::b33 * k3 ), H_calc_k4, t + Solver::a4 * t_step, savedStates );
    Sparse k5 = s.dgl_rungeFunction( rho + t_step * ( Solver::b41 * k1 + Solver::b42 * k2 + Solver::b43 * k3 + Solver::b44 * k4 ), H_calc_k5, t + Solver::a5 * t_step, savedStates );
    Sparse k6 = s.dgl_rungeFunction( rho + t_step * ( Solver::b51 * k1 + Solver::b52 * k2 + Solver::b53 * k3 + Solver::b54 * k4 + Solver::b55 * k5 ), H_calc_k6, t + Solver::a6 * t_step, savedStates );
    // Dichtematrix
    return rho + t_step * ( Solver::b61 * k1 + Solver::b63 * k3 + Solver::b64 * k4 + Solver::b65 * k5 + Solver::b66 * k6 );
}

std::pair<Sparse, double> ODESolver::iterateRungeKutta45( const Sparse &rho, System &s, const double t, const double t_step, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-6 ausrechnen
    Sparse H_calc_k1 = getHamilton( s, t );
    Sparse H_calc_k2 = getHamilton( s, t + Solver::a2 * t_step );
    Sparse H_calc_k3 = getHamilton( s, t + Solver::a3 * t_step );
    Sparse H_calc_k4 = getHamilton( s, t + Solver::a4 * t_step );
    Sparse H_calc_k5 = getHamilton( s, t + Solver::a5 * t_step );
    Sparse H_calc_k6 = getHamilton( s, t + Solver::a6 * t_step );
    // k1-6 ausrechnen
    Sparse k1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    Sparse k2 = s.dgl_rungeFunction( rho + t_step * Solver::b11 * k1, H_calc_k2, t + Solver::a2 * t_step, savedStates );
    Sparse k3 = s.dgl_rungeFunction( rho + t_step * ( Solver::b21 * k1 + Solver::b22 * k2 ), H_calc_k3, t + Solver::a3 * t_step, savedStates );
    Sparse k4 = s.dgl_rungeFunction( rho + t_step * ( Solver::b31 * k1 + Solver::b32 * k2 + Solver::b33 * k3 ), H_calc_k4, t + Solver::a4 * t_step, savedStates );
    Sparse k5 = s.dgl_rungeFunction( rho + t_step * ( Solver::b41 * k1 + Solver::b42 * k2 + Solver::b43 * k3 + Solver::b44 * k4 ), H_calc_k5, t + Solver::a5 * t_step, savedStates );
    Sparse k6 = s.dgl_rungeFunction( rho + t_step * ( Solver::b51 * k1 + Solver::b52 * k2 + Solver::b53 * k3 + Solver::b54 * k4 + Solver::b55 * k5 ), H_calc_k6, t + Solver::a6 * t_step, savedStates );

    Sparse drho = Solver::b61 * k1 + Solver::b63 * k3 + Solver::b64 * k4 + Solver::b65 * k5 + Solver::b66 * k6;
    Sparse ret = rho + t_step * drho;
    // Error
    Sparse H_calc_k7 = getHamilton( s, t + t_step );
    Sparse k7 = s.dgl_rungeFunction( ret, H_calc_k7, t + t_step, savedStates );
    dSparse errmat = ( drho - ( k1 * Solver::e1 + k3 * Solver::e3 + k4 * Solver::e4 + k5 * Solver::e5 + k6 * Solver::e6 + k7 * Solver::e7 ) ).cwiseAbs2();
    double err = errmat.sum() / drho.cwiseAbs2().sum();

    // Dichtematrix
    return std::make_pair( ret, err );
}

Sparse ODESolver::iterate( const Sparse &rho, System &s, const double t, const double t_step, std::vector<SaveState> &savedStates, const int dir ) {
    if ( s.parameters.numerics_rk_order == 4 )
        return iterateRungeKutta4( rho, s, t, t_step, savedStates );
    return iterateRungeKutta5( rho, s, t, t_step, savedStates );
}

bool ODESolver::calculate_runge_kutta( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output ) {
    if ( s.parameters.numerics_rk_order == 45 ) {
        calculate_runge_kutta_45( rho0, t_start, t_end, t_step_initial, rkTimer, progressbar, progressbar_name, s, output, do_output, s.parameters.numerics_rk_interpolate, s.parameters.numerics_rk_tol );
        return true;
    }

    Log::L3( "Setting up Runge-Kutta Solver...\n" );
    // Reserve Output Vector
    output.reserve( s.parameters.iterations_t_max + 1 );
    // Save initial value
    saveState( rho0, t_start, output );
    // Calculate Remaining
    Sparse rho = rho0;
    Log::L3( "Calculating Runge-Kutta Loop...\n" );
    for ( double t_t = t_start; t_t <= t_end; t_t += t_step_initial ) {
        // Runge-Kutta iteration
        rho = iterate( rho, s, t_t, t_step_initial, output );
        // Save Rho
        saveState( rho, t_t + t_step_initial, output );
        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
        }
    }
    Log::L3( "Done!\n" );
    return true;
}

bool ODESolver::calculate_runge_kutta_45( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output, bool interpolate, double tolerance ) {
    double t_step = t_step_initial;
    // Reserve Output Vector
    output.reserve( s.parameters.iterations_t_max + 1 );
    // Create temporary out vector:
    std::vector<SaveState> temp;

    // Save initial value
    saveState( rho0, t_start, temp );
    double t_t = t_start; // + s.parameters.numerics_rk_stepdelta; // + t_step;
    Sparse rho;
    while ( t_t <= t_end ) {
        // Runge-Kutta iteration
        auto rkret = iterateRungeKutta45( temp.back().mat, s, t_t, t_step, temp );
        double error = rkret.second;
        double dh = std::pow( tolerance / 2. / error, 0.25 );
        bool accept = false;
        if ( error < tolerance ) {
            accept = true;
        }
        double t_step_new = t_step * dh;
        if ( s.parameters.numerics_rk_usediscrete_timesteps ) {
            t_step_new = t_step;
            if ( dh < 1 and t_step_new - s.parameters.numerics_rk_stepdelta > 0 and t_step_new - s.parameters.numerics_rk_stepdelta > s.parameters.numerics_rk_stepmin ) {
                t_step_new -= s.parameters.numerics_rk_stepdelta * std::floor( 1.0 / dh );
            } else if ( dh > 1 and t_step_new + s.parameters.numerics_rk_stepdelta < s.parameters.numerics_rk_stepmax ) {
                t_step_new += s.parameters.numerics_rk_stepdelta * std::floor( dh );
            }
            //} else
            //    accept = true;
        }
        Log::L3( "[t = {}] - - - Local error: {} - dh = {}, current timestep is: {}, new timestep will be: {}\n", t_t, error, dh, t_step, t_step_new );
        if ( accept ) {
            t_t += t_step;
            saveState( rkret.first, t_t, temp );
            Log::L3( "[t = {}] - Accepdet step - Local error: {} - current timestep: {}, dh = {}\n", t_t, error, t_step, dh );
            // Progress and time output
            rkTimer.iterate();
            if ( do_output ) {
                Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
            }
        }
        t_step = t_step_new;
    }
    // Interpolate if neccessary
    if ( s.parameters.numerics_rk_interpolate ) {
        output.clear();
        // Do a very simple linear interpolation:
        int i = 1;
        for ( double t_t = t_start; t_t < t_end; t_t += t_step_initial ) {
            while ( i < temp.size() - 1 and t_t > temp[i].t ) {
                i++;
            }
            double first = i > 0 ? temp[i - 1].t : 0.0;
            double second = temp[i].t;
            double f = ( t_t - first ) / ( second - first );
            saveState( temp[i - 1].mat + f * ( temp[i].mat - temp[i - 1].mat ), t_t, output );
        }
    } else {
        output = temp;
    }
    return true;
}