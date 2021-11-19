#include "solver/solver_ode.h"
#include "solver/solver.h"

Sparse QDLC::Numerics::ODESolver::iterateRungeKutta4( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDLC::SaveState> &savedStates ) {
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

Sparse QDLC::Numerics::ODESolver::iterateRungeKutta5( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDLC::SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-6 ausrechnen
    Sparse H_calc_k1 = getHamilton( s, t );
    Sparse H_calc_k2 = getHamilton( s, t + RKCoefficients::a2 * t_step );
    Sparse H_calc_k3 = getHamilton( s, t + RKCoefficients::a3 * t_step );
    Sparse H_calc_k4 = getHamilton( s, t + RKCoefficients::a4 * t_step );
    Sparse H_calc_k5 = getHamilton( s, t + RKCoefficients::a5 * t_step );
    Sparse H_calc_k6 = getHamilton( s, t + RKCoefficients::a6 * t_step );
    // k1-6 ausrechnen
    Sparse k1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    Sparse k2 = s.dgl_rungeFunction( rho + t_step * RKCoefficients::b11 * k1, H_calc_k2, t + RKCoefficients::a2 * t_step, savedStates );
    Sparse k3 = s.dgl_rungeFunction( rho + t_step * ( RKCoefficients::b21 * k1 + RKCoefficients::b22 * k2 ), H_calc_k3, t + RKCoefficients::a3 * t_step, savedStates );
    Sparse k4 = s.dgl_rungeFunction( rho + t_step * ( RKCoefficients::b31 * k1 + RKCoefficients::b32 * k2 + RKCoefficients::b33 * k3 ), H_calc_k4, t + RKCoefficients::a4 * t_step, savedStates );
    Sparse k5 = s.dgl_rungeFunction( rho + t_step * ( RKCoefficients::b41 * k1 + RKCoefficients::b42 * k2 + RKCoefficients::b43 * k3 + RKCoefficients::b44 * k4 ), H_calc_k5, t + RKCoefficients::a5 * t_step, savedStates );
    Sparse k6 = s.dgl_rungeFunction( rho + t_step * ( RKCoefficients::b51 * k1 + RKCoefficients::b52 * k2 + RKCoefficients::b53 * k3 + RKCoefficients::b54 * k4 + RKCoefficients::b55 * k5 ), H_calc_k6, t + RKCoefficients::a6 * t_step, savedStates );
    // Dichtematrix
    return rho + t_step * ( RKCoefficients::b61 * k1 + RKCoefficients::b63 * k3 + RKCoefficients::b64 * k4 + RKCoefficients::b65 * k5 + RKCoefficients::b66 * k6 );
}

std::pair<Sparse, double> QDLC::Numerics::ODESolver::iterateRungeKutta45( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDLC::SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-6 ausrechnen
    Sparse H_calc_k1 = getHamilton( s, t );
    Sparse H_calc_k2 = getHamilton( s, t + RKCoefficients::a2 * t_step );
    Sparse H_calc_k3 = getHamilton( s, t + RKCoefficients::a3 * t_step );
    Sparse H_calc_k4 = getHamilton( s, t + RKCoefficients::a4 * t_step );
    Sparse H_calc_k5 = getHamilton( s, t + RKCoefficients::a5 * t_step );
    Sparse H_calc_k6 = getHamilton( s, t + RKCoefficients::a6 * t_step );
    Sparse H_calc_k7 = getHamilton( s, t + t_step );
    // k1-6 ausrechnen
    Sparse k1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    Sparse k2 = s.dgl_rungeFunction( rho + t_step * RKCoefficients::b11 * k1, H_calc_k2, t + RKCoefficients::a2 * t_step, savedStates );
    Sparse k3 = s.dgl_rungeFunction( rho + t_step * ( RKCoefficients::b21 * k1 + RKCoefficients::b22 * k2 ), H_calc_k3, t + RKCoefficients::a3 * t_step, savedStates );
    Sparse k4 = s.dgl_rungeFunction( rho + t_step * ( RKCoefficients::b31 * k1 + RKCoefficients::b32 * k2 + RKCoefficients::b33 * k3 ), H_calc_k4, t + RKCoefficients::a4 * t_step, savedStates );
    Sparse k5 = s.dgl_rungeFunction( rho + t_step * ( RKCoefficients::b41 * k1 + RKCoefficients::b42 * k2 + RKCoefficients::b43 * k3 + RKCoefficients::b44 * k4 ), H_calc_k5, t + RKCoefficients::a5 * t_step, savedStates );
    Sparse k6 = s.dgl_rungeFunction( rho + t_step * ( RKCoefficients::b51 * k1 + RKCoefficients::b52 * k2 + RKCoefficients::b53 * k3 + RKCoefficients::b54 * k4 + RKCoefficients::b55 * k5 ), H_calc_k6, t + RKCoefficients::a6 * t_step, savedStates );

    Sparse drho = RKCoefficients::b61 * k1 + RKCoefficients::b63 * k3 + RKCoefficients::b64 * k4 + RKCoefficients::b65 * k5 + RKCoefficients::b66 * k6;
    Sparse ret = rho + t_step * drho;
    // Error
    Sparse k7 = s.dgl_rungeFunction( ret, H_calc_k7, t + t_step, savedStates );
    dSparse errmat = ( drho - ( k1 * RKCoefficients::e1 + k3 * RKCoefficients::e3 + k4 * RKCoefficients::e4 + k5 * RKCoefficients::e5 + k6 * RKCoefficients::e6 + k7 * RKCoefficients::e7 ) ).cwiseAbs2();
    double err = errmat.sum() / drho.cwiseAbs2().sum();

    // Dichtematrix
    return std::make_pair( ret, err );
}

Sparse QDLC::Numerics::ODESolver::iterate( const Sparse &rho, System &s, const double t, const double t_step, std::vector<QDLC::SaveState> &savedStates, const int dir ) {
    if ( s.parameters.numerics_rk_order == 4 )
        return iterateRungeKutta4( rho, s, t, t_step, savedStates );
    return iterateRungeKutta5( rho, s, t, t_step, savedStates );
}

bool QDLC::Numerics::ODESolver::calculate_runge_kutta( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output ) {
    if ( s.parameters.numerics_rk_order == 45 ) {
        calculate_runge_kutta_45( rho0, t_start, t_end, rkTimer, progressbar, progressbar_name, s, output, do_output, s.parameters.numerics_rk_interpolate, s.parameters.numerics_rk_tol );
        return true;
    }
    size_t t_index = std::min<size_t>(size_t(std::lower_bound(s.parameters.grid_values.begin(), s.parameters.grid_values.end(), t_start) - s.parameters.grid_values.begin()),s.parameters.grid_values.size()-2);//s.parameters.grid_value_indices[t_start];
    double t_step_initial = s.parameters.grid_steps[t_index];
    //Log::L2("t_index = {}, t_step_initial = {}\n",t_index, t_step_initial);

    // Reserve Output Vector
    output.reserve( s.parameters.iterations_t_max + 1 );
    // Save initial value
    saveState( rho0, t_start, output );
    // Calculate Remaining
    Sparse rho = rho0;
    for ( double t_t = t_start; t_t <= t_end; t_t += t_step_initial ) {
        // Runge-Kutta iteration
        rho = iterate( rho, s, t_t, t_step_initial, output );
        // Save Rho
        saveState( rho, t_t + t_step_initial, output );
        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, rkTimer.getTotalIterationNumber() , s.parameters.iterations_t_max, progressbar_name );
        }
        // Adjust t_end until ground state is reached. we assume the ground state is the first entry of the DM
        if ( s.parameters.numerics_calculate_till_converged and t_t + t_step_initial > t_end and std::real( rho.coeff( 0, 0 ) ) < 0.999 ) {
            Log::L3( "[RKSOLVER] Adjusted Calculation end from {} to {}\n", t_end, t_end + 10.0 * t_step_initial );
            t_end += 10.0 * t_step_initial;
            s.parameters.iterations_t_max = (int)std::ceil( ( t_end - t_start ) / t_step_initial );
        }
        t_index = std::min<size_t>(t_index+1, s.parameters.grid_values.size()-2);
        t_step_initial = s.parameters.grid_steps[t_index];
    }
    if ( s.parameters.numerics_calculate_till_converged ) {
        s.parameters.numerics_calculate_till_converged = false;
        s.parameters.t_end = t_end;
        s.parameters.iterations_t_max = output.size();
        s.parameters.iterations_t_skip = std::max( 1.0, std::ceil( 1.0 * s.parameters.iterations_t_max / s.parameters.iterations_tau_resolution ) );
        //TODO: für den t_step vector und so solle hier s.parameters.adjust aufgerufen werden?!
        reset( s );
        Log::L1( "[RKSOLVER] Adjusted t_end to {}.\n", s.parameters.t_end );
    }
    //TODO: wenn step im grid < dt, dann interpolieren? kp. oder immer RK45? dann hätte man das problem gar nicht.
    // oder, wenn gridstep > step, dann rechne mit step aber speicher gridstep? aber dann was wenn gridstep/step nicht grade ist? keeeeekw
    // aktuell ist es so, das man genauso viele rhos speichern muss wie man am ende braucht, also die zwischen den skips liegen NICHT speichern, da die korrelationsfunktionen keine ahnung haben wie sie an den skip kommen.
    return true;
}

bool QDLC::Numerics::ODESolver::calculate_runge_kutta_45( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output, bool interpolate, double tolerance ) {
    size_t t_index = std::min<size_t>(size_t(std::lower_bound(s.parameters.grid_values.begin(), s.parameters.grid_values.end(), t_start) - s.parameters.grid_values.begin()),s.parameters.grid_values.size()-2);//s.parameters.grid_value_indices[t_start];
    double t_step_initial = s.parameters.t_step; //s.parameters.grid_steps[t_index];
    double t_step = t_step_initial;
    // Reserve Output Vector
    output.reserve( s.parameters.iterations_t_max + 1 );
    // Create temporary out vector:
    std::vector<QDLC::SaveState> temp;

    // Save initial value
    saveState( rho0, t_start, temp );
    double t_t = t_start; // + s.parameters.numerics_rk_stepdelta; // + t_step;
    while ( t_t <= t_end ) {
        // Runge-Kutta iteration
        auto rkret = iterateRungeKutta45( temp.back().mat, s, t_t, t_step, temp );
        double error = rkret.second;
        double dh = std::pow( tolerance / 2. / std::max( error, 1E-15 ), 0.25 );
        if ( std::isnan( dh ) )
            dh = 1.0;
        double t_step_new = t_step * dh; //FIXME: ohne diskrete schritte kann das hier null werden, check machen!
        bool accept = true;
        if ( error >= tolerance ) {
            accept = false;
        }
        if ( s.parameters.numerics_rk_usediscrete_timesteps ) {
            t_step_new = t_step;
            if ( dh < 1 ) {
                t_step_new -= s.parameters.numerics_rk_stepdelta * std::floor( 1.0 / dh );
            } else {
                t_step_new += s.parameters.numerics_rk_stepdelta * std::floor( dh );
            }
            if ( t_step_new < s.parameters.numerics_rk_stepmin ) {
                t_step_new = s.parameters.numerics_rk_stepmin;
                accept = true;
            } else if ( t_step_new > s.parameters.numerics_rk_stepmax ) {
                t_step_new = s.parameters.numerics_rk_stepmax;
                accept = true;
            }
        }
        Log::L3( "[RK45SOLVER{}] (t = {}) - Local error: {} - dh = {}, current timestep is: {}, new timestep will be: {}, accept current step = {}\n", omp_get_thread_num(), t_t, error, dh, t_step, t_step_new, accept );
        if ( accept ) {
            t_t += t_step;
            saveState( rkret.first, t_t, temp );
            Log::L3( "[RK45SOLVER{}] --> (t = {}) - Accepdet step - Local error: {} - current timestep: {}, dh = {}\n", omp_get_thread_num(), t_t, error, t_step, dh );
            // Progress and time output
            rkTimer.iterate();
            if ( do_output ) {
                Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, temp.size(), s.parameters.iterations_t_max, progressbar_name );
            }
            // Adjust t_end until ground state is reached.we assume the ground state is the first entry of the DM
            if ( s.parameters.numerics_calculate_till_converged and t_t + t_step > t_end and std::real( temp.back().mat.coeff( 0, 0 ) ) < 0.999 ) {
                t_end += 10.0 * t_step;
                Log::L3( "[RK45SOLVER{}] Adjusted Calculation end to {}\n", omp_get_thread_num(), t_end );
            }
        }
        t_step = t_step_new;
    }
    // Interpolate if neccessary
    if ( s.parameters.numerics_rk_interpolate ) {
        output.clear();
        output = Numerics::interpolate_curve( temp, t_start, t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false );
    } else {
        output = temp;
    }
    if ( s.parameters.numerics_calculate_till_converged ) {
        s.parameters.numerics_calculate_till_converged = false;
        s.parameters.t_end = t_end;
        s.parameters.iterations_t_max = output.size();
        s.parameters.iterations_t_skip = std::max( 1.0, std::ceil( 1.0 * s.parameters.iterations_t_max / s.parameters.iterations_tau_resolution ) );
        //TODO: für den t_step vector und so solle hier s.parameters.adjust aufgerufen werden?!
        reset( s );
        Log::L1( "[RK45SOLVER] Adjusted t_end to {}.\n", s.parameters.t_end );
    }
    if ( do_output )
        Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, rkTimer.getTotalIterationNumber(), rkTimer.getTotalIterationNumber(), progressbar_name );
    return true;
}