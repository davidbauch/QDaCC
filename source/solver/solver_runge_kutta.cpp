#include "solver/solver.h"

MatType ODESolver::iterateRungeKutta4( const MatType &rho, System &s, const double t, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-4 ausrechnen
    MatType H_calc_k1 = getHamilton( s, t );
    MatType H_calc_k23 = getHamilton( s, t + s.parameters.t_step * 0.5 );
    MatType H_calc_k4 = getHamilton( s, t + s.parameters.t_step );
    // k1-4 ausrechnen
    MatType rk1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    MatType rk2 = s.dgl_rungeFunction( rho + s.parameters.t_step * 0.5 * rk1, H_calc_k23, t + s.parameters.t_step * 0.5, savedStates );
    MatType rk3 = s.dgl_rungeFunction( rho + s.parameters.t_step * 0.5 * rk2, H_calc_k23, t + s.parameters.t_step * 0.5, savedStates );
    MatType rk4 = s.dgl_rungeFunction( rho + s.parameters.t_step * rk3, H_calc_k4, t + s.parameters.t_step, savedStates );
    // Dichtematrix
    return rho + s.parameters.t_step / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
}

MatType ODESolver::iterateRungeKutta5( const MatType &rho, System &s, const double t, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-6 ausrechnen
    MatType H_calc_k1 = getHamilton( s, t );
    MatType H_calc_k2 = getHamilton( s, t + a2 * s.parameters.t_step );
    MatType H_calc_k3 = getHamilton( s, t + a3 * s.parameters.t_step );
    MatType H_calc_k4 = getHamilton( s, t + a4 * s.parameters.t_step );
    MatType H_calc_k5 = getHamilton( s, t + a5 * s.parameters.t_step );
    MatType H_calc_k6 = getHamilton( s, t + a6 * s.parameters.t_step );
    // k1-6 ausrechnen
    MatType k1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    MatType k2 = s.dgl_rungeFunction( rho + s.parameters.t_step * b11 * k1, H_calc_k2, t + a2 * s.parameters.t_step, savedStates );
    MatType k3 = s.dgl_rungeFunction( rho + s.parameters.t_step * ( b21 * k1 + b22 * k2 ), H_calc_k3, t + a3 * s.parameters.t_step, savedStates );
    MatType k4 = s.dgl_rungeFunction( rho + s.parameters.t_step * ( b31 * k1 + b32 * k2 + b33 * k3 ), H_calc_k4, t + a4 * s.parameters.t_step, savedStates );
    MatType k5 = s.dgl_rungeFunction( rho + s.parameters.t_step * ( b41 * k1 + b42 * k2 + b43 * k3 + b44 * k4 ), H_calc_k5, t + a5 * s.parameters.t_step, savedStates );
    MatType k6 = s.dgl_rungeFunction( rho + s.parameters.t_step * ( b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4 + b55 * k5 ), H_calc_k6, t + a6 * s.parameters.t_step, savedStates );
    // Dichtematrix
    return rho + s.parameters.t_step * ( b61 * k1 + b63 * k3 + b64 * k4 + b65 * k5 + b66 * k6 );
}

MatType ODESolver::iterate( const MatType &rho, System &s, const double t, std::vector<SaveState> &savedStates, const int dir ) {
    int order = dir == DIR_T ? s.parameters.numerics_order_t : s.parameters.numerics_order_tau;
    if ( order == 4 ) {
        return iterateRungeKutta4( rho, s, t, savedStates );
    }
    return iterateRungeKutta5( rho, s, t, savedStates );
}


// Wenn timestamps leer ist (size = 0), dann einfach rk steps ausgeben
bool ODESolver::calculate_runge_kutta( MatType &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output ) {
    // Reserve Output Vector
    int t_iterations_max = (int)std::ceil( (t_end-t_start)/t_step_initial ) + 1;
    output.reserve( t_iterations_max );
    // Save initial value
    saveState( rho0, t_start, output );
    // Calculate Remaining
    MatType rho = rho0;
    for ( double t_t = t_start + t_step_initial; t_t <= t_end; t_t += t_step_initial ) {
        // Runge-Kutta iteration
        rho = iterate( rho, s, t_t, output );
        // Save Rho
        saveState( rho, t_t, output );
        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, t_iterations_max, progressbar_name );
        }
    }
    return true;
}