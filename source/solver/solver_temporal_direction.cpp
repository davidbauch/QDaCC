#include "solver/solver_ode.h"

bool ODESolver::calculate_t_direction( System &s ) {
    MatType rho = s.operatorMatrices.rho;
    
    Timer &rkTimer = createTimer( "RungeKutta-Main-Loop" );
    ProgressBar progressbar = ProgressBar( s.parameters.iterations_t_max, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    rkTimer.start();
    
    logs.level2( "Calculating t-direction from {} to {} at stepsize {}... ", s.parameters.t_start, s.parameters.t_end, s.parameters.t_step );

    // Calculate Time evolution on time vector timestamps.
    calculate_runge_kutta( rho, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );

    // Finalize
    outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, "T-Direction: ", PROGRESS_FORCE_OUTPUT );
    rkTimer.end();
    logs.level2( "Done! Saved {} states, cached {} hamilton matrices.\n", savedStates.size(), savedHamiltons.size() );
    
    
    // Interpolate Matrices?
    std::vector<SaveState> interpolate_savedstates;
    if ( s.parameters.numerics_interpolate_outputs ) {
        logs.level2("Calculating interpolated matrices for temporal properties...\n");
        interpolate_savedstates = Solver::calculate_smooth_curve( savedStates, s.parameters.t_start, s.parameters.t_end, std::max(s.parameters.iterations_t_max*5,2500), s.parameters.output_handlerstrings );
        logs.level2("Done!\n");
    } else {
        logs.level2("Using the non-interpolated matrices for temporal properties\n");
        interpolate_savedstates = savedStates;
    }

    // Calculate expectation values
    logs.level2("Calculating expectation values...\n");
    for ( long unsigned int i = 0; i < interpolate_savedstates.size(); i++) {
        s.traceValid( interpolate_savedstates.at(i).mat, interpolate_savedstates.at(i).t );
        s.expectationValues( interpolate_savedstates.at(i).mat, interpolate_savedstates.at(i).t, interpolate_savedstates );
    }
    logs.level2("Done!\n");
    return true;
}