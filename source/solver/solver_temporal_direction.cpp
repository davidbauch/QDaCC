#include "solver/solver_ode.h"

bool ODESolver::calculate_t_direction( System &s ) {
    Sparse rho = s.operatorMatrices.rho;

    Timer &rkTimer = Timers::create( "RungeKutta-Main-Loop" );
    ProgressBar progressbar = ProgressBar( s.parameters.iterations_t_max, 60, 0, BAR_VERTICAL, true, 0.1 );
    rkTimer.start();

    Log::L2( "Calculating t-direction from {} to {} at stepsize {}... ", s.parameters.t_start, s.parameters.t_end, s.parameters.t_step );

    // Calculate Time evolution on time vector timestamps.
    //TODO: if pathintegral then calculate_path_integral else calculate_runge_kutta
    if ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
        calculate_path_integral( rho, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );
    } else {
        calculate_runge_kutta( rho, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );
    }
    // Finalize
    Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, "T-Direction: ", PROGRESS_FORCE_OUTPUT );
    rkTimer.end();
    Log::L2( "Done! Saved {} states, cached {} hamilton matrices.\n", savedStates.size(), savedHamiltons.size() );

    // Interpolate Matrices?
    std::vector<SaveState> interpolate_savedstates;
    if ( s.parameters.numerics_interpolate_outputs ) {
        Log::L2( "Calculating interpolated matrices for temporal properties...\n" );
        interpolate_savedstates = Solver::calculate_smooth_curve( savedStates, s.parameters.t_start, s.parameters.t_end, std::max( s.parameters.iterations_t_max * 5, 2500 ), s.parameters.output_handlerstrings );
        Log::L2( "Done!\n" );
    } else {
        Log::L2( "Using the non-interpolated matrices for temporal properties\n" );
        interpolate_savedstates = savedStates;
    }

    // Calculate expectation values
    Timer &evalTimer = Timers::create( "Expectation-Value-Loop" );
    std::string pbname = "Expectation Values: ";
    ProgressBar progressbar2 = ProgressBar( interpolate_savedstates.size(), 60, 0, BAR_VERTICAL, true, 0.1 );
    evalTimer.start();
    Log::L2( "Calculating expectation values...\n" );
    for ( long unsigned int i = 0; i < interpolate_savedstates.size(); i++ ) {
        //if ( !s.traceValid( interpolate_savedstates.at( i ).mat, interpolate_savedstates.at( i ).t ) ) {
        //    return false;
        //}
        s.expectationValues( interpolate_savedstates.at( i ).mat, interpolate_savedstates.at( i ).t, interpolate_savedstates );
        evalTimer.iterate();
        Timers::outputProgress( s.parameters.output_handlerstrings, evalTimer, progressbar2, interpolate_savedstates.size(), pbname );
    }
    Timers::outputProgress( s.parameters.output_handlerstrings, evalTimer, progressbar2, interpolate_savedstates.size(), pbname, PROGRESS_FORCE_OUTPUT );
    evalTimer.end();
    Log::L2( "Done!\n" );

    return true;
}