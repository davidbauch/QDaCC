#include "solver/solver_ode.h"

bool QDLC::Numerics::ODESolver::calculate_t_direction( System &s ) {
    Sparse rho = s.operatorMatrices.rho;

    Timer &rkTimer = Timers::create( "RungeKutta-Main-Loop" );
    ProgressBar progressbar = ProgressBar();
    rkTimer.start();

    Log::L2( "[Solver] Calculating t-direction from {} to {} at stepsize {}...\n", s.parameters.t_start, s.parameters.t_end, s.parameters.t_step );

    // Calculate Time evolution on time vector timestamps.
    //TODO: if pathintegral then calculate_path_integral else calculate_runge_kutta
    if ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
        calculate_path_integral( rho, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );
    } else {
        calculate_runge_kutta( rho, s.parameters.t_start, s.parameters.t_end, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );
    }
    // Finalize
    rkTimer.end();
    Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), rkTimer.getTotalIterationNumber(), "T-Direction: ", Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[Solver] Done! Saved {} states.\n", savedStates.size() );
    Log::L2( "[Solver] Hamiltons: Attempts w/r: {}, Write: {}, Calc: {}, Read: {}. Done!\n", track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_calc, track_gethamilton_read );

    // Index Map:
    for (int i = 0; i < savedStates.size(); i++) {
        rho_index_map[getTimeAt(i)] = i;
    }
    // Calculate expectation values
    Timer &evalTimer = Timers::create( "Expectation-Value-Loop" );
    evalTimer.start();
    Log::L2( "[Solver] Calculating expectation values...\n" );
    //TODO: interpolation sollte ausschaltbar/wechselbar per parameter sien.
    s.expectationValues( ( s.parameters.input_correlation_resolution.count("Modified") ? Numerics::interpolate_curve(savedStates, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, s.parameters.numerics_maximum_threads, 3) : savedStates ), evalTimer );
    evalTimer.end();
    Log::L2( "[Solver] Done!\n" );

    return true;
}