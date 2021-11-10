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
        calculate_runge_kutta( rho, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );
    }
    // Finalize
    rkTimer.end();
    Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, rkTimer.getTotalIterationNumber(), rkTimer.getTotalIterationNumber(), "T-Direction: ", Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[Solver] Done! Saved {} states.\n", savedStates.size() );
    Log::L2( "[Solver] Hamiltons: Attempts w/r: {}, Write: {}, Calc: {}, Read: {}. Done!\n", track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_calc, track_gethamilton_read );

    // Interpolate Matrices? //TODO: das hier und RK45 interpolation ist doppelt. wenn rk45, dann hier interpolieren und am ende savedStates = interpolatedSavedStates
    std::vector<QDLC::SaveState> interpolate_savedstates;
    if ( s.parameters.numerics_interpolate_outputs ) {
        Log::L2( "[Solver] Calculating interpolated matrices for temporal properties...\n" );
        interpolate_savedstates = QDLC::Numerics::calculate_smooth_curve( savedStates, s.parameters.t_start, s.parameters.t_end, std::max( s.parameters.iterations_t_max * 5, 2500 ), s.parameters.output_handlerstrings );
        Log::L2( "[Solver] Done!\n" );
    } else {
        Log::L2( "[Solver] Using the non-interpolated matrices for temporal properties\n" );
        interpolate_savedstates = savedStates;
    }

    // Calculate expectation values
    Timer &evalTimer = Timers::create( "Expectation-Value-Loop" );
    evalTimer.start();
    Log::L2( "[Solver] Calculating expectation values...\n" );
    s.expectationValues( interpolate_savedstates, evalTimer );
    evalTimer.end();
    Log::L2( "[Solver] Done!\n" );

    return true;
}