#include "solver/solver_ode.h"

bool QDLC::Numerics::ODESolver::calculate_t_direction( System &s ) {
    Sparse rho = s.operatorMatrices.rho;

    Timer &rkTimer = Timers::create( "RungeKutta-Main-Loop" );
    ProgressBar progressbar = ProgressBar();
    rkTimer.start();

    Log::L2( "[Solver] Calculating t-direction from {} to {} at stepsize {}...\n", s.parameters.t_start, s.parameters.t_end, s.parameters.t_step );

    // Calculate Time evolution on time vector timestamps.
    if ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ) {
        calculate_path_integral( rho, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step_pathint, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );
    } else {
        // If using RK45 and calculating correlation functions, disable USING cached functions while enabling caching them, then reenable using the cached matrices for Tau direction
        if ( s.parameters.numerics_rk_order >= 45 and s.parameters.input_correlation.size() > 0 and s.parameters.numerics_use_saved_coefficients ) {
            s.parameters.numerics_force_caching = true;
            s.parameters.numerics_use_saved_coefficients = false;
        }
        calculate_runge_kutta( rho, s.parameters.t_start, s.parameters.t_end, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );
        // Reenabling if caching was disabled.
        if ( s.parameters.numerics_force_caching ) {
            s.parameters.numerics_force_caching = false;
            s.parameters.numerics_use_saved_coefficients = true;
        }
    }
    // Finalize
    rkTimer.end();
    Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), rkTimer.getTotalIterationNumber(), "T-Direction: ", Timers::PROGRESS_FORCE_OUTPUT );
    Log::L2( "[Solver] Done! Saved {} states.\n", savedStates.size() );
    Log::L2( "[Solver] Hamiltons: Attempts w/r: {}, Write: {}, Calc: {}, Read: {}. Done!\n", track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_calc, track_gethamilton_read );
    size_t sum = 0;
    std::for_each( s.savedCoefficients.begin(), s.savedCoefficients.end(), [&]( const std::pair<double, std::map<double, QDLC::SaveStateTau>> &m ) { sum += m.second.size(); } );
    Log::L2( "[Solver] Cached {} phonon matrices.\n", sum );

    // Interpolate Outputstates with spline interpolation
    auto output_states = s.parameters.numerics_interpolate_outputs ? Numerics::interpolate_curve( savedStates, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, s.parameters.numerics_maximum_threads, s.parameters.numerics_interpolate_method_time ) : savedStates;
    // Interpolate savedStates if RK45 was used or grid was modified
    if ( s.parameters.numerics_rk_order > 5 )
        savedStates = Numerics::interpolate_curve( savedStates, s.parameters.t_start, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau ); // Numerics::interpolate_curve(savedStates, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, s.parameters.numerics_maximum_threads, 0);

    // Index Map:
    for ( int i = 0; i < savedStates.size(); i++ ) {
        rho_index_map[getTimeAt( i )] = i;
    }

    // Calculate expectation values
    Timer &evalTimer = Timers::create( "Expectation-Value-Loop" );
    evalTimer.start();
    Log::L2( "[Solver] Calculating expectation values...\n" );
    s.expectationValues( output_states, evalTimer );
    evalTimer.end();
    Log::L2( "[Solver] Done!\n" );

    return true;
}