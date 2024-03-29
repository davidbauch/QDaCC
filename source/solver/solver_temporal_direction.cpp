#include "solver/solver_ode.h"

bool QDACC::Numerics::ODESolver::calculate_t_direction( System &s ) {
    MatrixMain rho = s.operatorMatrices.rho;

    Timer &rkTimer = Timers::create( "RungeKutta-Main-Loop" ).start();
    ProgressBar progressbar = ProgressBar();

    Log::L2( "[Solver] Calculating t-direction from {} to {} at stepsize {}...\n", s.parameters.t_start, s.parameters.t_end, s.parameters.t_step );

    // Calculate Time evolution on time vector timestamps.
    if ( s.parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) {
        calculate_path_integral( rho, s.parameters.t_start, s.parameters.t_end, rkTimer, progressbar, "T-Direction: ", s, savedStates, true );
        if (s.parameters.numerics_pathintegral_set_couplings_zero)
            for (auto& el : s.operatorMatrices.phonon_group_index_to_coupling_value)
                el = 0;
    } else {
        // If using RK45 and calculating correlation functions, disable USING cached functions while enabling caching them, then reenable using the cached matrices for Tau direction
        if ( s.parameters.numerics_rk_order == 45 and not s.parameters.input_correlation.empty() and s.parameters.numerics_use_saved_coefficients ) {
            s.parameters.numerics_enable_saving_coefficients = true;
            s.parameters.numerics_use_saved_coefficients = false;
            Log::L2( "[Solver] Disabling usage of saved phonon matrix coefficients while forcing matrix caching.\n" );
        } else if ( s.parameters.numerics_rk_order != 45 ) {
            // Adjust Temporal Stuff and Grid if no RK45 is used:
            s.parameters.post_adjust_input();
        }
        calculate_runge_kutta( rho, s.parameters.t_start, s.parameters.t_end, rkTimer, progressbar, "T-Direction", s, savedStates, true );
        // Reenabling if caching was disabled.
        if ( s.parameters.numerics_enable_saving_coefficients ) {
            s.parameters.numerics_enable_saving_coefficients = false;
            s.parameters.numerics_use_saved_coefficients = true;
            Log::L2( "[Solver] Re-enabling usage of saved phonon matrix coefficients. Cached {} matrices.\n", s.savedCoefficients.size() );
        }
    }
    // Finalize
    rkTimer.end();
    Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), rkTimer.getTotalIterationNumber(), "T-Direction: ", Timers::PROGRESS_FORCE_OUTPUT );

    // Post-Temporal Direction parameter adjustments. Also generates the Grid.
    s.parameters.post_adjust_input();

    Log::L2( "[Solver] Saved {} states.\n", savedStates.size() );
    Log::L2( "[Solver] Hamiltons: Attempts w/r: {}, Write: {}, Calc: {}, Read: {}.\n", track_gethamilton_calcattempt, track_gethamilton_write, track_gethamilton_calc, track_gethamilton_read );
    size_t sum = 0;
    std::ranges::for_each( s.savedCoefficients.begin(), s.savedCoefficients.end(), [&]( const std::pair<double, std::map<double, QDACC::SaveStateTau>> &m ) { sum += m.second.size(); } );
    Log::L2( "[Solver] Cached {} phonon matrices.\n", sum );

    // Interpolate Outputstates with (spline) interpolation. We interpolate seperately from the grid interpolation to allow for different user options.
    auto output_states = s.parameters.numerics_interpolate_outputs ? Numerics::interpolate_curve( savedStates, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, s.parameters.numerics_maximum_primary_threads, s.parameters.numerics_interpolate_method_time ) : savedStates;
    // Interpolate savedStates to the grid
    savedStates = Numerics::interpolate_curve( savedStates, s.parameters.t_start, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau ); // Numerics::interpolate_curve(savedStates, s.parameters.t_start, s.parameters.t_end, s.parameters.t_step, s.parameters.numerics_maximum_primary_threads, 0);

    s.parameters.iterations_t_max = savedStates.size();

    // Index Map:
    Log::L2( "[Solver] Creating Index map for {} density matrices...\n", savedStates.size() );
    // Output Warning if size of cached rhos is too large, indicating interpolation is necessary.
    if (savedStates.size() > 20000)
        Log::L2("[Solver] WARNING Size of saved states is very large! It may be advisable to interpolate to lower dimensions.\n");
    for ( int i = 0; i < savedStates.size(); i++ ) {
        const auto current_time = savedStates.at( i ).t;
        rho_index_map[current_time] = i;
    }
    Log::L2( "[Solver] Saved {} t-index pairs.\n", rho_index_map.size() );

    // Calculate expectation values
    Timer &evalTimer = Timers::create( "Expectation-Value-Loop" );
    evalTimer.start();
    Log::L2( "[Solver] Calculating expectation values...\n" );
    s.calculate_expectation_values( output_states, evalTimer );
    evalTimer.end();

    // Switch main direction trigger
    s.parameters.numerics_main_direction_done = true;

    return true;
}