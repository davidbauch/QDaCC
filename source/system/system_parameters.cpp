#include "system/parameters.h"

namespace std {
double max( const Parameter &p1, const Parameter &p2 ) {
    return max( p2.get(), p1.get() );
}
double min( const Parameter &p1, const Parameter &p2 ) {
    return min( p2.get(), p1.get() );
}
}; // namespace std

Parameters::Parameters( const std::vector<std::string> &arguments ) {
    numerics_calculate_spectrum = 0;
    numerics_calculate_g2 = 0;
    numerics_use_interactionpicture = 0;
    numerics_use_rwa = 0;
    numerics_order_timetrafo = TIMETRANSFORMATION_MATRIXEXPONENTIAL;
    numerics_order_t = 4;
    numerics_order_tau = 4;
    output_advanced_log = 0;
    output_handlerstrings = 0;
    output_operators = 0;
    iterations_t_skip = 1;
    iterations_tau_resolution = 1;
    iterations_w_resolution = 1;
    scale_parameters = false;
    scale_value = 1E12;
    iterations_t_max = 1;
    output_full_dm = false;
    output_no_dm = false;

    // Parsing input:
    Timer &timer_parseInput = Timers::create( "Parsing parameters", true, false );
    Log::wrapInBar( "Conversion of input variables", Log::BAR_SIZE_FULL, Log::LEVEL_2, Log::BAR_0 );
    Log::L2( "\n" );
    Log::L2( "Parsing input variables... " );
    timer_parseInput.start();
    if ( !parseInput( arguments ) ) {
        Log::L2( "Parsing input variables failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    timer_parseInput.end();
    Log::L2( "successful. Elapsed time is {}ms\n", timer_parseInput.getWallTime( TIMER_MILLISECONDS ) );

    // Scaling inputs:
    if ( scale_parameters ) {
        Log::L2( "Rescaling parameters to {}... ", scale_value );
        scaleInputs( scale_value );
        Log::L2( "Done!\n" );
    }

    // Adjusting inputs:
    Timer &timer_adjustInput = Timers::create( "Adjusting parameters", true, false );
    Log::L2( "Adjusting input variables... " );
    timer_adjustInput.start();
    if ( !adjustInput() ) {
        Log::L2( "Adjusting input variables failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    Log::L2( "successful. Elapsed time is {}ms\n", timer_adjustInput.getWallTime( TIMER_MILLISECONDS ) );
    timer_adjustInput.end();
}

int Parameters::index_to_state( const char mode, const int state ) {
    if ( mode == 'E' )
        return state % 4;
    if ( mode == 'H' )
        return (int)std::floor( state / ( 4 * ( p_max_photon_number + 1 ) ) );
    if ( mode == 'V' )
        return ( (int)std::floor( state / 4 ) ) % ( p_max_photon_number + 1 );
    return -1;
}

bool Parameters::parseInput( const std::vector<std::string> &arguments ) {
    Parse_Parameters params;
    // Look for --time, if not found, standard values are used (t0 = 0, t1 = 1ns, deltaT = auto)
    t_start = get_parameter<double>( "--time", "tstart" );
    t_end = get_parameter<double>( "--time", "tend" );
    t_step = get_parameter<double>( "--time", "tstep" );

    // Look for --system, if not found, standard system is used
    p_omega_atomic_G_H = get_parameter<double>( "--system", "we" );
    p_omega_cavity_V = get_parameter<double>( "--system", "wcV" );
    p_omega_cavity_H = get_parameter<double>( "--system", "wcH" );
    p_omega_coupling = get_parameter<double>( "--system", "coupling" );
    p_omega_cavity_loss = get_parameter<double>( "--system", "kappa" );
    p_omega_pure_dephasing = get_parameter<double>( "--system", "gammapure" );
    p_omega_decay = get_parameter<double>( "--system", "gamma" );
    p_deltaE = get_parameter<double>( "--system", "deltaE" );
    p_biexciton_bindingenergy = get_parameter<double>( "--system", "excitonBindEnergy" );

    // Look for --chirp, if not found, standard system is used (no chirp, everything zero)
    chirp_t = get_parameter_vector<Parameter>( "--chirp", "chirpT" );
    chirp_y = get_parameter_vector<Parameter>( "--chirp", "chirpY" );
    chirp_ddt = get_parameter_vector<Parameter>( "--chirp", "chirpDDT" );
    chirp_type = get_parameter( "--chirp", "chirpType" );
    // Adjust length if they are not passed
    map_vector_to_standard( chirp_y, chirp_t, 0.0 );
    map_vector_to_standard( chirp_y, chirp_ddt, 0.0 );

    // Look for --pulse, if not found, standard system is used (no pulse, everything zero)
    pulse_center = get_parameter_vector<Parameter>( "--pulse", "pulseCenter" );
    pulse_amp = get_parameter_vector<Parameter>( "--pulse", "pulseAmp" );
    pulse_omega = get_parameter_vector<Parameter>( "--pulse", "pulseFreq" );
    pulse_sigma = get_parameter_vector<Parameter>( "--pulse", "pulseSigma" );
    pulse_type = get_parameter_vector( "--pulse", "pulseType" );
    pulse_pol = get_parameter_vector( "--pulse", "pulsePol" );
    pulse_omega_chirp = get_parameter_vector<Parameter>( "--pulse", "pulseChirp" );
    // Adjust length if they are not passed
    map_vector_to_standard( pulse_amp, pulse_center, pulse_center.at( 0 ) );
    map_vector_to_standard( pulse_amp, pulse_omega, pulse_omega.at( 0 ) );
    map_vector_to_standard( pulse_amp, pulse_sigma, pulse_sigma.at( 0 ) );
    map_vector_to_standard( pulse_amp, pulse_type, pulse_type.at( 0 ) );
    map_vector_to_standard( pulse_amp, pulse_pol, pulse_pol.at( 0 ) );
    map_vector_to_standard( pulse_amp, pulse_omega_chirp, pulse_omega_chirp.at( 0 ) );

    // Look for --dimensions, if not found, standard system is used (maxphotons = 0, starting state = |g,0>)
    p_max_photon_number = get_parameter<int>( "--maxPhotons" );
    p_initial_state_electronic = get_parameter( "--initState", "initElectronicState" ).front();
    p_initial_state_photon_h = get_parameter<int>( "--initState", "initHorizontalPhotons" );
    p_initial_state_photon_v = get_parameter<int>( "--initState", "initVerticalPhotons" );
    p_initial_state = instr( "ghvbc", p_initial_state_electronic ) + 4 * ( p_max_photon_number + 1 ) * p_initial_state_photon_h + 4 * p_initial_state_photon_v;

    // Look for --spectrum, if not found, no spectrum is evaluated
    iterations_tau_resolution = get_parameter<int>( "--spectrum", "specTauRes" );
    spectrum_frequency_center = get_parameter<double>( "--spectrum", "specCenter" );
    spectrum_frequency_range = get_parameter<double>( "--spectrum", "specRange" );
    iterations_w_resolution = get_parameter<int>( "--spectrum", "specWRes" );
    numerics_calculate_spectrum_H = get_parameter_passed( "-spectrum" ) || get_parameter_passed( "-spectrumH" );
    numerics_calculate_spectrum_V = get_parameter_passed( "-spectrum" ) || get_parameter_passed( "-spectrumV" );

    // Look for (-RK4), -RK5, (-RK4T), (-RK4Tau), -RK5T, -RK5Tau
    numerics_calculate_g2 = get_parameter_passed( "-g2" ) || get_parameter_passed( "-g2s" ) || get_parameter_passed( "-g2H" ) || get_parameter_passed( "-g2V" ) || get_parameter_passed( "-g2C" );
    numerics_calculate_g2_H = get_parameter_passed( "-g2H" ) || get_parameter_passed( "-g2" );
    numerics_calculate_g2_V = get_parameter_passed( "-g2V" ) || get_parameter_passed( "-g2" );
    numerics_calculate_g2_C = get_parameter_passed( "-g2C" ) || get_parameter_passed( "-g2" );
    numerics_order_t = ( get_parameter_passed( "-RK5" ) || get_parameter_passed( "-RK5T" ) ? 5 : 4 );
    numerics_order_tau = ( get_parameter_passed( "-RK5" ) || get_parameter_passed( "-RK5Tau" ) ? 5 : 4 );
    numerics_order_highest = numerics_order_t;
    if ( numerics_order_tau > numerics_order_highest )
        numerics_order_highest = numerics_order_tau;
    numerics_use_interactionpicture = get_parameter_passed( "-noInteractionpic" ) ? 0 : 1;
    numerics_use_rwa = get_parameter_passed( "-noRWA" ) ? 0 : 1;
    numerics_maximum_threads = get_parameter<int>( "--Threads" );
    if ( numerics_maximum_threads == -1 )
        numerics_maximum_threads = omp_get_max_threads();
    output_handlerstrings = get_parameter_passed( "-noHandler" ) ? 0 : 1;
    output_operators = get_parameter_passed( "-outputOp" ) ? 2 : ( get_parameter_passed( "-outputHamiltons" ) ? 1 : ( get_parameter_passed( "-outputOpStop" ) ? 3 : 0 ) );
    numerics_order_timetrafo = get_parameter_passed( "-timeTrafoMatrixExponential" ) ? TIMETRANSFORMATION_MATRIXEXPONENTIAL : TIMETRANSFORMATION_ANALYTICAL;
    startCoherent = get_parameter_passed( "-startCoherent" ) || ( p_initial_state_electronic.front() == 'c' );
    output_full_dm = get_parameter_passed( "-fullDM" );
    output_no_dm = get_parameter_passed( "-noDM" );
    scale_parameters = get_parameter_passed( "-scale" );
    numerics_use_saved_coefficients = !get_parameter_passed( "-disableMatrixCaching" );
    numerics_use_saved_hamiltons = !get_parameter_passed( "-disableHamiltonCaching" );
    numerics_phonons_maximum_threads = ( !numerics_use_saved_coefficients || !get_parameter_passed( "-disableMainProgramThreading" ) ) ? numerics_maximum_threads : 1;
    numerics_output_raman_population = get_parameter_passed( "-raman" );
    numerics_use_simplified_g2 = get_parameter_passed( "-g2s" );
    logfilecounter = get_parameter<int>( "--lfc" );
    numerics_calculate_timeresolution_indistinguishability = get_parameter_passed( "-timedepInd" );
    numerics_output_electronic_emission = get_parameter_passed( "-oElec" );
    numerics_stretch_correlation_grid = false; //FIXME: Doesnt work right now
    numerics_interpolate_outputs = get_parameter_passed( "-interpolate" );

    // Phonon Parameters
    p_phonon_alpha = get_parameter<double>( "--phonons", "phononalpha" );
    p_phonon_wcutoff = get_parameter<double>( "--phonons", "phononwcutoff" );
    p_phonon_tcutoff = get_parameter<double>( "--phonons", "phonontcutoff" );
    p_phonon_T = get_parameter<double>( "--phonons", "temperature" );
    numerics_phonon_approximation_order = get_parameter<int>( "--phonons", "phononorder" );
    numerics_phonon_approximation_markov1 = get_parameter_passed( "-noMarkov" ) ? 0 : 1; // First Markov
    output_coefficients = get_parameter_passed( "-phononcoeffs" ) ? 1 : 0;
    p_phonon_adjust = !get_parameter_passed( "-noPhononAdjust" );
    p_phonon_pure_dephasing = convertParam<double>( "1mueV" );

    kb = 1.3806488E-23;   // J/K, scaling needs to be for energy
    hbar = 1.0545718E-34; // J/s, scaling will be 1

    subfolder = arguments.back();
    return true;
}

bool Parameters::scaleInputs( const double scaling ) {
    // Adjust normal parameters: time is multiplid by scaling, frequency divided
    t_start.setScale( scaling, Parameter::SCALE_TIME );
    t_end.setScale( scaling, Parameter::SCALE_TIME );
    t_step.setScale( scaling, Parameter::SCALE_TIME );
    p_omega_atomic_G_H.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_atomic_G_V.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_atomic_H_B.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_atomic_V_B.setScale( scaling, Parameter::SCALE_ENERGY );
    p_deltaE.setScale( scaling, Parameter::SCALE_ENERGY );
    p_biexciton_bindingenergy.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_cavity_H.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_cavity_V.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_coupling.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_cavity_loss.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_pure_dephasing.setScale( scaling, Parameter::SCALE_ENERGY );
    p_omega_decay.setScale( scaling, Parameter::SCALE_ENERGY );
    // Adjust chirp and pulse
    for ( int i = 0; i < (int)chirp_t.size(); i++ ) {
        chirp_t.at( i ).setScale( scaling, Parameter::SCALE_TIME );
        chirp_y.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
        chirp_ddt.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
    }
    for ( int i = 0; i < (int)pulse_center.size(); i++ ) {
        pulse_center.at( i ).setScale( scaling, Parameter::SCALE_TIME );
        if ( pulse_type.at( i ).compare( "gauss_pi" ) != 0 )
            pulse_amp.at( i ) = scaleVariable( pulse_amp.at( i ), 1.0 / scaling );
        pulse_omega.at( i ).setScale( scaling, Parameter::SCALE_ENERGY );
        pulse_sigma.at( i ).setScale( scaling, Parameter::SCALE_TIME );
        pulse_omega_chirp.at( i ).setScale( scaling * scaling, Parameter::SCALE_ENERGY );
    }
    // Adjusting spectrum
    spectrum_frequency_center.setScale( scaling, Parameter::SCALE_ENERGY );
    spectrum_frequency_range.setScale( scaling, Parameter::SCALE_ENERGY );
    // Phonons
    p_phonon_wcutoff.setScale( scaling, Parameter::SCALE_ENERGY );
    p_phonon_tcutoff.setScale( scaling, Parameter::SCALE_TIME );
    p_phonon_alpha.setScale( scaling * scaling, Parameter::SCALE_ENERGY );
    p_phonon_pure_dephasing.setScale( scaling, Parameter::SCALE_ENERGY );
    kb.setScale( scaling, Parameter::SCALE_ENERGY );
    return true;
}

double Parameters::getIdealTimestep() {
    if ( numerics_use_rwa )
        return 2. / 8. * M_PI / std::max( std::max( init_rabifrequenz, max_rabifrequenz ), p_omega_coupling + p_omega_cavity_loss + p_omega_decay + p_omega_pure_dephasing + ( vec_max( pulse_amp ) > 0 ? std::abs( p_omega_atomic_G_H - vec_max( pulse_omega ) ) : 0 ) );
    if ( !numerics_use_rwa )
        return 2. / 8. * M_PI / std::max( std::max( init_rabifrequenz, max_rabifrequenz ), p_omega_atomic_G_H + p_omega_cavity_H );
    return 1E-13;
}

bool Parameters::adjustInput() {
    // Calculate/Recalculate some parameters:
    // Adjust pulse area if pulse_type is "gauss_pi"
    for ( int i = 0; i < (int)pulse_amp.size(); i++ )
        if ( pulse_type.at( i ).compare( "gauss_pi" ) == 0 ) {
            //if ( pulse_amp.at( i ) < 500 )
            pulse_amp.at( i ) = pulse_amp.at( i ) * M_PI / ( std::sqrt( 2.0 * M_PI * pulse_sigma.at( i ) * std::sqrt( std::pow( pulse_omega_chirp.at( i ) / pulse_sigma.at( i ), 2.0 ) + std::pow( pulse_sigma.at( i ), 2.0 ) ) ) ) / 2.0; //https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.241306
            pulse_type.at( i ) = "gauss";
        }

    // Calculating remaining atomic frequencies depending on delta E and biexciton binding energy.
    p_omega_atomic_G_V = p_omega_atomic_G_H + p_deltaE / 2.0;
    p_omega_atomic_B = 2.0 * p_omega_atomic_G_H - p_biexciton_bindingenergy;
    p_omega_atomic_G_H = p_omega_atomic_G_H - p_deltaE / 2.0;
    p_omega_atomic_H_B = p_omega_atomic_B - p_omega_atomic_G_H;
    p_omega_atomic_V_B = p_omega_atomic_B - p_omega_atomic_G_V;

    // Calculate Rabi frequencies:
    init_rabifrequenz_G_H = rabiFrequency( p_omega_atomic_G_H - p_omega_cavity_H, p_omega_coupling, index_to_state( 'H', p_initial_state ) );
    init_rabifrequenz_G_V = rabiFrequency( p_omega_atomic_G_V - p_omega_cavity_V, p_omega_coupling, index_to_state( 'V', p_initial_state ) );
    init_rabifrequenz_H_B = rabiFrequency( p_omega_atomic_H_B - p_omega_cavity_H, p_omega_coupling, index_to_state( 'H', p_initial_state ) );
    init_rabifrequenz_V_B = rabiFrequency( p_omega_atomic_V_B - p_omega_cavity_V, p_omega_coupling, index_to_state( 'V', p_initial_state ) );
    max_rabifrequenz_G_H = rabiFrequency( p_omega_atomic_G_H - p_omega_cavity_H, p_omega_coupling, index_to_state( 'H', p_initial_state + 1 ) );
    max_rabifrequenz_G_V = rabiFrequency( p_omega_atomic_G_V - p_omega_cavity_V, p_omega_coupling, index_to_state( 'V', p_initial_state + 1 ) );
    max_rabifrequenz_H_B = rabiFrequency( p_omega_atomic_H_B - p_omega_cavity_H, p_omega_coupling, index_to_state( 'H', p_initial_state + 1 ) );
    max_rabifrequenz_V_B = rabiFrequency( p_omega_atomic_V_B - p_omega_cavity_V, p_omega_coupling, index_to_state( 'V', p_initial_state + 1 ) );
    init_rabifrequenz = vec_max<double>( {init_rabifrequenz_G_H, init_rabifrequenz_G_V, init_rabifrequenz_H_B, init_rabifrequenz_V_B} );
    max_rabifrequenz = vec_max<double>( {max_rabifrequenz_G_H, max_rabifrequenz_G_V, max_rabifrequenz_H_B, max_rabifrequenz_V_B} );

    // Calculate minimum step necessary to resolve Rabi-oscillation if step=-1
    if ( t_step == -1 ) {
        t_step = std::min( scaleVariable( 1E-13, scale_value ), getIdealTimestep() );
        t_step = std::max( std::numeric_limits<double>::epsilon(), t_step );
    }

    // Calculate the maximum dimensions for operator matrices (max states)
    maxStates = 4 * ( p_max_photon_number + 1 ) * ( p_max_photon_number + 1 ); // 4 Electronic states x N+1 Photonic states * 2 for H,V modes

    // Calculate stuff for RK
    iterations_t_max = (int)std::ceil( ( t_end - t_start ) / t_step );
    iterations_t_skip = std::max( 1.0, std::ceil( iterations_t_max / iterations_tau_resolution ) );
    iterations_wtau_skip = iterations_t_skip; //1;//iterations_t_skip;

    // Mandatory: rescale chirp ddt into chirp/ps
    if ( chirp_type.compare( "sine" ) != 0 )
        for ( long unsigned i = 0; i < chirp_ddt.size(); i++ )
            chirp_ddt.at( i ) = chirp_ddt.at( i ) * 1E12;
    // Calculate values for chirp:
    chirp_total = vec_max( chirp_y );

    // Initial and maximum system detuning (not taking into account custom chirps)
    init_detuning_G_H = ( p_omega_atomic_G_H - p_omega_cavity_H );
    init_detuning_G_V = ( p_omega_atomic_G_V - p_omega_cavity_V );
    init_detuning_H_B = ( p_omega_atomic_H_B - p_omega_cavity_H );
    init_detuning_V_B = ( p_omega_atomic_V_B - p_omega_cavity_V );
    max_detuning_G_H = ( (double)( init_detuning_G_H + chirp_total ) > (double)init_detuning_G_H ) ? double( init_detuning_G_H + chirp_total ) : (double)init_detuning_G_H;
    max_detuning_G_V = ( (double)( init_detuning_G_V + chirp_total ) > (double)init_detuning_G_V ) ? double( init_detuning_G_V + chirp_total ) : (double)init_detuning_G_V;
    max_detuning_H_B = ( (double)( init_detuning_H_B + chirp_total ) > (double)init_detuning_H_B ) ? double( init_detuning_H_B + chirp_total ) : (double)init_detuning_H_B;
    max_detuning_V_B = ( (double)( init_detuning_V_B + chirp_total ) > (double)init_detuning_V_B ) ? double( init_detuning_V_B + chirp_total ) : (double)init_detuning_V_B;
    init_detuning = vec_max<double>( {init_detuning_G_H, init_detuning_G_V, init_detuning_H_B, init_detuning_V_B} );
    max_detuning = vec_max<double>( {max_detuning_G_H, max_detuning_G_V, max_detuning_H_B, max_detuning_V_B} );

    // Adjust/calculate frequency range for spectrum
    if ( spectrum_frequency_center == -1 )
        spectrum_frequency_center = p_omega_cavity_H;
    if ( spectrum_frequency_range == -1 )
        spectrum_frequency_range = ( std::abs( max_detuning * 1.5 ) + ( p_omega_coupling + p_omega_cavity_loss + p_omega_decay + p_omega_pure_dephasing ) * 3.0 );

    // No phonon adjust if pathintegral is chosen
    if ( numerics_phonon_approximation_order == 5 ) {
        p_phonon_adjust = false;
    }
    // Calculate phonon stuff
    p_phonon_b = 1.0;
    if ( p_phonon_T >= 0 ) {
        double integral = 0;
        double stepsize = 0.01 * p_phonon_wcutoff;
        for ( double w = stepsize; w < 10 * p_phonon_wcutoff; w += stepsize ) {
            integral += stepsize * ( p_phonon_alpha * w * std::exp( -w * w / 2.0 / p_phonon_wcutoff / p_phonon_wcutoff ) / std::tanh( hbar * w / 2.0 / kb / p_phonon_T ) );
        }
        p_phonon_b = std::exp( -0.5 * integral );
        if ( p_phonon_adjust ) {
            p_omega_pure_dephasing = p_phonon_pure_dephasing * p_phonon_T;
            p_omega_decay = p_omega_decay * p_phonon_b * p_phonon_b;
        }
    }
    numerics_calculate_spectrum = numerics_calculate_spectrum_H || numerics_calculate_spectrum_V;
    numerics_saved_coefficients_max_size = (int)( ( t_end - t_start ) / t_step * 2.0 * ( p_phonon_tcutoff / t_step ) ) + 10;
    trace.reserve( iterations_t_max + 5 );

    numerics_saved_coefficients_cutoff = ( numerics_calculate_spectrum || numerics_calculate_g2 ) ? 0 : ( p_phonon_tcutoff / t_step ) * 5;

    return true;
}

void Parameters::log( const std::vector<std::string> &info ) {
    Log::wrapInBar( "System Parameters" );
    Log::L1( "Version: {} ({})\n\n", GLOBAL_PROGRAM_VERSION, GLOBAL_PROGRAM_LASTCHANGE );

    Log::wrapInBar( "Base Energies", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Biexciton binding energy: {:.8e} Hz - {:.8} meV\n", p_biexciton_bindingenergy, p_biexciton_bindingenergy.getSI( Parameter::UNIT_ENERGY_MEV ) );
    Log::L1( "Exciton fine structure splitting: {:.8e} Hz - {:.8f} mueV\n", p_deltaE, p_deltaE.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Transition Energy H - G: {:.8e} Hz - {:.8} eV - {:.8f} nm\n", p_omega_atomic_G_H, p_omega_atomic_G_H.getSI( Parameter::UNIT_ENERGY_EV ), p_omega_atomic_G_H.getSI( Parameter::UNIT_WAVELENGTH_NM ) );
    Log::L1( "Transition Energy V - G: {:.8e} Hz - {:.8} eV - {:.8f} nm\n", p_omega_atomic_G_V, p_omega_atomic_G_V.getSI( Parameter::UNIT_ENERGY_EV ), p_omega_atomic_G_V.getSI( Parameter::UNIT_WAVELENGTH_NM ) );
    Log::L1( "Transition Energy B - H: {:.8e} Hz - {:.8} eV - {:.8f} nm\n", p_omega_atomic_H_B, p_omega_atomic_H_B.getSI( Parameter::UNIT_ENERGY_EV ), p_omega_atomic_H_B.getSI( Parameter::UNIT_WAVELENGTH_NM ) );
    Log::L1( "Transition Energy B - V: {:.8e} Hz - {:.8} eV - {:.8f} nm\n", p_omega_atomic_V_B, p_omega_atomic_V_B.getSI( Parameter::UNIT_ENERGY_EV ), p_omega_atomic_V_B.getSI( Parameter::UNIT_WAVELENGTH_NM ) );
    Log::L1( "Biexciton Resonance: {:.8e} Hz - {:.8} eV - {:.8f} nm\n", p_omega_atomic_B, p_omega_atomic_B.getSI( Parameter::UNIT_ENERGY_EV ), p_omega_atomic_B.getSI( Parameter::UNIT_WAVELENGTH_NM ) );
    Log::L1( "Two Photon Resonance: {:.8e} Hz - {:.8} eV - {:.8f} nm\n", p_omega_atomic_B / 2.0, p_omega_atomic_B.getSI( Parameter::UNIT_ENERGY_EV ) / 2.0, p_omega_atomic_B.getSI( Parameter::UNIT_WAVELENGTH_NM ) / 2.0 );
    Log::L1( "Cavity Energy H: {:.8e} Hz - {:.8} eV - {:.8f} nm\n", p_omega_cavity_H, p_omega_cavity_H.getSI( Parameter::UNIT_ENERGY_EV ), p_omega_cavity_H.getSI( Parameter::UNIT_WAVELENGTH_NM ) );
    Log::L1( "Cavity Energy V: {:.8e} Hz - {:.8} eV - {:.8f} nm\n\n", p_omega_cavity_V, p_omega_cavity_V.getSI( Parameter::UNIT_ENERGY_EV ), p_omega_cavity_V.getSI( Parameter::UNIT_WAVELENGTH_NM ) );

    Log::wrapInBar( "Coupling Parameters", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Coupling strengh g: {:.8e} Hz - {:.8} mueV\n", p_omega_coupling, p_omega_coupling.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Photon loss rate k: {:.8e} Hz - {:.8} mueV - Q = {:.2f}\n", p_omega_cavity_loss, p_omega_cavity_loss.getSI( Parameter::UNIT_ENERGY_MUEV ), p_omega_cavity_H / p_omega_cavity_loss );
    Log::L1( "Atomic dephasing rate gamma_pure: {:.8e} Hz - {:.8} mueV\n", p_omega_pure_dephasing, p_omega_pure_dephasing.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "RAD rate gamma: {:.8e} Hz - {:.8} mueV\n\n", p_omega_decay, p_omega_decay.getSI( Parameter::UNIT_ENERGY_MUEV ) );

    Log::wrapInBar( "Initial System Parameters", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Initial state rho0 = |{0}><{0}| with maximum number of {1} photons\n", info.at( p_initial_state ), p_max_photon_number );
    Log::L1( "Initial detuning X_H - G - H = {:.8e} Hz - {:.8} mueV, Maximum detuning X_H - G - H = {:.8e} Hz - {:.8} mueV\n", init_detuning_G_H, init_detuning_G_H.getSI( Parameter::UNIT_ENERGY_MUEV ), max_detuning_G_H, max_detuning_G_H.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Initial detuning X_V - G - V = {:.8e} Hz - {:.8} mueV, Maximum detuning X_V - G - V = {:.8e} Hz - {:.8} mueV\n", init_detuning_G_V, init_detuning_G_V.getSI( Parameter::UNIT_ENERGY_MUEV ), max_detuning_G_V, max_detuning_G_V.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Initial detuning B - X_H - H = {:.8e} Hz - {:.8} mueV, Maximum detuning B - X_H - H = {:.8e} Hz - {:.8} mueV\n", init_detuning_H_B, init_detuning_H_B.getSI( Parameter::UNIT_ENERGY_MUEV ), max_detuning_H_B, max_detuning_H_B.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Initial detuning B - X_V - V = {:.8e} Hz - {:.8} mueV, Maximum detuning B - X_V - V = {:.8e} Hz - {:.8} mueV\n", init_detuning_V_B, init_detuning_V_B.getSI( Parameter::UNIT_ENERGY_MUEV ), max_detuning_V_B, max_detuning_V_B.getSI( Parameter::UNIT_ENERGY_MUEV ) );

    Log::L1( "Initial Rabi Frequency X_H - G = {:.8e} Hz - {:.8} mueV, Maximum Rabi Frequency X_H - B = {:.8e} Hz - {:.8} mueV\n", init_rabifrequenz_G_H, init_rabifrequenz_G_H.getSI( Parameter::UNIT_ENERGY_MUEV ), max_rabifrequenz_G_H, max_rabifrequenz_G_H.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Initial Rabi Frequency X_V - G = {:.8e} Hz - {:.8} mueV, Maximum Rabi Frequency X_V - G = {:.8e} Hz - {:.8} mueV\n", init_rabifrequenz_G_V, init_rabifrequenz_G_V.getSI( Parameter::UNIT_ENERGY_MUEV ), max_rabifrequenz_G_V, max_rabifrequenz_G_V.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Initial Rabi Frequency B - X_H = {:.8e} Hz - {:.8} mueV, Maximum Rabi Frequency B - X_H = {:.8e} Hz - {:.8} mueV\n", init_rabifrequenz_H_B, init_rabifrequenz_H_B.getSI( Parameter::UNIT_ENERGY_MUEV ), max_rabifrequenz_H_B, max_rabifrequenz_H_B.getSI( Parameter::UNIT_ENERGY_MUEV ) );
    Log::L1( "Initial Rabi Frequency B - X_V = {:.8e} Hz - {:.8} mueV, Maximum Rabi Frequency B - X_V = {:.8e} Hz - {:.8} mueV\n\n", init_rabifrequenz_V_B, init_rabifrequenz_V_B.getSI( Parameter::UNIT_ENERGY_MUEV ), max_rabifrequenz_V_B, max_rabifrequenz_V_B.getSI( Parameter::UNIT_ENERGY_MUEV ) );

    Log::wrapInBar( "Pulse", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( pulse_amp.size() > 0 && pulse_center.at( 0 ) != -1 && pulse_amp.at( 0 ) != 0 ) {
        for ( int i = 0; i < (int)pulse_amp.size(); i++ ) {
            Log::L1( "Exiting system at t_{0} = {1:.8e}\nAmplitude_{0} {2:.10e} ({3:.8} meV, {4:.2f} pi)\nFrequency_{0} {5:.10e} ({6:.8}eV)\nFWHM_{0} {7:.8e}\nPulseChirp_{0}: {8:.8e}\n", i, pulse_center.at( i ), pulse_amp.at( i ), pulse_amp.at( i ).getSI( Parameter::UNIT_ENERGY_MEV ), pulse_amp.at( i ) / ( M_PI / ( std::sqrt( 2.0 * M_PI ) * pulse_sigma.at( i ) ) / 2.0 ), pulse_omega.at( i ), pulse_omega.at( i ).getSI( Parameter::UNIT_ENERGY_MEV ), pulse_sigma.at( i ) * ( 2 * std::sqrt( 2 * std::log( 2 ) ) ), pulse_omega_chirp.at( i ) );
            Log::L1( "Used pulse_type_{0} - {1}\nUsed pulse_pol_{0} on mode {2}\n", i, pulse_type.at( i ), pulse_pol.at( i ) );
        }
        Log::L1( "\n\n" );
    } else
        Log::L1( "Not using pulse to exite system\n\n" );
    Log::wrapInBar( "Chirp", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( chirp_total != 0 ) {
        double total = 0;
        int current = 0;
        for ( int i = 0; i < (int)chirp_t.size() - 1; i++ ) {
            if ( chirp_y.at( i + 1 ) - chirp_y.at( i ) != 0.0 ) {
                //Log::inBar( "Chirp " + toStr( current++ ) );
                Log::L1( "Chirp {0:} between t0 = {1:.8e} ps\nChirp {0:} t1 = {2:.8e} ps\nTotal Chirp {0:}: {3:.8} mueV\n-> average rate Chirp {0:}: {4:.8} mueV/ps\n", current, chirp_t.at( i ), chirp_t.at( i + 1 ), chirp_y.at( i + 1 ).getSI( Parameter::UNIT_ENERGY_MUEV ) - chirp_y.at( i ).getSI( Parameter::UNIT_ENERGY_MUEV ), ( chirp_y.at( i + 1 ).getSI( Parameter::UNIT_ENERGY_MUEV ) - chirp_y.at( i ).getSI( Parameter::UNIT_ENERGY_MUEV ) ) / ( chirp_t.at( i + 1 ).getSI( Parameter::UNIT_TIME_PS ) - chirp_t.at( i ).getSI( Parameter::UNIT_TIME_PS ) ) );
                total += chirp_y.at( i + 1 ).getSI( Parameter::UNIT_ENERGY_MUEV ) - chirp_y.at( i ).getSI( Parameter::UNIT_ENERGY_MUEV );
                current++;
            }
        }
        if ( chirp_type.compare( "none" ) != 0 )
            Log::L1( "\nChirp of type '" + chirp_type + "' is used!\n" );
        Log::L1( "Total Chirp = {:.8} mueV\n\n", total );
    } else
        Log::L1( "Not using chirp\n\n" );

    Log::wrapInBar( "Photon Statistics" );
    Log::L1( "\n" );
    Log::wrapInBar( "G1 Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Anticipated tau-grid resolution is {}x{} resulting in {} skips per timestep\n", iterations_tau_resolution, iterations_tau_resolution, iterations_t_skip );
    if ( numerics_calculate_spectrum_H || numerics_calculate_spectrum_V ) {
        Log::L1( "Calcluating spectrum for modes {}\n", ( numerics_calculate_spectrum_H && numerics_calculate_spectrum_V ? "cavities H and V" : fmt::format( "cavity {}", numerics_calculate_spectrum_H ? "H" : "V" ) ) );
        Log::L1( "Center Frequency: {:.8e} Hz - {:.8} eV\n", spectrum_frequency_center, spectrum_frequency_center.getSI( Parameter::UNIT_ENERGY_EV ) );
        Log::L1( "Frequency Range: +/- {:.8e} Hz - +/- {:.8} mueV\n", spectrum_frequency_range, spectrum_frequency_range.getSI( Parameter::UNIT_ENERGY_MEV ) );
        Log::L1( "Maximum w-vector resolution is {}\n", iterations_w_resolution );
    } else {
        Log::L1( "Not calculating spectrum\n" );
    }
    Log::L1( "\n" );
    Log::wrapInBar( "G2 Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Calculate advanced statistics? - {}\n", numerics_calculate_g2 ? fmt::format( "YES{}", ( numerics_use_simplified_g2 ? " (simplified)" : "" ) ) : "NO" );
    Log::L1( "Calculated full indistinguishability? - {}\n", numerics_calculate_timeresolution_indistinguishability ? "YES" : "NO" );
    Log::L1( "\n" );

    Log::wrapInBar( "Phonons", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    if ( p_phonon_T >= 0 ) {
        std::vector<std::string> approximations = {"Transformation integral via d/dt chi = -i/hbar*[H,chi] + d*chi/dt onto interaction picture chi(t-tau)", "Transformation Matrix U(t,tau)=exp(-i/hbar*H_DQ_L(t)*tau) onto interaction picture chi(t-tau)", "No Transformation, only interaction picture chi(t-tau)", "Analytical Lindblad formalism", "Mixed", "Path Integral"};
        Log::L1( "Temperature = {} k\nCutoff energy = {} meV\nCutoff Time = {} ps\nAlpha = {}\n<B> = {}\nFirst Markov approximation used? (rho(t) = rho(t-tau)) - {}\nTransformation approximation used: {} - {}\n\n", p_phonon_T, p_phonon_wcutoff.getSI( Parameter::UNIT_ENERGY_MEV ), p_phonon_tcutoff * 1E12, p_phonon_alpha, p_phonon_b, ( numerics_phonon_approximation_markov1 == 1 ? "Yes" : "No" ), numerics_phonon_approximation_order, approximations.at( numerics_phonon_approximation_order ) );
    } else {
        Log::L1( "Not using phonons\n\n" );
    }

    Log::wrapInBar( "Numerical Parameters" );
    Log::L1( "\n" );
    Log::wrapInBar( "Time", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Timeborder start: {:.8e} s - {:.2f} ps\n", t_start, t_start * 1E12 );
    Log::L1( "Timeborder end: {:.8e} s - {:.2f} ps\n", t_end, t_end * 1E12 );
    Log::L1( "Timeborder delta: {:.8e} s - {:.2f} fs \n", t_step, t_step * 1E15 );
    Log::L1( "Ideal time delta for this calculation: {:.8e} s - {:.2f} fs, minimum possible: {:.8e} s - {:.2f} fs\n", getIdealTimestep(), getIdealTimestep() * 1E15, std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::epsilon() * 1E15 );
    int works = 1;
    if ( ( init_rabifrequenz != 0.0 ) && ( 3. * t_step > 2. * M_PI / init_rabifrequenz ) )
        works = 0;
    else if ( max_rabifrequenz != 0.0 && 3. * t_step > 2. * M_PI / max_rabifrequenz )
        works = 0;
    if ( !works ) {
        fmt::print( "{} WARNING: Step may be too small to resolve predicted oscillation: dT needed vs dT: {:.10e} < {:.10e}\n", PREFIX_WARNING, 2. / 3. * M_PI / std::max( init_rabifrequenz, max_rabifrequenz ), t_step );
        Log::L1( "WARNING: Step may be too small to resolve predicted oscillation: \n- delta T needed: {:.10e} \n- delta T used: {:.10e}\n", 2. / 3. * M_PI / std::max( init_rabifrequenz, max_rabifrequenz ), t_step );
    }
    Log::L1( "Time iterations (main loop) = {}\n\n", iterations_t_max );

    Log::wrapInBar( "Settings", Log::BAR_SIZE_HALF, Log::LEVEL_1, Log::BAR_1 );
    Log::L1( "Order of Runge-Kutta used: Time: RK{}, Spectrum: RK{}\n", numerics_order_t, numerics_order_tau );
    Log::L1( "Use rotating wave approximation (RWA)? - {}\n", ( ( numerics_use_rwa == 1 ) ? "YES" : "NO" ) );
    Log::L1( "Use interaction picture for calculations? - {}\n", ( ( numerics_use_interactionpicture == 1 ) ? "YES" : "NO" ) );
    Log::L1( "Time Transformation used? - {}\n", ( ( numerics_order_timetrafo == TIMETRANSFORMATION_ANALYTICAL ) ? "Analytic" : "Matrix Exponential" ) );
    Log::L1( "Threads used for primary calculations - {}\nThreads used for Secondary calculations - {}\n", numerics_phonons_maximum_threads, numerics_maximum_threads );
    Log::L1( "Used scaling for parameters? - {}\n", ( scale_parameters ? std::to_string( scale_value ) : "no" ) );
    if ( p_phonon_T )
        Log::L1( "Cache Phonon Coefficient Matrices? - {}\n", ( numerics_use_saved_coefficients ? fmt::format( "Yes (maximum {} matrices saved)", ( numerics_saved_coefficients_cutoff > 0 ) ? numerics_saved_coefficients_cutoff : numerics_saved_coefficients_max_size ) : "No" ) );
    if ( numerics_interpolate_outputs )
        Log::L1( "WARNING: Temporal outputs are interpolated!\n" );
    Log::L1( "\n\n" );
    if ( logfilecounter >= 0 )
        Log::L1( "Logfile ident number: {}\n\n", logfilecounter );
    Log::wrapInBar( "Program Log:", Log::BAR_SIZE_FULL, Log::LEVEL_2 );
    Log::L1( "\n" );
    Log::L2( "OutputHandlerStrings: {}\n", output_handlerstrings );
}

void Parameters::help() {
    fmt::print( "--help, -help, -h\tThis screen\n" );
    fmt::print( "--time [start] [end] [step]\n\t--tstart [start]\n\t--tend [end]\n\t--tstep [step]\n" );
    fmt::print( "--system [p_omega_atomic] [p_omega_cavity_H] [p_omega_cavity_V] [p_omega_coupling] [kappa] [gammapure] [gamma] [deltaE] [bexcitonbinding] else standard values (1.366eV,1.366eV,1.366eV,66mueV,66mueV,3mueV,1mueV,0,3meV) are used\n\t--we [1.366eV] First (H-polarized) exciton energy\n\t--wcH [1.366eV] Cavity Energy H\n\t--wcV [1.366eV] Cavity Energy V\n\t--coupling [66mueV] QD-Lightfield coupling. Same for H and V.\n\t--kappa [66mueV] Cavity Decay\n\t--gamma [1mueV] Radiative Decay\n\t--gammapure [3mueV] Pure Dephasing\n\t--deltaE [0eV] H and V energy difference\n\t--excitonBindEnergy [3meV] Biexciton binding energy\n" );
    fmt::print( "--chirp ['[Array Time]'] ['[Array Y]'] ['[Array d/dt]'] [type]\n\t--chirpT ['[Array Time]']\n\t--chirpY ['[Array Y]']\n\t--chirpDDT ['[Array d/dt]']\n\t--chirpType [type] where type = monotone, hermite, linear, spline\n" );
    fmt::print( "--pulse [Center] [Amplitude] [Frequency] [Sigma] [Type]\n\t-pulse for standard pulse\n\t--pulseCenter [Center]\n\t--pulseAmp [Amplitude]\n\t--pulseFreq [Frequency]\n\t--pulseSigma [Sigma]\n\t--pulseType [Type] where Type = cw, gauss, gauss_pi\n" );
    fmt::print( "--dimensions [maximum Photons] [Initial state]\n\t--maxPhotons [maximum Photons]\n\t--initState [Initial state], has to be smaller than (2*n+1)\n" );
    fmt::print( "--spectrum [Tau Resolution] [Center] [Range] [Omega Resolution] enables spectrum\n\t-spectrum enables spectrum centered at cavity\n\t--specTauRes [Grid resolution (int)] standard is 1000\n\t--specCenter [Center]\n\t--specRange [Range]\n\t--specWRes [w vector resolution (int)] standard is 1000\n" );
    fmt::print( "--phonons [Alpha] [W Cutoff] [t Cutoff ] [Temperature] Enables phonons with custom settings\n\t--temperature [temperature] Enables phonons with standard values and set temperature\n\t-phonons Enables phonons with standard values at T=3k\n\t--phononorder Sets the order of approximation to use\n\t-noMarkov disables first markov approximation\n\t-phononcoeffs Enables output of additional phonon coefficients\n" );
    fmt::print( "-g2 enables calculation of full G2 statistics\n\t-g2s enables calculation of simple G2 statistics\n-RK5 enables Runge Kutta of order 5 for T and Tau direction\n\t-RK5T enables Runge Kutta of order 5 for T direction\n\t-RK5Tau enables Runge Kutta of order 5 for Tau direction\n" );
    fmt::print( "-noInteractionpic disables Interaction picture - enabled by default\n-noRWA disables rotating wave approximation - enabled by default\n-timeTrafoMatrixExponential enables Time Transformation via Matrix exponential - disabled by default\n-startCoherent enables starting with a coherent state. Starting state is then ground state with alpha = initState\n-fullDM enables full output of densitymatrix including all offdiagonal terms.\n" );
    fmt::print( "--Threads [number] number of threads to use for both AKF and Spectrum integral calculation\n" );
    fmt::print( "Additional commands:\n\t-advLog Enables advanced logging\n\t-noHandler disables handler strings and enables loadbar output (for console)\n\t-output_operators, -outputHamiltons, -outputOperatorsStop Enables output of matrices (requires -advLog)\n\t-noRaman disables slow Raman calculations" );
}