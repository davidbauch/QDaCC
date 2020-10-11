#include "solver/solver_ode.h"

// Description: Calculates the Eberly-Wódkiewicz spectrum. Uses precalculated values from akf_mat
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false
bool ODESolver::calculate_spectrum( System &s, const Sparse &op_creator, const Sparse &op_annihilator, std::string fileOutputName, const int cache_index ) {
    // Send system command to change to single core mainprogram, because this memberfunction is already using multithreading
    s.command( Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
    // Cache matrix. Saves complex double entries of G1(t,tau)
    Dense akf_mat = Dense::Zero( dim, dim );
    // Calculate G1(t,tau) with given operator matrices
    calculate_g1( s, op_creator, op_annihilator, akf_mat, "spectrum" );
    // Create Timer and Progressbar for the spectrum loop
    Timer &timer = Timers::create( "Spectrum-Loop" );
    int totalIterations = getIterationNumberSpectrum( s );
    ProgressBar progressbar = ProgressBar( totalIterations );
    timer.start();
    Log::L2( "Calculating spectrum... Calculating frequencies... " );
    //Calculate frequencies:
    std::vector<double> spectrum_frequency_w;
    std::vector<Scalar> out;
    for ( int w = 0; w < s.parameters.iterations_w_resolution; w++ ) {
        spectrum_frequency_w.push_back( s.parameters.spectrum_frequency_center - ( s.parameters.spectrum_frequency_range ) + w / ( (double)s.parameters.iterations_w_resolution ) * ( 2. * ( s.parameters.spectrum_frequency_range ) ) );
        out.push_back( 0 );
    }
    Log::L2( "Done, calculating fourier transform via direct integral... " );
    // Calculate main fourier transform integral
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int spec_w = 0; spec_w < s.parameters.iterations_w_resolution; spec_w++ ) {
        std::vector<Scalar> expfunc;
        expfunc.reserve( dim / s.parameters.iterations_wtau_skip );
        for ( int spec_tau = 0; spec_tau < dim; spec_tau += s.parameters.iterations_wtau_skip ) {
            expfunc.emplace_back( std::exp( -1i * spectrum_frequency_w.at( spec_w ) * (double)(spec_tau)*s.parameters.t_step ) );
        }
        for ( long unsigned int i = 0; i < dim; i += s.parameters.iterations_t_skip ) {
            double t_t = ( (double)i ) * s.parameters.t_step;
            int dim2 = (int)std::floor( ( s.parameters.t_end - s.parameters.t_start - t_t ) / ( s.parameters.t_step ) );
            for ( int j = 0; j < dim2; j += s.parameters.iterations_wtau_skip ) {
                out.at( spec_w ) += expfunc.at( j / s.parameters.iterations_wtau_skip ) * akf_mat( i, j ) / s.parameters.t_step / s.parameters.t_step / ( (double)s.parameters.iterations_t_skip * s.parameters.iterations_wtau_skip );
            }
            Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "Spectrum: " );
        }
        timer.iterate();
    }
    // Final output and timer end
    Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "Spectrum", PROGRESS_FORCE_OUTPUT );
    timer.end();
    // Save output
    Log::L2( "Done, saving to {}... ", fileOutputName );
    std::string filepath = s.parameters.subfolder + fileOutputName;
    FILE *spectrumfile = std::fopen( filepath.c_str(), "w" );
    if ( !spectrumfile ) {
        Log::L2( "Failed to open outputfile for spectrum!\n" );
        return false;
    }
    for ( int spec_w = 0; spec_w < s.parameters.iterations_w_resolution; spec_w++ ) {
        fmt::print( spectrumfile, "{0:.15e}\t{1:.15e}\n", spectrum_frequency_w[spec_w], real( out[spec_w] * s.parameters.t_step * s.parameters.t_step * (double)( s.parameters.iterations_t_skip * s.parameters.iterations_t_skip ) ) );
    }
    std::fclose( spectrumfile );
    Log::L2( "Done!\n" );
    // Check if we need to save calculated g-functions
    if ( s.parameters.numerics_calculate_g2 ) {
        Log::L2( "Saving g1 to temporary matrix {}", cache_index );
        if ( cache_index == 1 )
            cache1 = akf_mat;
        else if ( cache_index == 2 )
            cache2 = akf_mat;
        Log::L2( " succesfull\n" );
    }
    // Sucessfull
    return true;
}