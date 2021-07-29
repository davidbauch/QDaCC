#include "solver/solver_ode.h"

// Description: Calculates the Eberly-Wódkiewicz spectrum. Uses precalculated values from akf_mat
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false

bool ODESolver::calculate_spectrum( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator, double frequency_center, double frequency_range, int resolution ) {
    // Send system command to change to single core mainprogram, because this memberfunction is already using multithreading
    s.command( Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
    // Calculate G1(t,tau) with given operator matrices
    std::string s_g1 = "G1-" + s_op_creator + "-" + s_op_annihilator;
    auto [op_creator, op_annihilator] = calculate_g1( s, s_op_creator, s_op_annihilator, s_g1 );
    Dense &akf_mat = cache[s_g1];

    // Create Timer and Progressbar for the spectrum loop
    Timer &timer = Timers::create( "Spectrum (" + s_g1 + ")" );
    int totalIterations = resolution; //getIterationNumberSpectrum( s );
    ProgressBar progressbar = ProgressBar( resolution );
    timer.start();
    Log::L2( "Calculating spectrum... Calculating frequencies...\n" );

    //Calculate frequencies:
    std::vector<Scalar> spectrum_frequency_w;
    std::vector<Scalar> out;
    for ( int w = 0; w < resolution; w++ ) {
        spectrum_frequency_w.push_back( frequency_center - ( frequency_range ) + w / ( (double)resolution ) * ( 2. * ( frequency_range ) ) );
        out.push_back( 0 );
    }
    Log::L2( "Done, calculating fourier transform via direct integral...\n" );
    Log::L2( "Size = {} x {}\n", akf_mat.rows(), akf_mat.cols() );
    // Calculate main fourier transform integral
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int spec_w = 0; spec_w < resolution; spec_w++ ) {
        std::vector<Scalar> expfunc;
        expfunc.reserve( dim / s.parameters.iterations_wtau_skip );
        for ( int spec_tau = 0; spec_tau < dim; spec_tau += s.parameters.iterations_wtau_skip ) {
            expfunc.emplace_back( std::exp( -1i * spectrum_frequency_w.at( spec_w ) * (double)(spec_tau)*s.parameters.t_step ) );
        }
        for ( long unsigned int i = 0; i < dim; i += s.parameters.iterations_t_skip ) {
            double t_t = ( (double)i ) * s.parameters.t_step;
            int dim2 = (int)std::floor( ( s.parameters.t_end - s.parameters.t_start - t_t ) / ( s.parameters.t_step ) );
            for ( int j = 0; j < dim2; j += s.parameters.iterations_wtau_skip ) {
                out.at( spec_w ) += expfunc.at( j / s.parameters.iterations_wtau_skip ) * akf_mat( i, j ); // / s.parameters.t_step / s.parameters.t_step / ( (double)s.parameters.iterations_t_skip * s.parameters.iterations_wtau_skip );
            }
        }
        Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "Spectrum (" + s_g1 + "): " );
        out.at( spec_w ) = std::real( out.at( spec_w ) );
        timer.iterate();
    }
    // Final output and timer end
    Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, totalIterations, "Spectrum (" + s_g1 + ")", PROGRESS_FORCE_OUTPUT );
    timer.end();
    // Save output
    to_output["Spectrum_frequency"][s_g1] = spectrum_frequency_w;
    to_output["Spectrum"][s_g1] = out;

    // Sucessfull
    return true;
}