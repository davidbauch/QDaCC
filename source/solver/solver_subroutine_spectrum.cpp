#include "solver/solver_ode.h"

// Description: Calculates the Eberly-Wódkiewicz spectrum. Uses precalculated values from akf_mat
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false

bool QDACC::Numerics::ODESolver::calculate_spectrum( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator, double frequency_center, double frequency_range, int resolution, int order, bool normalize, std::string s_g ) {
    // Set Number of Phonon cores to 1 because this memberfunction is already using multithreading
    // s.parameters.numerics_maximum_secondary_threads = 1;
    // Calculate G1/2(t,tau) with given operator matrices
    if ( s_g.size() == 0 ) {
        if ( order == 1 ) {
            s_g = get_operators_purpose( { s_op_creator, s_op_annihilator } );
            calculate_g1( s, s_op_creator, s_op_annihilator, s_g );
        } else {
            s_g = get_operators_purpose( { s_op_creator, s_op_creator, s_op_annihilator, s_op_annihilator } );
            calculate_g2( s, s_op_creator, s_op_creator, s_op_annihilator, s_op_annihilator, s_g );
        }
    }
    auto &akf_mat = cache[s_g];

    // Create Timer and Progressbar for the spectrum loop
    Timer &timer = Timers::create( "Spectrum (" + s_g + ")" );
    int totalIterations = resolution;
    auto progressbar = ProgressBar();
    timer.start();
    Log::L2( "[Spectrum] Calculating spectrum... Calculating frequencies...\n" );

    // Calculate frequencies:
    std::vector<Scalar> spectrum_frequency_w;
    std::vector<Scalar> out;
    for ( int w = 0; w < resolution; w++ ) {
        spectrum_frequency_w.push_back( frequency_center - frequency_range + w / ( (double)resolution ) * ( 2. * frequency_range ) );
        out.push_back( 0 );
    }
    Log::L2( "[Spectrum] Done, calculating fourier transform via direct integral...\n" );
    Log::L2( "[Spectrum] Matrix Size = {0} x {0}\n", akf_mat.dim() );
    // Calculate main fourier transform integral
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( int spec_w = 0; spec_w < resolution; spec_w++ ) {
        for ( long unsigned int i = 0; i < akf_mat.dim(); i++ ) {
            double dt = akf_mat.dt( i );
            for ( int j = 0; j < akf_mat.dim() - i; j++ ) { // -i because of triangular grid
                double tau = akf_mat.tau( j, i );
                double dtau = akf_mat.dtau( j, i );
                out.at( spec_w ) += std::exp( -1.0i * spectrum_frequency_w.at( spec_w ) * tau ) * akf_mat( i, j ) * dtau * dt;
            }
        }
        Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), totalIterations, "Spectrum (" + s_g + "): " );
        out.at( spec_w ) = std::real( out.at( spec_w ) );
        timer.iterate();
    }
    // Normalize
    if ( normalize ) {
        Scalar vec_min = QDACC::Misc::vec_filter( out, []( const Scalar &a, const Scalar &b ) { return std::real( a ) < std::real( b ); } );
        Scalar vec_max = QDACC::Misc::vec_filter( out, []( const Scalar &a, const Scalar &b ) { return std::real( a ) > std::real( b ); } );
        std::ranges::for_each( out, [&]( Scalar &val ) { val = ( val - vec_min ) / ( vec_max - vec_min ); } );
    }
    // Final output and timer end
    timer.end();
    Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), totalIterations, "Spectrum (" + s_g + ")", Timers::PROGRESS_FORCE_OUTPUT );
    // Save output
    to_output["Spectrum_frequency"][s_g] = spectrum_frequency_w;
    to_output["Spectrum"][s_g] = out;

    // Sucessfull
    return true;
}