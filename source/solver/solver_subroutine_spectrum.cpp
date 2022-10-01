#include "solver/solver_ode.h"

// Description: Calculates the Eberly-Wódkiewicz spectrum. Uses precalculated values from akf_mat
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false

bool QDLC::Numerics::ODESolver::calculate_spectrum( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator, double frequency_center, double frequency_range, int resolution, int order, bool normalize, std::string s_g ) {
    // Set Number of Phonon cores to 1 because this memberfunction is already using multithreading
    s.parameters.numerics_maximum_secondary_threads = 1;
    // Calculate G1/2(t,tau) with given operator matrices
    if ( s_g.size() == 0 ) {
        s_g = order == 1 ? get_operators_purpose( { s_op_creator, s_op_annihilator } ) : get_operators_purpose( { s_op_creator, s_op_creator, s_op_annihilator, s_op_annihilator } );
        if ( order == 1 )
            calculate_g1( s, s_op_creator, s_op_annihilator, s_g );
        else
            calculate_g2( s, s_op_creator, s_op_annihilator, s_op_creator, s_op_annihilator, s_g );
    }
    auto &akf_mat = cache[s_g];
    auto &akf_mat_time = cache[s_g + "_time"];

    // Create Timer and Progressbar for the spectrum loop
    Timer &timer = Timers::create( "Spectrum (" + s_g + ")" );
    int totalIterations = resolution;
    ProgressBar progressbar = ProgressBar();
    timer.start();
    Log::L2( "Calculating spectrum... Calculating frequencies...\n" );

    // Calculate frequencies:
    std::vector<Scalar> spectrum_frequency_w;
    std::vector<Scalar> out;
    for ( int w = 0; w < resolution; w++ ) {
        spectrum_frequency_w.push_back( frequency_center - ( frequency_range ) + w / ( (double)resolution ) * ( 2. * ( frequency_range ) ) );
        out.push_back( 0 );
    }
    Log::L2( "Done, calculating fourier transform via direct integral...\n" );
    double t_step = ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? s.parameters.t_step_pathint : s.parameters.t_step );
    Log::L2( "Size = {} x {}, using dt = {}\n", akf_mat.rows(), akf_mat.cols(), t_step );
    // Calculate main fourier transform integral
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( int spec_w = 0; spec_w < resolution; spec_w++ ) {
        for ( long unsigned int i = 0; i < akf_mat.rows(); i++ ) {
            double dt = Numerics::get_tdelta( akf_mat_time, 0, i );
            double t_t = std::real( akf_mat_time( i, 0 ) );
            for ( int j = 0; j < akf_mat.cols() - i; j++ ) { // -i because of triangular grid
                double tau = std::imag( akf_mat_time( i, j ) ) - std::real( akf_mat_time( i, j ) );
                double dtau = Numerics::get_taudelta( akf_mat_time, i, j );
                out.at( spec_w ) += std::exp( -1.0i * spectrum_frequency_w.at( spec_w ) * tau ) * akf_mat( i, j ) * dtau * dt;
            }
        }
        Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), totalIterations, "Spectrum (" + s_g + "): " );
        out.at( spec_w ) = std::real( out.at( spec_w ) );
        timer.iterate();
    }
    // Normalize
    if ( normalize ) {
        Scalar vec_min = QDLC::Misc::vec_filter( out, []( const Scalar &a, const Scalar &b ) { return std::real( a ) < std::real( b ); } );
        Scalar vec_max = QDLC::Misc::vec_filter( out, []( const Scalar &a, const Scalar &b ) { return std::real( a ) > std::real( b ); } );
        std::for_each( out.begin(), out.end(), [&]( Scalar &val ) { val = ( val - vec_min ) / ( vec_max - vec_min ); } );
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