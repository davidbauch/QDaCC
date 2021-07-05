#include "solver/solver.h"

std::vector<SaveState> Solver::calculate_smooth_curve( const std::vector<SaveState> &input, double t_start, double t_end, int num_of_points, bool output_handler ) {
    Log::L2( " : Setting up the smooth curve interpolator... \n" );
    int matrix_dimension = input.at( 0 ).mat.rows() * input.at( 0 ).mat.cols();
    Timer &timer = Timers::create( "Interpolator" );
    int maximum_iterations = matrix_dimension * 2 + num_of_points;
    ProgressBar progressbar = ProgressBar( maximum_iterations, 60, 0, BAR_VERTICAL, true, 0.1, {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"} );
    timer.start();

    // Stuff
    std::vector<SaveState> ret;
    ret.reserve( num_of_points );
    double delta_t = ( t_end - t_start ) / ( (double)num_of_points );

    Log::L2( " : : Copying initial SaveStates into seperate vectors...\n" );
    // Generate N^2 vectors from the initial density matrices
    long unsigned int input_length = input.size();
    std::vector<std::vector<double>> raw_real;
    std::vector<std::vector<double>> raw_imag;
    std::vector<double> time_raw;
    time_raw.reserve( input_length );
    for ( int k = 0; k < input_length; k++ ) {
        time_raw.emplace_back( input.at( k ).t );
    }
    for ( int i = 0; i < input.at( 0 ).mat.rows(); i++ ) {
        for ( int j = 0; j < input.at( 0 ).mat.cols(); j++ ) {
            std::vector<double> cur_real;
            std::vector<double> cur_imag;
            cur_real.reserve( input_length );
            cur_imag.reserve( input_length );
            for ( int k = 0; k < input_length; k++ ) {
                Dense mat = Dense( input.at( k ).mat );
                cur_real.emplace_back( std::real( mat( i, j ) ) );
                cur_imag.emplace_back( std::imag( mat( i, j ) ) );
            }
            raw_real.emplace_back( cur_real );
            raw_imag.emplace_back( cur_imag );
            // Time
            timer.iterate();
            Timers::outputProgress( output_handler, timer, progressbar, maximum_iterations, "Interpolator: " );
        }
    }

    Log::L2( " : : Interpolating vectors...\n" );
    // Interpolate each vector on its own
    std::vector<std::vector<double>> interpolated_real;
    std::vector<std::vector<double>> interpolated_imag;
    interpolated_real.reserve( matrix_dimension );
    interpolated_imag.reserve( matrix_dimension );
    Interpolant interpolant_real;
    Interpolant interpolant_imag;
    std::vector<double> time_out;
    time_out.reserve( input_length );
    for ( long unsigned int k = 0; k < num_of_points; k++ ) {
        time_out.emplace_back( t_start + k * delta_t );
    }
    for ( int i = 0; i < input.at( 0 ).mat.rows(); i++ ) {
        for ( int j = 0; j < input.at( 0 ).mat.cols(); j++ ) {
            int index = i * input.at( 0 ).mat.cols() + j;
            interpolant_real = Interpolant( time_raw, raw_real.at( index ), "monotone" );
            interpolant_imag = Interpolant( time_raw, raw_imag.at( index ), "monotone" );
            interpolated_real.emplace_back( interpolant_real.evaluate( time_out ) );
            interpolated_imag.emplace_back( interpolant_imag.evaluate( time_out ) );
            // Time
            timer.iterate();
            Timers::outputProgress( output_handler, timer, progressbar, maximum_iterations, "Interpolator: " );
        }
    }

    Log::L2( " : : Writing vectors back to matrices" );
    // Write new vectors back into matices
    maximum_iterations = num_of_points;
    for ( int k = 0; k < num_of_points; k++ ) {
        Dense cur( input.at( 0 ).mat.rows(), input.at( 0 ).mat.cols() );
        for ( int i = 0; i < input.at( 0 ).mat.rows(); i++ ) {
            for ( int j = 0; j < input.at( 0 ).mat.cols(); j++ ) {
                int index = i * input.at( 0 ).mat.cols() + j;
                cur( i, j ) = Scalar( interpolated_real.at( index ).at( k ), interpolated_imag.at( index ).at( k ) );
            }
        }
        ret.emplace_back( SaveState( cur.sparseView(), time_out.at( k ) ) );
        // Time
        timer.iterate();
        Timers::outputProgress( output_handler, timer, progressbar, maximum_iterations, "Interpolator: " );
    }
    Timers::outputProgress( output_handler, timer, progressbar, maximum_iterations, "Interpolator: ", PROGRESS_FORCE_OUTPUT );
    timer.end();
    Log::L2( " : Done!\n" );
    return ret;
}