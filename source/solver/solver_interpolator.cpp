// Interpolate std:vector linear/monotone
// Interpolate using G1/2 gridsettings for i'th index column
// --> no 2d interpolation needed then

#include "solver/solver.h"

std::vector<QDLC::SaveState> QDLC::Numerics::interpolate_curve( const std::vector<QDLC::SaveState> &input, double t_start, double t_end, const std::vector<double> &t_values, const std::vector<double> &t_steps, const std::map<double, size_t> &t_index, bool output_handler ) {
    size_t current_index = std::min<size_t>( size_t( std::lower_bound( t_values.begin(), t_values.end(), t_start ) - t_values.begin() ), t_values.size() - 2 ); //t_index.at(t_start);
    size_t num_of_points = (size_t)( t_values.size() - current_index );
    size_t max_index = t_values.size() - 2;
    std::vector<QDLC::SaveState> ret;
    //Log::L2("Interpolating from {} to {} to vector of size {} with timevector of size {}, anticipaed return suize is {}\n",t_start,t_end,input.size(), t_values.size(), t_values.size() - current_index);
    // Do a very simple linear interpolation. Linear and monotone interpolation should be changed by paramter
    size_t i = 1;
    for ( double t_t = t_start; t_t <= t_end; t_t += t_steps[current_index] ) {
        while ( i < input.size() - 1 and t_t > input[i].t ) {
            i++;
        }
        double first = i > 0 ? input[i - 1].t : 0.0;
        double second = input[i].t;
        double f = ( t_t - first ) / ( second - first );
        Sparse mat = input[i - 1].mat + f * ( input[i].mat - input[i - 1].mat );
        ret.push_back( { mat, t_t } );
        //current_index = std::min<size_t>( current_index + 1, t_values.size() - 2 );
        if ( current_index < max_index )
            current_index++;
    }
    //Log::L2("Done interpolating, returnign vector of size {}\n",ret.size());
    return ret;
}

std::vector<QDLC::SaveState> QDLC::Numerics::interpolate_curve( const std::vector<QDLC::SaveState> &input, double t_start, double t_end, double t_step, int threads, int order, bool output_handler ) {
    size_t num_of_points = (size_t)( ( t_end - t_start ) / t_step );
    std::vector<QDLC::SaveState> ret;
    Log::L2( "Interpolation order = {}\n", order );
    if ( order == 0 ) {
        // Linear Interpolation
        ret.reserve( num_of_points );

        // Do a very simple linear interpolation. Linear and monotone interpolation should be changed by paramter
        size_t i = 1;
        for ( double t_t = t_start; t_t < t_end; t_t += t_step ) {
            while ( i < input.size() - 1 and t_t > input[i].t ) {
                i++;
            }
            double first = i > 0 ? input[i - 1].t : 0.0;
            double second = input[i].t;
            double f = ( t_t - first ) / ( second - first );
            Sparse mat = input[i - 1].mat + f * ( input[i].mat - input[i - 1].mat );
            ret.push_back( { mat, t_t } );
        }
    } else if ( order == 1 ) {
        // Quintic Hermite Interpolation
        Log::L2( "Using Quintic Hermite Interpolation...\n" );
        // Derivatives
        std::vector<Sparse> first_derivative, second_derivative;
        Log::L2( "Calculating first derivatives...\n" );
        for ( int i = 0; i < input.size() - 1; i++ ) {
            first_derivative.emplace_back( ( input[i + 1].mat - input[i].mat ) / ( input[i + 1].t - input[i].t ) );
        }
        first_derivative.emplace_back( first_derivative.back() );
        Log::L2( "Calculating second derivatives...\n" );
        for ( int i = 0; i < first_derivative.size() - 1; i++ ) {
            second_derivative.emplace_back( ( first_derivative[i + 1] - first_derivative[i] ) / ( input[i + 1].t - input[i].t ) );
        }
        second_derivative.emplace_back( second_derivative.back() );
        Log::L2( "Sizes: {} {} {}\n", input.size(), first_derivative.size(), second_derivative.size() );

        // Interpolate
        size_t i = 1;
        ret.push_back( { input[0].mat, t_start } );
        for ( double t_t = t_start; t_t < t_end; t_t += t_step ) {
            while ( i < input.size() - 1 and t_t > input[i].t ) {
                i++;
            }
            i--;
            if ( input[i].t == t_t ) {
                ret.push_back( { input[i].mat, t_t } );
                continue;
            }
            double x = t_t;
            double xmx0 = x - input[i].t;
            double x1mx0 = input[i + 1].t - input[i].t;
            double xmx1 = x - input[i + 1].t;
            const Sparse &fx0 = input[i].mat;
            const Sparse &fx1 = input[i + 1].mat;
            const Sparse &fdx0 = first_derivative[i];
            const Sparse &fdx1 = first_derivative[i + 1];
            const Sparse &fddx0 = second_derivative[i];
            const Sparse &fddx1 = second_derivative[i + 1];

            double xmx0_x1mx0 = xmx0 / x1mx0;
            double xmx1_x1mx0 = xmx1 / x1mx0;

            auto d1 = ( fx1 - fx0 - fdx0 * x1mx0 - 0.5 * fddx0 * x1mx0 * x1mx0 ) * xmx0_x1mx0 * xmx0_x1mx0 * xmx0_x1mx0;
            auto d2 = ( 3.0 * fx0 - 3.0 * fx1 + ( 2 * fdx0 + fdx1 ) * x1mx0 + 0.5 * fddx0 * x1mx0 * x1mx0 ) * xmx0_x1mx0 * xmx0_x1mx0 * xmx0_x1mx0 * xmx1_x1mx0;
            auto d3 = ( 6.0 * fx1 - 6.0 * fx0 - 3.0 * ( fdx0 + fdx1 ) * x1mx0 + 0.5 * ( fddx1 - fddx0 ) * x1mx0 * x1mx0 ) * xmx0_x1mx0 * xmx0_x1mx0 * xmx0_x1mx0 * xmx1_x1mx0 * xmx1_x1mx0;
            Sparse mat = fx0 + fdx0 * xmx0 + 0.5 * fddx0 * xmx0 * xmx0 + d1 + d2 + d3;
            ret.push_back( { mat, t_t } );
        }
    } else {
        // Cubic monotone with library
        // Generate N^2 vectors from the initial density matrices
        Log::L2( "Using Monotone Hermite Interpolation\n" );
        Timer &interpolateTimer = Timers::create( "Monotone Spline Interpolation" );
        ProgressBar progressbar = ProgressBar();
        interpolateTimer.start();

        size_t matrix_dimension = input.front().mat.rows();
        long unsigned int input_length = input.size();
        std::vector<std::vector<double>> raw_real( matrix_dimension * matrix_dimension );
        std::vector<std::vector<double>> raw_imag( matrix_dimension * matrix_dimension );
        std::vector<double> time_raw;
        time_raw.reserve( input_length );
        for ( int k = 0; k < input_length; k++ ) {
            time_raw.emplace_back( input.at( k ).t );
        }
#pragma omp parallel for num_threads( threads )
        for ( int i = 0; i < matrix_dimension; i++ ) {
            for ( int j = 0; j < matrix_dimension; j++ ) {
                std::vector<double> cur_real;
                std::vector<double> cur_imag;
                cur_real.reserve( input_length );
                cur_imag.reserve( input_length );
                for ( int k = 0; k < input_length; k++ ) {
                    Dense mat = Dense( input.at( k ).mat );
                    cur_real.emplace_back( std::real( mat( i, j ) ) );
                    cur_imag.emplace_back( std::imag( mat( i, j ) ) );
                }
                raw_real[i * matrix_dimension + j] = cur_real;
                raw_imag[i * matrix_dimension + j] = cur_imag;
            }
            interpolateTimer.iterate();
            Timers::outputProgress( interpolateTimer, progressbar, i, matrix_dimension, "Copying DM Elements into single Vector" );
        }

        // Interpolate each vector on its own
        std::vector<std::vector<double>> interpolated_real( matrix_dimension * matrix_dimension );
        std::vector<std::vector<double>> interpolated_imag( matrix_dimension * matrix_dimension );
        std::vector<double> time_out;
        time_out.reserve( input_length );
        for ( long unsigned int k = 0; k < num_of_points; k++ ) {
            time_out.emplace_back( t_start + k * t_step );
        }
#pragma omp parallel for num_threads( threads )
        for ( int i = 0; i < matrix_dimension; i++ ) {
            for ( int j = 0; j < matrix_dimension; j++ ) {
                int index = i * matrix_dimension + j;
                Interpolant interpolant_real = Interpolant( time_raw, raw_real.at( index ), "monotone" );
                Interpolant interpolant_imag = Interpolant( time_raw, raw_imag.at( index ), "monotone" );
                interpolated_real[i * matrix_dimension + j] = interpolant_real.evaluate( time_out );
                interpolated_imag[i * matrix_dimension + j] = interpolant_imag.evaluate( time_out );
            }
            interpolateTimer.iterate();
            Timers::outputProgress( interpolateTimer, progressbar, i, matrix_dimension, "Interpolating Vectors" );
        }
        // Write new vectors back into matices
        for ( int k = 0; k < num_of_points; k++ ) {
            Dense cur( matrix_dimension, matrix_dimension );
#pragma omp parallel for num_threads( threads )
            for ( int i = 0; i < matrix_dimension; i++ ) {
                for ( int j = 0; j < matrix_dimension; j++ ) {
                    int index = i * matrix_dimension + j;
                    cur( i, j ) = Scalar( interpolated_real[index][k], interpolated_imag[index][k] );
                }
            }
            ret.push_back( { cur.sparseView(), time_out[k] } );
            ret.back().mat.makeCompressed();
            interpolateTimer.iterate();
            Timers::outputProgress( interpolateTimer, progressbar, k, num_of_points, "Writing Vectors back into DM" );
        }
        interpolateTimer.end();
        Timers::outputProgress( interpolateTimer, progressbar, interpolateTimer.getTotalIterationNumber(), interpolateTimer.getTotalIterationNumber(), "Monotone Spline Interpolation", true );
    }
    Log::L3( "Done Interpolating!\n" );
    return ret;
}
