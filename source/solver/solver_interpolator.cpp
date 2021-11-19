// Interpolate std:vector linear/monotone
// Interpolate using G1/2 gridsettings for i'th index column
// --> no 2d interpolation needed then

#include "solver/solver.h"

//TODO: alles interpolieren der DM (aus RK45 oder pathint) hier machen! art der interpolation (linear, cubisch, monotone, ...) per parameter übergeben!
std::vector<QDLC::SaveState> QDLC::Numerics::interpolate_curve( const std::vector<QDLC::SaveState> &input, double t_start, double t_end, const std::vector<double> &t_values, const std::vector<double> &t_steps, const std::map<double,size_t>& t_index, bool output_handler ) {
    size_t current_index = std::min<size_t>(size_t(std::lower_bound(t_values.begin(), t_values.end(), t_start) - t_values.begin()),t_values.size()-2);//t_index.at(t_start);
    size_t num_of_points = (size_t)( t_values.size() - current_index );
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
        current_index = std::min<size_t>( current_index + 1, t_values.size() - 2 );
    }
    //Log::L2("Done interpolating, returnign vector of size {}\n",ret.size());
    return ret;
}
std::vector<QDLC::SaveState> QDLC::Numerics::interpolate_curve( const std::vector<QDLC::SaveState> &input, double t_start, double t_end, double t_step, bool output_handler ) {
    size_t num_of_points = (size_t)( ( t_end - t_start ) / t_step );
    std::vector<QDLC::SaveState> ret;

    if ( true ) {
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
    } else if ( false ) {
        // Cubic interpolation
        // Derivatives
        std::vector<Sparse> dys, ms;
        std::vector<double> dxs;
        for ( size_t i = 0; i < input.size() - 1; i++ ) {
            dxs.emplace_back( input[i + 1].t - input[i].t );
            dys.emplace_back( input[i + 1].mat - input[i].mat );
            ms.emplace_back( ( input[i + 1].mat - input[i].mat ) / ( input[i + 1].t - input[i].t ) );
        }
        // Degree-1 coefficients
        std::vector<Sparse> c1s;
        c1s.reserve( dxs.size() - 1 );
        c1s.emplace_back( ms.front() );
        for ( size_t i = 0; i < dxs.size() - 1; i++ ) {
            Sparse m = ms[i];
            Sparse mNext = ms[i + 1];
            double dx = dxs[i];
            double dxNext = dxs[i + 1];
            double common = dx + dxNext;
            Sparse map_real = ( m * mNext ).unaryExpr( []( Scalar x ) { return std::real( x ) <= 0 ? Scalar( 0.0 ) : Scalar( 1.0 ); } );
            Sparse map_imag = ( m * mNext ).unaryExpr( []( Scalar x ) { return std::imag( x ) <= 0 ? Scalar( 0.0 ) : Scalar( 1.0 ); } );
            Sparse mat = 3.0 * common * ( ( common + dxNext ) * m.cwiseInverse() + ( common + dx ) * mNext.cwiseInverse() ).cwiseInverse();
            c1s.emplace_back( mat.real().cwiseProduct( map_real ) + 1.0i * mat.imag().cwiseProduct( map_imag ) );
        }
        // Degree-2 and Degree-3 coefficients
        std::vector<Sparse> c2s, c3s;
        c2s.reserve( c1s.size() );
        c3s.reserve( c1s.size() );
        for ( size_t i = 0; i < c1s.size() - 1; i++ ) {
            double invDx = 1.0 / dxs[i];
            Sparse common = c1s[i] + c1s[i + 1] - 2.0 * ms[i];
            c2s.emplace_back( ( ms[i] - c1s[i] - common ) * invDx );
            c3s.emplace_back( common * invDx * invDx );
        }

        //Log::L2("Sizes: mat: {}, c1s: {}, c2s: {}, c3s: {}\n",input.size(),c1s.size(),c2s.size(),c3s.size());
        // Interpolate
        size_t i;
        for ( double t_t = t_start; t_t < t_end; t_t += t_step ) {
            while ( i < input.size() - 1 and t_t > input[i].t ) {
                i++;
            }
            i = std::max<int>( 0, i - 1 );
            if ( i >= c3s.size() )
                break;
            if ( input[i].t == t_t ) {
                ret.push_back( { input[i].mat, t_t } );
                continue;
            }
            double diff = t_t - input[i].t;
            Sparse mat = input[i].mat + c1s[i] * diff + c2s[i] * diff * diff + c3s[i] * diff * diff * diff;
            ret.push_back( { mat, t_t } );
        }
    } else {
        // Cubic monotone with library
        // Generate N^2 vectors from the initial density matrices
        size_t matrix_dimension = input.front().mat.rows();
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
            }
        }

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
            time_out.emplace_back( t_start + k * t_step );
        }
        for ( int i = 0; i < input.at( 0 ).mat.rows(); i++ ) {
            for ( int j = 0; j < input.at( 0 ).mat.cols(); j++ ) {
                int index = i * input.at( 0 ).mat.cols() + j;
                interpolant_real = Interpolant( time_raw, raw_real.at( index ), "monotone" );
                interpolant_imag = Interpolant( time_raw, raw_imag.at( index ), "monotone" );
                interpolated_real.emplace_back( interpolant_real.evaluate( time_out ) );
                interpolated_imag.emplace_back( interpolant_imag.evaluate( time_out ) );
            }
        }
        // Write new vectors back into matices
        for ( int k = 0; k < num_of_points; k++ ) {
            Dense cur( input.at( 0 ).mat.rows(), input.at( 0 ).mat.cols() );
            for ( int i = 0; i < input.at( 0 ).mat.rows(); i++ ) {
                for ( int j = 0; j < input.at( 0 ).mat.cols(); j++ ) {
                    int index = i * input.at( 0 ).mat.cols() + j;
                    cur( i, j ) = Scalar( interpolated_real.at( index ).at( k ), interpolated_imag.at( index ).at( k ) );
                }
            }
            ret.emplace_back( QDLC::SaveState( cur.sparseView(), time_out.at( k ) ) );
            ret.back().mat.makeCompressed();
        }
    }
    Log::L3( "Done Interpolating!\n" );
    return ret;
}
