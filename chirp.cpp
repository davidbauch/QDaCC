#include "chirp.h"

Chirp::Chirp( Chirp::Inputs &_inputs ) : inputs( _inputs ) {
    counter_evaluated = 0;
    counter_returned = 0;
    logs.level2( "Creating Chirp with {} points of type {}... ", inputs.t.size(), inputs.type );
    int n = (int)( ( inputs.t_end - inputs.t_start ) / inputs.t_step * 2.0 + 5 );
    chirparray.reserve( n );
    timearray.reserve( n );
    //if (p.numerics_order_t == 5 || p.numerics_order_tau == 5)
    //    steps = { 0, 1./5.*p.t_step, 3./10.*p.t_step, 1./2.*p.t_step, 4./5.*p.t_step, 8./9.*p.t_step };
    //else
    steps = {0, 0.5 * inputs.t_step};
    if (inputs.type.compare("sine") != 0) {
        logs.level2( "Done initializing class, creating interpolant... " );
        interpolant = Interpolant( inputs.t, inputs.y, inputs.ddt, inputs.type );
        logs.level2( "Done creating interpolant, generating chirp... " );
    } else {
        logs.level2( "No interpolant class used." );
        inputs.isSineChirp = true;
    }
    generate();
    logs.level2( "Done!\n" );
}

double Chirp::evaluate( double t ) {
    if (inputs.isSineChirp) {
        double ret = 0;
        for (long unsigned int i = 0; i < inputs.y.size(); i++) {
            ret += inputs.y.at(i)*std::sin(6.283*(t-inputs.t.at(i))*inputs.ddt.at(i));
        }
        return ret;
    }
    return interpolant.evaluate( t );
}
double Chirp::evaluate_derivative( double t ) {
    return ( evaluate( t ) - evaluate( t - inputs.t_step ) ) / inputs.t_step;
}
double Chirp::evaluate_integral( double t ) {
    if ( chirparray_integral.size() == 0 || t == 0.0 )
        return evaluate( t ) * inputs.t_step;
    return chirparray_integral.at( std::max( (int)std::floor( t / inputs.t_step ), (int)chirparray_integral.size()-1 ) ) + evaluate( t ) * inputs.t_step;
}

void Chirp::generate() {
    for ( double t1 = inputs.t_start; t1 < inputs.t_end + 10 * inputs.t_step; t1 += inputs.t_step / 2.0 ) {
        chirparray.push_back( evaluate( t1 ) );
        chirparray_derivative.push_back( evaluate_derivative( t1 ) );
        chirparray_integral.push_back( evaluate_integral( t1 ) );
        timearray.push_back( t1 );
    }
    chirparray.shrink_to_fit();
    timearray.shrink_to_fit();
    size = chirparray.size();
    logs.level2( "chirparray.size() = {}... ", size );
}

void Chirp::fileOutput( std::string filepath ) {
    FILE *chirpfile = std::fopen( filepath.c_str(), "w" );
    if ( !chirpfile ) {
        logs.level2( "Failed to open outputfile for chirp!\n" );
        return;
    }
    fmt::print( chirpfile, "Time\tChirp\tDerivative\tIntegral\n" );
    for ( long unsigned int i = 0; i < chirparray.size(); i++ ) {
        fmt::print( chirpfile, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", timearray.at( i ), chirparray.at( i ), chirparray_derivative.at( i ), chirparray_integral.at( i ) );
    }
    std::fclose( chirpfile );
}

void Chirp::Inputs::add( double _t, double _y, double _ddt ) {
    t.emplace_back( _t );
    y.emplace_back( _y );
    ddt.emplace_back( _ddt );
}

void Chirp::Inputs::add( std::vector<double> &_t, std::vector<double> &_y, std::vector<double> &_ddt ) {
    if ( !( _t.size() == _y.size() && _t.size() == _ddt.size() ) ) {
        logs.level2( "Input arrays don't have the same length! No Vectors are created, initializing interpolant will fail!\n" );
        return;
    }
    for ( int i = 0; i < (int)_t.size(); i++ ) {
        t.emplace_back( _t.at( i ) );
        y.emplace_back( _y.at( i ) );
        ddt.emplace_back( _ddt.at( i ) );
    }
}

double Chirp::get( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate( t );
    }
    int i = std::max( 0.0, std::floor( t / inputs.t_step - 1 ) ) * steps.size();
    while ( timearray.at( i ) < t ) {
        i++;
    }
    if ( std::abs( timearray.at( i ) - t ) > 1E-15 ) {
        counter_evaluated++;
        return evaluate( t );
    }
    if ( i < 0 || i >= size ) {
        counter_evaluated++;
        return evaluate( t );
    }
    counter_returned++;
    return chirparray.at( i );
}
double Chirp::derivative( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_derivative( t );
    }
    int i = std::max( 0.0, std::floor( t / inputs.t_step - 1 ) ) * steps.size();
    while ( timearray.at( i ) < t ) {
        i++;
    }
    if ( std::abs( timearray.at( i ) - t ) > 1E-15 ) {
        counter_evaluated++;
        return evaluate_derivative( t );
    }
    if ( i < 0 || i >= size ) {
        counter_evaluated++;
        return evaluate_derivative( t );
    }
    counter_returned++;
    return chirparray_derivative.at( i );
}
double Chirp::integral( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_integral( t );
    }
    int i = std::max( 0.0, std::floor( t / inputs.t_step - 1 ) ) * steps.size();
    while ( timearray.at( i ) < t ) {
        i++;
    }
    if ( std::abs( timearray.at( i ) - t ) > 1E-15 ) {
        counter_evaluated++;
        return evaluate_integral( t );
    }
    if ( i < 0 || i >= size ) {
        counter_evaluated++;
        return evaluate_integral( t );
    }
    counter_returned++;
    return chirparray_integral.at( i );
}

void Chirp::log() {
    logs.level2( "Chirp evaluations/returns: {}/{}\n", counter_evaluated, counter_returned );
}