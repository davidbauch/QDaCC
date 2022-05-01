#include "chirp.h"
#include "solver/solver.h"

Chirp::Chirp( Chirp::Inputs &_inputs ) : inputs( _inputs ) {
    counter_evaluated = 0;
    counter_returned = 0;
    Log::L2( "Creating Chirp with {} points of type {}...\n", inputs.t.size(), inputs.type );
    int n = (int)( ( inputs.t_end - inputs.t_start ) / inputs.t_step * 2.0 + 5 );
    if ( inputs.order > 5 )
        steps = { QDLC::Numerics::RKCoefficients::a1 * inputs.t_step, QDLC::Numerics::RKCoefficients::a2 * inputs.t_step, QDLC::Numerics::RKCoefficients::a3 * inputs.t_step, QDLC::Numerics::RKCoefficients::a4 * inputs.t_step, QDLC::Numerics::RKCoefficients::a5 * inputs.t_step };
    else
        steps = { 0, 0.5 * inputs.t_step };
    if ( inputs.type.compare( "sine" ) != 0 ) {
        Log::L2( "Done initializing class, creating interpolant...\n" );
        interpolant = Interpolant( inputs.t, inputs.y, inputs.ddt, inputs.type );
        Log::L2( "Done creating interpolant, generating chirp...\n" );
    } else {
        Log::L2( "No interpolant class used, using sine chirp instead." );
        inputs.isSineChirp = true;
    }
    generate();
}

double Chirp::evaluate( double t ) {
    if ( inputs.isSineChirp ) {
        double ret = 0;
        for ( long unsigned int i = 0; i < inputs.y.size(); i++ ) {
            ret += inputs.y.at( i ) * std::sin( 6.283 * t * inputs.ddt.at( i ) + 6.283 * inputs.t.at( i ) );
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
    return integral( t - inputs.t_step ) + evaluate( t ) * inputs.t_step;
}

void Chirp::generate() {
    for ( double t1 = inputs.t_start; t1 < inputs.t_end + inputs.t_step * steps.size(); t1 += inputs.t_step ) {
        for ( int i = 0; i < (int)steps.size(); i++ ) {
            double t = t1 + steps[i];
            chirparray[t] = get( t );
        }
        chirparray_derivative[t1] = derivative( t1 );
        chirparray_integral[t1] = 0.0; //integral( t1 ); // FIXME: segmentation fault, just integrade properly.
    }
    size = chirparray.size();
    Log::L2( "chirparray.size() = {}... ", size );
}

void Chirp::fileOutput( std::string filepath ) {
    FILE *chirpfile = std::fopen( filepath.c_str(), "w" );
    if ( !chirpfile ) {
        Log::L2( "Failed to open outputfile for chirp!\n" );
        return;
    }
    fmt::print( chirpfile, "Time\tChirp\tDerivative\tIntegral\n" );
    for ( double t = inputs.t_start; t < inputs.t_end + inputs.t_step * steps.size(); t += inputs.t_step ) {
        fmt::print( chirpfile, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t, chirparray[t], chirparray_derivative[t], chirparray_integral[t] );
    }
    std::fclose( chirpfile );
}

void Chirp::Inputs::add( double _t, double _y, double _ddt ) {
    t.emplace_back( _t );
    y.emplace_back( _y );
    ddt.emplace_back( _ddt );
}

void Chirp::Inputs::add( std::vector<Parameter> &_t, std::vector<Parameter> &_y, std::vector<Parameter> &_ddt ) {
    if ( !( _t.size() == _y.size() && _t.size() == _ddt.size() ) ) {
        Log::L2( "Input arrays don't have the same length! No Vectors are created, initializing interpolant will fail!\n" );
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
    if ( chirparray.count( t ) == 0 ) {
        counter_evaluated++;
        double val = evaluate( t );
#pragma omp critical
        chirparray[t] = val;
        return val;
    }
    counter_returned++;
    return chirparray[t];
}
double Chirp::derivative( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_derivative( t );
    }
    if ( chirparray_derivative.count( t ) == 0 ) {
        counter_evaluated++;
        double val = evaluate_derivative( t );
#pragma omp critical
        chirparray_derivative[t] = val;
        return val;
    }
    counter_returned++;
    return chirparray_derivative[t];
}
double Chirp::integral( double t, bool force_evaluate ) {
    if ( force_evaluate ) {
        counter_evaluated++;
        return evaluate_integral( t );
    }
    if ( chirparray_integral.count( t ) == 0 ) {
        counter_evaluated++;
        double val = evaluate_integral( t );
#pragma omp critical
        chirparray_integral[t] = val;
        return val;
    }
    counter_returned++;
    return chirparray_integral[t];
}

void Chirp::log() {
    Log::L2( "Chirp evaluations/returns: {}/{}\n", counter_evaluated, counter_returned );
}