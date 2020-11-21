#include "system/system.h"

Scalar System::dgl_phonons_kernel( const double t ) {
    Scalar integral = 0;
    double stepsize = 0.01 * parameters.p_phonon_wcutoff;
    for ( double w = stepsize; w < 10 * parameters.p_phonon_wcutoff; w += stepsize ) {
        if ( t == 0.0 ) {
            integral += stepsize * ( parameters.p_phonon_alpha * w * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff ) * ( 1.0 / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) * ( ( 1 - std::cos( w * parameters.t_step ) ) + 1i * std::sin( w * parameters.t_step ) - 1i * w * parameters.t_step ) ) );
        } else {
            integral += stepsize * 2.0 * ( parameters.p_phonon_alpha * w * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff ) * ( 1 - std::cos( w * parameters.t_step ) ) * ( 1.0 / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) * ( std::cos( w * t ) - 1i * std::sin( w * t ) ) ) );
        }
    }
    return integral;
}

Scalar System::dgl_phonon_S_function( const int t_delta, const int i_n, const int j_n, const int i_nd, const int j_nd ) {
    Scalar result = 0;
    // If i or j = Groundstate, return 0
    if ( i_n != 0 && i_n == i_nd ) {
        result -= phi_vector.at( t_delta );
    }
    if ( j_n != 0 && j_n == j_nd ) {
        result -= std::conj( phi_vector.at( t_delta ) );
    }
    if ( i_n != 0 && i_n == j_nd ) {
        result += std::conj( phi_vector.at( t_delta ) );
    }
    if ( j_n != 0 && i_nd == j_n ) {
        result += phi_vector.at( t_delta );
    }
    return result;
}

void System::initialize_path_integral_functions() {
    Log::L2( "Initializing Path-Integral functions...\n" );
    // kernel in phi-vector schreiben
    int tau_max = 2 * std::ceil( parameters.p_phonon_tcutoff / parameters.t_step );
    tau_max = 25 > tau_max ? 25 : tau_max;
    std::vector<Scalar> SS( tau_max );
    phi_vector.reserve( tau_max );
    Log::L2( "Initializing Kernel Memory functions...\n" );
    for ( double tau = 0.0; tau <= parameters.t_step * tau_max; tau += parameters.t_step ) {
        phi_vector.emplace_back( dgl_phonons_kernel( tau ) );
    }
    Log::L2( "Initializing exp(S) function...\n" );
    // S vektor schreiben
    //Sparse S = dgl_phonon_S_function( 0 );
    //phonon_S.emplace_back( S );
    SS.emplace_back( phi_vector.at( 0 ) ); // + std::conj( phi_vector.at( 0 ) ) );
    for ( int l = 0; l < tau_max; l++ ) {
        //Sparse S = phonon_S.back();
        Scalar cSS = SS.back();
        for ( int ld = 0; ld < l; ld++ ) {
            //S += dgl_phonon_S_function( ld );
            cSS += phi_vector.at( ld ); // + std::conj( phi_vector.at( ld ) );
        }
        //phonon_S.emplace_back( S );
        SS.emplace_back( cSS );
    }
    Log::L2( "Outputting phonon functions to phonons.txt from phi_vector({})...\n", phi_vector.size() );
    // Output Phonon Functions
    FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
    fmt::print( fp_phonons, "t\treal(K(t))\timag(K(t))\treal(S_ij)\timag(S_ij)\n" );
    for ( double t = parameters.t_start; t < parameters.p_phonon_tcutoff; t += parameters.t_step ) {
        int i = std::floor( t / parameters.t_step );
        fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\n", t, std::real( phi_vector.at( i ) ), std::imag( phi_vector.at( i ) ), std::real( SS.at( i ) ), std::imag( SS.at( i ) ) );
    }
    std::fclose( fp_phonons );
    Log::L2( "Done...\n" );
}