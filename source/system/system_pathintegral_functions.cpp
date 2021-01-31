#include "system/system.h"

Scalar System::dgl_phonons_kernel( const double t, const double t_step ) {
    Scalar integral = 0;
    double stepsize = 1.0E-3 * parameters.p_phonon_wcutoff;
    double eV7 = convertParam<double>( "7.0eV" );
    double eV35 = -convertParam<double>( "3.5eV" );
    double v_c = 5110.0;
    double a_e = 5E-9;       //3E-9;
    double a_h = 0.87 * a_e; //a_e / 1.15;
    double rho = 5370.0;
    for ( double w = stepsize; w < 10.0 * parameters.p_phonon_wcutoff; w += stepsize ) {
        double J = w * parameters.p_phonon_alpha * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff );
        //double J = w * parameters.hbar * std::pow( eV7 * std::exp( -w * w * a_e * a_e / ( 4. * v_c * v_c ) ) - eV35 * std::exp( -w * w * a_h * a_h / ( 4. * v_c * v_c ) ), 2. ) / ( 4. * 3.1415 * 3.1415 * rho * std::pow( v_c, 5. ) );
        if ( t != 0.0 ) {
            integral += stepsize * 2.0 * J * ( 1.0 - std::cos( w * t_step ) ) * ( std::cos( w * t ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) - 1.i * std::sin( w * t ) );
        } else {
            integral += stepsize * J * ( ( 1.0 - std::cos( w * t_step ) ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) + 1.i * std::sin( w * t_step ) ); //- 1.i * w * t_step );
        }
    }
    return integral;
}

Scalar System::dgl_phonon_S_function( const int t_delta, const int i_n, const int j_n, const int i_nd, const int j_nd ) {
    return -phi_vector[t_delta] * operatorMatrices.phononCouplingFactor( i_n, i_nd ) - std::conj( phi_vector[t_delta] ) * operatorMatrices.phononCouplingFactor( j_n, j_nd ) + std::conj( phi_vector[t_delta] ) * operatorMatrices.phononCouplingFactor( i_n, j_nd ) + phi_vector[t_delta] * operatorMatrices.phononCouplingFactor( j_n, i_nd );
    // If i or j = Groundstate, return 0
    //Scalar result = 0;
    //if ( i_n == i_nd ) {
    //    result -= phi_vector[t_delta] * operatorMatrices.phononCouplingFactor[i_n] * operatorMatrices.phononCouplingFactor[i_nd];
    //}
    //if ( j_n == j_nd ) {
    //    result -= std::conj( phi_vector[t_delta] ) * operatorMatrices.phononCouplingFactor[j_n] * operatorMatrices.phononCouplingFactor[j_nd];
    //}
    //if ( i_n == j_nd ) {
    //    result += std::conj( phi_vector[t_delta] ) * operatorMatrices.phononCouplingFactor[i_n] * operatorMatrices.phononCouplingFactor[j_nd];
    //}
    //if ( i_nd == j_n ) {
    //    result += phi_vector[t_delta] * operatorMatrices.phononCouplingFactor[j_n] * operatorMatrices.phononCouplingFactor[i_nd];
    //}
    //return result;
}

void System::initialize_path_integral_functions() {
    Log::L2( "Initializing Path-Integral functions...\n" );
    // kernel in phi-vector schreiben
    int tau_max = parameters.p_phonon_nc + 1;
    phi_vector.reserve( tau_max );
    Log::L2( "Initializing Kernel Memory functions...\n" );
    for ( double tau = 0.0; tau < parameters.t_step * tau_max; tau += parameters.t_step ) {
        phi_vector.emplace_back( dgl_phonons_kernel( tau, parameters.t_step ) );
    }

    Log::L2( "Outputting phonon functions to phonons.txt from phi_vector({})...\n", phi_vector.size() );
    // Output Phonon Functions
    FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
    fmt::print( fp_phonons, "t\tabs(K(t))\treal(K(t))\timag(K(t))\treal(S_ij)\timag(S_ij)\n" );
    for ( double t = parameters.t_start; t < parameters.t_step * tau_max; t += parameters.t_step ) {
        int i = std::floor( t / ( parameters.t_step ) );
        fmt::print( fp_phonons, "{}\t{}\t{}\t{}\n", t, std::abs( phi_vector.at( i ) ), std::real( phi_vector.at( i ) ), std::imag( phi_vector.at( i ) ) );
    }
    // Lets output more than the 4-8 elements usually used
    //std::vector<Scalar> phi_vector_o;
    //for ( double tau = 0.0; tau <= parameters.t_step * tau_max; tau += parameters.t_step / 10.0 ) {
    //    phi_vector_o.emplace_back( dgl_phonons_kernel( tau, parameters.t_step / 10.0 ) );
    //}
    //// Output Phonon Functions
    //FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
    //fmt::print( fp_phonons, "t\tabs(K(t))\treal(K(t))\timag(K(t))\treal(S_ij)\timag(S_ij)\n" );
    //for ( double t = parameters.t_start; t < parameters.p_phonon_tcutoff; t += parameters.t_step / 10.0 ) {
    //    int i = std::floor( t / ( parameters.t_step / 10.0 ) );
    //    fmt::print( fp_phonons, "{}\t{}\t{}\t{}\n", t, std::abs( phi_vector_o.at( i ) ), std::real( phi_vector_o.at( i ) ), std::imag( phi_vector_o.at( i ) ) );
    //}
    std::fclose( fp_phonons );
    Log::L2( "Done...\n" );
}