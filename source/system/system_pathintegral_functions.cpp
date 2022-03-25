#include "system/system.h"

Scalar System::dgl_phonons_kernel( const double t, const double t_step ) {
    Scalar integral = 0;
    double stepsize = 1.0E-4 * parameters.p_phonon_wcutoff;
    for ( double w = 1E-10; w < 10.0 * parameters.p_phonon_wcutoff; w += stepsize ) {
        double J;
        if ( parameters.p_phonon_qd_ae == 0.0 )
            J = parameters.p_phonon_alpha * w * std::exp( -w * w / 2.0 / parameters.p_phonon_wcutoff / parameters.p_phonon_wcutoff );
        else
            J = w * parameters.hbar * std::pow( parameters.p_phonon_qd_de * std::exp( -w * w * parameters.p_phonon_qd_ae * parameters.p_phonon_qd_ae / ( 4. * parameters.p_phonon_qd_cs * parameters.p_phonon_qd_cs ) ) - parameters.p_phonon_qd_dh * std::exp( -w * w * parameters.p_phonon_qd_ae / parameters.p_phonon_qd_ratio * parameters.p_phonon_qd_ae / parameters.p_phonon_qd_ratio / ( 4. * parameters.p_phonon_qd_cs * parameters.p_phonon_qd_cs ) ), 2. ) / ( 4. * 3.1415 * 3.1415 * parameters.p_phonon_qd_rho * std::pow( parameters.p_phonon_qd_cs, 5. ) );
        // double J = w * parameters.hbar * std::pow( eV7 * std::exp( -w * w * a_e * a_e / ( 4. * v_c * v_c ) ) - eV35 * std::exp( -w * w * a_h * a_h / ( 4. * v_c * v_c ) ), 2. ) / ( 4. * 3.1415 * 3.1415 * rho * std::pow( v_c, 5. ) );
        if ( t < t_step / 2.0 ) {
            integral += stepsize * J * ( ( 1.0 - std::cos( w * t_step ) ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) + 1.0i * std::sin( w * t_step ) ); // - 1.i * w * t_step );
        } else {
            integral += stepsize * 2.0 * J * ( 1.0 - std::cos( w * t_step ) ) * ( std::cos( w * t ) / std::tanh( parameters.hbar * w / 2.0 / parameters.kb / parameters.p_phonon_T ) - 1.0i * std::sin( w * t ) );
        }
    }
    return integral;
}

Scalar System::dgl_phonon_S_function( const int t_delta, const int i_n, const int j_n, const int i_nd, const int j_nd ) {
    Scalar result = 0.0;
    if ( i_n == i_nd )
        result -= phi_vector_int[t_delta] * operatorMatrices.phonon_coupling_index_value[i_n];
    if ( j_n == j_nd )
        result -= std::conj( phi_vector_int[t_delta] ) * operatorMatrices.phonon_coupling_index_value[j_n];
    if ( i_n == j_nd )
        result += std::conj( phi_vector_int[t_delta] ) * operatorMatrices.phonon_coupling_index_value[i_n];
    if ( j_n == i_nd )
        result += phi_vector_int[t_delta] * operatorMatrices.phonon_coupling_index_value[j_n];
    return result;
}

void System::initialize_path_integral_functions() {
    Log::L2( "[System-Path-Integral] Initializing Path-Integral functions...\n" );
    // kernel in phi-vector schreiben
    int tau_max = parameters.p_phonon_nc + 1;
    Log::L2( "[System-Path-Integral] Initializing Kernel Memory functions...\n" );
    phi_vector_int.clear();

    // If required, determine tau_max and t_step_pathint automatically
    if ( parameters.t_step_pathint < 0 ) {
        double tau = 0.0;
        double last = 1.0;
        double first = std::abs( dgl_phonons_kernel( 0.0, parameters.numerics_subiterator_stepsize ) );
        while ( true ) {
            double current = std::abs( dgl_phonons_kernel( tau, parameters.numerics_subiterator_stepsize ) );
            if ( std::abs( 1.0 - current / last ) < 1E-2 or ( last != 1.0 and std::abs( current / first ) < 1E-3 ) ) {
                parameters.t_step_pathint = tau / ( 1.0 * parameters.p_phonon_nc );
                Log::L2( "[System-Path-Integral] Path Integral t-cutoff was automatically determined to t_cutoff = {}, resulting in a pathintegral timestep of {}\n", tau, parameters.t_step_pathint );
                break;
            }
            last = current;
            tau += parameters.numerics_subiterator_stepsize;
        }
        parameters.adjustInput();
    }

    // Precalculate Phi
    for ( double tau = 0.0; tau < parameters.t_step_pathint * tau_max; tau += parameters.t_step_pathint ) {
        Scalar kernel = dgl_phonons_kernel( tau, parameters.t_step_pathint );
        phi_vector[tau] = kernel;
        phi_vector_int.emplace_back( kernel );
    }

    Log::L2( "[PathIntegral] Outputting phonon functions to phonons.txt from phi_vector({})...\n", phi_vector.size() );
    // Output Phonon Functions
    FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
    fmt::print( fp_phonons, "t\tabs(K(t))\treal(K(t))\timag(K(t))\tabs(K(t))\treal(K(t))\timag(K(t))\n" );
    // for ( double t = 0.0; t < parameters.t_step_pathint * tau_max; t += parameters.t_step_pathint ) {
    //     int i = std::floor( t / ( parameters.t_step_pathint ) );
    //     auto ampfactor = parameters.t_step_pathint / parameters.t_step;
    //     auto phi = dgl_phonons_kernel( t, parameters.t_step ) * ampfactor;
    //     fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", t, std::abs( phi ), std::real( phi ), std::imag( phi ), std::abs( phi_vector[t] ), std::real( phi_vector[t] ), std::imag( phi_vector[t] ) );
    //     for ( double dt = t + parameters.t_step; dt < t + parameters.t_step_pathint; dt += parameters.t_step ) {
    //         auto phi = dgl_phonons_kernel( dt, parameters.t_step ) * ampfactor;
    //         fmt::print( fp_phonons, "{}\t{}\t{}\t{}\t \t \t \n", dt, std::abs( phi ), std::real( phi ), std::imag( phi ) );
    //     }
    // }
    for ( int i = 0; i < phi_vector_int.size(); i++ ) {
        fmt::print( fp_phonons, "{}\t{}\t{}\t{}\n", parameters.t_step_pathint * i, std::abs( phi_vector_int[i] ), std::real( phi_vector_int[i] ), std::imag( phi_vector_int[i] ) );
    }
    // Lets output more than the 4-8 elements usually used
    // std::vector<Scalar> phi_vector_o;
    // for ( double tau = 0.0; tau <= parameters.t_step * tau_max; tau += parameters.t_step / 10.0 ) {
    //    phi_vector_o.emplace_back( dgl_phonons_kernel( tau, parameters.t_step / 10.0 ) );
    //}
    //// Output Phonon Functions
    // FILE *fp_phonons = std::fopen( ( parameters.subfolder + "phonons.txt" ).c_str(), "w" );
    // fmt::print( fp_phonons, "t\tabs(K(t))\treal(K(t))\timag(K(t))\treal(S_ij)\timag(S_ij)\n" );
    // for ( double t = parameters.t_start; t < parameters.p_phonon_tcutoff; t += parameters.t_step / 10.0 ) {
    //     int i = std::floor( t / ( parameters.t_step / 10.0 ) );
    //     fmt::print( fp_phonons, "{}\t{}\t{}\t{}\n", t, std::abs( phi_vector_o.at( i ) ), std::real( phi_vector_o.at( i ) ), std::imag( phi_vector_o.at( i ) ) );
    // }
    std::fclose( fp_phonons );
    Log::L2( "[PathIntegral] Done...\n" );
}