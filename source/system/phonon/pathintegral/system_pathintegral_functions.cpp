#include "system/system.h"

using namespace QDACC;

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
                parameters.t_step_pathint = tau / ( 1.0 * tau_max );
                Log::L2( "[System-Path-Integral] Path Integral t-cutoff was automatically determined to t_cutoff = {}, resulting in a pathintegral timestep of {}\n", tau, parameters.t_step_pathint );
                break;
            }
            last = current;
            tau += parameters.numerics_subiterator_stepsize;
        }
        parameters.post_adjust_grids();
    }

    // Precalculate Phi
    for ( double tau = 0.0; tau < parameters.t_step_pathint * tau_max; tau += parameters.t_step_pathint ) {
        Scalar kernel = dgl_phonons_kernel( tau, parameters.t_step_pathint );
        phi_vector[tau] = kernel;
        phi_vector_int.emplace_back( kernel );
    }

    Log::L2( "[PathIntegral] Outputting phonon functions to phonons.txt from phi_vector({})...\n", phi_vector.size() );
    // Output Phonon Functions
    if ( parameters.output_dict.contains( "PIkernel" ) or parameters.output_dict.contains( "phononcoefficients" ) ) {
        auto &file = FileOutput::add_file( "phonon_PI" );
        file << std::format( "t\tabs(K(t))\treal(K(t))\timag(K(t))\tabs(K(t))\treal(K(t))\timag(K(t))\n" );
        for ( int i = 0; i < phi_vector_int.size(); i++ ) {
            file << std::format( "{}\t{}\t{}\t{}\n", parameters.t_step_pathint * i, std::abs( phi_vector_int[i] ), std::real( phi_vector_int[i] ), std::imag( phi_vector_int[i] ) );
        }
        file.close();
    }
    // Lets output more than the 4-8 elements usually used
    // std::vector<Scalar> phi_vector_o;
    // for ( double tau = 0.0; tau <= parameters.t_step * tau_max; tau += parameters.t_step / 10.0 ) {
    //    phi_vector_o.emplace_back( dgl_phonons_kernel( tau, parameters.t_step / 10.0 ) );
    //}
    //// Output Phonon Functions
    // FILE *fp_phonons = std::fopen( ( parameters.working_directory + "phonons.txt" ).c_str(), "w" );
    // fmt::print( fp_phonons, "t\tabs(K(t))\treal(K(t))\timag(K(t))\treal(S_ij)\timag(S_ij)\n" );
    // for ( double t = parameters.t_start; t < parameters.p_phonon_tcutoff; t += parameters.t_step / 10.0 ) {
    //     int i = std::floor( t / ( parameters.t_step / 10.0 ) );
    //     fmt::print( fp_phonons, "{}\t{}\t{}\t{}\n", t, std::abs( phi_vector_o.at( i ) ), std::real( phi_vector_o.at( i ) ), std::imag( phi_vector_o.at( i ) ) );
    // }
    Log::L2( "[PathIntegral] Done...\n" );
}