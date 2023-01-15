#include "solver/solver_ode.h"

void QDLC::Numerics::ODESolver::calculate_hamilton_eigenvalues( System &s, const int power ) {
    bool did_title = false;
    const std::string filename = "eigenvalues" + ( power > 1 ? "_" + std::to_string( power ) : "" );
    Log::L2( "[Solver] Calculating Hamilton Eigenvalues using H{}(t), saving to {}.txt\n", power > 1 ? fmt::format( "^{}", power ) : "", filename );
    auto &file = FileOutput::add_file( filename );
    for ( auto &tup : savedStates ) {
        auto &t = tup.t;
        Dense hamilton = Dense( s.dgl_get_hamilton( t ) );
        if ( power > 1 )
            hamilton = hamilton.pow( power ).eval();
        // Eigen::EigenSolver<Dense> eigensolver( hamilton );
        auto eigs = hamilton.eigenvalues() * QDLC::Math::ev_conversion;
        if ( not did_title ) {
            file << "Time";
            for ( int i = 0; i < eigs.size(); i++ )
                file << fmt::format( "\tReal({})", i );
            for ( int i = 0; i < eigs.size(); i++ )
                file << fmt::format( "\tImag({})", i );
            file << "\n";
            did_title = true;
        }
        file << fmt::format( "{:.8e}", t );
        for ( int i = 0; i < eigs.size(); i++ )
            file << fmt::format( "\t{:.10e}", eigs.real()( i ) );
        for ( int i = 0; i < eigs.size(); i++ )
            file << fmt::format( "\t{:.10e}", eigs.imag()( i ) );
        file << "\n";
    }
}

bool QDLC::Numerics::ODESolver::output_numerical_data( System &s ) {
    // Output Numerical RK Error
    if ( s.parameters.output_dict.contains( "rkerror" ) ) { // Chain...
        Log::L2( "[Solver] Outputting Numerical RK45 error...\n" );
        auto &file = FileOutput::add_file( "numerical" );
        // Header
        if ( not s.parameters.numerics_interpolate_outputs )
            file << "t\tRK_Error_Acceped\t";
        file << "t\tRK_Error\n";
        for ( int i = 0; i < rk_error.size(); i++ ) {
            if ( not s.parameters.numerics_interpolate_outputs and i < rk_error_accepted.size() ) {
                auto &[t_t, error] = rk_error_accepted[i];
                file << fmt::format( "{}\t{}\t", t_t, error );
            } else {
                file << "NaN\tNaN\t";
            }
            auto &[t_t, error, t_step, tries] = rk_error[i];
            file << fmt::format( "{}\t{}\t{}\t{}\n", t_t, error, t_step, tries );
        }
    }
    // TODO: putput list als dict, dann if "eigenvalues" in outputdict
    if ( s.parameters.output_dict.contains( "eigenvalues" ) ) {
        QDLC::Numerics::ODESolver::calculate_hamilton_eigenvalues( s );
    }
    if ( s.parameters.output_dict.contains( "eigenvalues2" ) ) {
        QDLC::Numerics::ODESolver::calculate_hamilton_eigenvalues( s, 2 );
    }
    return true;
}