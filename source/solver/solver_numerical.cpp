#include "solver/solver_ode.h"

void QDLC::Numerics::ODESolver::calculate_hamilton_eigenvalues( System &s ) {
    auto output_format = Eigen::IOFormat( -1, 0, "\t", ", ", "", "" );
    for ( auto &tup : savedStates ) {
        auto &t = tup.t;
        auto hamilton = s.dgl_get_hamilton( t );
        // Eigen::EigenSolver<Dense> eigensolver( hamilton );
        auto eigs = Dense( hamilton ).eigenvalues().real() * QDLC::Math::ev_conversion;
        fmt::print( s.fileoutput.fp_eigenvalues, "{}\t{}\n", t, eigs.format( output_format ) );
    }
}

bool QDLC::Numerics::ODESolver::output_numerical_data( System &s ) {
    // Output Numerical RK Error
    if ( s.parameters.numerics_output_rkerror ) { // Chain...
        Log::L2( "[Solver] Outputting Numerical RK45 error...\n" );
        // Header
        if ( not s.parameters.numerics_interpolate_outputs )
            fmt::print( s.fileoutput.fp_numerical, "t\tRK_Error_Acceped\t" );
        fmt::print( s.fileoutput.fp_numerical, "t\tRK_Error\n" );
        for ( int i = 0; i < rk_error.size(); i++ ) {
            if ( not s.parameters.numerics_interpolate_outputs and i < rk_error_accepted.size() ) {
                auto &[t_t, error] = rk_error_accepted[i];
                fmt::print( s.fileoutput.fp_numerical, "{}\t{}\t", t_t, error );
            } else {
                fmt::print( s.fileoutput.fp_numerical, "NaN\tNaN\t" );
            }
            auto &[t_t, error, t_step, tries] = rk_error[i];
            fmt::print( s.fileoutput.fp_numerical, "{}\t{}\t{}\t{}\n", t_t, error, t_step, tries );
        }
    }
    if ( s.parameters.output_eigenvalues ) {
        QDLC::Numerics::ODESolver::calculate_hamilton_eigenvalues( s );
    }
    return true;
}