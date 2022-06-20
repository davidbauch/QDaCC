#include "solver/solver_ode.h"

void QDLC::Numerics::ODESolver::calculate_hamilton_eigenvalues( System &s ) {
    bool did_title = false;
    for ( auto &tup : savedStates ) {
        auto &t = tup.t;
        auto hamilton = s.dgl_get_hamilton( t );
        // Eigen::EigenSolver<Dense> eigensolver( hamilton );
        auto eigs = Dense( hamilton ).eigenvalues() * QDLC::Math::ev_conversion;
        if ( not did_title ) {
            fmt::print( s.fileoutput.fp_eigenvalues, "Time" );
            for ( int i = 0; i < eigs.size(); i++ )
                fmt::print( s.fileoutput.fp_eigenvalues, "\tReal({})", i );
            for ( int i = 0; i < eigs.size(); i++ )
                fmt::print( s.fileoutput.fp_eigenvalues, "\tImag({})", i );
            fmt::print( s.fileoutput.fp_eigenvalues, "\n" );
            did_title = true;
        }
        fmt::print( s.fileoutput.fp_eigenvalues, "{:.8e}", t );
        for ( int i = 0; i < eigs.size(); i++ )
            fmt::print( s.fileoutput.fp_eigenvalues, "\t{:.10e}", eigs.real()( i ) );
        for ( int i = 0; i < eigs.size(); i++ )
            fmt::print( s.fileoutput.fp_eigenvalues, "\t{:.10e}", eigs.imag()( i ) );
        fmt::print( s.fileoutput.fp_eigenvalues, "\n" );
    }
}

bool QDLC::Numerics::ODESolver::output_numerical_data( System &s ) {
    // Output Numerical RK Error
    if ( s.parameters.numerics_output_rkerror ) { // Chain...
        LOG2( "[Solver] Outputting Numerical RK45 error...\n" );
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