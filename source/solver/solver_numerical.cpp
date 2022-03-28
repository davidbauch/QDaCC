#include "solver/solver_ode.h"

/**
 * @brief Outputs numerical data to file
 *
 * @param s System
 * @return true
 */
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
    Log::L2( "[Solver] Done!\n" );
    return true;
}