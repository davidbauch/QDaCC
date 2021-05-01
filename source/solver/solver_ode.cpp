#include "solver/solver_ode.h"

ODESolver::ODESolver( System &s ) {
    Log::L2( "Creating ODESolver Class... " );
    reset( s );
    savedStates.clear();
    savedStates.reserve( dim );
    Log::L2( "Done!\n" );
}

Sparse ODESolver::getHamilton( System &s, const double t, bool use_saved_hamiltons ) {
    if ( s.parameters.numerics_use_saved_hamiltons ) { ////|| use_saved_hamiltons ) {
        if ( savedHamiltons.size() == 0 || t > savedHamiltons.back().t ) {
            track_gethamilton_calcattempt++;
            // Calculate new Hamilton for t, save it, return it
            saveHamilton( s.dgl_getHamilton( t ), t );
            track_gethamilton_write++;
            return savedHamiltons.back().mat;
        } else {
            long unsigned int approx = std::floor( t / s.parameters.t_step / 2.0 );
            if ( !( approx >= savedHamiltons.size() || approx < 0 ) ) { //Fixed: > --> >=
                // Look for corresponding Hamilton. t-direction has to be calculated first if used with multiple threads
                while ( t > savedHamiltons.at( approx ).t )
                    approx++;
                while ( t < savedHamiltons.at( approx ).t )
                    approx--;
                if ( approx < (int)savedHamiltons.size() ) {
                    track_gethamilton_read++;
                    return savedHamiltons.at( approx ).mat;
                }
            }
        }
    }
    track_gethamilton_calc++;
    return s.dgl_getHamilton( t );
}

void ODESolver::saveState( const Sparse &mat, const double t, std::vector<SaveState> &savedStates ) {
    savedStates.emplace_back( mat, t );
}

void ODESolver::saveHamilton( const Sparse &mat, const double t ) {
    savedHamiltons.emplace_back( mat, t );
}

bool ODESolver::queueNow( System &s, int &curIt ) {
    bool queue = s.parameters.numerics_calculate_spectrum_H || s.parameters.numerics_calculate_spectrum_V || s.parameters.numerics_calculate_g2;
    if ( queue && curIt % s.parameters.iterations_t_skip == 0 ) {
        curIt = 1;
        return true;
    }
    curIt++;
    return false;
}

int ODESolver::reset( System &s ) {
    track_gethamilton_read = 0;
    track_gethamilton_write = 0;
    track_gethamilton_calc = 0;
    track_gethamilton_calcattempt = 0;
    dim = (int)std::ceil( ( s.parameters.t_end - s.parameters.t_start ) / s.parameters.t_step ) + 10; //(int)( s.parameters.iterations_t_max / s.parameters.iterations_t_skip ) + 10;
    return dim;
}

int ODESolver::getIterationNumberTau( System &s ) {
    int num = 0;
    // Tau Direction Iteration steps
    for ( int i = 0; i < (int)savedStates.size(); i += s.parameters.iterations_t_skip ) {
        double t_t = getTimeAt( i );
        for ( double t_tau = t_t + s.parameters.t_step; t_tau < s.parameters.t_end; t_tau += s.parameters.t_step ) { // t + +s.parameters.t_step
            num++;
        }
    }
    return num;
}

// Description: Function to calculate the number of iterations used for spectru calculations
// Type: ODESolver private function
// @param s: [&System] Class providing set of system functions
// @return: [int] Number of spectrum iterations

int ODESolver::getIterationNumberSpectrum( System &s ) {
    int num = 0;
    // Spectrum steps
    for ( int spec_w = 0; spec_w < s.parameters.iterations_w_resolution; spec_w++ ) {
        num++;
    }
    return num;
}

// Description: Function to extract the time corresponding to a certain iteration of saved states
// Type: ODESolver private function
// @param i: [int] Iteration number
// @return: [double] Time corresponding to iteration number i

double ODESolver::getTimeAt( int i ) {
    return savedStates.at( i ).t;
}

// Description: Function to extract the matrix corresponding to a certain iteration of saved states
// Type: ODESolver private function
// @param i: [int] Iteration number
// @return: [Sparse] (Density-) Matrix corresponding to iteration number i

Sparse ODESolver::getRhoAt( int i ) {
    return savedStates.at( i ).mat;
}