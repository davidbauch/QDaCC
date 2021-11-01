#include "solver/solver_ode.h"

ODESolver::ODESolver( System &s ) {
    Log::L2( "Creating ODESolver Class... " );
    reset( s );
    track_gethamilton_read = 0;
    track_gethamilton_write = 0;
    track_gethamilton_calc = 0;
    track_gethamilton_calcattempt = 0;
    savedStates.clear();
    savedStates.reserve( dim );
    Log::L2( "Done!\n" );
}

Sparse ODESolver::getHamilton( System &s, const double t, bool use_saved_hamiltons ) {
    if ( s.parameters.numerics_use_saved_hamiltons ) {
        if ( savedHamiltons.count( t ) == 0 ) {
#pragma omp critical
            saveHamilton( s.dgl_getHamilton( t ), t );
            track_gethamilton_write++;
            track_gethamilton_calc++;
        }
        track_gethamilton_read++;
        return savedHamiltons[t];
    }
    track_gethamilton_calc++;
    return s.dgl_getHamilton( t );
}

void ODESolver::saveState( const Sparse &mat, const double t, std::vector<SaveState> &savedStates ) {
    savedStates.emplace_back( mat, t );
}

void ODESolver::saveHamilton( const Sparse &mat, const double t ) {
    savedHamiltons[t] = mat;
}

int ODESolver::reset( System &s ) {
    //TODO: remove dim
    dim = (int)std::ceil( ( s.parameters.t_end - s.parameters.t_start ) / ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? s.parameters.t_step_pathint : s.parameters.t_step ) ) + 10; //(int)( s.parameters.iterations_t_max / s.parameters.iterations_t_skip ) + 10;
    return dim;
}

int ODESolver::getIterationNumberTau( System &s ) {
    int num = 0;
    // Tau Direction Iteration steps
    for ( int i = 0; i < (int)savedStates.size(); i += s.parameters.iterations_t_skip ) {
        double t_t = getTimeAt( i );
        for ( double t_tau = t_t + ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? s.parameters.t_step_pathint : s.parameters.t_step ); t_tau < s.parameters.t_end; t_tau += ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? s.parameters.t_step_pathint : s.parameters.t_step ) ) { // t + +s.parameters.t_step
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

// Helperfunctions for g1/2 function string-matrix switch
std::tuple<std::string, std::string> ODESolver::get_operator_strings( const std::string &operators ) {
    std::string s_creator = "";
    std::string s_annihilator = "";
    for ( auto split_s_op : QDLC::String::splitline( operators, '+' ) ) {
        if ( s_creator.size() > 0 ) {
            s_creator += "+";
            s_annihilator += "+";
        }
        if ( std::isupper( split_s_op.front() ) ) {
            s_annihilator += split_s_op;
            std::reverse( split_s_op.begin(), split_s_op.end() );
            s_creator += split_s_op;
        } else {
            s_creator += split_s_op + "bd";
            s_annihilator += split_s_op + "b";
        }
    }
    return std::make_tuple( s_creator, s_annihilator );
}

std::string ODESolver::get_operators_purpose( const std::vector<std::string> &operators, int order ) {
    std::string ret = "G" + std::to_string( order );
    for ( auto &el : operators )
        ret += "-" + el;
    return ret;
}

std::tuple<Sparse, Sparse> ODESolver::get_operators_matrices( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator ) {
    Sparse op_creator = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_annihilator = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    for ( auto &split_s_op_creator : QDLC::String::splitline( s_op_creator, '+' ) )
        op_creator += s.operatorMatrices.el_transitions.count( split_s_op_creator ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_creator].hilbert : s.operatorMatrices.ph_transitions[split_s_op_creator].hilbert;
    for ( auto &split_s_op_annihilator : QDLC::String::splitline( s_op_annihilator, '+' ) )
        op_annihilator += s.operatorMatrices.el_transitions.count( split_s_op_annihilator ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_annihilator].hilbert : s.operatorMatrices.ph_transitions[split_s_op_annihilator].hilbert;
    return std::make_tuple( op_creator, op_annihilator );
    //return std::make_tuple( op_creator, op_creator.adjoint() );
}