#include "solver/solver_ode.h"

QDLC::Numerics::ODESolver::ODESolver( System &s ) {
    Log::L2( "[System] Creating ODESolver Class...\n" );
    track_gethamilton_read = 0;
    track_gethamilton_write = 0;
    track_gethamilton_calc = 0;
    track_gethamilton_calcattempt = 0;
    auto size = savedStates.size();
    savedStates.clear();
    savedStates.reserve( size );
}

// TODO: interpolation für tau direction, dann wieder cachen in t wie mit phononen
Sparse QDLC::Numerics::ODESolver::getHamilton( System &s, const double t, bool use_saved_hamiltons ) {
    if ( s.parameters.numerics_use_saved_hamiltons ) {
        if ( savedHamiltons.count( t ) == 0 ) {
#pragma omp critical
            save_hamilton( s.dgl_get_hamilton( t ), t );
            track_gethamilton_write++;
            track_gethamilton_calc++;
        }
        track_gethamilton_read++;
        return savedHamiltons[t];
    }
    track_gethamilton_calc++;
    return s.dgl_get_hamilton( t );
}

void QDLC::Numerics::ODESolver::saveState( const Sparse &mat, const double t, std::vector<QDLC::SaveState> &savedStates ) {
    savedStates.emplace_back( mat, t );
}

void QDLC::Numerics::ODESolver::save_hamilton( const Sparse &mat, const double t ) {
    savedHamiltons[t] = mat;
}

double &QDLC::Numerics::ODESolver::get_time_at( int i ) {
    return savedStates.at( i ).t;
}

Sparse &QDLC::Numerics::ODESolver::get_rho_at( int i ) {
    return savedStates.at( i ).mat;
}

std::tuple<std::string, std::string> QDLC::Numerics::ODESolver::get_operator_strings( const std::string &operators ) {
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

std::string QDLC::Numerics::ODESolver::get_operators_purpose( const std::vector<std::string> &operators, int order ) {
    std::string ret = "G" + std::to_string( order );
    for ( auto &el : operators )
        ret += "-" + el;
    return ret;
}

std::tuple<Sparse, Sparse> QDLC::Numerics::ODESolver::get_operators_matrices( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator ) {
    Sparse op_creator = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_annihilator = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    for ( auto &split_s_op_creator : QDLC::String::splitline( s_op_creator, '+' ) )
        op_creator += s.operatorMatrices.el_transitions.count( split_s_op_creator ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_creator].hilbert : s.operatorMatrices.ph_transitions[split_s_op_creator].hilbert;
    for ( auto &split_s_op_annihilator : QDLC::String::splitline( s_op_annihilator, '+' ) )
        op_annihilator += s.operatorMatrices.el_transitions.count( split_s_op_annihilator ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_annihilator].hilbert : s.operatorMatrices.ph_transitions[split_s_op_annihilator].hilbert;
    return std::make_tuple( op_creator, op_annihilator );
    // return std::make_tuple( op_creator, op_creator.adjoint() );
}