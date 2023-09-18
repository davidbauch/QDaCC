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
    Log::L2("[System] Allocated memory for {} states.\n",savedStates.size());
}

// TODO: interpolation für tau direction, dann wieder cachen in t wie mit phononen
Sparse QDLC::Numerics::ODESolver::getHamilton( System &s, const double t ) {
    // Don't use saved Hamiltons, instead just calculate new Operator
    if ( not s.parameters.numerics_use_saved_hamiltons ) {
        track_gethamilton_calc++;
        return s.dgl_get_hamilton( t );
    }
    // Use saved Hamiltons
    // If main direction is not done, calculate and cache Operator if it doesnt exist, then return Operator
    if ( not s.parameters.numerics_main_direction_done ) {
        if ( not savedHamiltons.contains( t ) ) {
            save_hamilton( s.dgl_get_hamilton( t ), t );
            track_gethamilton_write++;
            track_gethamilton_calc++;
        }
        track_gethamilton_read++;
        return savedHamiltons[t];
    }
    // Main direction is done, find operator or return interpolated operator
    if ( savedHamiltons.contains( t ) ) {
        track_gethamilton_read++;
        return savedHamiltons[t];
    }
    // For the edge cases just calculate a new Operator.
    if ( t <= s.parameters.t_start + s.parameters.t_step or t >= s.parameters.t_end - s.parameters.t_step ) {
        track_gethamilton_calc++;
        return s.dgl_get_hamilton( t );
    }
    // Hacky way to find [min,<t>,max].
    const auto &greater_or_equal_than_t = savedHamiltons.lower_bound( t );
    const auto &smaller_than_t = std::prev( greater_or_equal_than_t );
    //  Interpolate
    track_gethamilton_read++;
    return QDLC::Math::lerp( smaller_than_t->second, greater_or_equal_than_t->second, ( t - smaller_than_t->first ) / ( greater_or_equal_than_t->first - smaller_than_t->first ) );
}

void QDLC::Numerics::ODESolver::saveState( const Sparse &mat, const double t, std::vector<QDLC::SaveState> &savedStates ) {
    savedStates.emplace_back( mat, t );
}

void QDLC::Numerics::ODESolver::save_hamilton( const Sparse &mat, const double t ) {
#pragma omp critical
    savedHamiltons[t] = mat;
}

double &QDLC::Numerics::ODESolver::get_time_at( int i ) {
    return savedStates.at( i ).t;
}

Sparse &QDLC::Numerics::ODESolver::get_rho_at( int i ) {
    return savedStates.at( i ).mat;
}

std::tuple<std::string, std::string> QDLC::Numerics::ODESolver::get_operator_strings( System &s, const std::string &operators ) {
    std::string s_creator = "";
    std::string s_annihilator = "";
    for ( auto split_s_op : QDLC::String::splitline( operators, '+' ) ) {
        if ( s_creator.size() > 0 ) {
            s_creator += "+";
            s_annihilator += "+";
        }
        if ( s.operatorMatrices.el_transitions.contains( split_s_op ) ) {
            s_annihilator += split_s_op;
            s_creator += s.operatorMatrices.el_transitions[split_s_op].name_transposed;
        } else {
            s_creator += split_s_op + "bd";
            s_annihilator += split_s_op + "b";
        }
    }
    return std::make_tuple( s_creator, s_annihilator );
}

std::string QDLC::Numerics::ODESolver::get_operators_purpose( const std::vector<std::string> &operators ) {
    const int order = operators.size() / 2;
    std::string ret = "G" + std::to_string( order );
    for ( auto &el : operators )
        ret += "-" + el;
    return ret;
}

Sparse QDLC::Numerics::ODESolver::get_operators_matrix( System &s, const std::string &s_op ) {
    if ( s_op == "internal_identitymatrix" ) {
        return s.operatorMatrices.identity;
    }
    Sparse ret( s.parameters.maxStates, s.parameters.maxStates );
    for ( auto &split_s_op : QDLC::String::splitline( s_op, '+' ) )
        ret += s.operatorMatrices.el_transitions.contains( split_s_op ) ? s.operatorMatrices.el_transitions[split_s_op].hilbert : s.operatorMatrices.ph_transitions[split_s_op].hilbert;
    return ret;
}
std::tuple<Sparse, Sparse> QDLC::Numerics::ODESolver::get_operators_matrices( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator ) {
    const auto op_creator = get_operators_matrix( s, s_op_creator );
    const auto op_annihilator = get_operators_matrix( s, s_op_annihilator );
    return std::make_tuple( op_creator, op_annihilator );
    // return std::make_tuple( op_creator, op_creator.adjoint() );
}