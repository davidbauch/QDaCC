#include "solver/solver_ode.h"

/**
 * These functions calculate the n-Dimensional G^(n) correlation function.
 *
 * The n-th order correlation function is defined as:
 * `G^(n)(t1,...,tn) = < a_1^d(t1) * a_2^d(t2) * ... * a_n^d(tn) * a_n+1(tn) * a_n+2(tn-1) * ... * a_2n(t) >`
 * `= Tr( rho * a_1^d(t1) * a_2^d(t2) * ... * a_n^d(tn) * a_n+1(tn) * a_n+2(tn-1) * ... * a_2n(t) )`
 *
 * The most used correlation functions are G^(1), G^(2) and G^(3). Using the unitary
 * time transformation `U(t0,t1) = exp(-i int_t0^t1 H(t) dt)` and the Quantum
 * Regression Theorem, we can rewrite the correlation functions as:
 *
 * `G^(1)(t,tau) = < a_1^d(t) * a_2(t+tau) >`
 * `             = Tr (rho * a_1^d(t) * a_2(t+tau) )`
 * `             = Tr( [a_2 * rho(t)]->tau * a_1^d )`
 * `G^(2)(t,tau) = < a_1^d(t) * a_2^d(t+tau) * a_3(t+tau) * a_4(t) >`
 * `             = Tr( rho * a_1^d(t) * a_2^d(t+tau) * a_3(t+tau) * a_4(t))`
 * `             = Tr( [a_4 * rho(t) * a_1^d]->tau * a_2^d * a_3 )`
 * `G^(3)(t,tau1,tau2) = < a_1^d(t) * a_2^d(t+tau1) * a_3^d(t+tau2) * a_4(t+tau2) * a_5(t+tau1) * a_6(t) >`
 * `                   = Tr( rho * a_1^d(t) * a_2^d(t+tau1) * a_3^d(t+tau2) * a_4(t+tau2) * a_5(t+tau1) * a_6(t) )`
 * `                   = Tr ( [ [ a_6 * rho(t) * a_1^d ]->tau1 a_2^d * a_5 ]->tau2 * a_3^d * a_4 ) )`
 * ...
 * `G^(n)(t,tau1,...,taun-1) = < a_1^d(t) * a_2^d(t+tau1) * ... * a_n^d(t+tau_n-1) * a_n+1(t+tau_n-1) * ... * a_2n(t) >`
 * `                         = Tr( rho * a_1^d(t) * a_2^d(t+tau1) * ... * a_n^d(t+tau_n-1) * a_n+1(t+tau_n-1) * ... * a_2n(t) )`
 * `                         = Tr( [ [ [ [ [ a_2n * rho(t) * a_1^d ]->tau1 * a_2^d * a_2n-1 ]->tau2 ] * a_3^d * a_2n-2 ]->tau3 ... * a_n-1^d a_n+2 ]->taun-1 * a_n^d * a_n+1 )`
 *
 * Effectively, we are chaining regular G(2) correlation calculations together.
 * With the MultidimensionalCacheMatrix, we can, in theory, calculate any G^(n)
 * as long as we have enough memory to save the resulting N^n matrices.
 */

void QDACC::Numerics::ODESolver::calculate_g1( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[CorrelationFunction] G1(tau) for {} already exists.\n", purpose );
        return;
    }
    calculate_g1( s, std::vector<std::string>{ s_op_i }, s_op_j, { purpose } );
}

void QDACC::Numerics::ODESolver::calculate_g1( System &s, const std::vector<std::string> &s_op_i, const std::string &s_op_j, const std::vector<std::string> &purposes ) {
    const auto size = s_op_i.size();
    calculate_g2( s, "internal_identitymatrix", std::vector<std::string>{ size, "internal_identitymatrix" }, s_op_i, s_op_j, purposes );
}

void QDACC::Numerics::ODESolver::calculate_g2( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[CorrelationFunction] G2(tau) for {} already exists.\n", purpose );
        return;
    }
    calculate_g2( s, s_op_i, std::vector<std::string>{ s_op_j }, { s_op_k }, s_op_l, { purpose } );
}

void QDACC::Numerics::ODESolver::calculate_g3( System &s, const std::string &s_op_i, const std::string &s_op_j, const std::string &s_op_k, const std::string &s_op_l, const std::string &s_op_m, const std::string &s_op_n, const std::string &purpose ) {
    if ( cache.contains( purpose ) ) {
        Log::L2( "[CorrelationFunction] G3(tau) for {} already exists.\n", purpose );
        return;
    }
    calculate_g3( s, s_op_i, std::vector<std::string>{ s_op_j }, { s_op_k }, {s_op_l}, {s_op_m}, s_op_n, { purpose } );
}