#include "solver/solver_ode.h"

bool QDACC::Numerics::ODESolver::calculate_indistinguishability( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator ) {
    // Set Number of Phonon cores to 1 because this memberfunction is already using multithreading
    s.parameters.numerics_maximum_secondary_threads = 1;
    // Progress
    ProgressBar progressbar = ProgressBar();

    // Calculate G2(t,tau) with given operator matrices
    std::string s_g1 = get_operators_purpose( { s_op_creator, s_op_annihilator } );
    std::string s_g2 = get_operators_purpose( { s_op_creator, s_op_creator, s_op_annihilator, s_op_annihilator } );

    // Get Sparse Operator Matrices
    const auto &[op_creator, op_annihilator] = get_operators_matrices( s, s_op_creator, s_op_annihilator );

    // Calculate G-Functions if neccessary
    calculate_g1( s, s_op_creator, s_op_annihilator, s_g1 );
    calculate_g2( s, s_op_creator, s_op_creator, s_op_annihilator, s_op_annihilator, s_g2 );

    auto &akf_mat_g1 = cache[s_g1];
    auto &akf_mat_g2 = cache[s_g2];

    int pbsize = akf_mat_g1.dim();

    Sparse M1 = op_creator * op_annihilator;

    std::string fout = s_op_creator + "-" + s_op_annihilator;
    Timer &timer = Timers::create( "Indistinguishability " + fout ).start();

    size_t T = std::min<size_t>( pbsize, savedStates.size() );

    std::vector<Scalar> top( T, 0 ), bottom( T, 0 ), topv( T, 0 ), outp( T, 0 ), outpv( T, 0 );
    std::vector<Scalar> time;
    for ( int i = 0; i < T; i++ ) // This could also be a row extraction using Eigen
        time.emplace_back( akf_mat_g1.t( i ) );

#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( int upper_limit = 0; upper_limit < T; upper_limit++ ) {
        for ( int i = 0; i <= upper_limit; i++ ) {
            double t_t = akf_mat_g1.t( i );
            auto rho = get_rho_at( rho_index_map[t_t] );
            int j = upper_limit - i;
            double t_tau = akf_mat_g1.t( i + j );
            auto rho_tau = get_rho_at( rho_index_map[t_tau] );
            Scalar gpop = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M1, t_tau );
            Scalar gbot = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator, t_t );
            double dt = akf_mat_g1.dt( i );
            double dtau = akf_mat_g1.dtau( j, i );
            top[upper_limit] += ( gpop + akf_mat_g2.get( i, j ) - akf_mat_g1.get( i, j ) * std::conj( akf_mat_g1.get( i, j ) ) ) * dt * dtau;
            bottom[upper_limit] += ( 2.0 * gpop - gbot * std::conj( gbot ) ) * dt * dtau;
            topv[upper_limit] += std::pow( std::abs( akf_mat_g1.get( i, j ) ), 2.0 ) * dt * dtau;
        }
        Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), pbsize, "Indistinguishability (Simplified) (" + fout + "): " );
        timer.iterate();
    }

    Scalar topsum = 0, bottomsum = 0, topsumv = 0, bottomsumv = 0;

    for ( int i = 0; i < top.size(); i++ ) {
        double dt = akf_mat_g1.dt( i );
        double t_t = akf_mat_g1.t( i );
        topsum += top[i];
        bottomsum += bottom[i];
        topsumv += topv[i];
        bottomsumv += s.dgl_expectationvalue<Sparse, Scalar>( get_rho_at( rho_index_map[t_t] ), M1, t_t ) * dt;
        outpv[i] = 2.0 * topsumv / std::pow( bottomsumv, 2.0 );
        outp[i] = 1.0 - std::abs( topsum / bottomsum );
    }
    // Final output and timer end
    timer.end();
    Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), pbsize, "Indistinguishability (" + fout + "): ", Timers::PROGRESS_FORCE_OUTPUT );

    // Add to Fileoutput:
    if ( to_output["Indist"].size() == 0 )
        add_to_output( "Indist", "Time", time, to_output );
    add_to_output( "Indist", fout, outp, to_output );
    add_to_output( "Visibility", fout, outpv, to_output );

    Log::L1( "Final Indistinguishability: {} {}\n", std::real( outp.back() ), fout );
    Log::L1( "Final Visibility: {} {}\n", std::real( outpv.back() ), fout );

    return true;
}