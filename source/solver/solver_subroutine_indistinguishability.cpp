#include "solver/solver_ode.h"

bool QDLC::Numerics::ODESolver::calculate_indistinguishability( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator ) {
    // Set Number of Phonon cores to 1 because this memberfunction is already using multithreading
    s.parameters.numerics_phonons_maximum_threads = 1;
    // Progress
    ProgressBar progressbar = ProgressBar();

    // Calculate G2(t,tau) with given operator matrices
    std::string s_g1 = get_operators_purpose( { s_op_creator, s_op_annihilator }, 1 );
    std::string s_g2 = get_operators_purpose( { s_op_creator, s_op_annihilator, s_op_creator, s_op_annihilator }, 2 );

    auto [op_creator, op_annihilator] = calculate_g1( s, s_op_creator, s_op_annihilator, s_g1 );
    calculate_g2( s, s_op_creator, s_op_annihilator, s_op_creator, s_op_annihilator, s_g2 );

    auto &akf_mat_g1 = cache[s_g1];
    auto &akf_mat_g2 = cache[s_g2];
    auto &akf_mat_g1_time = cache[s_g1 + "_time"];
    auto &akf_mat_g2_time = cache[s_g2 + "_time"];

    int pbsize = akf_mat_g1.rows();

    Sparse M1 = op_creator * op_annihilator;

    std::string fout = s_op_creator + "-" + s_op_annihilator;
    Timer &timer = Timers::create( "Indistinguishability " + fout );
    timer.start();

    size_t T = std::min<size_t>( akf_mat_g1.rows(), savedStates.size() );

    std::vector<Scalar> top, bottom, topv, outp, outpv, time;

    for ( int i = 0; i < T; i++ ) {
        outp.emplace_back( 0 );
        outpv.emplace_back( 0 );
        top.emplace_back( 0 );
        bottom.emplace_back( 0 );
        topv.emplace_back( 0 );
        time.emplace_back( std::real( akf_mat_g1_time( i, 0 ) ) );
    }

    /**
     * @brief Iteratively extracts the next "section" of a triangular integral, e.g
     *
     *  :
     *  + :
     *  : + :
     *  : : + :
     *  : : : + :
     *  : : : : + :
     *  : : : : : + :
     * Where the "+" row is between current_iteration and upper_limit.
     *
     */
    // for ( int current_iteration = 0; current_iteration < T; current_iteration++ ) {
    //     Scalar result = 0.0;
    //     for ( int t = current_iteration; t < upper_limit; t++ ) {
    //         for ( int tau = 0; tau < upper_limit - current_iteration; tau++ ) {
    //             top[t] += function( i, j );
    //         }
    //     }
    // }

#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int upper_limit = 0; upper_limit < T; upper_limit++ ) {
        for ( int i = 0; i <= upper_limit; i++ ) {
            double t_t = std::real( akf_mat_g1_time( i, 0 ) );
            auto rho = getRhoAt( rho_index_map[t_t] );
            // for ( int j = 0; j < T - i and i + j < T; j++ ) {
            int j = upper_limit - i;
            double t_tau = std::real( akf_mat_g1_time( i + j, 0 ) ); // Important: t+tau (i+j)!
            auto rho_tau = getRhoAt( rho_index_map[t_tau] );
            Scalar gpop = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M1, t_tau );
            Scalar gbot = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator, t_t );
            double dt = Numerics::get_tdelta( akf_mat_g1_time, 0, i );
            double dtau = Numerics::get_taudelta( akf_mat_g1_time, i, j ); // Numerics::get_tdelta( akf_mat_g1_time, 0, i + j );
            top[upper_limit] += ( gpop + akf_mat_g2( i, j ) - akf_mat_g1( i, j ) * std::conj( akf_mat_g1( i, j ) ) ) * dt * dtau;
            bottom[upper_limit] += ( 2.0 * gpop - gbot * std::conj( gbot ) ) * dt * dtau;
            topv[upper_limit] += std::pow( std::abs( akf_mat_g1( i, j ) ), 2.0 ) * dt * dtau;
        }
        Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), pbsize, "Indistinguishability (Simplified) (" + fout + "): " );
        timer.iterate();
    }

    Scalar topsum = 0;
    Scalar bottomsum = 0;
    Scalar topsumv = 0;
    Scalar bottomsumv = 0;
    for ( int i = 0; i < top.size(); i++ ) {
        double dt = Numerics::get_tdelta( akf_mat_g1_time, 0, i );
        double t_t = std::real( akf_mat_g1_time( i, 0 ) );
        topsum += top[i];
        bottomsum += bottom[i];
        topsumv += topv[i];
        bottomsumv += s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( rho_index_map[t_t] ), M1, t_t ) * dt;
        outpv[i] = 2.0 * topsumv / std::pow( bottomsumv, 2.0 );
        outp[i] = 1.0 - std::abs( topsum / bottomsum );
    }
    // Final output and timer end
    timer.end();
    Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), pbsize, "Indistinguishability (" + fout + "): ", Timers::PROGRESS_FORCE_OUTPUT );

    // Add to Fileoutput:
    if ( to_output["Indist"].size() == 0 )
        to_output["Indist"]["Time"] = time;
    to_output["Indist"][fout] = outp;
    to_output["Visibility"][fout] = outpv;
    Log::L1( "Final Indistinguishability: {} {}\n", std::real( outp.back() ), fout );
    Log::L1( "Final Visibility: {} {}\n", std::real( outpv.back() ), fout );
    return true;
}