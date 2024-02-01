#include "solver/solver_ode.h"

/**
 * @brief Calculate the indistinguishability of two operators
 * 
 * @param s System to calculate the indistinguishability for
 * @param s_op_creator Operator name of the creator operator
 * @param s_op_annihilator Operator name of the annihilator operator
 *
 * The indistinguishability is calculated by the following formula:
 * \mathcal{I}_{i} = 1-p_{c,i} = 1-\frac{\int_0^{t_\text{max}}\int_0^{{t_\text{max}}-t} 2G_{\text{HOM},i}^{(2)}(t,t')dt'dt}{\int_0^{t_\text{max}}\int_0^{{t_\text{max}}-t}\left( 2G_{\text{pop},i}^{(2)}(t,t')- \abs{\ev{\hat{b}_i(t+t')} \ev{\hat{b}_i^\dagger(t)}}^2 \right)dt'dt}
 * I = 1-p_c, with p_c as the counts of coincidences and I as the indistinguishability
 * p_c = int_0^t_max int_0^(t_max-t) 2G_HOM(t,t')dt'dt / int_0^t_max int_0^(t_max-t) (2G_pop(t,t') - |<b(t+t')b^+(t)>|^2)dt'dt
 * G_HOM = 0.5* (G_pop(t,t') + G2(t,t') - |G1(t,t')|^2), G_pop = <b^+b>(t) <b^+b>(t+t')
 * 
 * The Integrals int_0^t_max int_0^(t_max-t) are calculated by a triangular integral.
 * The triangular integral looks as follows in an equidistant grid:
 * 
 * | # - - - - - - - - - |
 * | # # - - - - - - - - |
 * | # # # - - - - - - - |
 * | o # # # - - - - - - |
 * | # o # # # - - - - - |
 * | # # o # # # - - - - |
 * | # # # o # # # - - - |
 * | # # # # o # # # - - |
 * | # # # # # o # # # - |
 * | # # # # # # o # # # |
 * 
 * Where '#' are integrated values.
 * 
 * The integral is iterated by advancing the upper limit (o) of the integral by one step and adding the previous integral value to the current one.
 * On thes lower diagonal line, t will vary from 0 to upper_limit while t_tau = t_index + tau_index = upper_limit will always be constant.
 */

bool QDACC::Numerics::ODESolver::calculate_indistinguishability( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator ) {
    // Set Number of Phonon cores to 1 because this memberfunction is already using multithreading
    s.parameters.numerics_maximum_secondary_threads = 1;
    // Progress
    ProgressBar progressbar = ProgressBar();

    // Calculate G2(t,tau) with given operator matrices
    std::string s_g1 = get_operators_purpose( { s_op_creator, s_op_annihilator } );
    std::string s_g2 = get_operators_purpose( { s_op_creator, s_op_creator, s_op_annihilator, s_op_annihilator } );

    // Get Operator Matrices
    const auto &[op_creator, op_annihilator] = get_operators_matrices( s, s_op_creator, s_op_annihilator );

    // Calculate G-Functions if neccessary
    calculate_g1( s, s_op_creator, s_op_annihilator, s_g1 );
    calculate_g2( s, s_op_creator, s_op_creator, s_op_annihilator, s_op_annihilator, s_g2 );

    auto &akf_mat_g1 = cache[s_g1];
    auto &akf_mat_g2 = cache[s_g2];

    int pbsize = akf_mat_g1.dim();

    MatrixMain M1 = op_creator * op_annihilator;

    std::string fout = s_op_creator + "-" + s_op_annihilator;
    Timer &timer = Timers::create( "Indistinguishability " + fout ).start();

    size_t maximum_time = std::min<size_t>( pbsize, savedStates.size() );

    std::vector<Scalar> top( maximum_time, 0 ), bottom( maximum_time, 0 ), top_vis( maximum_time, 0 ), outp( maximum_time, 0 ), outpv( maximum_time, 0 );
    std::vector<Scalar> time;
    for ( size_t t = 0; t < maximum_time; t++ ) {
        double t_t = s.getTimeOf( 0, { t, 0 } );
        time.emplace_back( t_t );
    }
    // Do a Triangular Integral over G2(t,tau) and G1(t,tau)
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t upper_limit = 0; upper_limit < maximum_time; upper_limit++ ) {
        const size_t t_tau = s.getTimeOf( 0, { upper_limit, 0 } ); // = t_t + s.getTimeOf( 1, { t, tau } )
        const auto &current_tau_state = savedStates.at( upper_limit ); // = savedStates.at( rho_index_map[t_tau] );
        auto &rho_tau = current_tau_state.mat;
        for ( size_t t = 0; t <= upper_limit; t++ ) {
            // Current Time and Density Matrix
            double t_t = s.getTimeOf( 0, { t, 0 } );
            const auto t_t_index = rho_index_map[t_t];
            const auto &current_state = savedStates.at( t_t_index );
            auto &rho = current_state.mat;
            // t_tau = t+t' is always equal to upper_limit because of the triangular integral, but the tau index is not.            
            size_t tau = upper_limit - t;
            //double t_tau = t_t + s.getTimeOf( 1, { t, tau } );
            //const auto t_tau_index = rho_index_map[t_tau];
            //const auto &current_tau_state = savedStates.at( t_tau_index );
            //auto &rho_tau = current_tau_state.mat;
            // Delta Times
            double dt = s.getDeltaTimeOf( 0, { t, tau } );
            double dtau = s.getDeltaTimeOf( 1, { t, tau } );
            // Population Normalization constants
            Scalar gpop = s.dgl_expectationvalue<MatrixMain>( rho, M1, t_t ) * s.dgl_expectationvalue<MatrixMain>( rho_tau, M1, t_tau );
            Scalar gbot = s.dgl_expectationvalue<MatrixMain>( rho_tau, op_annihilator, t_tau ) * s.dgl_expectationvalue<MatrixMain>( rho, op_creator, t_t );
            // Required summands for indistinguishability
            const auto g1 = akf_mat_g1.get( t, tau );
            const auto g2 = akf_mat_g2.get( t, tau );
            const auto g1_abs2 = std::pow( std::abs( g1 ), 2.0 );
            const auto gbot_abs2 = std::pow( std::abs( gbot ), 2.0 );
            // Upper and Lower Integral for Indistinguishability
            top[upper_limit] += ( gpop + g2 - g1_abs2 ) * dt * dtau;
            bottom[upper_limit] += ( 2.0 * gpop - gbot_abs2 ) * dt * dtau;
            top_vis[upper_limit] += g1_abs2 * dt * dtau;
        }
        Timers::outputProgress( timer, progressbar, timer.getTotalIterationNumber(), pbsize, "Indistinguishability (Simplified) (" + fout + "): " );
        timer.iterate();
    }

    Scalar topsum = 0, bottomsum = 0, topsumv = 0, bottomsumv = 0;

    for ( size_t t = 0; t < top.size(); t++ ) {
        double dt = s.getDeltaTimeOf( 0, { t, 0 } );
        double t_t = s.getTimeOf( 0, { t, 0 } );
        topsum += top[t];
        bottomsum += bottom[t];
        topsumv += top_vis[t];
        const auto t_t_index = rho_index_map[t_t];
        const auto &current_state = savedStates.at( t_t_index );
        bottomsumv += s.dgl_expectationvalue<MatrixMain>( current_state.mat, M1, t_t ) * dt;
        outpv[t] = 2.0 * topsumv / std::pow( bottomsumv, 2.0 );
        outp[t] = 1.0 - std::abs( topsum / bottomsum );
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