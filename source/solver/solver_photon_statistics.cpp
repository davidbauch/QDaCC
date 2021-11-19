#include "solver/solver_ode.h"
#include <cmath>
#include <complex>

// Description: Calculates G2 Indistinguishability
bool QDLC::Numerics::ODESolver::calculate_indistinguishability( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator ) {
    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
    s.command( QDLC::Numerics::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
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

    std::vector<Scalar> outp, outpv;
    std::vector<Scalar> time;

    Sparse M1 = op_creator * op_annihilator;

    std::string fout = s_op_creator + "-" + s_op_annihilator;
    Timer &timer = Timers::create( "Indistinguishability " + fout );
    timer.start();

    size_t T = std::min<size_t>( akf_mat_g1.rows(), savedStates.size() ); // Make sure G1 and G2 have same rows/cols/dts?

    std::vector<Scalar> top, bottom, topv;
    for ( int i = 0; i < T; i++ ) {
        outp.emplace_back( 0 );
        outpv.emplace_back( 0 );
        top.emplace_back( 0 );
        bottom.emplace_back( 0 );
        topv.emplace_back( 0 );
        time.emplace_back( std::real( akf_mat_g1_time( i, 0 ) ) );
    }
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < T; i++ ) {
        double t_t = std::real( akf_mat_g1_time( i, 0 ) );
        auto rho = getRhoAt( rho_index_map[t_t] );
        for ( int j = 0; j < T - i and i + j < T; j++ ) {
            double t_tau = std::real( akf_mat_g1_time( i + j, 0 ) );
            auto rho_tau = getRhoAt( rho_index_map[t_tau] );
            Scalar gpop = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M1, t_tau );
            Scalar gbot = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator, t_t );
            double dt = Numerics::get_tdelta( akf_mat_g1_time, j, i );
            double dtau = Numerics::get_taudelta( akf_mat_g1_time, i, j );
            top[i] += ( gpop + akf_mat_g2( i, j ) - akf_mat_g1( i, j ) * std::conj( akf_mat_g1( i, j ) ) ) * dt * dtau;
            bottom[i] += ( 2.0 * gpop - gbot * std::conj( gbot ) ) * dt * dtau;
            topv[i] += std::pow( std::abs( akf_mat_g1( i, j ) ), 2.0 ) * dt * dtau;
        }
        Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, timer.getTotalIterationNumber(), pbsize, "Indistinguishability (Simplified) (" + fout + "): " );
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
    Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, timer.getTotalIterationNumber(), pbsize, "Indistinguishability (" + fout + "): ", Timers::PROGRESS_FORCE_OUTPUT );

    // Add to Fileoutput:
    if ( to_output["Indist"].size() == 0 )
        to_output["Indist"]["Time"] = time;
    to_output["Indist"][fout] = outp;
    to_output["Visibility"][fout] = outpv;
    Log::L1( "Final Indistinguishability: {} {}\n", std::real( outp.back() ), fout );
    Log::L1( "Final Visibility: {} {}\n", std::real( outpv.back() ), fout );
    return true;
}

// Description: Calculates Concurrence
bool QDLC::Numerics::ODESolver::calculate_concurrence( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2 ) {
    Log::L2( "[Concurrence] Conc for modes {} {} and {} {}\n", s_op_creator_1, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_2 );

    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
    s.command( QDLC::Numerics::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
    // Progress
    ProgressBar progressbar = ProgressBar();

    std::string fout = s_op_creator_1 + "-" + s_op_annihilator_1 + "-" + s_op_creator_2 + "-" + s_op_annihilator_2;
    // Calculate G2(t,tau) with given operator matrices
    std::string s_g2_1111 = get_operators_purpose( { s_op_creator_1, s_op_annihilator_1, s_op_creator_1, s_op_annihilator_1 }, 2 );
    std::string s_g2_1122 = get_operators_purpose( { s_op_creator_1, s_op_annihilator_2, s_op_creator_1, s_op_annihilator_2 }, 2 );
    std::string s_g2_1212 = get_operators_purpose( { s_op_creator_1, s_op_annihilator_1, s_op_creator_2, s_op_annihilator_2 }, 2 );
    std::string s_g2_1221 = get_operators_purpose( { s_op_creator_1, s_op_annihilator_2, s_op_creator_2, s_op_annihilator_1 }, 2 );
    std::string s_g2_2121 = get_operators_purpose( { s_op_creator_2, s_op_annihilator_2, s_op_creator_1, s_op_annihilator_1 }, 2 );
    std::string s_g2_2112 = get_operators_purpose( { s_op_creator_2, s_op_annihilator_1, s_op_creator_1, s_op_annihilator_2 }, 2 );
    std::string s_g2_2211 = get_operators_purpose( { s_op_creator_2, s_op_annihilator_1, s_op_creator_2, s_op_annihilator_1 }, 2 );
    std::string s_g2_2222 = get_operators_purpose( { s_op_creator_2, s_op_annihilator_2, s_op_creator_2, s_op_annihilator_2 }, 2 );

    calculate_g2( s, s_op_creator_1, s_op_annihilator_1, s_op_creator_1, s_op_annihilator_1, s_g2_1111 );
    calculate_g2( s, s_op_creator_1, s_op_annihilator_2, s_op_creator_1, s_op_annihilator_2, s_g2_1122 );
    calculate_g2( s, s_op_creator_2, s_op_annihilator_2, s_op_creator_1, s_op_annihilator_1, s_g2_2121 );
    auto [op_creator_2, op_annihilator_1, op_creator_1, op_annihilator_2] = calculate_g2( s, s_op_creator_2, s_op_annihilator_1, s_op_creator_1, s_op_annihilator_2, s_g2_2112 );
    calculate_g2( s, s_op_creator_1, s_op_annihilator_2, s_op_creator_2, s_op_annihilator_1, s_g2_1221 );
    calculate_g2( s, s_op_creator_2, s_op_annihilator_2, s_op_creator_2, s_op_annihilator_2, s_g2_2222 );
    cache[s_g2_2211] = cache[s_g2_1122].conjugate();
    cache[s_g2_2211+"_time"] = cache[s_g2_1122+"_time"];
    cache[s_g2_1212] = cache[s_g2_2121].conjugate();
    cache[s_g2_1212+"_time"] = cache[s_g2_2121+"_time"];

    int pbsize = 2 * cache[s_g2_1212].rows();

    std::map<std::string, std::vector<Scalar>> rho;
    std::map<std::string, std::vector<Scalar>> rho_g2zero;
    std::map<std::string, Sparse> matmap_g2zero = { { s_g2_1111, op_creator_1 * op_creator_1 * op_annihilator_1 * op_annihilator_1 },
                                                    { s_g2_1122, op_creator_1 * op_creator_1 * op_annihilator_2 * op_annihilator_2 },
                                                    { s_g2_1212, op_creator_1 * op_creator_2 * op_annihilator_1 * op_annihilator_2 },
                                                    { s_g2_1221, op_creator_1 * op_creator_2 * op_annihilator_2 * op_annihilator_1 },
                                                    { s_g2_2121, op_creator_2 * op_creator_1 * op_annihilator_2 * op_annihilator_1 },
                                                    { s_g2_2112, op_creator_2 * op_creator_1 * op_annihilator_1 * op_annihilator_2 },
                                                    { s_g2_2211, op_creator_2 * op_creator_2 * op_annihilator_1 * op_annihilator_1 },
                                                    { s_g2_2222, op_creator_2 * op_creator_2 * op_annihilator_2 * op_annihilator_2 } };

    Timer &timer_c = Timers::create( "Concurrence (" + fout + ")" );
    timer_c.start();

    auto T = std::min<size_t>( cache[s_g2_1111].rows(), savedStates.size() );
    auto &mat_time = cache[s_g2_1111 + "_time"];

    for ( long unsigned int t = 0; t < T; t++ ) {
        for ( auto &mode : { s_g2_1111, s_g2_1122, s_g2_1212, s_g2_1221, s_g2_2121, s_g2_2112, s_g2_2211, s_g2_2222 } ) {
            rho[mode].emplace_back( 0 );
            rho_g2zero[mode].emplace_back( 0 );
        }
    }
    for ( auto &mode : { s_g2_1111, s_g2_1122, s_g2_1212, s_g2_1221, s_g2_2121, s_g2_2112, s_g2_2211, s_g2_2222 } ) {
        auto &gmat_time = cache[mode + "_time"];
        bool didout = false;
        rho[mode][0] = cache[mode]( 0, 0 );
        rho_g2zero[mode][0] = s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( 0 ), matmap_g2zero[mode], getTimeAt( 0 ) ) * Numerics::get_tdelta( gmat_time, 0, 0 );
        for ( long unsigned int i = 1; i < T; i++ ) {
            double dt = Numerics::get_tdelta( gmat_time, 0, i );
            rho[mode][i] = rho[mode][i - 1];
            for ( int tau = 0; tau < T - i; tau++ ) {
                double dtau = Numerics::get_taudelta( gmat_time, i, tau );
                rho[mode][i] += cache[mode]( i, tau ) * dt * dtau;
            }
            double t_t = std::real( mat_time( i, 0 ) ); //Note:: all correlation functions have to have the same times. cache[mode+"_time"] else.
            rho_g2zero[mode][i] = rho_g2zero[mode][i - 1] + s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( rho_index_map[t_t] ), matmap_g2zero[mode], t_t ) * dt;

            if ( not didout ) {
                didout = true;
                timer_c.iterate();
                Timers::outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Concurrence (" + fout + "): " );
            }
        }
    }
    // Calculate EigenValues
    std::vector<Scalar> output, output_g2zero;
    std::vector<Scalar> time;
    std::vector<Dense> twophotonmatrix, twophotonmatrix_g2zero;
    for ( long unsigned int i = 0; i < T; i++ ) {
        output.emplace_back( 0 );
        output_g2zero.emplace_back( 0 );
        time.emplace_back( 0 );
        twophotonmatrix.emplace_back( Dense::Zero( 4, 4 ) );
        twophotonmatrix_g2zero.emplace_back( Dense::Zero( 4, 4 ) );
    }
    Dense spinflip = Dense::Zero( 4, 4 );
    spinflip( 0, 3 ) = -1;
    spinflip( 1, 2 ) = 1;
    spinflip( 2, 1 ) = 1;
    spinflip( 3, 0 ) = -1;
    Log::L3( "[Concurrence] Spinflip Matrix: {}\n", spinflip );
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_threads )
    for ( long unsigned int k = 0; k < T; k++ ) { // cache[s_g2_1111].rows() statt T?
        //Log::L3( "Creating 2 photon matrix\n" );
        Dense rho_2phot = Dense::Zero( 4, 4 );
        Dense rho_2phot_g2zero = Dense::Zero( 4, 4 );
        rho_2phot( 0, 0 ) = rho[s_g2_1111][k];
        rho_2phot( 3, 3 ) = rho[s_g2_2222][k];
        rho_2phot( 0, 3 ) = rho[s_g2_2211][k];
        rho_2phot( 3, 0 ) = rho[s_g2_1122][k];
        rho_2phot( 1, 2 ) = rho[s_g2_2121][k];
        rho_2phot( 2, 1 ) = rho[s_g2_1212][k];
        rho_2phot( 1, 1 ) = std::abs( rho[s_g2_2112][k] ) != 0 ? rho[s_g2_2112][k] : 1E-100;
        rho_2phot( 2, 2 ) = std::abs( rho[s_g2_1221][k] ) != 0 ? rho[s_g2_1221][k] : 1E-100;
        rho_2phot = rho_2phot / rho_2phot.trace();
        rho_2phot_g2zero( 0, 0 ) = rho_g2zero[s_g2_1111][k];
        rho_2phot_g2zero( 3, 3 ) = rho_g2zero[s_g2_2222][k];
        rho_2phot_g2zero( 0, 3 ) = rho_g2zero[s_g2_2211][k];
        rho_2phot_g2zero( 3, 0 ) = rho_g2zero[s_g2_1122][k];
        rho_2phot_g2zero( 1, 2 ) = rho_g2zero[s_g2_2121][k];
        rho_2phot_g2zero( 2, 1 ) = rho_g2zero[s_g2_1212][k];
        rho_2phot_g2zero( 1, 1 ) = std::abs( rho_g2zero[s_g2_2112][k] ) != 0 ? rho_g2zero[s_g2_2112][k] : 1E-100;
        rho_2phot_g2zero( 2, 2 ) = std::abs( rho_g2zero[s_g2_1221][k] ) != 0 ? rho_g2zero[s_g2_1221][k] : 1E-100;
        rho_2phot_g2zero = rho_2phot_g2zero / rho_2phot_g2zero.trace();
        //Log::L3( "Normalizing 2 photon matrix\n" );
        //if ( std::abs( rho_2phot.trace() ) != 0 ) {
        //Log::L3( "Rho_2phot = {}\n", rho_2phot );
        //Log::L3( "Calculating sqrt(rho)\n" );
        //Eigen::MatrixPower<Dense> Mpow( rho_2phot );
        Dense sqrtrho2phot = rho_2phot.sqrt();               //Mpow( 0.5 );
        Dense sqrtrho2phot_g2zero = rho_2phot_g2zero.sqrt(); //Mpow( 0.5 );
        //Log::L3( "Calculating R\n" );
        Dense R = sqrtrho2phot * spinflip * rho_2phot * spinflip * sqrtrho2phot;
        Dense R_g2zero = sqrtrho2phot_g2zero * spinflip * rho_2phot_g2zero * spinflip * sqrtrho2phot_g2zero;
        //Log::L3( "R = {}\n", R );
        //Log::L3( "Calculating sqrt(R)\n" );
        //Eigen::MatrixPower<Dense> SMPow( R );
        auto R5 = R.sqrt();               //SMPow( 0.5 );
        auto R5_g2zero = R_g2zero.sqrt(); //SMPow( 0.5 );
        //Log::L3( "Calculating Eigenvalues\n" );
        auto eigenvalues = R5.eigenvalues();
        auto eigenvalues_g2zero = R5_g2zero.eigenvalues();
        // Sometimes, the numerical method for eigenvalue evaluation yields crap (first element really big, rest shifted), so we check this here:
        if ( QDLC::Math::abs2( eigenvalues( 3 ) ) > 5.0 ) {
            eigenvalues( 3 ) = eigenvalues( 2 );
            eigenvalues( 2 ) = eigenvalues( 1 );
            eigenvalues( 1 ) = eigenvalues( 0 );
            eigenvalues( 0 ) = 0.0;
        }
        if ( QDLC::Math::abs2( eigenvalues_g2zero( 3 ) ) > 5.0 ) {
            eigenvalues_g2zero( 3 ) = eigenvalues_g2zero( 2 );
            eigenvalues_g2zero( 2 ) = eigenvalues_g2zero( 1 );
            eigenvalues_g2zero( 1 ) = eigenvalues_g2zero( 0 );
            eigenvalues_g2zero( 0 ) = 0.0;
        }
        //Log::L1( "rho2phot = {}\n\nsqrtrho2phot = {}\n\nR = {}\n\nRS = {}\nEigenvalues at t = {} are {}\n", rho_2phot, sqrtrho2phot, R, R5, getTimeAt( i ), eigenvalues );
        auto conc = eigenvalues( 3 ) - eigenvalues( 2 ) - eigenvalues( 1 ) - eigenvalues( 0 );
        //Log::L1("Eigenvalues {}: C = {} - {} - {} - {}\n",k,eigenvalues(3),eigenvalues(2),eigenvalues(1),eigenvalues(0));
        auto conc_g2zero = eigenvalues_g2zero( 3 ) - eigenvalues_g2zero( 2 ) - eigenvalues_g2zero( 1 ) - eigenvalues_g2zero( 0 );
        output.at( k ) = conc;
        output_g2zero.at( k ) = conc_g2zero;
        time.at( k ) = std::real( mat_time( k, 0 ) );
        twophotonmatrix.at( k ) = rho_2phot;
        twophotonmatrix_g2zero.at( k ) = rho_2phot_g2zero;
        timer_c.iterate();
    }
    // Final output and timer end
    timer_c.end();
    Timers::outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Concurrence (" + fout + ")", Timers::PROGRESS_FORCE_OUTPUT );
    // Add to Fileoutput:
    if ( to_output["Conc"].size() == 0 )
        to_output["Conc"]["Time"] = time;
    to_output["Conc"][fout] = output;
    to_output["Conc_g2zero"][fout] = output_g2zero;
    to_output_m["TwoPMat"][fout] = twophotonmatrix;
    to_output_m["TwoPMat"][fout + "_g2zero"] = twophotonmatrix_g2zero;
    Log::L1( "Final Concurrence: {:.10f} ({:.10f} with g2(0)) {}\n", std::real( output.back() ), std::real( output_g2zero.back() ), fout );
    return true;
}

bool QDLC::Numerics::ODESolver::calculate_wigner( System &s, const std::string &s_mode, const double x, const double y, const int resolution, const int skips ) {
    // Partially trace the mode to wigner:
    std::vector<Dense> reduced_rho;
    std::vector<Scalar> time;
    int base = s.operatorMatrices.el_states.count( s_mode ) != 0 ? s.operatorMatrices.el_states[s_mode].base : s.operatorMatrices.ph_states[s_mode].base;
    Log::L2( "[Wigner] Calculating Wigner function for mode {}/{}\n", s_mode, base );
    for ( int i = 0; i < savedStates.size(); i += skips ) {
        reduced_rho.emplace_back( s.partialTrace( getRhoAt( i ), base ) );
        time.emplace_back( getTimeAt( i ) );
    }

    double g = std::sqrt( 2 );
    // Calculate Wigner functions over time:
    // Calculate Meshgrid
    auto [X_mat, Y_mat] = QDLC::Matrix::meshgrid( -x, -y, x, y, resolution );
    //Dense A = 0.5 * g * ( X_mat + 1.0i * Y_mat );
    Dense A = ( X_mat + 1.0i * Y_mat );

    std::vector<Dense> wigner( reduced_rho.size(), Dense::Zero( A.rows(), A.cols() ) );
    Log::L2( "[Wigner] First Matrix for Wigner:\n{}\n", reduced_rho.front() );

    ProgressBar progressbar = ProgressBar();
    Timer &timer_w = Timers::create( "Wigner (" + s_mode + ")" );
    timer_w.start();
    // Iterate
    //Calculates the Wigner function at coordinate alpha
    auto Wigner = [&]( Dense &DM, Scalar alpha ) {
        double Wval = 0.0;
        for ( int m = 0; m < DM.rows(); m++ ) {
            if ( QDLC::Math::abs2( DM( m, m ) ) > 0.0 ) {
                Wval += std::real( DM( m, m ) * std::pow( -1., (double)m ) * std::assoc_laguerre( m, 0, 4. * QDLC::Math::abs2( alpha ) ) );
            }
        }

        for ( int m = 0; m < DM.rows(); m++ ) {
            for ( int n = m + 1; n < DM.rows(); n++ ) {
                if ( QDLC::Math::abs2( DM( m, n ) ) > 0.0 ) {
                    Wval += 2. * std::real( DM( m, n ) * std::pow( -1., (double)m ) * std::pow( 2. * alpha, (double)( n - m ) ) * std::sqrt( QDLC::Math::factorial( m ) / QDLC::Math::factorial( n ) ) * std::assoc_laguerre( m, n - m, 4. * QDLC::Math::abs2( alpha ) ) );
                }
            }
        }
        return ( 1 / ( 2. * 3.1415 ) * std::exp( -2. * QDLC::Math::abs2( alpha ) ) * Wval );
    };
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < reduced_rho.size(); i++ ) {
        Dense W = Dense::Zero( A.rows(), A.cols() );
        for ( int k = 0; k < A.rows(); k++ ) {
            for ( int l = 0; l < A.cols(); l++ ) {
                W( k, l ) = Wigner( reduced_rho[i], A( k, l ) );
            }
        }
        wigner.at( i ) = W;
        timer_w.iterate();
        Timers::outputProgress( s.parameters.output_handlerstrings, timer_w, progressbar, timer_w.getTotalIterationNumber(), reduced_rho.size(), "Wigner (" + s_mode + "): " );
    }
    timer_w.end();
    Timers::outputProgress( s.parameters.output_handlerstrings, timer_w, progressbar, timer_w.getTotalIterationNumber(), reduced_rho.size(), "Wigner (" + s_mode + "): ", Timers::PROGRESS_FORCE_OUTPUT );
    // Add to Fileoutput:
    if ( to_output["Wigner"].size() == 0 )
        to_output["Wigner"]["Time"] = time;
    to_output_m["Wigner"][s_mode] = wigner;
    to_output_m["Wigner"]["rho_" + s_mode] = reduced_rho;
    return true;
}

// Description: Calculates and outputs G1 and G2 photon statistics
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator(annihilator)_X: set of operator matrices and their conjugated forms
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false

bool QDLC::Numerics::ODESolver::calculate_advanced_photon_statistics( System &s ) {
    // Calculate Spectra
    auto &spectrum_s = s.parameters.input_correlation["Spectrum"];
    for ( int i = 0; i < spectrum_s.string_v["Modes"].size(); i++ ) {
        const auto &[s_creator, s_annihilator] = get_operator_strings( spectrum_s.string_v["Modes"][i] );
        calculate_spectrum( s, s_creator, s_annihilator, spectrum_s.numerical_v["Center"][i], spectrum_s.numerical_v["Range"][i], spectrum_s.numerical_v["resW"][i] );
    }
    // Calculate Indist
    auto &indist_s = s.parameters.input_correlation["Indist"];
    for ( int i = 0; i < indist_s.string_v["Modes"].size(); i++ ) {
        const auto &[s_creator, s_annihilator] = get_operator_strings( indist_s.string_v["Modes"][i] );
        calculate_indistinguishability( s, s_creator, s_annihilator );
    }
    // Calculate Conc
    auto &conc_s = s.parameters.input_correlation["Conc"];
    for ( auto &modes : conc_s.string_v["Modes"] ) {
        std::vector<std::string> s_creator, s_annihilator;
        for ( auto &mode : QDLC::String::splitline( modes, '-' ) ) {
            const auto &[ss_creator, ss_annihilator] = get_operator_strings( mode );
            s_creator.emplace_back( ss_creator );
            s_annihilator.emplace_back( ss_annihilator );
        }
        calculate_concurrence( s, s_creator[0], s_annihilator[0], s_creator[1], s_annihilator[1] );
    }
    // Calculate G1/G2 functions
    auto &gs_s = s.parameters.input_correlation["GFunc"];
    for ( int i = 0; i < gs_s.string_v["Modes"].size(); i++ ) {
        double t_step = ( s.parameters.numerics_phonon_approximation_order == PHONON_PATH_INTEGRAL ? s.parameters.t_step_pathint : s.parameters.t_step );
        /// TODO : in funktion
        auto modes = gs_s.string_v["Modes"][i];
        int order = std::abs( gs_s.numerical_v["Order"][i] );
        const auto &[s_creator, s_annihilator] = get_operator_strings( modes );
        std::string purpose = order == 1 ? get_operators_purpose( { s_creator, s_annihilator }, 1 ) : get_operators_purpose( { s_creator, s_annihilator, s_creator, s_annihilator }, 2 );
        Sparse creator, annihilator;
        if ( order == 1 ) {
            auto [a, b] = calculate_g1( s, s_creator, s_annihilator, purpose );
            creator = std::move( a );
            annihilator = std::move( b );
        } else {
            auto [a, b, _discard1, _discard2] = calculate_g2( s, s_creator, s_annihilator, s_creator, s_annihilator, purpose );
            creator = std::move( a );
            annihilator = std::move( b );
        }
        // Directly output corresponding matrix here so G1/2 functions calculated by other function calls are not output if they are not demanded.
        auto &gmat = cache[purpose];
        auto &gmat_time = cache[purpose + "_time"];
        // G2(t,tau)
        if ( gs_s.numerical_v["Integrated"][i] == 0 || gs_s.numerical_v["Integrated"][i] == 2 ) {
            Log::L2( "Saving G{} function matrix to {}_m.txt...\n", order, purpose );
            FILE *f_gfunc = std::fopen( ( s.parameters.subfolder + purpose + "_m.txt" ).c_str(), "w" );
            fmt::print( f_gfunc, "Time\tTau\tAbs\tReal\tImag\n" );
            for ( int k = 0; k < gmat.rows(); k++ ) {
                double t_t = std::real( gmat_time( k, 0 ) );
                for ( int l = 0; l < gmat.cols(); l++ ) {
                    double t_tau = std::imag( gmat_time( k, l ) ) - std::real( gmat_time( k, l ) );
                    fmt::print( f_gfunc, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t_t, t_tau, std::abs( gmat( k, l ) ), std::real( gmat( k, l ) ), std::imag( gmat( k, l ) ) );
                }
                fmt::print( f_gfunc, "\n" );
            }
            std::fclose( f_gfunc );
        }
        // G2(0)
        std::vector<Scalar> topv, g2ofzero;
        for ( int i = 0; i < std::min<int>( gmat.cols(), savedStates.size() ); i++ ) {
            g2ofzero.emplace_back( 0 );
            topv.emplace_back( 0 );
        }
#pragma omp parallel for schedule( dynamic ) num_threads( s.parameters.numerics_maximum_threads )
        for ( int i = 0; i < std::min<int>( gmat.cols(), savedStates.size() ); i++ ) {
            for ( int j = 0; j < gmat.cols() - i; j++ ) {
                topv[i] += gmat( i, j );
            }
        }
        Scalar topsumv = 0;
        Scalar bottomsumv = 0;
        for ( int k = 0; k < topv.size(); k++ ) {
            double t_t = std::real( gmat_time( k, 0 ) );
            int t = rho_index_map[t_t];
            topsumv += topv[k];
            bottomsumv += s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( t ), creator * annihilator, t_t );
            g2ofzero[k] = 2.0 * topsumv / std::pow( bottomsumv, 2.0 );
        }
        // G2(t,0) and G2(tau)
        if ( gs_s.numerical_v["Integrated"][i] == 1 || gs_s.numerical_v["Integrated"][i] == 2 ) {
            Log::L2( "[PhotonStatistics] Saving G{} integrated function to {}.txt...\n", order, purpose );
            FILE *f_gfunc = std::fopen( ( s.parameters.subfolder + purpose + ".txt" ).c_str(), "w" );
            fmt::print( f_gfunc, "Time\tAbs(g{0}(tau))\tReal(g{0}(tau))\tImag(g{0}(tau))\tAbs(g{0}(t,0))\tReal(g{0}(t,0))\tImag(g{0}(t,0))\tAbs(g{0}(0))\tReal(g{0}(0))\tImag(g{0}(0))\n", order );
            for ( int l = 0; l < topv.size(); l++ ) { //gmat.cols()
                Scalar g2oftau = 0;
                for ( int k = 0; k < gmat.rows(); k++ ) {
                    g2oftau += gmat( k, l ) * Numerics::get_tdelta( gmat_time, l, k );
                }
                double t_tau = std::imag( gmat_time( 0, l ) );
                size_t tau_index = rho_index_map[t_tau];
                Scalar g2oft = s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( tau_index ), creator * creator * annihilator * annihilator, t_tau ); // / std::pow( s.dgl_expectationvalue<Sparse, Scalar>( getRhoAt( l ), creator * annihilator, getTimeAt( l ) ), 2.0 );
                fmt::print( f_gfunc, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", t_tau, std::abs( g2oftau ), std::real( g2oftau ), std::imag( g2oftau ), std::abs( g2oft ), std::real( g2oft ), std::imag( g2oft ), std::abs( g2ofzero[l] ), std::real( g2ofzero[l] ), std::imag( g2ofzero[l] ) );
            }
            std::fclose( f_gfunc );
        }
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    // Calculate Conc
    auto &wigner_s = s.parameters.input_correlation["Wigner"];
    for ( int i = 0; i < wigner_s.string_v["Modes"].size(); i++ ) {
        calculate_wigner( s, wigner_s.string_v["Modes"][i], wigner_s.numerical_v["X"][i], wigner_s.numerical_v["Y"][i], wigner_s.numerical_v["Res"][i], wigner_s.numerical_v["Skip"][i] );
    }

    // Output Spectra and Rest in seperate Files
    for ( auto &[mode, data] : to_output["Spectrum"] ) {
        Log::L2( "[PhotonStatistics] Saving Emission Spectrum to spectrum_" + mode + ".txt...\n" );
        FILE *f_spectrum = std::fopen( ( s.parameters.subfolder + "spectrum_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_spectrum, "Omega\t{}\n", mode );
        for ( int i = 0; i < to_output["Spectrum"][mode].size(); i++ ) {
            fmt::print( f_spectrum, "{:.8e}\t{:.8e}\n", std::real( to_output["Spectrum_frequency"][mode][i] ), std::real( to_output["Spectrum"][mode][i] ) );
        }
        std::fclose( f_spectrum );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output["Indist"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Indistinguishability and Visibility to indist_" + mode + ".txt...\n" );
        FILE *f_indist = std::fopen( ( s.parameters.subfolder + "indist_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_indist, "Time\tIndist_{}\tVisibility_{}\n", mode, mode );
        for ( int i = 0; i < to_output["Indist"][mode].size(); i++ ) {
            fmt::print( f_indist, "{:.8e}\t{:.8e}\t{:.8e}\n", std::real( to_output["Indist"]["Time"][i] ), std::real( to_output["Indist"][mode][i] ), std::real( to_output["Visibility"][mode][i] ) );
        }
        std::fclose( f_indist );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output["Conc"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Concurrence to conc_" + mode + ".txt...\n" );
        FILE *f_indist = std::fopen( ( s.parameters.subfolder + "conc_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_indist, "Time\t{}\t{}(g2(0))\n", mode, mode );
        for ( int i = 0; i < to_output["Conc"][mode].size(); i++ ) {
            fmt::print( f_indist, "{:.8e}\t{:.8e}\t{:.8e}\n", std::real( to_output["Conc"]["Time"][i] ), std::real( to_output["Conc"][mode][i] ), std::real( to_output["Conc_g2zero"][mode][i] ) );
        }
        std::fclose( f_indist );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output_m["TwoPMat"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Two-Photon Matrix to twopmat_" + mode + ".txt...\n" );
        FILE *f_twophot = std::fopen( ( s.parameters.subfolder + "twopmat_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_twophot, "Time\t" );
        std::vector<std::string> modes = { "11", "12", "21", "22" };
        for ( int k = 0; k < 4; k++ ) {
            for ( int l = 0; l < 4; l++ ) {
                fmt::print( f_twophot, "Re({}{})\t", modes[k], modes[l] );
            }
        }
        for ( int k = 0; k < 4; k++ ) {
            for ( int l = 0; l < 4; l++ ) {
                fmt::print( f_twophot, "Im({}{})\t", modes[k], modes[l] );
            }
        }
        fmt::print( f_twophot, "\n" );
        for ( int i = 0; i < to_output_m["TwoPMat"][mode].size(); i++ ) {
            fmt::print( f_twophot, "{:.8e}\t", std::real( to_output["Conc"]["Time"][i] ) );
            for ( int k = 0; k < 4; k++ ) {
                for ( int l = 0; l < 4; l++ ) {
                    fmt::print( f_twophot, "{:.8e}\t", std::real( to_output_m["TwoPMat"][mode][i]( k, l ) ) );
                }
                for ( int k = 0; k < 4; k++ ) {
                    for ( int l = 0; l < 4; l++ ) {
                        fmt::print( f_twophot, "{:.8e}\t", std::imag( to_output_m["TwoPMat"][mode][i]( k, l ) ) );
                    }
                }
            }
            fmt::print( f_twophot, "\n" );
        }
        std::fclose( f_twophot );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    for ( auto &[mode, data] : to_output_m["Wigner"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "[PhotonStatistics] Saving Wigner function to wigner_" + mode + ".txt...\n" );
        FILE *f_wigner = std::fopen( ( s.parameters.subfolder + "wigner_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_wigner, "Time\t{}\n", mode );
        for ( int i = 0; i < data.size(); i++ ) {
            fmt::print( f_wigner, "{:.8e}\t", std::real( to_output["Wigner"]["Time"][i] ) );
            auto &currentwigner = data[i];
            for ( int k = 0; k < currentwigner.rows(); k++ ) {
                for ( int l = 0; l < currentwigner.cols(); l++ ) {
                    fmt::print( f_wigner, "{:.8e}\t", std::real( currentwigner( k, l ) ) );
                }
            }
            fmt::print( f_wigner, "\n" );
        }
        std::fclose( f_wigner );
        Log::L2( "[PhotonStatistics] Done!\n" );
    }
    return true;
}