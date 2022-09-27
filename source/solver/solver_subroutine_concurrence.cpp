#include "solver/solver_ode.h"

// Description: Calculates Concurrence
bool QDLC::Numerics::ODESolver::calculate_concurrence( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2 ) {
    Log::L2( "[Concurrence] Conc for modes {} {} and {} {}\n", s_op_creator_1, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_2 );

    // Set Number of Phonon cores to 1 because this memberfunction is already using multithreading
    s.parameters.numerics_maximum_secondary_threads = 1;
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

    // Get Sparse Operator Matrices
    const auto &[op_creator_1, op_annihilator_1] = get_operators_matrices( s, s_op_creator_1, s_op_annihilator_1 );
    const auto &[op_creator_2, op_annihilator_2] = get_operators_matrices( s, s_op_creator_2, s_op_annihilator_2 );

    // Calculate G2's if neccessary
    // calculate_g2( s, s_op_creator_1, s_op_annihilator_1, s_op_creator_1, s_op_annihilator_1, s_g2_1111 );
    // calculate_g2( s, s_op_creator_1, s_op_annihilator_2, s_op_creator_1, s_op_annihilator_2, s_g2_1122 );
    // calculate_g2( s, s_op_creator_2, s_op_annihilator_2, s_op_creator_1, s_op_annihilator_1, s_g2_2121 );
    // calculate_g2( s, s_op_creator_2, s_op_annihilator_1, s_op_creator_1, s_op_annihilator_2, s_g2_2112 );
    // calculate_g2( s, s_op_creator_1, s_op_annihilator_2, s_op_creator_2, s_op_annihilator_1, s_g2_1221 );
    // calculate_g2( s, s_op_creator_2, s_op_annihilator_2, s_op_creator_2, s_op_annihilator_2, s_g2_2222 );
    // Make sure to only evaluate the tau-direction when absolutely neccessary
    calculate_g2( s, s_op_creator_1, { s_op_annihilator_1, s_op_annihilator_2 }, { s_op_creator_1, s_op_creator_2 }, s_op_annihilator_1, { s_g2_1111, s_g2_1221 } );
    calculate_g2( s, s_op_creator_2, { s_op_annihilator_1, s_op_annihilator_2 }, { s_op_creator_1, s_op_creator_2 }, s_op_annihilator_2, { s_g2_2112, s_g2_2222 } );
    calculate_g2( s, s_op_creator_1, { s_op_annihilator_2, s_op_annihilator_1 }, { s_op_creator_1, s_op_creator_2 }, s_op_annihilator_2, { s_g2_1122, s_g2_1212 } );
    calculate_g2( s, s_op_creator_2, { s_op_annihilator_2, s_op_annihilator_1 }, { s_op_creator_1, s_op_creator_2 }, s_op_annihilator_1, { s_g2_2121, s_g2_2211 } );
    // cache[s_g2_2211] = cache[s_g2_1122].conjugate();
    // cache[s_g2_2211 + "_time"] = cache[s_g2_1122 + "_time"];
    // cache[s_g2_1212] = cache[s_g2_2121].conjugate();
    // cache[s_g2_1212 + "_time"] = cache[s_g2_2121 + "_time"];

    // Note: This will probably be either removed completely, or implemented correctly.
    if ( s.parameters.input_correlation["Conc"].numerical_v["Center"].size() > 0 ) {
        cache["concurrence_total_" + fout] = cache[s_g2_1111] + cache[s_g2_1122] + cache[s_g2_1212] + cache[s_g2_1221] + cache[s_g2_2121] + cache[s_g2_2112] + cache[s_g2_2211] + cache[s_g2_2222];
        cache["concurrence_total_" + fout + "_time"] = cache[s_g2_1212 + "_time"];
        calculate_spectrum( s, "none", "none", s.parameters.input_correlation["Conc"].numerical_v["Center"].front(), s.parameters.input_correlation["Conc"].numerical_v["Range"].front(), s.parameters.input_correlation["Conc"].numerical_v["resW"].front(), 2, false, "concurrence_total_" + fout );
        for ( auto mode : { s_g2_1111, s_g2_1122, s_g2_1212, s_g2_1221, s_g2_2121, s_g2_2112, s_g2_2211, s_g2_2222 } ) {
            calculate_spectrum( s, "none", "none", s.parameters.input_correlation["Conc"].numerical_v["Center"].front(), s.parameters.input_correlation["Conc"].numerical_v["Range"].front(), s.parameters.input_correlation["Conc"].numerical_v["resW"].front(), 2, false, mode );
        }
    }

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

    for ( size_t t = 0; t < T; t++ ) {
        for ( auto mode : { s_g2_1111, s_g2_1122, s_g2_1212, s_g2_1221, s_g2_2121, s_g2_2112, s_g2_2211, s_g2_2222 } ) {
            rho[mode].emplace_back( 0 );
            rho_g2zero[mode].emplace_back( 0 );
        }
    }
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( auto mode : { s_g2_1111, s_g2_1122, s_g2_1212, s_g2_1221, s_g2_2121, s_g2_2112, s_g2_2211, s_g2_2222 } ) {
        auto &gmat_time = cache[mode + "_time"];
        for ( int upper_limit = 0; upper_limit < T; upper_limit++ ) {
            // G2(t,tau)
            rho[mode][upper_limit] = upper_limit > 0 ? rho[mode][upper_limit - 1] : 0.0;
            for ( int i = 0; i <= upper_limit; i++ ) {
                double dt = Numerics::get_tdelta( gmat_time, 0, i );
                int tau = upper_limit - i;
                // for ( int tau = 0; tau < T - i; tau++ ) {
                double dtau = Numerics::get_taudelta( gmat_time, i, tau ); // Numerics::get_tdelta( gmat_time, 0, i + tau );
                rho[mode][upper_limit] += cache[mode]( i, tau ) * dt * dtau;
            }
            // G2(t,0)
            double dt = Numerics::get_tdelta( gmat_time, 0, upper_limit );
            double t_t = std::real( gmat_time( upper_limit, 0 ) ); // Note: all correlation functions have to have the same times. cache[mode+"_time"] else.
            rho_g2zero[mode][upper_limit] = upper_limit > 0 ? rho_g2zero[mode][upper_limit - 1] + s.dgl_expectationvalue<Sparse, Scalar>( get_rho_at( rho_index_map[t_t] ), matmap_g2zero[mode], t_t ) * dt : s.dgl_expectationvalue<Sparse, Scalar>( get_rho_at( 0 ), matmap_g2zero[mode], get_time_at( 0 ) ) * dt;
            if ( mode == s_g2_1111 ) {
                timer_c.iterate();
                Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Concurrence (" + fout + "): " );
            }
        }
    }
    // Calculate EigenValues
    std::vector<Scalar> output( T, 0 );
    std::vector<Scalar> output_g2zero( T, 0 );
    std::vector<Scalar> output_simple( T, 0 );
    std::vector<Scalar> output_g2zero_simple( T, 0 );
    std::vector<Scalar> output_fidelity( T, 0 );
    std::vector<Scalar> output_fidelity_g2zero( T, 0 );
    std::vector<Scalar> time( T, 0 );
    std::vector<Dense> twophotonmatrix( T, Dense::Zero( 4, 4 ) );
    std::vector<Dense> twophotonmatrix_g2zero( T, Dense::Zero( 4, 4 ) );

    Dense spinflip = Dense::Zero( 4, 4 );
    spinflip( 0, 3 ) = -1;
    spinflip( 1, 2 ) = 1;
    spinflip( 2, 1 ) = 1;
    spinflip( 3, 0 ) = -1;
    Log::L2( "[Concurrence] Spinflip Matrix: {}\n", spinflip );
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t k = 0; k < T; k++ ) { // cache[s_g2_1111].rows() statt T?
        // Log::L3( "Creating 2 photon matrix\n" );
        Dense rho_2phot = Dense::Zero( 4, 4 );
        Dense rho_2phot_g2zero = Dense::Zero( 4, 4 );
        rho_2phot( 0, 0 ) = rho[s_g2_1111][k];
        rho_2phot( 3, 3 ) = rho[s_g2_2222][k]; // std::cout << "Setting rho(0,0) to " << rho_2phot( 0, 0 ) << ", Setting rho(3,3) to " << rho_2phot( 3, 3 ) << std::endl;
        rho_2phot( 0, 3 ) = rho[s_g2_2211][k];
        rho_2phot( 3, 0 ) = rho[s_g2_1122][k];
        rho_2phot( 1, 2 ) = rho[s_g2_2121][k];
        rho_2phot( 2, 1 ) = rho[s_g2_1212][k];
        rho_2phot( 1, 1 ) = std::abs( rho[s_g2_2112][k] ) != 0 ? rho[s_g2_2112][k] : 1E-100;
        rho_2phot( 2, 2 ) = std::abs( rho[s_g2_1221][k] ) != 0 ? rho[s_g2_1221][k] : 1E-100;
        rho_2phot = rho_2phot / std::abs( rho_2phot.trace() );
        rho_2phot_g2zero( 0, 0 ) = rho_g2zero[s_g2_1111][k];
        rho_2phot_g2zero( 3, 3 ) = rho_g2zero[s_g2_2222][k];
        rho_2phot_g2zero( 0, 3 ) = rho_g2zero[s_g2_2211][k];
        rho_2phot_g2zero( 3, 0 ) = rho_g2zero[s_g2_1122][k];
        rho_2phot_g2zero( 1, 2 ) = rho_g2zero[s_g2_2121][k];
        rho_2phot_g2zero( 2, 1 ) = rho_g2zero[s_g2_1212][k];
        rho_2phot_g2zero( 1, 1 ) = std::abs( rho_g2zero[s_g2_2112][k] ) != 0 ? rho_g2zero[s_g2_2112][k] : 1E-100;
        rho_2phot_g2zero( 2, 2 ) = std::abs( rho_g2zero[s_g2_1221][k] ) != 0 ? rho_g2zero[s_g2_1221][k] : 1E-100;
        rho_2phot_g2zero = rho_2phot_g2zero / std::abs( rho_2phot_g2zero.trace() );
        // Log::L3( "Normalizing 2 photon matrix\n" );
        // if ( std::abs( rho_2phot.trace() ) != 0 ) {
        // Log::L3( "Rho_2phot = {}\n", rho_2phot );
        // Log::L3( "Calculating sqrt(rho)\n" );
        // Eigen::MatrixPower<Dense> Mpow( rho_2phot );
        Dense sqrtrho2phot = rho_2phot.sqrt();               // Mpow( 0.5 );
        Dense sqrtrho2phot_g2zero = rho_2phot_g2zero.sqrt(); // Mpow( 0.5 );
        // Log::L3( "Calculating R\n" );
        Dense R = sqrtrho2phot * spinflip * rho_2phot * spinflip * sqrtrho2phot;
        Dense R_g2zero = sqrtrho2phot_g2zero * spinflip * rho_2phot_g2zero * spinflip * sqrtrho2phot_g2zero;
        // Log::L3( "R = {}\n", R );
        // Log::L3( "Calculating sqrt(R)\n" );
        // Eigen::MatrixPower<Dense> SMPow( R );
        Dense R5 = R.sqrt();               // SMPow( 0.5 );
        Dense R5_g2zero = R_g2zero.sqrt(); // SMPow( 0.5 );
        // Log::L3( "Calculating Eigenvalues\n" );
        // auto eigenvalues = R5.eigenvalues();
        // auto eigenvalues_g2zero = R5_g2zero.eigenvalues();
        Eigen::SelfAdjointEigenSolver<Dense> eigensolver( R5 );
        auto eigenvalues = eigensolver.eigenvalues();
        Eigen::SelfAdjointEigenSolver<Dense> eigensolver_g2zero( R5_g2zero );
        auto eigenvalues_g2zero = eigensolver_g2zero.eigenvalues();

        // Sometimes, the numerical method for eigenvalue evaluation yields crap (first element really big, rest shifted), so we check this here:
        // if ( QDLC::Math::abs2( eigenvalues( 3 ) ) > 50.0 ) {
        //    eigenvalues( 3 ) = eigenvalues( 2 );
        //    eigenvalues( 2 ) = eigenvalues( 1 );
        //    eigenvalues( 1 ) = eigenvalues( 0 );
        //    eigenvalues( 0 ) = 0.0;
        //}
        // if ( QDLC::Math::abs2( eigenvalues_g2zero( 3 ) ) > 50.0 ) {
        //    eigenvalues_g2zero( 3 ) = eigenvalues_g2zero( 2 );
        //    eigenvalues_g2zero( 2 ) = eigenvalues_g2zero( 1 );
        //    eigenvalues_g2zero( 1 ) = eigenvalues_g2zero( 0 );
        //    eigenvalues_g2zero( 0 ) = 0.0;
        //}
        // Log::L1( "rho2phot = {}\n\nsqrtrho2phot = {}\n\nR = {}\n\nRS = {}\nEigenvalues at t = {} are {}\n", rho_2phot, sqrtrho2phot, R, R5, get_time_at( i ), eigenvalues );
        auto conc = eigenvalues( 3 ) - eigenvalues( 2 ) - eigenvalues( 1 ) - eigenvalues( 0 );
        double fidelity = std::pow( std::real( R5.trace() ), 2.0 );
        // Log::L2( "Eigenvalues {} (size of vec: {}): C = {} - {} - {} - {}\n", k, eigenvalues.size(), eigenvalues( 3 ), eigenvalues( 2 ), eigenvalues( 1 ), eigenvalues( 0 ) );
        auto conc_g2zero = eigenvalues_g2zero( 3 ) - eigenvalues_g2zero( 2 ) - eigenvalues_g2zero( 1 ) - eigenvalues_g2zero( 0 );
        double fidelity_g2zero = std::pow( std::real( R5_g2zero.trace() ), 2.0 );
        output.at( k ) = conc;
        output_simple.at( k ) = 2.0 * std::abs( rho_2phot( 3, 0 ) / rho_2phot.trace() );
        output_fidelity.at( k ) = fidelity;
        output_g2zero.at( k ) = conc_g2zero;
        output_g2zero_simple.at( k ) = 2.0 * std::abs( rho_2phot_g2zero( 3, 0 ) / rho_2phot_g2zero.trace() );
        output_fidelity_g2zero.at( k ) = fidelity_g2zero;
        time.at( k ) = std::real( mat_time( k, 0 ) );
        // std::cout << "Rho(3,0) = "<<rho_2phot( 3, 0 )<<", rho.trace() = "<<rho_2phot.trace()<<", Rho before saving :\n" << rho_2phot.format(Eigen::IOFormat( 4, 0, ", ", "\n", "[", "]" )) << std::endl;
        // std::cout << "Rho_g20(3,0) = "<<rho_2phot_g2zero( 3, 0 )<<", rho_g20.trace() = "<<rho_2phot_g2zero.trace()<<", Rho_g20 before saving :\n" << rho_2phot_g2zero.format(Eigen::IOFormat( 4, 0, ", ", "\n", "[", "]" )) << std::endl;
        twophotonmatrix.at( k ) = rho_2phot;
        twophotonmatrix_g2zero.at( k ) = rho_2phot_g2zero;
        timer_c.iterate();
    }
    // Final output and timer end
    timer_c.end();
    Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Concurrence (" + fout + ")", Timers::PROGRESS_FORCE_OUTPUT );
    // Add to Fileoutput:
    if ( to_output["Conc"].size() == 0 )
        to_output["Conc"]["Time"] = time;
    to_output["Conc"][fout] = output;
    to_output["Conc_simple"][fout] = output_simple;
    to_output["Conc_fidelity"][fout] = output_fidelity;
    to_output["Conc_g2zero"][fout] = output_g2zero;
    to_output["Conc_g2zero_simple"][fout] = output_g2zero_simple;
    to_output["Conc_g2zero_fidelity"][fout] = output_fidelity_g2zero;
    to_output_m["TwoPMat"][fout] = twophotonmatrix;
    to_output_m["TwoPMat"][fout + "_g2zero"] = twophotonmatrix_g2zero;
    Log::L1( "Final Concurrence: {:.10f} ({:.10f} simple) {}\n", std::real( output.back() ), std::real( output_simple.back() ), fout );
    Log::L1( "Final Concurrence with g2(0): {:.10f} ({:.10f} simple) {}\n", std::real( output_g2zero.back() ), std::real( output_g2zero_simple.back() ), fout );
    return true;
}