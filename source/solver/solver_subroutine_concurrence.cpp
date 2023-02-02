#include "solver/solver_ode.h"
#include "solver/solver_analytical_eigenvalues.h"
// Description: Calculates Concurrence

static Dense _generate_spinflip() {
    Dense spinflip = Dense::Zero(4,4);
    spinflip << 0, 0, 0, -1,
        0, 0, 1, 0,
        0, 1, 0, 0,
        -1, 0, 0, 0;
    return spinflip;
}

static Dense _rho_to_tpdm( const size_t i, const std::map<std::string, std::vector<Scalar>> &rho,
                          const std::string &a11, const std::string &a12, const std::string &a13, const std::string &a14,
                          const std::string &a21, const std::string &a22, const std::string &a23, const std::string &a24,
                          const std::string &a31, const std::string &a32, const std::string &a33, const std::string &a34,
                          const std::string &a41, const std::string &a42, const std::string &a43, const std::string &a44 ) {
    Dense ret = Dense::Zero(4,4);
    ret << rho.at( a11 )[i], rho.at( a12 )[i], rho.at( a13 )[i], rho.at( a14 )[i],
        rho.at( a21 )[i], rho.at( a22 )[i], rho.at( a23 )[i], rho.at( a24 )[i],
        rho.at( a31 )[i], rho.at( a32 )[i], rho.at( a33 )[i], rho.at( a34 )[i],
        rho.at( a41 )[i], rho.at( a42 )[i], rho.at( a43 )[i], rho.at( a44 )[i];
    return ret / std::abs( ret.trace() );
}

// TODO: Remodel, simplify, reuse
// TODO: approximations: outers, outers+inners, all, use_conjugated
bool QDLC::Numerics::ODESolver::calculate_concurrence( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2 ) {
    Log::L2( "[Concurrence] Conc for modes {} {} and {} {}\n", s_op_creator_1, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_2 );

    // Set Number of Phonon cores to 1 because this memberfunction is already using multithreading
    s.parameters.numerics_maximum_secondary_threads = 1;
    // Progress
    ProgressBar progressbar = ProgressBar();

    std::string fout = s_op_creator_1 + "-" + s_op_annihilator_1 + "-" + s_op_creator_2 + "-" + s_op_annihilator_2;

    // Calculate G2(t,tau) with given operator matrices
    // TODO: Make this a function
    std::string s_g2_1111 = get_operators_purpose( { s_op_creator_1, s_op_creator_1, s_op_annihilator_1, s_op_annihilator_1 } );
    std::string s_g2_1122 = get_operators_purpose( { s_op_creator_1, s_op_creator_1, s_op_annihilator_2, s_op_annihilator_2 } );
    std::string s_g2_1212 = get_operators_purpose( { s_op_creator_1, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_2 } );
    std::string s_g2_1221 = get_operators_purpose( { s_op_creator_1, s_op_creator_2, s_op_annihilator_2, s_op_annihilator_1 } );
    std::string s_g2_2121 = get_operators_purpose( { s_op_creator_2, s_op_creator_1, s_op_annihilator_2, s_op_annihilator_1 } );
    std::string s_g2_2112 = get_operators_purpose( { s_op_creator_2, s_op_creator_1, s_op_annihilator_1, s_op_annihilator_2 } );
    std::string s_g2_2211 = get_operators_purpose( { s_op_creator_2, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_1 } );
    std::string s_g2_2222 = get_operators_purpose( { s_op_creator_2, s_op_creator_2, s_op_annihilator_2, s_op_annihilator_2 } );

    const std::vector modes = { s_g2_1111, s_g2_1122, s_g2_1212, s_g2_1221, s_g2_2121, s_g2_2112, s_g2_2211, s_g2_2222 };


    // Calculate G2's if neccessary
    // Make sure to only evaluate the tau-direction when absolutely neccessary
    // calculate_g2( s, s_op_creator_1, s_op_creator_1,s_op_annihilator_1,  s_op_annihilator_1, { s_g2_1111 } );
    // calculate_g2( s, s_op_creator_1, s_op_creator_2,s_op_annihilator_2,  s_op_annihilator_1, { s_g2_1221 } );
    // calculate_g2( s, s_op_creator_2, s_op_creator_1,s_op_annihilator_1,  s_op_annihilator_2, { s_g2_2112 } );
    // calculate_g2( s, s_op_creator_2, s_op_creator_2,s_op_annihilator_2,  s_op_annihilator_2, { s_g2_2222 } );
    // calculate_g2( s, s_op_creator_2, s_op_creator_2,s_op_annihilator_1,  s_op_annihilator_1, { s_g2_2211 } );
    // calculate_g2( s, s_op_creator_1, s_op_creator_2,s_op_annihilator_1,  s_op_annihilator_2, { s_g2_1212 } );
    
    calculate_g2( s, s_op_creator_1, { s_op_creator_1, s_op_creator_2 }, { s_op_annihilator_1, s_op_annihilator_2 }, s_op_annihilator_1, { s_g2_1111, s_g2_1221 } );
    calculate_g2( s, s_op_creator_2, { s_op_creator_1, s_op_creator_2 }, { s_op_annihilator_1, s_op_annihilator_2 }, s_op_annihilator_2, { s_g2_2112, s_g2_2222 } );
    calculate_g2( s, s_op_creator_1, { s_op_creator_1, s_op_creator_2 }, { s_op_annihilator_2, s_op_annihilator_1 }, s_op_annihilator_2, { s_g2_1122, s_g2_1212 } );
    calculate_g2( s, s_op_creator_2, { s_op_creator_1, s_op_creator_2 }, { s_op_annihilator_2, s_op_annihilator_1 }, s_op_annihilator_1, { s_g2_2121, s_g2_2211 } );
    
    // this is like in seidelamann
    //calculate_g2( s, s_op_creator_1, { s_op_creator_1, s_op_creator_2 }, { s_op_annihilator_1, s_op_annihilator_2 }, s_op_annihilator_1, { s_g2_1111, s_g2_1212 } );
    //calculate_g2( s, s_op_creator_2, { s_op_creator_1, s_op_creator_2 }, { s_op_annihilator_1, s_op_annihilator_2 }, s_op_annihilator_2, { s_g2_2121, s_g2_2222 } );
    //calculate_g2( s, s_op_creator_1, { s_op_creator_1, s_op_creator_2 }, { s_op_annihilator_2, s_op_annihilator_1 }, s_op_annihilator_2, { s_g2_1122, s_g2_1221 } );
    //calculate_g2( s, s_op_creator_2, { s_op_creator_1, s_op_creator_2 }, { s_op_annihilator_2, s_op_annihilator_1 }, s_op_annihilator_1, { s_g2_2112, s_g2_2211 } );
    
    // TODONEXT: das hier triggerable --concFullEval (ALLE elemente der 2PM ausrehcnen!), eigenvalues ausgeben triggerable via --output.
    //  TODO: Make this triggerable with user flag. Assuming 2211 = 1122* is correct, but when evaluating it numerically, integration errors contribute to a significant degree.
    // cache[s_g2_2211] = cache[s_g2_1122].conjugate( s_g2_2211 );
    // cache[s_g2_1212] = cache[s_g2_2121].conjugate( s_g2_1212 );

    // Note: This will probably be either removed completely, or implemented correctly.
    if ( s.parameters.input_correlation["Conc"].property_set["Center"].size() > 0 ) {
        const Dense combined = cache[s_g2_1111].get() + cache[s_g2_1122].get() + cache[s_g2_1212].get() + cache[s_g2_1221].get() + cache[s_g2_2121].get() + cache[s_g2_2112].get() + cache[s_g2_2211].get() + cache[s_g2_2222].get();
        const std::string combined_name = "concurrence_total_" + fout;
        cache[combined_name] = CacheMatrix( combined, cache[s_g2_1111].get_time(), combined_name );
        calculate_spectrum( s, "none", "none", s.parameters.input_correlation["Conc"].property_set["Center"].front(), s.parameters.input_correlation["Conc"].property_set["Range"].front(), (int)s.parameters.input_correlation["Conc"].property_set["resW"].front(), 2, false, combined_name );
        for ( const auto &mode : modes ) {
            calculate_spectrum( s, "none", "none", s.parameters.input_correlation["Conc"].property_set["Center"].front(), s.parameters.input_correlation["Conc"].property_set["Range"].front(), (int)s.parameters.input_correlation["Conc"].property_set["resW"].front(), 2, false, mode );
        }
    }

    // Cache Progressbar Size
    int pbsize = 2 * cache[s_g2_1212].dim();

    // Get Sparse Operator Matrices
    const auto &[op_creator_1, op_annihilator_1] = get_operators_matrices( s, s_op_creator_1, s_op_annihilator_1 );
    const auto &[op_creator_2, op_annihilator_2] = get_operators_matrices( s, s_op_creator_2, s_op_annihilator_2 );
    // Generate multiplied operator matrices for G2(0) evaluation
    std::map<std::string, Sparse> matmap_g2zero = { { s_g2_1111, op_creator_1 * op_creator_1 * op_annihilator_1 * op_annihilator_1 },
                                                    { s_g2_1122, op_creator_1 * op_creator_1 * op_annihilator_2 * op_annihilator_2 },
                                                    { s_g2_1212, op_creator_1 * op_creator_2 * op_annihilator_1 * op_annihilator_2 },
                                                    { s_g2_1221, op_creator_1 * op_creator_2 * op_annihilator_2 * op_annihilator_1 },
                                                    { s_g2_2121, op_creator_2 * op_creator_1 * op_annihilator_2 * op_annihilator_1 },
                                                    { s_g2_2112, op_creator_2 * op_creator_1 * op_annihilator_1 * op_annihilator_2 },
                                                    { s_g2_2211, op_creator_2 * op_creator_2 * op_annihilator_1 * op_annihilator_1 },
                                                    { s_g2_2222, op_creator_2 * op_creator_2 * op_annihilator_2 * op_annihilator_2 } };

    Timer &timer_c = Timers::create( "Concurrence (" + fout + ")" ).start();

    // Maximum Time / Dimension of the cache matrices
    auto T = std::min<size_t>( cache[s_g2_1111].dim(), savedStates.size() );

    // Generate DM cache matrices
    std::map<std::string, std::vector<Scalar>> rho;
    std::map<std::string, std::vector<Scalar>> rho_g2zero;
    // And fill them with zeros
    for ( const auto mode : modes ) {
        rho[mode] = std::vector<Scalar>( T, 0 );
        rho_g2zero[mode] = std::vector<Scalar>( T, 0 );
    }
    rho["zero"] = std::vector<Scalar>( T, 1E-100 );
    rho_g2zero["zero"] = std::vector<Scalar>( T, 1E-100 );

#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( auto mode : modes ) {
        auto &gmat = cache[mode];
        // Integrate from 0 to t=T. Running index is t = t
        for ( int t = 0; t < T; t++ ) {
            // G2(t,tau)
            // Lower Triangle Integral -> Move to seperate function!
            rho[mode][t] = t > 0 ? rho[mode][t - 1] : 0.0;
            for ( int i = 0; i <= t; i++ ) {
                double dt = gmat.dt( i );
                // Tau index is
                int tau = t - i;
                // for ( int tau = 0; tau < T - i; tau++ ) {
                double dtau = gmat.dtau( tau, i ); // FIXED: dtau instead of dt!
                rho[mode][t] += cache[mode]( i, tau ) * dt * dtau;
            }
            // G2(t,0) - No Triangular Integral
            double dt = gmat.dt( t );
            double t_t = gmat.t( t ); // std::real( gmat_time( t, 0 ) ); // Note: all correlation functions have to have the same times. cache[mode+"_time"] else.
            rho_g2zero[mode][t] = t > 0 ? rho_g2zero[mode][t - 1] + s.dgl_expectationvalue<Sparse, Scalar>( get_rho_at( rho_index_map[t_t] ), matmap_g2zero[mode], t_t ) * dt : s.dgl_expectationvalue<Sparse, Scalar>( get_rho_at( 0 ), matmap_g2zero[mode], get_time_at( 0 ) ) * dt;
            if ( mode == s_g2_1111 ) {
                timer_c.iterate();
                Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Concurrence (" + fout + "): " );
            }
        }
    }

    // Generate Cache Vectors
    to_output["Conc"]["Time"] = std::vector<Scalar>( T, 0 );
    to_output["Conc"][fout] = std::vector<Scalar>( T, 0 );
    to_output["Conc_simple"][fout] = std::vector<Scalar>( T, 0 );
    to_output["Conc_analytical"][fout] = std::vector<Scalar>( T, 0 );
    to_output["Conc_fidelity"][fout] = std::vector<Scalar>( T, 0 );
    to_output["Conc_g2zero"][fout] = std::vector<Scalar>( T, 0 );
    to_output["Conc_g2zero_simple"][fout] = std::vector<Scalar>( T, 0 );
    to_output["Conc_g2zero_fidelity"][fout] = std::vector<Scalar>( T, 0 );
    to_output_m["TwoPMat"][fout] = std::vector<Dense>( T, Dense::Zero( 4, 4 ) );
    to_output_m["TwoPMat"][fout + "_g2zero"] = std::vector<Dense>( T, Dense::Zero( 4, 4 ) );
    to_output_m["ConcEV"][fout + "_EV"] = std::vector<Dense>( T, Vector::Zero( 4 ) );

    // Generate Referenced Shortcuts
    auto &time = to_output["Conc"]["Time"];
    auto &output = to_output["Conc"][fout];
    auto &output_simple = to_output["Conc_simple"][fout];
    auto &output_analytical = to_output["Conc_analytical"][fout];
    auto &output_fidelity = to_output["Conc_fidelity"][fout];
    auto &output_g2zero = to_output["Conc_g2zero"][fout];
    auto &output_g2zero_simple = to_output["Conc_g2zero_simple"][fout];
    auto &output_fidelity_g2zero = to_output["Conc_g2zero_simple"][fout];
    auto &twophotonmatrix = to_output_m["TwoPMat"][fout];
    auto &twophotonmatrix_g2zero = to_output_m["TwoPMat"][fout + "_g2zero"];
    auto &oeigenvalues = to_output_m["ConcEV"][fout + "_EV"];

    // Generate const spinflip matrix
    const Dense spinflip = _generate_spinflip();

    // Calculate two-photon densitymatrices and calculate concurrence
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t k = 0; k < T; k++ ) {
        // for rho in rho, rho_g2zero
        auto rho_2phot = _rho_to_tpdm( k, rho,
                                 s_g2_1111, "zero", "zero", s_g2_2211,
                                 "zero", s_g2_1212, s_g2_2112, "zero",
                                 "zero", s_g2_1221, s_g2_2121, "zero",
                                 s_g2_1122, "zero", "zero", s_g2_2222 );

        auto rho_2phot_g2zero = _rho_to_tpdm( k, rho_g2zero,
                                 s_g2_1111, "zero", "zero", s_g2_2211,
                                 "zero", s_g2_1212, s_g2_2112, "zero",
                                 "zero", s_g2_1221, s_g2_2121, "zero",
                                 s_g2_1122, "zero", "zero", s_g2_2222 );

        // TODO: two modes: this and seidelmann (Different Types of Photon Entanglement from a Constantly Driven Quantum Emitter Inside a Cavity)
        Dense sqrtrho2phot = rho_2phot.sqrt();               // Mpow( 0.5 );
        Dense sqrtrho2phot_g2zero = rho_2phot_g2zero.sqrt(); // Mpow( 0.5 );

        Dense R = sqrtrho2phot * spinflip * rho_2phot * spinflip * sqrtrho2phot;
        Dense R_g2zero = sqrtrho2phot_g2zero * spinflip * rho_2phot_g2zero * spinflip * sqrtrho2phot_g2zero;
        ;
        Dense R5 = R.sqrt();               // SMPow( 0.5 );
        Dense R5_g2zero = R_g2zero.sqrt(); // SMPow( 0.5 );

        Eigen::SelfAdjointEigenSolver<Dense> eigensolver( R5 );
        auto eigenvalues = eigensolver.eigenvalues();
        oeigenvalues.at( k ) = eigenvalues;

        Eigen::SelfAdjointEigenSolver<Dense> eigensolver_g2zero( R5_g2zero );
        auto eigenvalues_g2zero = eigensolver_g2zero.eigenvalues();
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
        time.at( k ) = cache[s_g2_1111].t( k );
        const auto eigenvalues_analytical = analytical_eigenvalues( rho_2phot / rho_2phot.trace() );
        output_analytical.at( k ) = eigenvalues_analytical( 3 ) - eigenvalues_analytical( 2 ) - eigenvalues_analytical( 1 ) - eigenvalues_analytical( 0 );
        // std::cout << "Eigenvalues: " << eigenvalues_analytical.format( Eigen::IOFormat( 4, 0, ", ", " ", "[", "]" ) ) << std::endl;
        //  std::cout << "Rho(3,0) = "<<rho_2phot( 3, 0 )<<", rho.trace() = "<<rho_2phot.trace()<<", Rho before saving :\n" << rho_2phot.format(Eigen::IOFormat( 4, 0, ", ", "\n", "[", "]" )) << std::endl;
        //  std::cout << "Rho_g20(3,0) = "<<rho_2phot_g2zero( 3, 0 )<<", rho_g20.trace() = "<<rho_2phot_g2zero.trace()<<", Rho_g20 before saving :\n" << rho_2phot_g2zero.format(Eigen::IOFormat( 4, 0, ", ", "\n", "[", "]" )) << std::endl;
        twophotonmatrix.at( k ) = rho_2phot;
        twophotonmatrix_g2zero.at( k ) = rho_2phot_g2zero;
        timer_c.iterate();
    }
    // Final output and timer end
    timer_c.end();
    Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Concurrence (" + fout + ")", Timers::PROGRESS_FORCE_OUTPUT );
    Log::L1( "Final Concurrence: {:.10f} ({:.10f} simple) {}\n", std::real( output.back() ), std::real( output_simple.back() ), fout );
    Log::L1( "Final Concurrence with g2(0): {:.10f} ({:.10f} simple) {}\n", std::real( output_g2zero.back() ), std::real( output_g2zero_simple.back() ), fout );
    return true;
}