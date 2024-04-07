#include "solver/solver_ode.h"
#include "solver/solver_analytical_eigenvalues.h"
// Description: Calculates Concurrence

/**
 * @brief Generating the Pauli matrices sigma_+ * sigma_+^t
 * The Spinmatrix
 * | ..  i |
 * | -i .. |
 * is tensored with itself to
 * | .. .. .. -1 |
 * | .. ..  1 .. |
 * | ..  1 .. .. |
 * | -1 .. .. .. |
 *
 * @return Dense
 */
Dense _generate_spinflip() {
    Dense spinflip = Dense::Zero( 4, 4 );
    spinflip << 0, 0, 0, -1,
        0, 0, 1, 0,
        0, 1, 0, 0,
        -1, 0, 0, 0;
    return spinflip;
}

/**
 * @brief Calculates the fidelity matrix F = sqrt(rho) * spinflip * rho.conjugated * spinflip * sqrt(rho)
 * according to https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.80.2245
 * @param rho Input two photon density matrix
 * @param spinflip Spinflip matrix
 * @return Dense
 */
Dense _fidelity_matrix_wootters( const Dense &rho, const Dense &spinflip ) {
    const Dense sqrt = rho.sqrt();
    const Dense R = sqrt * spinflip * rho.conjugate() * spinflip * sqrt;
    const Dense R5 = R.sqrt();
    return R5;
}

/**
 * @brief Calculates the fidelity matrix F = rho * spinflip * rho.conjugated * spinflip
 * according to https://onlinelibrary.wiley.com/doi/full/10.1002/qute.202000108
 * @param rho Input two photon density matrix
 * @param spinflip Spinflip matrix
 * @return Dense
 */
Dense _fidelity_matrix_seidelmann( const Dense &rho, const Dense &spinflip ) {
    const Dense M = rho * spinflip * rho.conjugate() * spinflip;
    return M;
}

/**
 * @brief Calculates the concurrence eigenvalues from the fidelity matrix
 * @param fidelity_matrix Fidelity matrix
 * @return Dense
 */
Dense _concurrence_eigenvalues( const Dense &fidelity_matrix ) {
    // check if fidelity_matrix is approximately real:
    const double threshold = 1E-14;
    const bool is_approximately_real = fidelity_matrix.real().isApprox( fidelity_matrix, threshold );
    Vector eigenvalues;
    if ( is_approximately_real ) {
        Log::L2( "[Concurrence] Fidelity matrix is approximately real. Using General EigenSolver class.\n" );
        dDense fidelity_matrix_real = fidelity_matrix.real();
        Eigen::EigenSolver<dDense> eigensolver( fidelity_matrix_real );
        eigenvalues = eigensolver.eigenvalues().eval();
    } else {
        Eigen::ComplexEigenSolver<Dense> eigensolver( fidelity_matrix );
        eigenvalues = eigensolver.eigenvalues().eval();
    }
    // The eigenvalues have to be real, positive and sorted. To ensure this, we take the absolute value and sort them in descending order.
    std::ranges::sort( eigenvalues, []( const auto &a, const auto &b ) { return std::real( a ) > std::real( b ); } ); // std::greater<double>()
    return eigenvalues;
}

bool QDACC::Numerics::ODESolver::calculate_concurrence( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2, const std::string& method ) {
    Log::L2( "[Concurrence] Conc for modes {} {} and {} {}\n", s_op_creator_1, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_2 );

    // Progress
    auto progressbar = ProgressBar();

    // Generate the same name that is used in the TPM function. Maybe switch to a key based parameter instead.
    std::string fout = s_op_creator_1 + "-" + s_op_annihilator_1 + "-" + s_op_creator_2 + "-" + s_op_annihilator_2;
    // This is most certainly never the case, but we may as well make sure.
    if ( not to_output_m["TwoPMat"].contains( fout ) ) {
        Log::L2( "[Concurrence] Required TPM is not in cache!\n" );
        return false;
    }

    // Cache Progressbar Size
    auto cache_size = to_output_m["TwoPMat"][fout].size();
    int pbsize = 2 * cache_size;
    Timer &timer_c = Timers::create( "Concurrence (" + fout + ")" ).start();

    // Maximum Time / Dimension of the cache matrices
    auto maximum_time = std::min<size_t>( cache_size, savedStates.size() );

    // Generate Referenced Shortcuts
    auto &output = add_to_output( "Conc", fout, std::vector<Scalar>( maximum_time, 0 ), to_output );
    auto &output_simple = add_to_output( "Conc_simple", fout, std::vector<Scalar>( maximum_time, 0 ), to_output );
    auto &output_analytical = add_to_output( "Conc_analytical", fout, std::vector<Scalar>( maximum_time, 0 ), to_output );
    auto &output_fidelity = add_to_output( "Conc_fidelity", fout, std::vector<Scalar>( maximum_time, 0 ), to_output );
    auto &output_g2zero = add_to_output( "Conc_g2zero", fout, std::vector<Scalar>( maximum_time, 0 ), to_output );
    auto &output_g2zero_simple = add_to_output( "Conc_g2zero_simple", fout, std::vector<Scalar>( maximum_time, 0 ), to_output );
    auto &output_fidelity_g2zero = add_to_output( "Conc_g2zero_fidelity", fout, std::vector<Scalar>( maximum_time, 0 ), to_output );
    auto &output_eigenvalues = add_to_output( "ConcEV", fout + "_EV", std::vector<Dense>( maximum_time, Vector::Zero( 4 ) ), to_output_m );
    // Get the TPM References
    auto &twophotonmatrix = to_output_m["TwoPMat"][fout];
    auto &twophotonmatrix_g2zero = to_output_m["TwoPMat"][fout + "_g2zero"];
    // Generate spinflip matrix
    const Dense spinflip = _generate_spinflip();

    const int fidelity_method = method == "wootters" ? 0 : 1; // 0 Wootters, 1 Seidelmann
    Log::L2( "[Concurrence] Using {} method.\n", fidelity_method == 0 ? "Wootters" : "Seidelmann" );

    // Calculate two-photon densitymatrices and calculate concurrence
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t k = 0; k < maximum_time; k++ ) {
        // Neumann-Iterated TPDM
        auto rho_2phot = twophotonmatrix[k];
        // Add little bit of noise to avoid numerical issues
        rho_2phot.array() += 1E-200;
        // Normalize
        rho_2phot /= rho_2phot.trace();

        const auto fidelity_matrix = fidelity_method ? _fidelity_matrix_seidelmann( rho_2phot, spinflip ) : _fidelity_matrix_wootters( rho_2phot, spinflip );
        auto eigenvalues = _concurrence_eigenvalues( fidelity_matrix );
        if ( fidelity_method )
            eigenvalues = eigenvalues.cwiseSqrt().eval();
        auto conc = eigenvalues( 0 ) - eigenvalues( 1 ) - eigenvalues( 2 ) - eigenvalues( 3 );
        double fidelity = std::pow( std::real( fidelity_matrix.trace() ), 2.0 );

        // G2(0) TPDM
        auto rho_2phot_g2zero = twophotonmatrix_g2zero[k];
        // Add little bit of noise to avoid numerical issues
        rho_2phot_g2zero.array() += 1E-200;
        // Normalize
        rho_2phot_g2zero /= rho_2phot_g2zero.trace();

        const auto fidelity_matrix_g2zero = fidelity_method ? _fidelity_matrix_seidelmann( rho_2phot_g2zero, spinflip ) : _fidelity_matrix_wootters( rho_2phot_g2zero, spinflip );
        auto eigenvalues_g2zero = _concurrence_eigenvalues( fidelity_matrix_g2zero );
        if ( fidelity_method )
            eigenvalues_g2zero = eigenvalues_g2zero.cwiseSqrt().eval();
        auto conc_g2zero = eigenvalues_g2zero( 0 ) - eigenvalues_g2zero( 1 ) - eigenvalues_g2zero( 2 ) - eigenvalues_g2zero( 3 );
        double fidelity_g2zero = std::pow( std::real( fidelity_matrix_g2zero.trace() ), 2.0 );

        // Analytical Eigenvalues.
        auto eigenvalues_analytical = analytical_eigenvalues( rho_2phot );
        std::ranges::sort( eigenvalues_analytical, []( const auto &a, const auto &b ) { return std::real( a ) > std::real( b ); } );
        const auto conc_analytical = eigenvalues_analytical( 0 ) - eigenvalues_analytical( 1 ) - eigenvalues_analytical( 2 ) - eigenvalues_analytical( 3 );

        // Cache to output arrays
        output_eigenvalues.at( k ) = eigenvalues;
        output.at( k ) = conc;
        output_simple.at( k ) = std::abs( rho_2phot( 0, 3 ) ) + std::abs( rho_2phot( 3, 0 ) );
        output_fidelity.at( k ) = fidelity;
        output_g2zero.at( k ) = conc_g2zero;
        output_g2zero_simple.at( k ) = std::abs( rho_2phot_g2zero( 0, 3 ) ) + std::abs( rho_2phot_g2zero( 3, 0 ) );
        output_fidelity_g2zero.at( k ) = fidelity_g2zero;
        output_analytical.at( k ) = conc_analytical;
        timer_c.iterate();
    }
    // Final output and timer end
    timer_c.end();
    Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Concurrence (" + fout + ")", Timers::PROGRESS_FORCE_OUTPUT );

    Log::L1( "Final Concurrence: {:.10f} ({:.10f} simple) {}\n", std::real( output.back() ), std::real( output_simple.back() ), fout );
    Log::L1( "Final Concurrence with g2(0): {:.10f} ({:.10f} simple) {}\n", std::real( output_g2zero.back() ), std::real( output_g2zero_simple.back() ), fout );

    return true;
}