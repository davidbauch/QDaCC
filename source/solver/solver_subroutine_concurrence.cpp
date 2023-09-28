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
static Dense _generate_spinflip() {
    Dense spinflip = Dense::Zero( 4, 4 );
    spinflip << 0, 0, 0, -1,
        0, 0, 1, 0,
        0, 1, 0, 0,
        -1, 0, 0, 0;
    return spinflip;
}

/**
 * @brief Constructs the Density Matrix rho from rho_ij(t) element vectors
 *
 * @param i Current vector index
 * @param rho Map of rho matrix elements. Each entry contains a vector rho_ij(t), where t is the vector index
 * @return Dense
 */
static Dense _rho_to_tpdm( const size_t i, const std::map<std::string, std::vector<Scalar>> &rho ) {
    Dense ret = Dense::Zero( 4, 4 );
    ret << rho.at( "1111" )[i], rho.at( "1211" )[i], rho.at( "2111" )[i], rho.at( "2211" )[i],
        rho.at( "1112" )[i], rho.at( "1212" )[i], rho.at( "2112" )[i], rho.at( "2212" )[i],
        rho.at( "1121" )[i], rho.at( "1221" )[i], rho.at( "2121" )[i], rho.at( "2221" )[i],
        rho.at( "1122" )[i], rho.at( "1222" )[i], rho.at( "2122" )[i], rho.at( "2222" )[i];
    // Add a tiny number that is otherwise out of reach of the TPM elements to avoid real zeros.
    ret.array() += 1E-200;
    // Norm by trace.
    return ret / ret.trace();
}

/**
 * @brief Generates any of the available creator-annihilator permutations of the two photon density matrix rho
 *
 * @tparam input_type Can be any type
 * @param input Input Mode. E.g. '1111' for rho_11, '1122' for rho_14, etc.
 * @param a1,a2,b1,b2 Creator-Annihilator operators
 * @return std::vector<input_type>
 */
template <typename input_type>
static std::vector<input_type> _get_fourway_permutation( const std::string &input, const input_type &a1, const input_type &a2, const input_type &b1, const input_type &b2 ) {
    return { input[0] == '1' ? a1 : a2, input[1] == '1' ? a1 : a2, input[2] == '1' ? b1 : b2, input[3] == '1' ? b1 : b2 };
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
        Log::L2("[Concurrence] Fidelity matrix is approximately real. Using General EigenSolver class.\n");
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

/**
 * @brief
 * The concurrence can be evaluated in several priority modes:
 *           HH        HV        VH        VV                         HH        HV        VH        VV
 *       ----------------------------------------                  -------------------------------------
 *  HH   |  HHHH      HVHH      VHHH      VVHH  |             HH   |  1         3         3         1  |
 *       |                                      |                  |                                   |
 *  HV   |  HHHV      HVHV      VHHV      VVHV  |   ==\       HV   |  3         2         2         3  |
 *       |                                      |   == )           |                                   |
 *  VH   |  HHVH      HVVH      VHVH      VVVH  |   ==/       VH   |  3         2         2         3  |
 *       |                                      |                  |                                   |
 *  VV   |  HHVV      HVVV      VHVV      VVVV  |             VV   |  1         3         3         1  |
 *       ----------------------------------------                  -------------------------------------
 * where for a given priority, all elements with a higher priority are set to zero.
 *
 * @param s
 * @param s_op_creator_1
 * @param s_op_annihilator_1
 * @param s_op_creator_2
 * @param s_op_annihilator_2
 */
bool QDACC::Numerics::ODESolver::calculate_concurrence( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2, const int matrix_priority_evaluation, const double spec_center, const double spec_range, const double spec_res ) {
    Log::L2( "[Concurrence] Conc for modes {} {} and {} {} with priority {}\n", s_op_creator_1, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_2, matrix_priority_evaluation );

    // Set Number of Phonon cores to 1 because this memberfunction is already using multithreading
    s.parameters.numerics_maximum_secondary_threads = 1;
    // Progress
    auto progressbar = ProgressBar();

    std::string fout = s_op_creator_1 + "-" + s_op_annihilator_1 + "-" + s_op_creator_2 + "-" + s_op_annihilator_2;
    // Calculate G2(t,tau) with given operator matrices
    std::vector modes = { "1111", "1122", "2211", "2222" };
    if ( matrix_priority_evaluation == 2 )
        modes.insert( modes.end(), { "1212", "1221", "2112", "2121" } );
    if ( matrix_priority_evaluation == 3 )
        modes.insert( modes.end(), { "1211", "2122", "2111", "1222", "1112", "2221", "1121", "2212" } );

    std::map<std::string, std::string> mode_purpose;
    for ( const auto &mode : modes ) {
        const auto permutation = _get_fourway_permutation( mode, s_op_creator_1, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_2 );
        mode_purpose[mode] = get_operators_purpose( permutation );
    }

    // Get Sparse Operator Matrices
    const auto &[op_creator_1, op_annihilator_1] = get_operators_matrices( s, s_op_creator_1, s_op_annihilator_1 );
    const auto &[op_creator_2, op_annihilator_2] = get_operators_matrices( s, s_op_creator_2, s_op_annihilator_2 );
    // Generate multiplied operator matrices for G2(0) evaluation
    std::map<std::string, Sparse> mode_matrix;
    for ( const auto &mode : modes ) {
        const auto permutation = _get_fourway_permutation( mode, op_creator_1, op_creator_2, op_annihilator_1, op_annihilator_2 );
        mode_matrix[mode] = permutation[0] * permutation[1] * permutation[2] * permutation[3];
    }

    // Calculate G2's if neccessary.
    // Make sure to only evaluate the tau-direction when absolutely neccessary
    for ( const auto &mode : { "1111", "1122", "2211", "2222" } ) {
        std::vector<std::string> current_purpose;
        std::vector<std::string> current_creator_2;
        std::vector<std::string> current_annihilator_1;
        // This could be done much smoother using the c++23 views and pipe features of ranges, but at this point the compiler is not yet ready for it
        // Hence, this is very inefficient because we skip a lot of modes, but its only done once so it doesnt really matter
        for ( const auto &other : modes ) {
            if ( other[0] != mode[0] or other[3] != mode[3] )
                continue;
            current_purpose.push_back( mode_purpose[other] );
            // Add all creator and annihilators to vetor
            current_creator_2.push_back( other[1] == '1' ? s_op_creator_1 : s_op_creator_2 );
            current_annihilator_1.push_back( other[2] == '1' ? s_op_annihilator_1 : s_op_annihilator_2 );
        }
        // Calculate G2 function
        auto current_creator_1 = mode[0] == '1' ? s_op_creator_1 : s_op_creator_2;
        auto current_annihilator_2 = mode[3] == '1' ? s_op_annihilator_1 : s_op_annihilator_2;
        // TODO: MPI this.
        Log::L2( "[Concurrence] Queued calculation of {}\n", std::accumulate( current_purpose.begin(), current_purpose.end(), std::string{}, []( const auto &a, const auto &b ) { return a + b + ", "; } ) );
        calculate_g2( s, current_creator_1, current_creator_2, current_annihilator_1, current_annihilator_2, current_purpose );
    }

    // TODONEXT: das hier triggerable --concFullEval (ALLE elemente der 2PM ausrehcnen!), eigenvalues ausgeben triggerable via --output.
    //  TODO: Make this triggerable with user flag. Assuming 2211 = 1122* is correct, but when evaluating it numerically, integration errors contribute to a significant degree.
    // cache[s_g2_2211] = cache[s_g2_1122].conjugate( s_g2_2211 );
    // cache[s_g2_1212] = cache[s_g2_2121].conjugate( s_g2_1212 );
    // Note: This will probably be either removed completely, or implemented correctly.
    if ( spec_range > 0 ) {
        const auto matrix_dim = cache[mode_purpose.at( "1111" )].dim();
        // Don't use accumulate here, because the Dense Matrix in cache is hiding behind the .get() call.
        Dense combined = Dense::Zero( matrix_dim, matrix_dim );
        for ( const auto &mode : modes ) {
            combined += cache[mode_purpose[mode]].get();
        }
        const std::string combined_name = "concurrence_total_" + fout;
        cache[combined_name] = CacheMatrix( combined, cache[mode_purpose.at( "1111" )].get_time(), combined_name );
        calculate_spectrum( s, "none", "none", spec_center, spec_range, (int)spec_res, 2, false, combined_name );
        for ( const auto &[n, mode] : mode_purpose ) {
            calculate_spectrum( s, "none", "none", spec_center, spec_range, (int)spec_res, 2, false, mode );
        }
    }

    // Cache Progressbar Size
    int pbsize = 2 * cache[mode_purpose.at( "1111" )].dim();

    Timer &timer_c = Timers::create( "Concurrence (" + fout + ")" ).start();

    // Maximum Time / Dimension of the cache matrices
    auto T = std::min<size_t>( cache[mode_purpose.at( "1111" )].dim(), savedStates.size() );

    // Generate DM cache matrices
    using tpdm_t = std::map<std::string, std::vector<Scalar>>;
    tpdm_t rho;
    tpdm_t rho_g2zero;
    // And fill them with zeros
    for ( const auto &mode : { "1111", "1211", "2111", "2211", "1112", "1212", "2112", "2212", "1121", "1221", "2121", "2221", "1122", "1222", "2122", "2222" } ) {
        rho[mode] = std::vector<Scalar>( T, 0 );
        rho_g2zero[mode] = std::vector<Scalar>( T, 0 );
    }

    // #pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( const auto &[mode, purpose] : mode_purpose ) {
        const auto &gmat = cache.at( purpose );
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
                rho[mode][t] += gmat( i, tau ) * dt * dtau;
            }
            // G2(t,0) - No Triangular Integral
            double dt = gmat.dt( t );
            double t_t = gmat.t( t ); // std::real( gmat_time( t, 0 ) ); // Note: all correlation functions have to have the same times. cache[purpose+"_time"] else.
            rho_g2zero[mode][t] = t > 0 ? rho_g2zero[mode][t - 1] + s.dgl_expectationvalue<Sparse, Scalar>( get_rho_at( rho_index_map[t_t] ), mode_matrix[mode], t_t ) * dt : s.dgl_expectationvalue<Sparse, Scalar>( get_rho_at( 0 ), mode_matrix[mode], get_time_at( 0 ) ) * dt;
            if ( mode == "1111" ) {
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
    auto &output_fidelity_g2zero = to_output["Conc_g2zero_fidelity"][fout];
    auto &twophotonmatrix = to_output_m["TwoPMat"][fout];
    auto &twophotonmatrix_g2zero = to_output_m["TwoPMat"][fout + "_g2zero"];
    auto &output_eigenvalues = to_output_m["ConcEV"][fout + "_EV"];

    // Generate spinflip matrix
    const Dense spinflip = _generate_spinflip();

    const int fidelity_method = 0; // 0 Wootters, 1 Seidelmann
    // Calculate two-photon densitymatrices and calculate concurrence
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t k = 0; k < T; k++ ) {
        // Neumann-Iterated TPDM
        const auto rho_2phot = _rho_to_tpdm( k, rho );
        const auto fidelity_matrix = fidelity_method ? _fidelity_matrix_seidelmann( rho_2phot, spinflip ) : _fidelity_matrix_wootters( rho_2phot, spinflip );
        auto eigenvalues = _concurrence_eigenvalues( fidelity_matrix );
        if ( fidelity_method )
            eigenvalues = eigenvalues.cwiseSqrt().eval();
        auto conc = eigenvalues( 0 ) - eigenvalues( 1 ) - eigenvalues( 2 ) - eigenvalues( 3 );
        double fidelity = std::pow( std::real( fidelity_matrix.trace() ), 2.0 );

        // G2(0) TPDM
        const auto rho_2phot_g2zero = _rho_to_tpdm( k, rho_g2zero );
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
        time.at( k ) = cache[mode_purpose.at( "1111" )].t( k );
        output_analytical.at( k ) = conc_analytical;
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