#include "solver/solver_ode.h"

// Description: Calculates the Two Photon Density Matrix from a set of static G2 correlation functions

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
    return ret;
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
 * @brief
 * The Two Photon Matrix can be evaluated in several priority modes:
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
bool QDACC::Numerics::ODESolver::calculate_static_two_photon_matrix( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2, const int matrix_priority_evaluation, const double spec_center, const double spec_range, const double spec_res ) {
    Log::L2( "[TPM] TPM for modes {} {} and {} {} with priority {}\n", s_op_creator_1, s_op_creator_2, s_op_annihilator_1, s_op_annihilator_2, matrix_priority_evaluation );

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

    // Get MatrixMain Operator Matrices
    const auto &[op_creator_1, op_annihilator_1] = get_operators_matrices( s, s_op_creator_1, s_op_annihilator_1 );
    const auto &[op_creator_2, op_annihilator_2] = get_operators_matrices( s, s_op_creator_2, s_op_annihilator_2 );
    // Generate multiplied operator matrices for G2(0) evaluation
    std::map<std::string, MatrixMain> mode_matrix;
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
        Log::L2( "[TPM] Queued calculation of {}\n", std::accumulate( current_purpose.begin(), current_purpose.end(), std::string{}, []( const auto &a, const auto &b ) { return a + b + ", "; } ) );
        calculate_g2( s, current_creator_1, current_creator_2, current_annihilator_1, current_annihilator_2, current_purpose );
    }

    // TODONEXT: das hier triggerable --concFullEval (ALLE elemente der 2PM ausrehcnen!), eigenvalues ausgeben triggerable via --output.
    //  TODO: Make this triggerable with user flag. Assuming 2211 = 1122* is correct, but when evaluating it numerically, integration errors contribute to a significant degree.
    // cache[s_g2_2211] = cache[s_g2_1122].conjugate( s_g2_2211 );
    // cache[s_g2_1212] = cache[s_g2_2121].conjugate( s_g2_1212 );
    // Note: This will probably be either removed completely, or implemented correctly.
    auto &first_matrix = cache[mode_purpose.at( "1111" )];
    if ( spec_range > 0 ) {
        const auto matrix_dim = first_matrix.dim();

        const std::string combined_name = "tpm_total_" + fout;
        cache[combined_name] = MultidimensionalCacheMatrix( { matrix_dim, matrix_dim }, combined_name );
        auto &combined = cache[combined_name];

        // Sum Matrix and set Time Matrix
        for ( const auto &mode : modes )
            combined.addMatrix( cache[mode_purpose[mode]] );

        calculate_spectrum( s, "none", "none", spec_center, spec_range, (int)spec_res, 2, false, combined_name );
        for ( const auto &[n, mode] : mode_purpose ) {
            calculate_spectrum( s, "none", "none", spec_center, spec_range, (int)spec_res, 2, false, mode );
        }
    }

    // Maximum Time / Dimension of the cache matrices
    auto maximum_time = std::min<size_t>( first_matrix.dim(), savedStates.size() );

    // Generate DM cache matrices
    using tpdm_t = std::map<std::string, std::vector<Scalar>>;
    tpdm_t rho;
    tpdm_t rho_g2zero;
    // And fill them with zeros
    for ( const auto &mode : { "1111", "1211", "2111", "2211", "1112", "1212", "2112", "2212", "1121", "1221", "2121", "2221", "1122", "1222", "2122", "2222" } ) {
        rho[mode] = std::vector<Scalar>( maximum_time, 0 );
        rho_g2zero[mode] = std::vector<Scalar>( maximum_time, 0 );
    }

    int pbsize = 2 * maximum_time;
    Timer &timer_c = Timers::create( "Two Photon Matrix (" + fout + ")" ).start();

//#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( const auto &[mode, purpose] : mode_purpose ) {
        const auto &gmat = cache.at( purpose );
        // Triangular Integral over G2(t,tau)
        for ( size_t upper_limit = 0; upper_limit < maximum_time; upper_limit++ ) {
            // G2(upper_limit,tau)
            // Lower Triangle Integral -> Move to seperate function!
            rho[mode][upper_limit] = upper_limit > 0 ? rho[mode][upper_limit - 1] : 0.0;
            for ( size_t t = 0; t <= upper_limit; t++ ) {
                // Tau index
                size_t tau = upper_limit - t;
                double dt = s.getDeltaTimeOf( 0, { t, tau } );
                double dtau = s.getDeltaTimeOf( 1, { t, tau } );
                rho[mode][upper_limit] += gmat.get( t, tau ) * dt * dtau;
            }
            // G2(upper_limit,0) - No Triangular Integral
            double dt = s.getDeltaTimeOf( 0, { upper_limit, 0 } );
            double t_t = s.getTimeOf( 0, { upper_limit, 0 } ); // std::real( gmat_time( upper_limit, 0 ) ); // Note: all correlation functions have to have the same times. cache[purpose+"_time"] else.
            const auto time_index = rho_index_map[t_t];
            const auto &current_state = savedStates.at( time_index );
            rho_g2zero[mode][upper_limit] = upper_limit > 0 ? rho_g2zero[mode][upper_limit - 1] + s.dgl_expectationvalue<MatrixMain>( current_state.mat, mode_matrix[mode], t_t ) * dt : s.dgl_expectationvalue<MatrixMain>( savedStates.at( 0 ).mat, mode_matrix[mode], savedStates.at( 0 ).t ) * dt;
            if ( mode == "1111" ) {
                timer_c.iterate();
                Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Two Photon Matrix (" + fout + "): " );
            }
        }
    }

    // Generate Referenced Shortcuts
    auto &time = add_to_output( "TwoPMat", "Time", std::vector<Scalar>( maximum_time, 0 ), to_output, true /* Overwrite Existing */ );
    auto &twophotonmatrix = add_to_output( "TwoPMat", fout, std::vector<Dense>( maximum_time, Dense::Zero( 4, 4 ) ), to_output_m );
    auto &twophotonmatrix_g2zero = add_to_output( "TwoPMat", fout + "_g2zero", std::vector<Dense>( maximum_time, Dense::Zero( 4, 4 ) ), to_output_m );

    // Calculate two-photon densitymatrices
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t k = 0; k < maximum_time; k++ ) {
        // Neumann-Iterated TPDM
        const auto rho_2phot = _rho_to_tpdm( k, rho );

        // G2(0) TPDM
        const auto rho_2phot_g2zero = _rho_to_tpdm( k, rho_g2zero );

        // Cache to output arrays
        time.at( k ) = s.getTimeOf( 0, { k, 0 } );
        twophotonmatrix.at( k ) = rho_2phot;
        twophotonmatrix_g2zero.at( k ) = rho_2phot_g2zero;
        timer_c.iterate();
    }
    // Final output and timer end
    timer_c.end();
    Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Two Photon Matrix (" + fout + ")", Timers::PROGRESS_FORCE_OUTPUT );

    return true;
}