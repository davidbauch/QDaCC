#include "solver/solver_ode.h"

// Description: Calculates the Two Photon Density Matrix from a set of time bin G2 correlation functions
// This functions assumes the G2 correlation functions to be calculated already and stored in the cache.
// If the functions do not exist, this function will return.

/**
 * @brief Constructs the Density Matrix rho from rho_ij(t) element vectors
 *
 * @param i Current vector index
 * @param rho Map of rho matrix elements. Each entry contains a vector rho_ij(t), where t is the vector index
 * @return Dense
 */
static Dense _rho_to_tpdm( const size_t i, const std::map<std::string, std::vector<Scalar>> &rho ) {
    Dense ret = Dense::Zero( 4, 4 );
    ret << rho.at( "EEEE" )[i], rho.at( "ELEE" )[i], rho.at( "LEEE" )[i], rho.at( "LLEE" )[i],
        rho.at( "EEEL" )[i], rho.at( "ELEL" )[i], rho.at( "LEEL" )[i], rho.at( "LLEL" )[i],
        rho.at( "EELE" )[i], rho.at( "ELLE" )[i], rho.at( "LELE" )[i], rho.at( "LLLE" )[i],
        rho.at( "EELL" )[i], rho.at( "ELLL" )[i], rho.at( "LELL" )[i], rho.at( "LLLL" )[i];
    // Add a tiny number that is otherwise out of reach of the TPM elements to avoid real zeros.
    ret.array() += 1E-200;
    // Norm by trace.
    return ret / ret.trace();
}

/**
 * @brief
 * The Two Photon Matrix 
 *           EE        EL        LE        LL    
 *       ----------------------------------------
 *  EE   |  EEEE      ELEE      LEEE      LLEE  |
 *       |                                      |
 *  EL   |  EEEL      ELEL      LEEL      LLEL  |
 *       |                                      |
 *  LE   |  EELE      ELLE      LELE      LLLE  |
 *       |                                      |
 *  LL   |  EELL      ELLL      LELL      LLLL  |
 *       ----------------------------------------
 * is always fully evaluated. No approximations are made.
 *
 * @param s
 * @param purpose
 */
bool QDACC::Numerics::ODESolver::calculate_timebin_two_photon_matrix( System &s, const std::string& purpose ) {
    Log::L2( "[TPM] TPM for purpose: {}\n", purpose );

    // Progress
    auto progressbar = ProgressBar();

    // Calculate G2's if neccessary. If one of the G2s already exists, we assume that all of them do.
    if (not cache.contains(purpose + "_EEEE")) {
        Log::L1( "[TPM] Required G2 do not exist.\n" );
        return false;
    }

    auto& first_matrix = cache.at( purpose + "_EEEE" );

    // Maximum Time / Dimension of the cache matrices
    auto maximum_time = std::min<size_t>( first_matrix.dim(), savedStates.size() );

    // Generate DM cache matrices
    using tpdm_t = std::map<std::string, std::vector<Scalar>>;
    tpdm_t rho;
    // And fill them with zeros
    for ( const auto &mode : { "EEEE", "EEEL", "EELE", "EELL", "ELEE", "ELEL", "ELLE", "ELLL", "LEEE", "LEEL", "LELE", "LELL", "LLEE", "LLEL", "LLLE", "LLLL" } ) {
        rho[mode] = std::vector<Scalar>( maximum_time, 0 );
    }

    int pbsize = 2 * maximum_time;
    Timer &timer_c = Timers::create( "Two Photon Matrix (" + purpose + ")" ).start();

//#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( const auto &mode : { "EEEE", "EEEL", "EELE", "EELL", "ELEE", "ELEL", "ELLE", "ELLL", "LEEE", "LEEL", "LELE", "LELL", "LLEE", "LLEL", "LLLE", "LLLL" } ) {
        const auto &gmat = cache.at( purpose + "_" + mode );
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
            if ( mode == "EEEE" ) {
                timer_c.iterate();
                Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Two Photon Matrix (" + purpose + "): " );
            }
        }
    }

    // Generate Referenced Shortcuts
    auto &time = add_to_output( "TwoPMat", "Time", std::vector<Scalar>( maximum_time, 0 ), to_output, false /* Overwrite Existing */ );
    auto &twophotonmatrix = add_to_output( "TwoPMat", purpose, std::vector<Dense>( maximum_time, Dense::Zero( 4, 4 ) ), to_output_m );

    // Calculate two-photon densitymatrices
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t k = 0; k < maximum_time; k++ ) {
        // Neumann-Iterated TPDM
        const auto rho_2phot = _rho_to_tpdm( k, rho );

        // Cache to output arrays
        time.at( k ) = s.getTimeOf( 0, { k, 0 } );
        twophotonmatrix.at( k ) = rho_2phot;
        timer_c.iterate();
    }
    // Final output and timer end
    timer_c.end();
    Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Two Photon Matrix (" + purpose + ")", Timers::PROGRESS_FORCE_OUTPUT );

    return true;
}