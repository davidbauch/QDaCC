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
    Dense ret = Dense::Zero( 8, 8 );
    ret << rho.at( "EEEEEE" )[i], rho.at( "EEEEEL" )[i], rho.at( "EEEELE" )[i], rho.at( "EEEELL" )[i], rho.at( "EEELLE" )[i], rho.at( "EEELEE" )[i], rho.at( "EEELEL" )[i], rho.at( "EEELLL" )[i], rho.at( "EELEEE" )[i],
        rho.at( "EELEEL" )[i], rho.at( "EELELE" )[i], rho.at( "EELELL" )[i], rho.at( "EELLLE" )[i], rho.at( "EELLEE" )[i], rho.at( "EELLEL" )[i], rho.at( "EELLLL" )[i], rho.at( "ELEEEE" )[i], rho.at( "ELEEEL" )[i],
        rho.at( "ELEELE" )[i], rho.at( "ELEELL" )[i], rho.at( "ELELLE" )[i], rho.at( "ELELEE" )[i], rho.at( "ELELEL" )[i], rho.at( "ELELLL" )[i], rho.at( "ELLEEE" )[i], rho.at( "ELLEEL" )[i], rho.at( "ELLELE" )[i],
        rho.at( "ELLELL" )[i], rho.at( "ELLLLE" )[i], rho.at( "ELLLEE" )[i], rho.at( "ELLLEL" )[i], rho.at( "ELLLLL" )[i], rho.at( "LLEEEE" )[i], rho.at( "LLEEEL" )[i], rho.at( "LLEELE" )[i], rho.at( "LLEELL" )[i],
        rho.at( "LLELLE" )[i], rho.at( "LLELEE" )[i], rho.at( "LLELEL" )[i], rho.at( "LLELLL" )[i], rho.at( "LEEEEE" )[i], rho.at( "LEEEEL" )[i], rho.at( "LEEELE" )[i], rho.at( "LEEELL" )[i], rho.at( "LEELLE" )[i],
        rho.at( "LEELEE" )[i], rho.at( "LEELEL" )[i], rho.at( "LEELLL" )[i], rho.at( "LELEEE" )[i], rho.at( "LELEEL" )[i], rho.at( "LELELE" )[i], rho.at( "LELELL" )[i], rho.at( "LELLLE" )[i], rho.at( "LELLEE" )[i],
        rho.at( "LELLEL" )[i], rho.at( "LELLLL" )[i], rho.at( "LLLEEE" )[i], rho.at( "LLLEEL" )[i], rho.at( "LLLELE" )[i], rho.at( "LLLELL" )[i], rho.at( "LLLLLE" )[i], rho.at( "LLLLEE" )[i], rho.at( "LLLLEL" )[i],
        rho.at( "LLLLLL" )[i];

    // Add a tiny number that is otherwise out of reach of the TPM elements to avoid real zeros.
    ret.array() += 1E-200;
    // Norm by trace.
    return ret / ret.trace();
}

/**
 * @brief
 * The Two Photon Matrix 
 *              EEE    EEL    ELE    LEE    ELL    LEL    LLE    LLL
 *          ---------------------------------------------------------------
 *   EEE    |   EEEEEE EEEEEL EEEELE EEELEE EEEELL EEELEL EEELLE EEELLL   |
 *   EEL    |   EELEEE EELEEL EELELE EELLEE EELELL EELLEL EELLLE EELLLL   |
 *   ELE    |   ELEEEE ELEEEL ELEELE ELELEE ELEELL ELELEL ELELLE ELELLL   |
 *   LEE    |   LEEEEE LEEEEL LEEELE LEELEE LEEELL LEELEL LEELLE LEELLL   | 
 *   ELL    |   ELLEEE ELLEEL ELLELE ELLLEE ELLELL ELLLEL ELLLLE ELLLLL   | 
 *   LEL    |   LELEEE LELEEL LELELE LELLEE LELELL LELLEL LELLLE LELLLL   | 
 *   LLE    |   LLEEEE LLEEEL LLEELE LLELEE LLEELL LLELEL LLELLE LLELLL   | 
 *   LLL    |   LLLEEE LLLEEL LLLELE LLLLEE LLLELL LLLLEL LLLLLE LLLLLL   |
 *          ---------------------------------------------------------------
 * is always fully evaluated. No approximations are made.
 *
 * @param s
 * @param purpose
 */
bool QDACC::Numerics::ODESolver::calculate_timebin_three_photon_matrix( System &s, const std::string &purpose ) {
    Log::L2( "[TPM] TPM for purpose: {}\n", purpose );

    // Progress
    auto progressbar = ProgressBar();

    // Calculate G2's if neccessary. If one of the G2s already exists, we assume that all of them do.
    if ( not cache.contains( purpose + "_EEEEEE" ) ) {
        Log::L1( "[TPM] Required G3 do not exist.\n" );
        return false;
    }

    auto &first_matrix = cache.at( purpose + "_EEEEEE" );

    // Maximum Time / Dimension of the cache matrices
    auto maximum_time = std::min<size_t>( first_matrix.dim(), savedStates.size() );

    // Generate DM cache matrices
    using tpdm_t = std::map<std::string, std::vector<Scalar>>;
    tpdm_t rho;
    // And fill them with zeros
    for ( const auto &mode : {
              "EEEEEE", "EEEEEL", "EEEELE", "EEEELL", "EEELLE", "EEELEE", "EEELEL", "EEELLL", "EELEEE", "EELEEL", "EELELE", "EELELL", "EELLLE", "EELLEE", "EELLEL", "EELLLL",
              "ELEEEE", "ELEEEL", "ELEELE", "ELEELL", "ELELLE", "ELELEE", "ELELEL", "ELELLL", "ELLEEE", "ELLEEL", "ELLELE", "ELLELL", "ELLLLE", "ELLLEE", "ELLLEL", "ELLLLL",
              "LLEEEE", "LLEEEL", "LLEELE", "LLEELL", "LLELLE", "LLELEE", "LLELEL", "LLELLL", "LEEEEE", "LEEEEL", "LEEELE", "LEEELL", "LEELLE", "LEELEE", "LEELEL", "LEELLL",
              "LELEEE", "LELEEL", "LELELE", "LELELL", "LELLLE", "LELLEE", "LELLEL", "LELLLL", "LLLEEE", "LLLEEL", "LLLELE", "LLLELL", "LLLLLE", "LLLLEE", "LLLLEL", "LLLLLL",
          } ) {
        rho[mode] = std::vector<Scalar>( maximum_time, 0 );
    }

    int pbsize = 2 * maximum_time;
    Timer &timer_c = Timers::create( "Three Photon Matrix (" + purpose + ")" ).start();

    //#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( const auto &mode : { "EEEEEE", "EEEEEL", "EEEELE", "EEEELL", "EEELLE", "EEELEE", "EEELEL", "EEELLL", "EELEEE", "EELEEL", "EELELE", "EELELL", "EELLLE", "EELLEE", "EELLEL", "EELLLL",
                               "ELEEEE", "ELEEEL", "ELEELE", "ELEELL", "ELELLE", "ELELEE", "ELELEL", "ELELLL", "ELLEEE", "ELLEEL", "ELLELE", "ELLELL", "ELLLLE", "ELLLEE", "ELLLEL", "ELLLLL",
                               "LLEEEE", "LLEEEL", "LLEELE", "LLEELL", "LLELLE", "LLELEE", "LLELEL", "LLELLL", "LEEEEE", "LEEEEL", "LEEELE", "LEEELL", "LEELLE", "LEELEE", "LEELEL", "LEELLL",
                               "LELEEE", "LELEEL", "LELELE", "LELELL", "LELLLE", "LELLEE", "LELLEL", "LELLLL", "LLLEEE", "LLLEEL", "LLLELE", "LLLELL", "LLLLLE", "LLLLEE", "LLLLEL", "LLLLLL" } ) {
        const auto &gmat = cache.at( purpose + "_" + mode );
        // Triangular Integral over G3(t,tau)
        for ( size_t upper_limit = 0; upper_limit < maximum_time; upper_limit++ ) {
            // G3(upper_limit,tau1,tau2)
            rho[mode][upper_limit] = upper_limit > 0 ? rho[mode][upper_limit - 1] : 0.0;
            for ( size_t t = 0; t <= upper_limit; t++ ) {
                // Tau index
                size_t tau1 = upper_limit - t;
                for ( size_t t2 = 0; t2 <= tau1; t2++ ) {
                    size_t tau2 = tau1 - t2;
                    double dt = s.getDeltaTimeOf( 0, { t, tau1, tau2 } );
                    double dtau = s.getDeltaTimeOf( 1, { t, tau1, tau2 } ) * s.getDeltaTimeOf( 2, { t, tau1, tau2 } );
                    rho[mode][upper_limit] += gmat.get( {t, tau1, tau2} ) * dt * dtau;
                }
            }
            // G2(upper_limit,0) - No Triangular Integral
            double dt = s.getDeltaTimeOf( 0, { upper_limit, 0, 0 } );
            double t_t = s.getTimeOf( 0, { upper_limit, 0, 0 } );
            const auto time_index = rho_index_map[t_t];
            const auto &current_state = savedStates.at( time_index );
            if ( mode == "EEEEEE" ) {
                timer_c.iterate();
                Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Two Photon Matrix (" + purpose + "): " );
            }
        }
    }

    // Generate Referenced Shortcuts
    auto &time = add_to_output( "TwoPMat", "Time", std::vector<Scalar>( maximum_time, 0 ), to_output, false /* Overwrite Existing */ );
    auto &twophotonmatrix = add_to_output( "TwoPMat", purpose, std::vector<Dense>( maximum_time, Dense::Zero( 8, 8 ) ), to_output_m );

    // Calculate two-photon densitymatrices
#pragma omp parallel for schedule( dynamic ) shared( timer_c ) num_threads( s.parameters.numerics_maximum_primary_threads )
    for ( size_t k = 0; k < maximum_time; k++ ) {
        // Neumann-Iterated TPDM
        const auto rho_3phot = _rho_to_tpdm( k, rho );

        // Cache to output arrays
        time.at( k ) = s.getTimeOf( 0, { k, 0, 0 } );
        twophotonmatrix.at( k ) = rho_3phot;
        timer_c.iterate();
    }
    // Final output and timer end
    timer_c.end();
    Timers::outputProgress( timer_c, progressbar, timer_c.getTotalIterationNumber(), pbsize, "Three Photon Matrix (" + purpose + ")", Timers::PROGRESS_FORCE_OUTPUT );

    return true;
}