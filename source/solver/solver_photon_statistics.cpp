#include "solver/solver_ode.h"

// Description: Calculates G2 photon statistics
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator(annihilator)_X: set of operator matrices and their conjugated forms
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false

bool ODESolver::calculate_advanced_photon_statistics( System &s, const Sparse &op_creator_1, const Sparse &op_annihilator_1, const Sparse &op_creator_2, const Sparse &op_annihilator_2, std::string fileOutputName ) {
    //bool calculate_full_g2_of_tau = !s.parameters.numerics_use_simplified_g2;
    bool output_full_g2 = false;
    bool output_full_ind = s.parameters.numerics_calculate_timeresolution_indistinguishability;
    bool calculate_concurrence_with_g2_of_zero = s.parameters.numerics_use_simplified_g2;
    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
    s.command( Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
    int pbsize = (int)savedStates.size() / s.parameters.iterations_t_skip;
    // Progress
    ProgressBar progressbar = ProgressBar( pbsize );
    // Cache matrix. Saves complex double entries of G1(t,tau). Will be created even if they are not needed
    Dense akf_mat_g2_11 = Dense::Zero( dim, dim );
    Dense akf_mat_g2_22 = Dense::Zero( dim, dim );
    Dense akf_mat_g2_12 = Dense::Zero( dim, dim );
    // Calculate G2(t,tau) with given operator matrices
    if ( !calculate_concurrence_with_g2_of_zero ) {
        Log::L2( "Calculating 3 seperate G2(t,tau) matrices...\n" );
        if ( s.parameters.numerics_calculate_g2_H || s.parameters.numerics_calculate_g2_C )
            calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_1, op_annihilator_1, akf_mat_g2_11, "Mode 11" );
        if ( s.parameters.numerics_calculate_g2_V || s.parameters.numerics_calculate_g2_C )
            calculate_g2( s, op_creator_2, op_annihilator_2, op_creator_2, op_annihilator_2, akf_mat_g2_22, "Mode 22" );
        if ( s.parameters.numerics_calculate_g2_C )
            calculate_g2( s, op_creator_1, op_annihilator_2, op_creator_1, op_annihilator_2, akf_mat_g2_12, "Mode 12" );
    }
    // Modified Dimension Step. If Upscaling of G1/G2 is active, step for matrices are i instead of t_skip
    int mat_step = ( s.parameters.numerics_stretch_correlation_grid ? 1 : s.parameters.iterations_t_skip );
    Log::L2( "Using mat_step = {}\n", mat_step );
    // Modified T-Step
    double deltaT = s.parameters.t_step * ( 1.0 * mat_step );
    // Calculating G2(tau)
    Timer &timer_g2 = Timers::create( "G2(tau integral" );
    timer_g2.start();
    Log::L2( "Calculating G2(tau)... " );
    // Prepare output vector
    std::vector<Scalar> g2_11, g2_22, g2_12;
    std::vector<Scalar> g2_11_zero, g2_22_zero, g2_12_zero;
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        g2_11.push_back( 0 );
        g2_22.push_back( 0 );
        g2_12.push_back( 0 );
        g2_11_zero.push_back( 0 );
        g2_22_zero.push_back( 0 );
        g2_12_zero.push_back( 0 );
    }
    Sparse M1, M2, M3, rho;
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        int k = i / mat_step;
        if ( !calculate_concurrence_with_g2_of_zero ) {
            for ( int t = 0; t < dim; t++ ) {
                g2_11.at( k ) += akf_mat_g2_11( k, t ) * s.parameters.t_step;
                g2_22.at( k ) += akf_mat_g2_22( k, t ) * s.parameters.t_step;
                g2_12.at( k ) += akf_mat_g2_12( k, t ) * s.parameters.t_step;
            }
        }
        double t_t = getTimeAt( i );
        rho = getRhoAt( i );
        if ( s.parameters.numerics_calculate_g2_H ) {
            M1 = op_creator_1 * op_creator_1 * op_annihilator_1 * op_annihilator_1;
            M2 = op_creator_1 * op_annihilator_1;
            g2_11_zero.at( k ) = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ); // / std::pow( s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t_t ), 2 );
        }
        if ( s.parameters.numerics_calculate_g2_V ) {
            M1 = op_creator_2 * op_creator_2 * op_annihilator_2 * op_annihilator_2;
            M2 = op_creator_2 * op_annihilator_2;
            g2_22_zero.at( k ) = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ); // / std::pow( s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t_t ), 2 );
        }
        if ( s.parameters.numerics_calculate_g2_C || s.parameters.numerics_calculate_g2_V || s.parameters.numerics_calculate_g2_H ) {
            M1 = op_creator_1 * op_creator_1 * op_annihilator_2 * op_annihilator_2;
            M2 = op_creator_1 * op_annihilator_1;
            M3 = op_creator_2 * op_creator_2;
            g2_12_zero.at( k ) = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ); // / ( s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t_t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, M3, t_t ) );
        }
        Timers::outputProgress( s.parameters.output_handlerstrings, timer_g2, progressbar, pbsize, "G2(tau): " );
        timer_g2.iterate();
    }
    timer_g2.end();
    // Calculating Concurrence
    // Prepare Concurrence output vector
    std::vector<Scalar> rho_11, rho_22, rho_12;
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        rho_11.push_back( 0 );
        rho_22.push_back( 0 );
        rho_12.push_back( 0 );
    }
    if ( s.parameters.numerics_calculate_g2_C ) {
        Timer &timer_c = Timers::create( "Concurrence" );
        timer_c.start();
        Log::L2( "Calculating Concurrence... " );
        // Main integral
        Log::L2( "Concurrence integral timestep: {}\n", deltaT );
        //#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
        for ( long unsigned int T = 0; T < savedStates.size(); T += mat_step ) {
            int t = T / mat_step;
            for ( int i = 0; i < T; i += mat_step ) { ////+= s.parameters.iterations_t_skip
                ////int k = i / s.parameters.iterations_t_skip;
                // Calculate for t
                if ( !calculate_concurrence_with_g2_of_zero ) {
                    for ( int tau = 0; tau < t - i; tau++ ) {
                        rho_11.at( t ) += akf_mat_g2_11( i, tau ) * s.parameters.t_step * deltaT;
                        rho_22.at( t ) += akf_mat_g2_22( i, tau ) * s.parameters.t_step * deltaT;
                        rho_12.at( t ) += akf_mat_g2_12( i, tau ) * s.parameters.t_step * deltaT;
                    }
                } else {
                    int k = i / mat_step;
                    rho_11.at( t ) += g2_11_zero.at( k ) * deltaT;
                    rho_22.at( t ) += g2_22_zero.at( k ) * deltaT;
                    rho_12.at( t ) += g2_12_zero.at( k ) * deltaT;
                }
                Timers::outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, pbsize, "Concurrence: " );
            }
            timer_c.iterate();
        }
        // Final output and timer end
        Timers::outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, pbsize, "Concurrence", PROGRESS_FORCE_OUTPUT );
        timer_c.end();
    }
    // Calculating indistinguishability
    Dense akf_mat_g1_1 = Dense::Zero( dim, dim );
    Dense akf_mat_g1_2 = Dense::Zero( dim, dim );
    std::vector<double> p1, p2;
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        p1.push_back( 0 );
        p2.push_back( 0 );
    }
    if ( !calculate_concurrence_with_g2_of_zero ) {
        M1 = op_creator_1 * op_annihilator_1;
        M2 = op_creator_2 * op_annihilator_2;
        if ( !cache1.isZero() )
            akf_mat_g1_1 = cache1;
        else if ( s.parameters.numerics_calculate_g2_H ) //TODO: hier vieleicht weg, nur g1 benutzen wenn auch das spektrum ausgerechnet wurde, sonst g1=0
            calculate_g1( s, op_creator_1, op_annihilator_1, akf_mat_g1_1, "g1_1" );
        if ( !cache2.isZero() )
            akf_mat_g1_2 = cache2;
        else if ( s.parameters.numerics_calculate_g2_V )
            calculate_g1( s, op_creator_2, op_annihilator_2, akf_mat_g1_2, "g1_2" );
        Timer &timer = Timers::create( "Indistinguishability" );
        timer.start();
        long unsigned int tstart = output_full_ind ? 0 : savedStates.size() - mat_step;
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
        for ( long unsigned int T = tstart; T < savedStates.size(); T += mat_step ) {
            Scalar top_1 = 0;
            Scalar top_2 = 0;
            Scalar bottom_1 = 0;
            Scalar bottom_2 = 0;
            for ( int i = 0; i < T; i += mat_step ) {
                int k = i / mat_step;
                auto rho = getRhoAt( i );
                double t = getTimeAt( i );
                for ( int j = 0; j < T - i; j += 1 ) { //&& j + i < (int)savedStates.size() ////mat_step
                    int l = j / mat_step;
                    auto rho_tau = getRhoAt( i + j );
                    double t_tau = getTimeAt( i + j );
                    Scalar gpop_1 = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M1, t_tau );
                    Scalar gpop_2 = s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M2, t_tau );
                    Scalar gbot_1 = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator_1, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator_1, t );
                    Scalar gbot_2 = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator_2, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator_2, t );
                    top_1 += ( gpop_1 + akf_mat_g2_11( i, j ) - akf_mat_g1_1( i, j ) * std::conj( akf_mat_g1_1( i, j ) ) ); //k,j nicht k, l
                    bottom_1 += 2.0 * gpop_1 - gbot_1 * std::conj( gbot_1 );
                    //fmt::print( "what contributes: gpop {}, g2_11 {}, g1_1 {}, gbot {}, TOP {}, BOT {}\n", gpop_1, akf_mat_g2_11( k, l ), std::pow( std::abs( akf_mat_g1_1( k, l ) ), 2.0 ), std::pow( std::abs( gbot_1 ), 2.0 ), top_1, bottom_1 );
                    top_2 += ( gpop_2 + akf_mat_g2_22( i, j ) - akf_mat_g1_2( i, j ) * std::conj( akf_mat_g1_2( i, j ) ) ); //T,j nicht k,l
                    bottom_2 += 2.0 * gpop_2 - gbot_2 * std::conj( gbot_2 );
                }
                // If we dont output the full time resoluted indistinguishability, we instead output the time resoluted integral for T = T_max
                if ( !output_full_ind ) {
                    p1.at( i / mat_step ) = ( 1.0 - std::abs( top_1 / bottom_1 ) );
                    p2.at( i / mat_step ) = ( 1.0 - std::abs( top_2 / bottom_2 ) );
                    Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability (Simplified): " );
                    timer.iterate();
                }
            }
            p1.at( T / mat_step ) = ( 1.0 - std::real( top_1 / bottom_1 ) );
            p2.at( T / mat_step ) = ( 1.0 - std::real( top_2 / bottom_2 ) );
            Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability (Full): " );
            timer.iterate();
        }
        // Final output and timer end
        Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability: ", PROGRESS_FORCE_OUTPUT );
        timer.end();
    }
    // Save output
    Log::L2( "Done, saving to {}... ", fileOutputName );
    std::string filepath = s.parameters.subfolder + fileOutputName;
    FILE *concurrencefile = std::fopen( filepath.c_str(), "w" );
    if ( !concurrencefile ) {
        Log::L2( "Failed to open outputfile for advanced photon statistics!\n" );
        return false;
    }
    fmt::print( concurrencefile, "Time\tMatrixelement_11\tMatrixelement_22\tConcurrence\tG2_11(tau)\tG2_22(tau)\tG2_12(tau)\tG2(tau=0)_11\tG2(tau=0)_22\tG2(tau=0)_12\tI_1\tI_2\n" );
    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
        int k = i / mat_step;
        fmt::print( concurrencefile, "{:.5e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", getTimeAt( i ), std::abs( rho_11.at( k ) / ( rho_11.at( k ) + rho_22.at( k ) ) ), std::abs( rho_22.at( k ) / ( rho_11.at( k ) + rho_22.at( k ) ) ), std::abs( 2.0 * std::abs( rho_12.at( k ) ) / ( rho_11.at( k ) + rho_22.at( k ) ) ), std::real( g2_11.at( k ) ), std::real( g2_22.at( k ) ), std::real( g2_12.at( k ) ), std::real( g2_11_zero.at( k ) ), std::real( g2_22_zero.at( k ) ), std::real( g2_12_zero.at( k ) ), p1.at( k ), p2.at( k ) );
    }
    std::fclose( concurrencefile );
    Log::L2( "Done!\n" );
    Log::L1( "Final Concurrence: {}\n", 2.0 * std::abs( rho_12.back() / ( rho_11.back() + rho_22.back() ) ) );
    if ( !calculate_concurrence_with_g2_of_zero )
        Log::L1( "Final Indistinguishability: {} H, {} V\n", p1.back(), p2.back() );
    //if ()
    //Log::L1( "Final G2(0): {}\n", std::real( 2.0 * std::abs( rho_12.back() ) / ( rho_11.back() + rho_22.back() ) ) );

    if ( output_full_g2 ) {
        Log::L2( "Saving full, non integrated G2(t,tau) matrices... " );
        FILE *ganz = std::fopen( ( s.parameters.subfolder + "g2(ttau)_matrix.txt" ).c_str(), "w" );
        fmt::print( ganz, "t\ttau\tG2_11(t,tau)\tG2_22(t,tau)\tG2_12(t,tau)\n" );
        for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
            int k = i / mat_step;
            for ( long unsigned int j = 0; j < savedStates.size(); j += mat_step ) {
                int l = j / mat_step;
                fmt::print( ganz, "{0:.5e}\t{1:.5e}\t{2:.8e}\t{3:.8e}\t{4:.8e}\n", getTimeAt( i ), getTimeAt( j ), std::abs( akf_mat_g2_11( k, l ) ), std::abs( akf_mat_g2_22( k, l ) ), std::abs( akf_mat_g2_12( k, l ) ) );
            }
            fmt::print( ganz, "\n" );
        }
        std::fclose( ganz );
        Log::L2( "Done!\n" );
    }
    // Sucessfull
    return true;
}