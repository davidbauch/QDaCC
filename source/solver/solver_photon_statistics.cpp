#include "solver/solver_ode.h"

// Description: Calculates G2 Indistinguishability
bool ODESolver::calculate_indistinguishability( System &s, const std::string &s_op_creator, const std::string &s_op_annihilator ) {
    Sparse op_creator = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_annihilator = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    for ( auto &split_s_op_creator : splitline( s_op_creator, '+' ) )
        op_creator += s.operatorMatrices.el_transitions.count( split_s_op_creator ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_creator].hilbert : s.operatorMatrices.ph_transitions[split_s_op_creator].hilbert;
    for ( auto &split_s_op_annihilator : splitline( s_op_annihilator, '+' ) )
        op_annihilator += s.operatorMatrices.el_transitions.count( split_s_op_annihilator ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_annihilator].hilbert : s.operatorMatrices.ph_transitions[split_s_op_annihilator].hilbert;

    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
    s.command( Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
    // Progress
    int pbsize = (int)savedStates.size() / s.parameters.iterations_t_skip;
    ProgressBar progressbar = ProgressBar( pbsize );

    // Calculate G2(t,tau) with given operator matrices
    std::string s_g1 = "G1-" + s_op_creator + "-" + s_op_annihilator;
    std::string s_g2 = "G2-" + s_op_creator + "-" + s_op_annihilator + "-" + s_op_creator + "-" + s_op_annihilator;
    if ( cache.count( s_g1 ) == 0 )
        cache[s_g1] = calculate_g1( s, op_creator, op_annihilator, s_g1 );
    if ( cache.count( s_g2 ) == 0 )
        cache[s_g2] = calculate_g2( s, op_creator, op_annihilator, op_creator, op_annihilator, s_g2 );

    Dense &akf_mat_g1 = cache[s_g1];
    Dense &akf_mat_g2 = cache[s_g2];

    std::vector<Scalar> outp; //( std::ceil( savedStates.size() / mat_step ), 0 );
    std::vector<Scalar> time; //( std::ceil( savedStates.size() / mat_step ), 0 );

    Sparse M1 = op_creator * op_annihilator;

    std::string fout = s_op_creator + "-" + s_op_annihilator;
    Timer &timer = Timers::create( "Indistinguishability " + fout );
    timer.start();

    int mat_step = ( s.parameters.numerics_stretch_correlation_grid ? 1 : s.parameters.iterations_t_skip );
    Log::L2( "Using mat_step = {}\n", mat_step );

    std::vector<Scalar> top;
    std::vector<Scalar> bottom;
    auto T = savedStates.size();
    for ( int i = 0; i < T; i += mat_step ) {
        outp.emplace_back( 0 );
        top.emplace_back( 0 );
        bottom.emplace_back( 0 );
        time.emplace_back( getTimeAt( i ) );
    }
#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < T; i += mat_step ) {
        int k = i / mat_step;
        auto rho = getRhoAt( i );
        double t = getTimeAt( i );
        for ( int j = 0; j < T - i; j += 1 ) { //&& j + i < (int)savedStates.size() ////mat_step
            int l = j / mat_step;
            auto rho_tau = getRhoAt( i + j );
            double t_tau = getTimeAt( i + j );
            Scalar gpop = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M1, t_tau );
            Scalar gbot = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator, t );
            top[k] += ( gpop + akf_mat_g2( i, j ) - akf_mat_g1( i, j ) * std::conj( akf_mat_g1( i, j ) ) ); //k,j nicht k, l
            bottom[k] += 2.0 * gpop - gbot * std::conj( gbot );
        }
        Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability (Simplified) (" + fout + "): " );
        timer.iterate();
    }
    Scalar topsum = 0;
    Scalar bottomsum = 0;
    for ( int k = 0; k < top.size(); k++ ) {
        topsum += top[k];
        bottomsum += bottom[k];
        outp[k] = ( 1.0 - std::abs( topsum / bottomsum ) ); // .at( i / mat_step ) =
    }
    // Final output and timer end
    Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability (" + fout + "): ", PROGRESS_FORCE_OUTPUT );
    timer.end();

    // Add to Fileoutput:
    if ( to_output["Indist"].size() == 0 )
        to_output["Indist"]["Time"] = time;
    to_output["Indist"][fout] = outp;
    Log::L1( "Final Indistinguishability: {} {}\n", std::real( outp.back() ), fout );
    return true;
}

// Description: Calculates Concurrence
bool ODESolver::calculate_concurrence( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2 ) {
    // Find Matrices
    Sparse op_creator_1 = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_creator_2 = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_annihilator_1 = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    Sparse op_annihilator_2 = Sparse( s.parameters.maxStates, s.parameters.maxStates );
    for ( auto &split_s_op_creator_1 : splitline( s_op_creator_1, '+' ) )
        op_creator_1 += s.operatorMatrices.el_transitions.count( split_s_op_creator_1 ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_creator_1].hilbert : s.operatorMatrices.ph_transitions[split_s_op_creator_1].hilbert;
    for ( auto &split_s_op_creator_2 : splitline( s_op_creator_2, '+' ) )
        op_creator_2 += s.operatorMatrices.el_transitions.count( split_s_op_creator_2 ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_creator_2].hilbert : s.operatorMatrices.ph_transitions[split_s_op_creator_2].hilbert;
    for ( auto &split_s_op_annihilator_1 : splitline( s_op_annihilator_1, '+' ) )
        op_annihilator_1 += s.operatorMatrices.el_transitions.count( split_s_op_annihilator_1 ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_annihilator_1].hilbert : s.operatorMatrices.ph_transitions[split_s_op_annihilator_1].hilbert;
    for ( auto &split_s_op_annihilator_2 : splitline( s_op_annihilator_2, '+' ) )
        op_annihilator_2 += s.operatorMatrices.el_transitions.count( split_s_op_annihilator_2 ) != 0 ? s.operatorMatrices.el_transitions[split_s_op_annihilator_2].hilbert : s.operatorMatrices.ph_transitions[split_s_op_annihilator_2].hilbert;

    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
    s.command( Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
    // Progress
    int pbsize = (int)savedStates.size() / s.parameters.iterations_t_skip;
    ProgressBar progressbar = ProgressBar( pbsize );

    std::string fout = s_op_creator_1 + "-" + s_op_annihilator_1 + "-" + s_op_creator_2 + "-" + s_op_annihilator_2;
    // Calculate G2(t,tau) with given operator matrices
    std::string s_g1_1 = "G1-" + s_op_creator_1 + "-" + s_op_annihilator_1;
    std::string s_g1_2 = "G1-" + s_op_creator_2 + "-" + s_op_annihilator_2;
    std::string s_g2_1111 = "G2-" + s_op_creator_1 + "-" + s_op_annihilator_1 + "-" + s_op_creator_1 + "-" + s_op_annihilator_1;
    std::string s_g2_1122 = "G2-" + s_op_creator_1 + "-" + s_op_annihilator_2 + "-" + s_op_creator_1 + "-" + s_op_annihilator_2;
    std::string s_g2_1212 = "G2-" + s_op_creator_1 + "-" + s_op_annihilator_1 + "-" + s_op_creator_2 + "-" + s_op_annihilator_2;
    std::string s_g2_1221 = "G2-" + s_op_creator_1 + "-" + s_op_annihilator_2 + "-" + s_op_creator_2 + "-" + s_op_annihilator_1;
    std::string s_g2_2121 = "G2-" + s_op_creator_2 + "-" + s_op_annihilator_2 + "-" + s_op_creator_1 + "-" + s_op_annihilator_1;
    std::string s_g2_2112 = "G2-" + s_op_creator_2 + "-" + s_op_annihilator_1 + "-" + s_op_creator_1 + "-" + s_op_annihilator_2;
    std::string s_g2_2211 = "G2-" + s_op_creator_2 + "-" + s_op_annihilator_1 + "-" + s_op_creator_2 + "-" + s_op_annihilator_1;
    std::string s_g2_2222 = "G2-" + s_op_creator_2 + "-" + s_op_annihilator_2 + "-" + s_op_creator_2 + "-" + s_op_annihilator_2;
    if ( cache.count( s_g1_1 ) == 0 )
        cache[s_g1_1] = calculate_g1( s, op_creator_1, op_annihilator_1, s_g1_1 );
    if ( cache.count( s_g1_2 ) == 0 )
        cache[s_g1_2] = calculate_g1( s, op_creator_2, op_annihilator_2, s_g1_2 );
    if ( cache.count( s_g2_1111 ) == 0 )
        cache[s_g2_1111] = calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_1, op_annihilator_1, s_g2_1111 );
    if ( cache.count( s_g2_1122 ) == 0 )
        cache[s_g2_1122] = calculate_g2( s, op_creator_1, op_annihilator_2, op_creator_1, op_annihilator_2, s_g2_1122 );
    if ( cache.count( s_g2_2121 ) == 0 )
        cache[s_g2_2121] = calculate_g2( s, op_creator_2, op_annihilator_2, op_creator_1, op_annihilator_1, s_g2_2121 );
    if ( cache.count( s_g2_2112 ) == 0 )
        cache[s_g2_2112] = calculate_g2( s, op_creator_2, op_annihilator_1, op_creator_1, op_annihilator_2, s_g2_2112 );
    if ( cache.count( s_g2_1221 ) == 0 )
        cache[s_g2_1221] = calculate_g2( s, op_creator_1, op_annihilator_2, op_creator_2, op_annihilator_1, s_g2_1221 );
    if ( cache.count( s_g2_2222 ) == 0 )
        cache[s_g2_2222] = calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_1, op_annihilator_1, s_g2_2222 );

    Dense &akf_mat_g1_1 = cache[s_g1_1];
    Dense &akf_mat_g1_2 = cache[s_g1_1];
    Dense &akf_mat_g2_1111 = cache[s_g2_1111];
    Dense &akf_mat_g2_1122 = cache[s_g2_1122];
    Dense &akf_mat_g2_1221 = cache[s_g2_1221];
    Dense &akf_mat_g2_2121 = cache[s_g2_2121];
    Dense &akf_mat_g2_2112 = cache[s_g2_2112];
    Dense &akf_mat_g2_2222 = cache[s_g2_2222];

    int mat_step = ( s.parameters.numerics_stretch_correlation_grid ? 1 : s.parameters.iterations_t_skip );
    Log::L2( "Using mat_step = {}\n", mat_step );

    std::map<std::string, std::vector<Scalar>> rho;
    //for ( auto &mode : { s_g2_1111, s_g2_1122, s_g2_1221, s_g2_2121, s_g2_2112, s_g2_2222, s_g2_2211, s_g2_1212 } ) {
    //    rho[mode] = std::vector<Scalar>( std::ceil( savedStates.size() / mat_step ), 0 );
    //}

    Timer &timer_c = Timers::create( "Concurrence (" + fout + ")" );
    timer_c.start();

    double deltaT = s.parameters.t_step * ( 1.0 * mat_step );
    Log::L2( "Calculating Concurrence... Integral Timestep: {}\n", deltaT );
    for ( long unsigned int T = 0; T < savedStates.size(); T += mat_step ) {
        int t = T / mat_step;
        rho[s_g2_1111].emplace_back( 0 );
        rho[s_g2_1122].emplace_back( 0 );
        rho[s_g2_1221].emplace_back( 0 );
        rho[s_g2_2121].emplace_back( 0 );
        rho[s_g2_2112].emplace_back( 0 );
        rho[s_g2_2222].emplace_back( 0 );
        rho[s_g2_2211].emplace_back( 0 );
        rho[s_g2_1212].emplace_back( 0 );
        for ( int i = 0; i < T; i += mat_step ) {
            for ( int tau = 0; tau < T - i; tau++ ) { // FIX: tau < t-i --> tau < T-i, da G2 mit t_step gespeichert wird! (i,T same scaling)
                rho[s_g2_1111].back() += akf_mat_g2_1111( i, tau ) * s.parameters.t_step * deltaT;
                rho[s_g2_1122].back() += akf_mat_g2_1122( i, tau ) * s.parameters.t_step * deltaT;
                rho[s_g2_1221].back() += akf_mat_g2_1221( i, tau ) * s.parameters.t_step * deltaT;
                rho[s_g2_2121].back() += akf_mat_g2_2121( i, tau ) * s.parameters.t_step * deltaT;
                rho[s_g2_2112].back() += akf_mat_g2_2112( i, tau ) * s.parameters.t_step * deltaT;
                rho[s_g2_2222].back() += akf_mat_g2_2222( i, tau ) * s.parameters.t_step * deltaT;
            }
            rho[s_g2_2211].back() = std::conj( rho[s_g2_1122].at( t ) );
            rho[s_g2_1212].back() = std::conj( rho[s_g2_2121].at( t ) );
        }
        timer_c.iterate();
        Timers::outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, pbsize, "Concurrence (" + fout + "): " );
    }
    // Calculate EigenValues
    std::vector<Scalar> output;
    std::vector<Scalar> time;
    std::vector<Dense> twophotonmatrix;

    Dense spinflip = Dense::Zero( 4, 4 );
    spinflip( 0, 3 ) = -1;
    spinflip( 1, 2 ) = 1;
    spinflip( 2, 1 ) = 1;
    spinflip( 3, 0 ) = -1;
    Log::L3( "Spinflip Matrix: {}\n", spinflip );
    Scalar conc = 0;
    for ( long unsigned int i = 0; i < rho[s_g2_1111].size(); i += mat_step ) {
        int k = i / mat_step;
        //Log::L3( "Creating 2 photon matrix\n" );
        Dense rho_2phot = Dense::Zero( 4, 4 );
        rho_2phot( 0, 0 ) = rho[s_g2_1111][k];
        rho_2phot( 3, 3 ) = rho[s_g2_2222][k];
        rho_2phot( 0, 3 ) = rho[s_g2_2211][k];
        rho_2phot( 3, 0 ) = rho[s_g2_1122][k];
        rho_2phot( 1, 2 ) = rho[s_g2_2121][k];
        rho_2phot( 2, 1 ) = rho[s_g2_1212][k];
        rho_2phot = rho_2phot / rho_2phot.trace();
        rho_2phot( 1, 1 ) = std::abs( rho[s_g2_2112][k] ) != 0 ? rho[s_g2_2112][k] : 1E-100;
        rho_2phot( 2, 2 ) = std::abs( rho[s_g2_1221][k] ) != 0 ? rho[s_g2_1221][k] : 1E-100;
        //Log::L3( "Normalizing 2 photon matrix\n" );
        //if ( std::abs( rho_2phot.trace() ) != 0 ) {
        //Log::L3( "Rho_2phot = {}\n", rho_2phot );
        //Log::L3( "Calculating sqrt(rho)\n" );
        //Eigen::MatrixPower<Dense> Mpow( rho_2phot );
        Dense sqrtrho2phot = rho_2phot.sqrt(); //Mpow( 0.5 );
        //Log::L3( "Calculating R\n" );
        Dense R = sqrtrho2phot * spinflip * rho_2phot * spinflip * sqrtrho2phot;
        //Log::L3( "R = {}\n", R );
        //Log::L3( "Calculating sqrt(R)\n" );
        //Eigen::MatrixPower<Dense> SMPow( R );
        auto R5 = R.sqrt(); //SMPow( 0.5 );
        //Log::L3( "Calculating Eigenvalues\n" );
        auto eigenvalues = R5.eigenvalues();
        //Log::L1( "rho2phot = {}\n\nsqrtrho2phot = {}\n\nR = {}\n\nRS = {}\nEigenvalues at t = {} are {}\n", rho_2phot, sqrtrho2phot, R, R5, getTimeAt( i ), eigenvalues );
        conc = eigenvalues( 3 ) - eigenvalues( 2 ) - eigenvalues( 1 ) - eigenvalues( 0 );
        output.emplace_back( conc );
        time.emplace_back( getTimeAt( i ) );
        twophotonmatrix.emplace_back( rho_2phot );
        //}
    }
    // Final output and timer end
    Timers::outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, pbsize, "Concurrence (" + fout + ")", PROGRESS_FORCE_OUTPUT );
    timer_c.end();
    // Add to Fileoutput:
    if ( to_output["Conc"].size() == 0 )
        to_output["Conc"]["Time"] = time;
    to_output["Conc"][fout] = output;
    to_output_m["TwoPMat"][fout] = twophotonmatrix;
    Log::L1( "Final Concurrence: {} {}\n", std::abs( output.back() ), fout );
    return true;
}

// Description: Calculates and outputs G1 and G2 photon statistics
// Type: ODESolver public function
// @param s: [&System] Class providing set of system functions
// @param op_creator(annihilator)_X: set of operator matrices and their conjugated forms
// @param fileOutputName: [std::string] Name of output file
// @return: [bool] True if calculations were sucessfull, else false

bool ODESolver::calculate_advanced_photon_statistics( System &s ) {
    // Calculate Spectra
    auto &spectrum_s = s.parameters.input_correlation["Spectrum"];
    for ( int i = 0; i < spectrum_s.string_v["Modes"].size(); i++ ) {
        std::string s_creator, s_annihilator;
        for ( auto split_s_op : splitline( spectrum_s.string_v["Modes"][i], '+' ) ) {
            if ( s_creator.size() > 0 ) {
                s_creator += "+";
                s_annihilator += "+";
            }
            if ( std::isupper( split_s_op.front() ) ) {
                s_annihilator += split_s_op;
                std::reverse( split_s_op.begin(), split_s_op.end() );
                s_creator += split_s_op;
            } else {
                s_creator += split_s_op + "bd";
                s_annihilator += split_s_op + "b";
            }
        }
        calculate_spectrum( s, s_creator, s_annihilator, spectrum_s.numerical_v["Center"][i], spectrum_s.numerical_v["Range"][i], spectrum_s.numerical_v["resW"][i] );
    }
    // Calculate Indist
    auto &indist_s = s.parameters.input_correlation["Indist"];
    for ( int i = 0; i < indist_s.string_v["Modes"].size(); i++ ) {
        std::string s_creator, s_annihilator;
        for ( auto split_s_op : splitline( indist_s.string_v["Modes"][i], '+' ) ) {
            if ( s_creator.size() > 0 ) {
                s_creator += "+";
                s_annihilator += "+";
            }
            if ( std::isupper( split_s_op.front() ) ) {
                s_annihilator += split_s_op;
                std::reverse( split_s_op.begin(), split_s_op.end() );
                s_creator += split_s_op;
            } else {
                s_creator += split_s_op + "bd";
                s_annihilator += split_s_op + "b";
            }
        }
        calculate_indistinguishability( s, s_creator, s_annihilator );
    }
    // Calculate Conc
    auto &conc_s = s.parameters.input_correlation["Conc"];
    for ( auto &modes : conc_s.string_v["Modes"] ) {
        std::vector<std::string> s_creator, s_annihilator;
        for ( auto &mode : splitline( modes, '-' ) ) {
            if ( std::isupper( mode.front() ) ) {
                s_annihilator.emplace_back( mode );
                auto ss_creator = mode;
                std::reverse( ss_creator.begin(), ss_creator.end() );
                s_creator.emplace_back( ss_creator );
            } else {
                s_creator.emplace_back( mode + "bd" );
                s_annihilator.emplace_back( mode + "b" );
            }
        }
        calculate_concurrence( s, s_creator[0], s_annihilator[0], s_creator[1], s_annihilator[1] );
    }
    // Calculate G1/G2 functions
    auto &gs_s = s.parameters.input_correlation["GFunc"];
    for ( int i = 0; i < gs_s.string_v["Modes"].size(); i++ ) {
        auto modes = gs_s.string_v["Modes"][i];
        int order = std::abs( gs_s.numerical_v["Order"][i] );
        std::string s_creator, s_annihilator;
        Sparse creator, annihilator;
        for ( auto &mode : splitline( modes, '-' ) ) {
            if ( std::isupper( mode.front() ) ) {
                s_annihilator = mode;
                auto s_creator = mode;
                std::reverse( s_creator.begin(), s_creator.end() );
                creator = s.operatorMatrices.el_transitions[s_creator].hilbert;
                annihilator = s.operatorMatrices.el_transitions[s_annihilator].hilbert;
            } else {
                s_creator = mode + "bd";
                s_annihilator = mode + "b";
                creator = s.operatorMatrices.el_transitions[s_creator].hilbert;
                annihilator = s.operatorMatrices.el_transitions[s_annihilator].hilbert;
            }
        }
        std::string purpose = order == 1 ? "G1-" + s_creator + "-" + s_annihilator : "G2-" + s_creator + "-" + s_annihilator + "-" + s_creator + "-" + s_annihilator;
        if ( cache.count( purpose ) == 0 )
            if ( order == 1 )
                calculate_g1( s, creator, annihilator, purpose );
            else
                calculate_g2( s, creator, annihilator, creator, annihilator, purpose );
        // Directly output corresponding matrix here so G1/2 functions calculated by other function calls are not output if they are not demanded.
        auto &gmat = cache[purpose];
        if ( gs_s.numerical_v["Integrated"][i] == 0 || gs_s.numerical_v["Integrated"][i] == 2 ) {
            Log::L2( "Saving G{} function matrix to {}_m.txt...\n", order, purpose );
            FILE *f_gfunc = std::fopen( ( s.parameters.subfolder + purpose + "_m.txt" ).c_str(), "w" );
            for ( int k = 0; k < gmat.rows(); k += s.parameters.iterations_t_skip ) {
                for ( int l = 0; l < gmat.cols(); l++ ) {
                    fmt::print( f_gfunc, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", 1.0 * k * s.parameters.t_step, 1.0 * l * s.parameters.t_step, std::abs( gmat( k, l ) ), std::real( gmat( k, l ) ), std::imag( gmat( k, l ) ) );
                }
                fmt::print( f_gfunc, "\n" );
            }
            std::fclose( f_gfunc );
        }
        if ( gs_s.numerical_v["Integrated"][i] == 1 || gs_s.numerical_v["Integrated"][i] == 2 ) {
            Log::L2( "Saving G{} integrated function to {}.txt...\n", order, purpose );
            FILE *f_gfunc = std::fopen( ( s.parameters.subfolder + purpose + ".txt" ).c_str(), "w" );
            fmt::print( f_gfunc, "Time\tAbs\tReal\tImag\n" );
            for ( int l = 0; l < gmat.cols(); l++ ) {
                Scalar integral = 0;
                for ( int k = 0; k < gmat.rows(); k += s.parameters.iterations_t_skip ) {
                    integral += gmat( k, l );
                }
                integral *= s.parameters.t_step * s.parameters.iterations_t_skip;
                fmt::print( f_gfunc, "{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", 1.0 * l * s.parameters.t_step, std::abs( integral ), std::real( integral ), std::imag( integral ) );
            }
            std::fclose( f_gfunc );
        }
        Log::L2( "Done!\n" );
    }
    // Output Spectra and Rest in seperate Files
    for ( auto &[mode, data] : to_output["Spectrum"] ) {
        Log::L2( "Saving Emission Spectrum to spectrum_" + mode + ".txt...\n" );
        FILE *f_spectrum = std::fopen( ( s.parameters.subfolder + "spectrum_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_spectrum, "Omega\t{}\n", mode );
        for ( int i = 0; i < to_output["Spectrum"][mode].size(); i++ ) {
            fmt::print( f_spectrum, "{:.8e}\t{:.8e}\n", std::real( to_output["Spectrum_frequency"][mode][i] ), std::real( to_output["Spectrum"][mode][i] ) );
        }
        std::fclose( f_spectrum );
        Log::L2( "Done!\n" );
    }
    for ( auto &[mode, data] : to_output["Indist"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "Saving Indistinguishability to indist_" + mode + ".txt...\n" );
        FILE *f_indist = std::fopen( ( s.parameters.subfolder + "indist_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_indist, "Time\t{}\n", mode );
        for ( int i = 0; i < to_output["Indist"][mode].size(); i++ ) {
            fmt::print( f_indist, "{:.8e}\t{:.8e}\n", std::real( to_output["Indist"]["Time"][i] ), std::real( to_output["Indist"][mode][i] ) );
        }
        std::fclose( f_indist );
        Log::L2( "Done!\n" );
    }
    for ( auto &[mode, data] : to_output["Conc"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "Saving Concurrence to conc_" + mode + ".txt...\n" );
        FILE *f_indist = std::fopen( ( s.parameters.subfolder + "conc_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_indist, "Time\t{}\n", mode );
        for ( int i = 0; i < to_output["Conc"][mode].size(); i++ ) {
            fmt::print( f_indist, "{:.8e}\t{:.8e}\n", std::real( to_output["Conc"]["Time"][i] ), std::real( to_output["Conc"][mode][i] ) );
        }
        std::fclose( f_indist );
        Log::L2( "Done!\n" );
    }
    for ( auto &[mode, data] : to_output_m["TwoPMat"] ) {
        if ( mode.compare( "Time" ) == 0 )
            continue;
        Log::L2( "Saving Two-Photon Matrix to twopmat_" + mode + ".txt...\n" );
        FILE *f_indist = std::fopen( ( s.parameters.subfolder + "twopmat_" + mode + ".txt" ).c_str(), "w" );
        fmt::print( f_indist, "Time\t{}\n", mode );
        for ( int i = 0; i < to_output_m["TwoPMat"][mode].size(); i++ ) {
            fmt::print( f_indist, "{:.8e}\t", std::real( to_output["Conc"]["Time"][i] ) );
            for ( int k = 0; k < 4; k++ ) {
                for ( int l = 0; l < 4; l++ ) {
                    fmt::print( f_indist, "{:.8e}\t", std::abs( to_output_m["TwoPMat"][mode][i]( k, l ) ) );
                }
            }
            fmt::print( f_indist, "\n" );
        }
        std::fclose( f_indist );
        Log::L2( "Done!\n" );
    }
    return true;
}
//bool ODESolver::calculate_advanced_photon_statistics( System &s, const std::string &s_op_creator_1, const std::string &s_op_annihilator_1, const std::string &s_op_creator_2, const std::string &s_op_annihilator_2, std::string fileOutputName ) {
//    // Find Matrices
//    Sparse op_creator_1 = s.operatorMatrices.el_transitions.find( s_op_creator_1 ) != 0 ? s.operatorMatrices.el_transitions[s_op_creator_1].hilbert : s.operatorMatrices.ph_transitions[s_op_creator_1].hilbert;
//    Sparse op_creator_2 = s.operatorMatrices.el_transitions.find( s_op_creator_2 ) != 0 ? s.operatorMatrices.el_transitions[s_op_creator_2].hilbert : s.operatorMatrices.ph_transitions[s_op_creator_2].hilbert;
//    Sparse op_annihilator_1 = s.operatorMatrices.el_transitions.find( s_op_annihilator_1 ) != 0 ? s.operatorMatrices.el_transitions[s_op_annihilator_1].hilbert : s.operatorMatrices.ph_transitions[s_op_annihilator_1].hilbert;
//    Sparse op_annihilator_2 = s.operatorMatrices.el_transitions.find( s_op_annihilator_2 ) != 0 ? s.operatorMatrices.el_transitions[s_op_annihilator_2].hilbert : s.operatorMatrices.ph_transitions[s_op_annihilator_2].hilbert;
//
//    //bool calculate_full_g2_of_tau = !s.parameters.numerics_use_simplified_g2;
//    bool output_full_g2 = false;
//    bool output_full_ind = s.parameters.numerics_calculate_timeresolution_indistinguishability;
//    bool calculate_concurrence_with_g2_of_zero = s.parameters.numerics_use_simplified_g2;
//    // Send system command to change to single core subprogram, because this memberfunction is already using multithreading
//    s.command( Solver::CHANGE_TO_SINGLETHREADED_MAINPROGRAM );
//    int pbsize = (int)savedStates.size() / s.parameters.iterations_t_skip;
//    // Progress
//    ProgressBar progressbar = ProgressBar( pbsize );
//    // Cache matrix. Saves complex double entries of G1(t,tau). Will be created even if they are not needed
//    Dense akf_mat_g2_HHHH = Dense::Zero( dim, dim ); // a
//    Dense akf_mat_g2_VVVV = Dense::Zero( dim, dim ); // b
//    Dense akf_mat_g2_HHVV = Dense::Zero( dim, dim ); // c-id
//    Dense akf_mat_g2_VHVH = Dense::Zero( dim, dim ); // e+if
//    Dense akf_mat_g2_VHHV = Dense::Zero( dim, dim ); // g
//    Dense akf_mat_g2_HVVH = Dense::Zero( dim, dim ); // h
//    //Dense akf_mat_g2_VVHH = Dense::Zero( dim, dim ); // c+id
//    //Dense akf_mat_g2_VVVH = Dense::Zero( dim, dim ); // 0
//    //Dense akf_mat_g2_VVHV = Dense::Zero( dim, dim ); // 0
//    //Dense akf_mat_g2_VHVV = Dense::Zero( dim, dim ); // 0
//    //Dense akf_mat_g2_VHHH = Dense::Zero( dim, dim ); // 0
//    //Dense akf_mat_g2_HVVV = Dense::Zero( dim, dim ); // 0
//    //Dense akf_mat_g2_HVHV = Dense::Zero( dim, dim ); // e-if
//    //Dense akf_mat_g2_HVHH = Dense::Zero( dim, dim ); // 0
//    //Dense akf_mat_g2_HHVH = Dense::Zero( dim, dim ); // 0
//    //Dense akf_mat_g2_HHHV = Dense::Zero( dim, dim ); // 0
//    //calculate_g2( s, op_creator_2, op_annihilator_1, op_creator_2, op_annihilator_1, akf_mat_g2_VVHH, "Mode VVHH" );
//    //calculate_g2( s, op_creator_2, op_annihilator_2, op_creator_2, op_annihilator_1, akf_mat_g2_VVVH, "Mode VVVH" );
//    //calculate_g2( s, op_creator_2, op_annihilator_1, op_creator_2, op_annihilator_2, akf_mat_g2_VVHV, "Mode VVHV" );
//    //calculate_g2( s, op_creator_2, op_annihilator_2, op_creator_1, op_annihilator_2, akf_mat_g2_VHVV, "Mode VHVV" );
//    //calculate_g2( s, op_creator_2, op_annihilator_1, op_creator_1, op_annihilator_1, akf_mat_g2_VHHH, "Mode VHHH" );
//    //calculate_g2( s, op_creator_1, op_annihilator_2, op_creator_2, op_annihilator_2, akf_mat_g2_HVVV, "Mode HVVV" );
//    //calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_2, op_annihilator_2, akf_mat_g2_HVHV, "Mode HVHV" );
//    //calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_2, op_annihilator_1, akf_mat_g2_HVHH, "Mode HVHH" );
//    //calculate_g2( s, op_creator_1, op_annihilator_2, op_creator_1, op_annihilator_1, akf_mat_g2_HHVH, "Mode HHVH" );
//    //calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_1, op_annihilator_2, akf_mat_g2_HHHV, "Mode HHHV" );
//    std::vector<Scalar> rho_VVHH, rho_VVVH, rho_VVHV, rho_VHVV, rho_VHVH, rho_VHHV, rho_VHHH, rho_HVVV, rho_HVVH, rho_HVHV, rho_HVHH, rho_HHVH, rho_HHHV;
//
//    // Calculate G2(t,tau) with given operator matrices
//    if ( !calculate_concurrence_with_g2_of_zero ) {
//        Log::L2( "Calculating up to 6 seperate G2(t,tau) matrices...\n" );
//        if ( s.parameters.numerics_calculate_g2_H || s.parameters.numerics_calculate_g2_C )
//            calculate_g2( s, op_creator_1, op_annihilator_1, op_creator_1, op_annihilator_1, akf_mat_g2_HHHH, "Mode HHHH" );
//        if ( s.parameters.numerics_calculate_g2_V || s.parameters.numerics_calculate_g2_C )
//            calculate_g2( s, op_creator_2, op_annihilator_2, op_creator_2, op_annihilator_2, akf_mat_g2_VVVV, "Mode VVVV" );
//        if ( s.parameters.numerics_calculate_g2_C ) {
//            calculate_g2( s, op_creator_1, op_annihilator_2, op_creator_1, op_annihilator_2, akf_mat_g2_HHVV, "Mode HHVV" );
//            calculate_g2( s, op_creator_2, op_annihilator_2, op_creator_1, op_annihilator_1, akf_mat_g2_VHVH, "Mode VHVH" );
//            calculate_g2( s, op_creator_2, op_annihilator_1, op_creator_1, op_annihilator_2, akf_mat_g2_VHHV, "Mode VHHV" );
//            calculate_g2( s, op_creator_1, op_annihilator_2, op_creator_2, op_annihilator_1, akf_mat_g2_HVVH, "Mode HVVH" );
//        }
//    }
//    // Modified Dimension Step. If Upscaling of G1/G2 is active, step for matrices are i instead of t_skip
//    int mat_step = ( s.parameters.numerics_stretch_correlation_grid ? 1 : s.parameters.iterations_t_skip );
//    Log::L2( "Using mat_step = {}\n", mat_step );
//    // Modified T-Step
//    double deltaT = s.parameters.t_step * ( 1.0 * mat_step );
//    // Calculating G2(tau)
//    Timer &timer_g2 = Timers::create( "G2(tau integral" );
//    timer_g2.start();
//    Log::L2( "Calculating G2(tau)... " );
//    // Prepare output vector
//    std::vector<Scalar> g2_11, g2_22, g2_12;
//    std::vector<Scalar> g2_11_zero, g2_22_zero, g2_12_zero;
//    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
//        g2_11.push_back( 0 );
//        g2_22.push_back( 0 );
//        g2_12.push_back( 0 );
//        g2_11_zero.push_back( 0 );
//        g2_22_zero.push_back( 0 );
//        g2_12_zero.push_back( 0 );
//    }
//    Sparse M1, M2, M3, rho;
//    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
//        int k = i / mat_step;
//        if ( !calculate_concurrence_with_g2_of_zero ) {
//            for ( int t = 0; t < dim; t++ ) {
//                g2_11.at( k ) += akf_mat_g2_HHHH( i, t ) * s.parameters.t_step; // FIX: akf(k,t) -> akf(i,t), da G2 mit t_step gespeichert wird!
//                g2_22.at( k ) += akf_mat_g2_VVVV( i, t ) * s.parameters.t_step; // FIX: akf(k,t) -> akf(i,t), da G2 mit t_step gespeichert wird!
//                g2_12.at( k ) += akf_mat_g2_HHVV( i, t ) * s.parameters.t_step; // FIX: akf(k,t) -> akf(i,t), da G2 mit t_step gespeichert wird!
//            }
//        }
//        double t_t = getTimeAt( i );
//        rho = getRhoAt( i );
//        if ( s.parameters.numerics_calculate_g2_H ) {
//            M1 = op_creator_1 * op_creator_1 * op_annihilator_1 * op_annihilator_1;
//            M2 = op_creator_1 * op_annihilator_1;
//            g2_11_zero.at( k ) = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ); // / std::pow( s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t_t ), 2 );
//        }
//        if ( s.parameters.numerics_calculate_g2_V ) {
//            M1 = op_creator_2 * op_creator_2 * op_annihilator_2 * op_annihilator_2;
//            M2 = op_creator_2 * op_annihilator_2;
//            g2_22_zero.at( k ) = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ); // / std::pow( s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t_t ), 2 );
//        }
//        if ( s.parameters.numerics_calculate_g2_C || s.parameters.numerics_calculate_g2_V || s.parameters.numerics_calculate_g2_H ) {
//            M1 = op_creator_1 * op_creator_1 * op_annihilator_2 * op_annihilator_2;
//            M2 = op_creator_1 * op_annihilator_1;
//            M3 = op_creator_2 * op_creator_2;
//            g2_12_zero.at( k ) = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t_t ); // / ( s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t_t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, M3, t_t ) );
//        }
//        Timers::outputProgress( s.parameters.output_handlerstrings, timer_g2, progressbar, pbsize, "G2(tau): " );
//        timer_g2.iterate();
//    }
//    timer_g2.end();
//
//    // Calculating Concurrence
//    // Prepare Concurrence output vector
//    std::vector<Scalar> rho_HHHH, rho_VVVV, rho_HHVV;
//    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
//        rho_HHHH.push_back( 0 );
//        rho_VVVV.push_back( 0 );
//        rho_HHVV.push_back( 0 );
//
//        rho_VVHH.push_back( 0 );
//        rho_VVVH.push_back( 0 );
//        rho_VVHV.push_back( 0 );
//        rho_VHVV.push_back( 0 );
//        rho_VHVH.push_back( 0 );
//        rho_VHHV.push_back( 0 );
//        rho_VHHH.push_back( 0 );
//        rho_HVVV.push_back( 0 );
//        rho_HVVH.push_back( 0 );
//        rho_HVHV.push_back( 0 );
//        rho_HVHH.push_back( 0 );
//        rho_HHVH.push_back( 0 );
//        rho_HHHV.push_back( 0 );
//    }
//    if ( s.parameters.numerics_calculate_g2_C ) {
//        Timer &timer_c = Timers::create( "Concurrence" );
//        timer_c.start();
//        Log::L2( "Calculating Concurrence... " );
//        // Main integral
//        Log::L2( "Concurrence integral timestep: {}\n", deltaT );
//        //#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
//        for ( long unsigned int T = 0; T < savedStates.size(); T += mat_step ) {
//            int t = T / mat_step;
//            for ( int i = 0; i < T; i += mat_step ) { ////+= s.parameters.iterations_t_skip
//                ////int k = i / s.parameters.iterations_t_skip;
//                // Calculate for t
//                if ( !calculate_concurrence_with_g2_of_zero ) {
//                    for ( int tau = 0; tau < T - i; tau++ ) { // FIX: tau < t-i --> tau < T-i, da G2 mit t_step gespeichert wird! (i,T same scaling)
//                        //rho_VVHH.at( t ) += akf_mat_g2_VVHH( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_VVVH.at( t ) += akf_mat_g2_VVVH( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_VVHV.at( t ) += akf_mat_g2_VVHV( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_VHVV.at( t ) += akf_mat_g2_VHVV( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_VHHH.at( t ) += akf_mat_g2_VHHH( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_HVVV.at( t ) += akf_mat_g2_HVVV( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_HVHV.at( t ) += akf_mat_g2_HVHV( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_HVHH.at( t ) += akf_mat_g2_HVHH( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_HHVH.at( t ) += akf_mat_g2_HHVH( i, tau ) * s.parameters.t_step * deltaT;
//                        //rho_HHHV.at( t ) += akf_mat_g2_HHHV( i, tau ) * s.parameters.t_step * deltaT;
//
//                        rho_VHVH.at( t ) += akf_mat_g2_VHVH( i, tau ) * s.parameters.t_step * deltaT;
//                        rho_VHHV.at( t ) += akf_mat_g2_VHHV( i, tau ) * s.parameters.t_step * deltaT;
//                        rho_HVVH.at( t ) += akf_mat_g2_HVVH( i, tau ) * s.parameters.t_step * deltaT;
//                        rho_HHHH.at( t ) += akf_mat_g2_HHHH( i, tau ) * s.parameters.t_step * deltaT;
//                        rho_VVVV.at( t ) += akf_mat_g2_VVVV( i, tau ) * s.parameters.t_step * deltaT;
//                        rho_HHVV.at( t ) += akf_mat_g2_HHVV( i, tau ) * s.parameters.t_step * deltaT;
//                    }
//                    rho_VVHH.at( t ) = std::conj( rho_HHVV.at( t ) );
//                    rho_HVHV.at( t ) = std::conj( rho_VHVH.at( t ) );
//                } else {
//                    int k = i / mat_step;
//                    rho_HHHH.at( t ) += g2_11_zero.at( k ) * deltaT;
//                    rho_VVVV.at( t ) += g2_22_zero.at( k ) * deltaT;
//                    rho_HHVV.at( t ) += g2_12_zero.at( k ) * deltaT;
//                }
//                Timers::outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, pbsize, "Concurrence: " );
//            }
//            timer_c.iterate();
//        }
//        // Final output and timer end
//        Timers::outputProgress( s.parameters.output_handlerstrings, timer_c, progressbar, pbsize, "Concurrence", PROGRESS_FORCE_OUTPUT );
//        timer_c.end();
//    }
//    // Calculating indistinguishability
//    Dense akf_mat_g1_1 = cache.count() Dense::Zero( dim, dim );
//    Dense akf_mat_g1_2 = Dense::Zero( dim, dim );
//    std::vector<double> p1, p2;
//    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
//        p1.push_back( 0 );
//        p2.push_back( 0 );
//    }
//    if ( !calculate_concurrence_with_g2_of_zero ) {
//        M1 = op_creator_1 * op_annihilator_1;
//        M2 = op_creator_2 * op_annihilator_2;
//        if ( !cache1.isZero() )
//            akf_mat_g1_1 = cache1;
//        else if ( s.parameters.numerics_calculate_g2_H ) //TODO: hier vieleicht weg, nur g1 benutzen wenn auch das spektrum ausgerechnet wurde, sonst g1=0
//            calculate_g1( s, op_creator_1, op_annihilator_1, akf_mat_g1_1, "g1_1" );
//        if ( !cache2.isZero() )
//            akf_mat_g1_2 = cache2;
//        else if ( s.parameters.numerics_calculate_g2_V )
//            calculate_g1( s, op_creator_2, op_annihilator_2, akf_mat_g1_2, "g1_2" );
//        Timer &timer = Timers::create( "Indistinguishability" );
//        timer.start();
//        long unsigned int tstart = output_full_ind ? 0 : savedStates.size() - mat_step;
//#pragma omp parallel for schedule( dynamic ) shared( timer ) num_threads( s.parameters.numerics_maximum_threads )
//        for ( long unsigned int T = tstart; T < savedStates.size(); T += mat_step ) {
//            Scalar top_1 = 0;
//            Scalar top_2 = 0;
//            Scalar bottom_1 = 0;
//            Scalar bottom_2 = 0;
//            for ( int i = 0; i < T; i += mat_step ) {
//                int k = i / mat_step;
//                auto rho = getRhoAt( i );
//                double t = getTimeAt( i );
//                for ( int j = 0; j < T - i; j += 1 ) { //&& j + i < (int)savedStates.size() ////mat_step
//                    int l = j / mat_step;
//                    auto rho_tau = getRhoAt( i + j );
//                    double t_tau = getTimeAt( i + j );
//                    Scalar gpop_1 = s.dgl_expectationvalue<Sparse, Scalar>( rho, M1, t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M1, t_tau );
//                    Scalar gpop_2 = s.dgl_expectationvalue<Sparse, Scalar>( rho, M2, t ) * s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, M2, t_tau );
//                    Scalar gbot_1 = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator_1, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator_1, t );
//                    Scalar gbot_2 = s.dgl_expectationvalue<Sparse, Scalar>( rho_tau, op_annihilator_2, t_tau ) * s.dgl_expectationvalue<Sparse, Scalar>( rho, op_creator_2, t );
//                    top_1 += ( gpop_1 + akf_mat_g2_HHHH( i, j ) - akf_mat_g1_1( i, j ) * std::conj( akf_mat_g1_1( i, j ) ) ); //k,j nicht k, l
//                    bottom_1 += 2.0 * gpop_1 - gbot_1 * std::conj( gbot_1 );
//                    //fmt::print( "what contributes: gpop {}, g2_11 {}, g1_1 {}, gbot {}, TOP {}, BOT {}\n", gpop_1, akf_mat_g2_HHHH( k, l ), std::pow( std::abs( akf_mat_g1_1( k, l ) ), 2.0 ), std::pow( std::abs( gbot_1 ), 2.0 ), top_1, bottom_1 );
//                    top_2 += ( gpop_2 + akf_mat_g2_VVVV( i, j ) - akf_mat_g1_2( i, j ) * std::conj( akf_mat_g1_2( i, j ) ) ); //T,j nicht k,l
//                    bottom_2 += 2.0 * gpop_2 - gbot_2 * std::conj( gbot_2 );
//                }
//                // If we dont output the full time resoluted indistinguishability, we instead output the time resoluted integral for T = T_max
//                if ( !output_full_ind ) {
//                    p1.at( i / mat_step ) = ( 1.0 - std::abs( top_1 / bottom_1 ) );
//                    p2.at( i / mat_step ) = ( 1.0 - std::abs( top_2 / bottom_2 ) );
//                    Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability (Simplified): " );
//                    timer.iterate();
//                }
//            }
//            p1.at( T / mat_step ) = ( 1.0 - std::real( top_1 / bottom_1 ) );
//            p2.at( T / mat_step ) = ( 1.0 - std::real( top_2 / bottom_2 ) );
//            Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability (Full): " );
//            timer.iterate();
//        }
//        // Final output and timer end
//        Timers::outputProgress( s.parameters.output_handlerstrings, timer, progressbar, pbsize, "Indistinguishability: ", PROGRESS_FORCE_OUTPUT );
//        timer.end();
//    }
//    // Save output
//    Log::L2( "Done, saving to {}... ", fileOutputName );
//
//    //REMOVE: Test
//    std::string filepathtwo = s.parameters.subfolder + fileOutputName + "_twophotondensitymatrix.txt";
//    FILE *twophotdensitymatrix = std::fopen( filepathtwo.c_str(), "w" );
//    fmt::print( twophotdensitymatrix, "Time\tVVVV\tVVVH\tVVHV\tVVHH\tVHVV\tVHVH\tVHHV\tVHHH\tHVVV\tHVVH\tHVHV\tHVHH\tHHVV\tHHVH\tHHHV\tHHHH\n" );
//
//    std::string filepath = s.parameters.subfolder + fileOutputName + ".txt";
//    FILE *concurrencefile = std::fopen( filepath.c_str(), "w" );
//    if ( !concurrencefile ) {
//        Log::L2( "Failed to open outputfile for advanced photon statistics!\n" );
//        return false;
//    }
//    fmt::print( concurrencefile, "Time\tMatrixelement_11\tMatrixelement_22\tConcurrence\tG2_11(tau)\tG2_22(tau)\tG2_12(tau)\tG2(tau=0)_11\tG2(tau=0)_22\tG2(tau=0)_12\tI_1\tI_2\n" );
//    Dense spinflip = Dense::Zero( 4, 4 );
//    spinflip( 0, 3 ) = -1;
//    spinflip( 1, 2 ) = 1;
//    spinflip( 2, 1 ) = 1;
//    spinflip( 3, 0 ) = -1;
//    Log::L3( "Spinflip Matrix: {}\n", spinflip );
//    Scalar conc = 0;
//    for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
//        int k = i / mat_step;
//        Log::L3( "Creating 2 photon matrix\n" );
//        Dense rho_2phot = Dense::Zero( 4, 4 );
//        rho_2phot( 0, 0 ) = rho_HHHH[k];
//        rho_2phot( 3, 3 ) = rho_VVVV[k];
//        rho_2phot( 0, 3 ) = rho_VVHH[k];
//        rho_2phot( 3, 0 ) = rho_HHVV[k];
//        rho_2phot( 1, 1 ) = rho_VHHV[k];
//        rho_2phot( 1, 2 ) = rho_VHVH[k];
//        rho_2phot( 2, 2 ) = rho_HVVH[k];
//        rho_2phot( 2, 1 ) = rho_HVHV[k];
//        Log::L3( "Normalizing 2 photon matrix\n" );
//        if ( std::abs( rho_2phot.trace() ) != 0 ) {
//            rho_2phot = rho_2phot / rho_2phot.trace();
//            Log::L3( "Rho_2phot = {}\n", rho_2phot );
//            Log::L3( "Calculating sqrt(rho)\n" );
//            //Eigen::MatrixPower<Dense> Mpow( rho_2phot );
//            Dense sqrtrho2phot = rho_2phot.sqrt(); //Mpow( 0.5 );
//            Log::L3( "Calculating R\n" );
//            Dense R = sqrtrho2phot * spinflip * rho_2phot * spinflip * sqrtrho2phot;
//            Log::L3( "R = {}\n", R );
//            Log::L3( "Calculating sqrt(R)\n" );
//            //Eigen::MatrixPower<Dense> SMPow( R );
//            R = R.sqrt(); //SMPow( 0.5 );
//            Log::L3( "Calculating Eigenvalues\n" );
//            auto eigenvalues = R.eigenvalues();
//            Log::L3( "Eigenvalues at t = {} are {}\n", getTimeAt( i ), eigenvalues );
//            conc = eigenvalues( 3 ) - eigenvalues( 2 ) - eigenvalues( 1 ) - eigenvalues( 0 );
//        }
//        fmt::print( concurrencefile, "{:.5e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\n", getTimeAt( i ), std::abs( rho_HHHH.at( k ) / ( rho_HHHH.at( k ) + rho_VVVV.at( k ) ) ), std::abs( rho_VVVV.at( k ) / ( rho_HHHH.at( k ) + rho_VVVV.at( k ) ) ), std::abs( conc ), std::real( g2_11.at( k ) ), std::real( g2_22.at( k ) ), std::real( g2_12.at( k ) ), std::real( g2_11_zero.at( k ) ), std::real( g2_22_zero.at( k ) ), std::real( g2_12_zero.at( k ) ), p1.at( k ), p2.at( k ) );
//        fmt::print( twophotdensitymatrix, "{:.5e}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", getTimeAt( i ), ( rho_VVVV[k] ), ( rho_VVVH[k] ), ( rho_VVHV[k] ), ( rho_VVHH[k] ), ( rho_VHVV[k] ), ( rho_VHVH[k] ), ( rho_VHHV[k] ), ( rho_VHHH[k] ), ( rho_HVVV[k] ), ( rho_HVVH[k] ), ( rho_HVHV[k] ), ( rho_HVHH[k] ), ( rho_HHVV[k] ), ( rho_HHVH[k] ), ( rho_HHHV[k] ), ( rho_HHHH[k] ) );
//    }
//    std::fclose( concurrencefile );
//
//    std::fclose( twophotdensitymatrix );
//
//    Log::L2( "Done!\n" );
//    Log::L1( "Final Concurrence: {} (non-normalized: {})\n", conc, 2.0 * std::abs( rho_HHVV.back() ) );
//    if ( !calculate_concurrence_with_g2_of_zero )
//        Log::L1( "Final Indistinguishability: {} H, {} V\n", p1.back(), p2.back() );
//    //if ()
//    //Log::L1( "Final G2(0): {}\n", std::real( 2.0 * std::abs( rho_HHVV.back() ) / ( rho_HHHH.back() + rho_VVVV.back() ) ) );
//
//    if ( output_full_g2 ) {
//        Log::L2( "Saving full, non integrated G2(t,tau) matrices... " );
//        FILE *ganz = std::fopen( ( s.parameters.subfolder + fileOutputName + "_g2(ttau)_matrix.txt" ).c_str(), "w" );
//        fmt::print( ganz, "t\ttau\tG2_11(t,tau)\tG2_22(t,tau)\tG2_12(t,tau)\n" );
//        for ( long unsigned int i = 0; i < savedStates.size(); i += mat_step ) {
//            int k = i / mat_step;
//            for ( long unsigned int j = 0; j < savedStates.size(); j += mat_step ) {
//                int l = j / mat_step;
//                fmt::print( ganz, "{0:.5e}\t{1:.5e}\t{2:.8e}\t{3:.8e}\t{4:.8e}\n", getTimeAt( i ), getTimeAt( j ), std::abs( akf_mat_g2_HHHH( k, l ) ), std::abs( akf_mat_g2_VVVV( k, l ) ), std::abs( akf_mat_g2_HHVV( k, l ) ) );
//            }
//            fmt::print( ganz, "\n" );
//        }
//        std::fclose( ganz );
//        Log::L2( "Done!\n" );
//    }
//    // Sucessfull
//    return true;
//}