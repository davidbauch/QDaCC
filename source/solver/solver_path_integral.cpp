#include "solver/solver_ode.h"

Sparse QDLC::Numerics::ODESolver::calculate_propagator_single( System &s, size_t tensor_dim, double t0, double t_step, int i, int j, std::vector<QDLC::SaveState> &output, const Sparse &one ) {
    Sparse projector = Sparse( tensor_dim, tensor_dim );
    projector.coeffRef( i, j ) = 1;
    Sparse M = iterate( projector, s, t0, t_step, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );

    Sparse map;
    if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
        map = QDLC::Matrix::sparse_projector( M );
    }
    for ( double tau = t_step; tau < s.parameters.t_step_pathint; tau += t_step ) {
        M = iterate( M, s, t0 + tau, t_step, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
    }
    if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
        return M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
    }
    M.makeCompressed();
    return M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
}

std::vector<std::vector<Sparse>> &QDLC::Numerics::ODESolver::calculate_propagator_vector( System &s, size_t tensor_dim, double t0, double t_step, std::vector<QDLC::SaveState> &output ) {
    if ( pathint_propagator.count( t0 ) > 0 ) {
        return pathint_propagator[t0];
    }
    Sparse one = Dense::Identity( tensor_dim, tensor_dim ).sparseView();
    std::vector<std::vector<Sparse>> ret( tensor_dim, { tensor_dim, Sparse( tensor_dim, tensor_dim ) } );
    // Calculate first by hand to ensure the Hamilton gets calculated correclty
    ret[0][0] = calculate_propagator_single( s, tensor_dim, t0, t_step, 0, 0, output, one ); //pathint_propagator[-1][0][0] );
// Calculate the remaining propagators
#pragma omp parallel for num_threads( s.parameters.numerics_phonons_maximum_threads )
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( i == 0 && j == 0 ) continue;
            ret[i][j] = calculate_propagator_single( s, tensor_dim, t0, t_step, i, j, output, one ); //pathint_propagator[-1][i][j] );
        }
    }
    // Only save propagator vector if correlation functions are calculated.
    if ( cache.size() > 0 ) {
        Log::L3( "[PathIntegral] Caching propagator vector for t = {}\n", t0 );
        pathint_propagator[t0] = ret;
        return pathint_propagator[t0];
    } else {
        pathint_propagator[-1] = ret;
    }
    return pathint_propagator[-1];
}

bool QDLC::Numerics::ODESolver::calculate_path_integral( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output ) {
    // Generate list of needed G1 and G2 functions. To save on the excessive RAM usage by caching the ADM for every timestep, we calculate the tau-direction for every t step here.
    std::map<std::string, std::vector<Sparse>> g12_settings;
    int g12_counter = 0;
    {
        // Order of Matrices is op1,op3,op2,op4 for g2, 1,1,op1,op2 for g1
        Sparse ident = Dense::Identity( rho0.rows(), rho0.cols() ).sparseView();
        auto &spectrum_s = s.parameters.input_correlation["Spectrum"];
        for ( int i = 0; i < spectrum_s.string_v["Modes"].size(); i++ ) {
            const auto &[s_creator, s_annihilator] = get_operator_strings( spectrum_s.string_v["Modes"][i] );
            std::string g1 = get_operators_purpose( { s_creator, s_annihilator }, 1 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings.count( g1 ) == 0 )
                g12_settings[g1] = { creator, annihilator };
        }
        // Calculate Indist
        auto &indist_s = s.parameters.input_correlation["Indist"];
        for ( int i = 0; i < indist_s.string_v["Modes"].size(); i++ ) {
            const auto &[s_creator, s_annihilator] = get_operator_strings( indist_s.string_v["Modes"][i] );
            std::string g1 = get_operators_purpose( { s_creator, s_annihilator }, 1 );
            std::string g2 = get_operators_purpose( { s_creator, s_annihilator, s_creator, s_annihilator }, 2 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings.count( g1 ) == 0 )
                g12_settings[g1] = { creator, annihilator };
            if ( g12_settings.count( g2 ) == 0 )
                g12_settings[g2] = { creator, annihilator, creator, annihilator };
        }
        // Calculate Conc
        auto &conc_s = s.parameters.input_correlation["Conc"];
        for ( auto &modes : conc_s.string_v["Modes"] ) {
            const auto mode = QDLC::String::splitline( modes, '-' );
            const auto &[s_creator_1, s_annihilator_1] = get_operator_strings( mode[0] );
            const auto &[s_creator_2, s_annihilator_2] = get_operator_strings( mode[1] );

            std::string g2_1111 = get_operators_purpose( { s_creator_1, s_annihilator_1, s_creator_1, s_annihilator_1 }, 2 );
            std::string g2_1212 = get_operators_purpose( { s_creator_1, s_annihilator_1, s_creator_2, s_annihilator_2 }, 2 );
            std::string g2_1221 = get_operators_purpose( { s_creator_1, s_annihilator_2, s_creator_2, s_annihilator_1 }, 2 );
            std::string g2_2112 = get_operators_purpose( { s_creator_2, s_annihilator_1, s_creator_1, s_annihilator_2 }, 2 );
            std::string g2_2211 = get_operators_purpose( { s_creator_2, s_annihilator_1, s_creator_2, s_annihilator_1 }, 2 );
            std::string g2_2222 = get_operators_purpose( { s_creator_2, s_annihilator_2, s_creator_2, s_annihilator_2 }, 2 );

            auto [creator_1, annihilator_1] = get_operators_matrices( s, s_creator_1, s_annihilator_1 );
            auto [creator_2, annihilator_2] = get_operators_matrices( s, s_creator_2, s_annihilator_2 );
            if ( g12_settings.count( g2_1111 ) == 0 )
                g12_settings[g2_1111] = { creator_1, annihilator_1, creator_1, annihilator_1 };
            if ( g12_settings.count( g2_1212 ) == 0 )
                g12_settings[g2_1212] = { creator_1, annihilator_1, creator_1, annihilator_1 };
            if ( g12_settings.count( g2_1221 ) == 0 )
                g12_settings[g2_1221] = { creator_1, annihilator_1, creator_1, annihilator_1 };
            if ( g12_settings.count( g2_2112 ) == 0 )
                g12_settings[g2_2112] = { creator_1, annihilator_1, creator_1, annihilator_1 };
            if ( g12_settings.count( g2_2211 ) == 0 )
                g12_settings[g2_2211] = { creator_1, annihilator_1, creator_1, annihilator_1 };
            if ( g12_settings.count( g2_2222 ) == 0 )
                g12_settings[g2_2222] = { creator_1, annihilator_1, creator_1, annihilator_1 };
        }
        // Calculate G1/G2 functions
        auto &gs_s = s.parameters.input_correlation["GFunc"];
        for ( int i = 0; i < gs_s.string_v["Modes"].size(); i++ ) {
            int order = std::abs( gs_s.numerical_v["Order"][i] );
            const auto &[s_creator, s_annihilator] = get_operator_strings( gs_s.string_v["Modes"][i] );
            std::string g = order == 1 ? get_operators_purpose( { s_creator, s_annihilator }, 1 ) : get_operators_purpose( { s_creator, s_annihilator, s_creator, s_annihilator }, 2 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings.count( g ) == 0 )
                if ( order == 1 )
                    g12_settings[g] = { creator, annihilator };
                else
                    g12_settings[g] = { creator, annihilator, creator, annihilator };
        }
        int matdim = std::min( int( std::floor( ( t_end - t_start ) / s.parameters.t_step_pathint ) / s.parameters.iterations_t_skip ) + 1, s.parameters.iterations_tau_resolution ) + 1;
        for ( auto &[purpose, matrices] : g12_settings ) {
            Log::L2( "[PathIntegral] Calculating G-Function with purpose {} in place with path integral.\n", purpose );
            cache[purpose] = Dense::Zero( matdim, matdim );
            cache[purpose + "_time"] = Dense::Zero( matdim, matdim );
        }
    }

    Log::L2( "[PathIntegral] Setting up Path-Integral Solver...\n" );
    output.reserve( s.parameters.iterations_t_max + 1 );
    int tensor_dim = rho0.rows();
    size_t total_progressbar_iterations = std::floor( ( s.parameters.t_end - s.parameters.t_start ) / s.parameters.t_step_pathint * ( 1 + 0.5 * g12_settings.size() * ( s.parameters.t_end - s.parameters.t_start ) / s.parameters.t_step_pathint ) );

    std::set<int> different_dimensions;
    for ( int i = 0; i < s.operatorMatrices.phononCouplingIndex.size(); i++ ) {
        different_dimensions.insert( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[i] : i );
    }

    pathint_tensor_dimensions = { tensor_dim };
    for ( int i = 1; i < s.parameters.p_phonon_nc; i++ )
        pathint_tensor_dimensions.emplace_back( different_dimensions.size() );

    std::stringstream adm_tensor_dimensions;
    std::copy( std::begin( pathint_tensor_dimensions ), std::end( pathint_tensor_dimensions ), std::ostream_iterator<int>( adm_tensor_dimensions, ", " ) );
    Log::L2( "[PathIntegral] ADM Tensor dimensions will reach: [{}]\n", adm_tensor_dimensions.str() );

    int adm_multithreading_cores = 1; //std::min<int>( s.parameters.numerics_phonons_maximum_threads, different_dimensions.size() );
    Log::L2( "[PathIntegral] Using #{} cpu cores for the ADM propagation.\n", adm_multithreading_cores );

    std::vector<std::vector<int>> adm_multithreaded_indices( different_dimensions.size() );
    for ( int i = 0; i < tensor_dim; i++ ) {
        adm_multithreaded_indices[s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[i] : i].emplace_back( i );
    }
    Log::L2( "[PathIntegral] CPU Threaded vectors are:\n" );
    for ( auto &v : adm_multithreaded_indices ) {
        Log::L2( "[PathIntegral] - {}\n", v );
    }

    Tensor adm_tensor = Tensor( adm_multithreading_cores );

    // First step is just rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );

    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( QDLC::Math::abs2( rho.coeff( i, j ) ) == 0 )
                continue;
            iVector ii = iVector( 1 );
            iVector jj = iVector( 1 );
            ii( 0 ) = i;
            jj( 0 ) = j;
            adm_tensor.addTriplet( ii, jj, rho.coeff( i, j ), 0 );
        }
    }
    adm_tensor.swapAndClearCache( tensor_dim, adm_multithreading_cores );

    Log::L3( "[PathIntegral] Filled Initial Tensor with elements: [" );
    for ( auto &[i0, a0] : adm_tensor.get() )
        for ( auto &[j0, a1] : a0 )
            for ( auto &[sparse_index_x, c0] : a1 )
                for ( auto &[sparse_index_y, value] : c0 ) {
                    Log::L3( "{} ({},{}), ", value, sparse_index_x, sparse_index_y );
                }
    Log::L3( "] ({} elements)\n", adm_tensor.nonZeros() );

    // Iterate Path integral for further time steps
    double t_start_new = t_start; // + s.parameters.t_step_pathint * ( s.parameters.p_phonon_nc - 1 );
    for ( double t_t = t_start_new; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Calculate Correlation functions:
        if ( g12_settings.size() > 0 and g12_counter % s.parameters.iterations_t_skip == 0 ) {
            //auto cached_adms_dimensions = adm_tensor.dimensions;
            for ( auto &[purpose, matrices] : g12_settings ) {
                Log::L3( "[PathIntegral] Calculating sub-rk for {}\n", purpose );
                //adm_tensor.dimensions = cached_adms_dimensions;
                //adm_tensor.rescale_dimensions();
                int order = matrices.size() == 4 ? 2 : 1;
                std::vector<QDLC::SaveState> temp;
                auto &gmat = cache[purpose];
                auto &timemat = cache[purpose + "_time"];
                calculate_path_integral_correlation( adm_tensor, rho, t_t, t_end, t_step_initial, rkTimer, progressbar, total_progressbar_iterations, purpose, s, temp, do_output, matrices, adm_multithreaded_indices, adm_multithreading_cores );
                for ( int32_t j = 0; j < std::min<int32_t>( temp.size(), gmat.rows() * s.parameters.iterations_t_skip ); j += s.parameters.iterations_t_skip ) {
                    double t_tau = temp.at( j ).t;
                    //Log::L3( "Filling gmat for t = {}, tau = {}\n", t_t, t_tau );
                    if ( order == 2 )
                        //gmat( i / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = ( ( matrices[OP2] * matrices[OP3] ).cwiseProduct( temp.at( j ).mat ) ).sum();
                        gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[2] * matrices[1], t_tau );
                    else 
                        //gmat( i / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = ( matrices[1].cwiseProduct( temp.at( j ).mat ) ).sum();
                        gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = std::conj( s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[1], t_tau ) ); // conj ammenakoyum?????
                        //gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[1], t_tau );
                    timemat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = Scalar( t_t, t_tau );
                    
                }
            }
            //adm_tensor.dimensions = cached_adms_dimensions;
            //adm_tensor.rescale_dimensions();
        }
        // Path-Integral iteration

        // Calculate Propagators for current time
        auto t0 = omp_get_wtime();
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, t_step_initial, output );
        t0 = ( omp_get_wtime() - t0 );

        double cur_min = 1;

        // Main Iteration loop
        Log::L3( "[PathIntegral] Current Tensor Size: {} elements at a total of {} MB\n", adm_tensor.nonZeros(), adm_tensor.size() / 1024. / 1024. );
        auto t1 = omp_get_wtime();
        double ts;
        double total_append_time = 0;
        double total_time = 0;

        bool addeddimension = std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) < ( s.parameters.p_phonon_nc - 2 );

#pragma omp parallel for collapse( 2 ) num_threads( adm_multithreading_cores )
        for ( auto &idx : adm_multithreaded_indices )
            for ( auto &jdx : adm_multithreaded_indices )
                for ( auto &i_n_m1 : idx )
                    for ( auto &j_n_m1 : jdx )
                        for ( int l = 0; l < propagator[i_n_m1][j_n_m1].outerSize(); ++l )
                            for ( Sparse::InnerIterator M( propagator[i_n_m1][j_n_m1], l ); M; ++M ) {
                                int i_n = M.row();
                                int j_n = M.col();
                                int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[i_n] : i_n;
                                int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[j_n] : j_n;
                                for ( auto &[sparse_index_x, inner] : adm_tensor.get()[i_n_m1][j_n_m1] )
                                    for ( auto &[sparse_index_y, value] : inner ) {
                                        if ( QDLC::Math::abs2( value ) == 0 ) continue;
                                        auto tt = omp_get_wtime();
                                        //Log::L3( "[PathIntegral] (T{}) handling ({} > {}),({} > {}) --> {}\n", omp_get_thread_num(), gi_n, sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), gj_n, sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), value );
                                        //for ( int l = 0; l < propagator[sparse_index_x( 0 )][sparse_index_y( 0 )].outerSize(); ++l )
                                        //    for ( Sparse::InnerIterator M( propagator[sparse_index_x( 0 )][sparse_index_y( 0 )], l ); M; ++M ) {

                                        //Log::L3( "[PathIntegral] --- Correlation Indices for tau = 0: i = {}, i' = {}, j = {}, j' = {}\n", gi_n, gi_n, gj_n, gj_n );
                                        Scalar phonon_s = s.dgl_phonon_S_function( 0, gi_n, gj_n, gi_n, gj_n );
                                        for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                                            int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_x( 0 )] : sparse_index_x( 0 ) ) : sparse_index_x( tau ) );
                                            int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_y( 0 )] : sparse_index_y( 0 ) ) : sparse_index_y( tau ) );
                                            phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                                            //Log::L3( "[PathIntegral] --- Correlation Indices for tau = {}: i = {}, i' = {}, j = {}, j' = {}\n", tau + 1, gi_n, gi_nd, gj_n, gj_nd );
                                        }

                                        Scalar val = M.value() * value * std::exp( phonon_s );
                                        double abs = QDLC::Math::abs2( val );
                                        auto appendtime = omp_get_wtime();
                                        // Add Element to Triplet list if they are diagonal elements or if they surpass the given threshold.
                                        if ( i_n == j_n || abs >= s.parameters.numerics_pathintegral_squared_threshold ) {
                                            // Add new indices to vectors:
                                            if ( addeddimension ) {
                                                iVector new_sparse_index_x = iVector::Zero( sparse_index_x.size() + 1 );
                                                iVector new_sparse_index_y = iVector::Zero( sparse_index_y.size() + 1 );
                                                for ( int i = 0; i < sparse_index_x.size(); i++ ) {
                                                    new_sparse_index_x( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_x( i )] : sparse_index_x( i ) );
                                                    new_sparse_index_y( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_y( i )] : sparse_index_y( i ) );
                                                }
                                                new_sparse_index_x( 0 ) = i_n;
                                                new_sparse_index_y( 0 ) = j_n;
                                                adm_tensor.addTriplet( new_sparse_index_x, new_sparse_index_y, val, omp_get_thread_num() );
                                            } else {
                                                cur_min = cur_min != 0.0 && cur_min < abs ? cur_min : abs;
                                                if ( s.parameters.numerics_pathint_partially_summed )
                                                    adm_tensor.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n, s.operatorMatrices.phononCouplingIndex[sparse_index_x( 0 )], s.operatorMatrices.phononCouplingIndex[sparse_index_y( 0 )] );
                                                else
                                                    adm_tensor.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n );
                                            }
                                        }
                                        total_append_time += omp_get_wtime() - appendtime;
                                        total_time += omp_get_wtime() - tt;
                                    }
                            }

        t1 = ( omp_get_wtime() - t1 );
        adm_tensor.swapAndClearCache( tensor_dim, adm_multithreading_cores );

        Log::L3( "[PathIntegral] Reducing ADM\n" );
        // Calculate the reduced density matrix by tracing over all past times for each entry
        auto t2 = omp_get_wtime();
        Dense newrho = Dense::Zero( tensor_dim, tensor_dim );

        for ( auto &[i0, a0] : adm_tensor.get() )
            for ( auto &[j0, a1] : a0 )
                for ( auto &[sparse_index_x, c0] : a1 )
                    for ( auto &[sparse_index_y, value] : c0 ) {
                        int i_n = sparse_index_x( 0 );
                        int j_n = sparse_index_y( 0 );
                        newrho( i_n, j_n ) += value;
                    }

        //rho = newrho.sparseView() / newrho.trace();
        rho = newrho.sparseView(); // / newrho.trace();
        t2 = ( omp_get_wtime() - t2 );
        Log::L3( "[PathIntegral] Iteration: {}, time taken: [ Propagator: {:.4f}s, ADM Advancing: {:.4f}s (Partial append time: {:.4f}\%), ADM Setting (Parallel): {:.4f}s, ADM Reduction: {:.4f}s ], Trace: {}, Elements: {}\n", t_t, t0, t1, 100.0 * total_append_time / total_time, ts, t2, s.getTrace<Scalar>( rho ), adm_tensor.nonZeros() );

        // Dynamic Cutoff
        if ( s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max > 0 ) {
            double ratio = (double)adm_tensor.nonZeros() / s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max;
            ratio = std::exp( ratio ) * std::pow( ratio, 5 );              // Scaled squared
            s.parameters.numerics_pathintegral_squared_threshold *= ratio; // Adjust Cutoff by squared ratio
            if ( s.parameters.numerics_pathintegral_squared_threshold < 1E-30 ) s.parameters.numerics_pathintegral_squared_threshold = 1E-30;
            if ( s.parameters.numerics_pathintegral_squared_threshold > 1E4 ) s.parameters.numerics_pathintegral_squared_threshold = 1E4;
            Log::L3( "[PathIntegral] Adjusted Cutoff Threshold to {} (x{})\n", std::sqrt( s.parameters.numerics_pathintegral_squared_threshold ), ratio );
        }

        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );
        // Progress and time output
        rkTimer.iterate();
        g12_counter++;
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
        }
    }

    // TODO MAYBE: interpolate density matrices to t_step instead of t_step_path -> indist, spectrum integral more smooth!
    // additionally, correlation functions do not need to concern about timestep

    Log::L3( "Done!\n" );
    return true;
}

bool QDLC::Numerics::ODESolver::calculate_path_integral_correlation( Tensor adm_correlation, Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, size_t total_progressbar_iterations, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output, const std::vector<Sparse> &matrices, const std::vector<std::vector<int>> &adm_multithreaded_indices, int adm_multithreading_cores ) {
    int order = matrices.size() == 4 ? 2 : 1;
    Log::L3( "- [PathIntegralCorrelation] Setting up Path-Integral Solver...\n" );
    int tensor_dim = rho0.rows();
    output.reserve( s.parameters.iterations_t_max + 1 );
    //saveState( rho0, t_start, output );
    bool modified = false;
    //saveState(rho0,t_start,output);

    for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Path-Integral iteration
        // Calculate Propagators for current time
        auto t0 = omp_get_wtime();
        auto propagator = calculate_propagator_vector( s, tensor_dim, t_t, t_step_initial, output );
        t0 = ( omp_get_wtime() - t0 );
        double cur_min = 1;

        // DAS hier funktioniert mit propagator als ref nicht, weil der dann bearbeitet wird!!
        if ( not modified ) {
            modified = true;
            for ( int i = 0; i < propagator.size(); i++ ) {
                for ( int j = 0; j < propagator[i].size(); j++ ) {
                    if ( order == 2 )
                        propagator[i][j] = s.dgl_timetrafo( matrices[3], t_start ) * propagator[i][j] * s.dgl_timetrafo( matrices[0], t_start );
                    else
                        propagator[i][j] = propagator[i][j] * s.dgl_timetrafo( matrices[0], t_start );
                }
            }
        }

        bool addeddimension = std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) < ( s.parameters.p_phonon_nc - 2 );

        // Main Iteration loop
        Log::L3( "- [PathIntegralCorrelation] Current Correlation Tensor Size: {} elements at a total of {} MB\r", adm_correlation.nonZeros(), adm_correlation.size() / 1024. / 1024. );
        auto t1 = omp_get_wtime();
        double total_append_time = 0;
        double total_time = 0;


#pragma omp parallel for collapse( 2 ) num_threads( adm_multithreading_cores )
        for ( auto &idx : adm_multithreaded_indices )
            for ( auto &jdx : adm_multithreaded_indices )
                for ( auto &i_n_m1 : idx )
                    for ( auto &j_n_m1 : jdx )
                        for ( int l = 0; l < propagator[i_n_m1][j_n_m1].outerSize(); ++l )
                            for ( Sparse::InnerIterator M( propagator[i_n_m1][j_n_m1], l ); M; ++M ) {
                                int i_n = M.row();
                                int j_n = M.col();
                                int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[i_n] : i_n;
                                int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[j_n] : j_n;
                                for ( auto &[sparse_index_x, inner] : adm_correlation.get()[i_n_m1][j_n_m1] )
                                    for ( auto &[sparse_index_y, value] : inner ) {
                                        if ( QDLC::Math::abs2( value ) == 0 ) continue;
                                        auto tt = omp_get_wtime();
                                        Scalar phonon_s = s.dgl_phonon_S_function( 0, gi_n, gj_n, gi_n, gj_n );
                                        for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                                            int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_x( 0 )] : sparse_index_x( 0 ) ) : sparse_index_x( tau ) );
                                            int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_y( 0 )] : sparse_index_y( 0 ) ) : sparse_index_y( tau ) );
                                            phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                                        }
                                        Scalar val = M.value() * value * std::exp( phonon_s );
                                        double abs = QDLC::Math::abs2( val );
                                        auto appendtime = omp_get_wtime();
                                        // Add Element to Triplet list if they are diagonal elements or if they surpass the given threshold.
                                        if ( i_n == j_n || abs >= s.parameters.numerics_pathintegral_squared_threshold ) {
                                            // Add new indices to vectors:
                                            if ( addeddimension ) {
                                                iVector new_sparse_index_x = iVector::Zero( sparse_index_x.size() + 1 );
                                                iVector new_sparse_index_y = iVector::Zero( sparse_index_y.size() + 1 );
                                                for ( int i = 0; i < sparse_index_x.size(); i++ ) {
                                                    new_sparse_index_x( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_x( i )] : sparse_index_x( i ) );
                                                    new_sparse_index_y( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_y( i )] : sparse_index_y( i ) );
                                                }
                                                new_sparse_index_x( 0 ) = i_n;
                                                new_sparse_index_y( 0 ) = j_n;
                                                adm_correlation.addTriplet( new_sparse_index_x, new_sparse_index_y, val, omp_get_thread_num() );
                                            } else {
                                                cur_min = cur_min != 0.0 && cur_min < abs ? cur_min : abs;
                                                if ( s.parameters.numerics_pathint_partially_summed )
                                                    adm_correlation.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n, s.operatorMatrices.phononCouplingIndex[sparse_index_x( 0 )], s.operatorMatrices.phononCouplingIndex[sparse_index_y( 0 )] );
                                                else
                                                    adm_correlation.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n );
                                            }
                                        }
                                        total_append_time += omp_get_wtime() - appendtime;
                                        total_time += omp_get_wtime() - tt;
                                    }
                            }

        t1 = ( omp_get_wtime() - t1 );

        auto ts = omp_get_wtime();
        adm_correlation.swapAndClearCache( tensor_dim, adm_multithreading_cores );
        ts = omp_get_wtime() - ts;

        Log::L3( "- [PathIntegralCorrelation] Reducing Correlation ADM\n" );
        // Calculate the reduced density matrix by tracing over all past times for each entry
        auto t2 = omp_get_wtime();
        Dense cache = Dense::Zero( tensor_dim, tensor_dim );
        for ( auto &[i0, a0] : adm_correlation.get() )
            for ( auto &[j0, a1] : a0 )
                for ( auto &[sparse_index_x, c0] : a1 )
                    for ( auto &[sparse_index_y, value] : c0 ) {
                        int i_n = sparse_index_x( 0 );
                        int j_n = sparse_index_y( 0 );
                        cache( i_n, j_n ) += value;
                    }
        auto rho = cache.sparseView();
        t2 = ( omp_get_wtime() - t2 );
        Log::L3( "- [PathIntegralCorrelation] Correlation Iteration: {}, time taken: [ Propagator: {:.4f}s, ADM Advancing: {:.4f}s (Partial append time: {:.4f}\%), ADM Setting: {:.4f}s, ADM Reduction: {:.4f}s ]\n", t_t, t0, t1, 100.0 * total_append_time / total_time, ts, t2 );

        // Dynamic Cutoff
        if ( s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max > 0 ) {
            double ratio = (double)adm_correlation.nonZeros() / s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max;
            ratio = std::exp( ratio ) * std::pow( ratio, 5 );              // Scaled squared
            s.parameters.numerics_pathintegral_squared_threshold *= ratio; // Adjust Cutoff by squared ratio
            if ( s.parameters.numerics_pathintegral_squared_threshold < 1E-30 ) s.parameters.numerics_pathintegral_squared_threshold = 1E-30;
            if ( s.parameters.numerics_pathintegral_squared_threshold > 1E4 ) s.parameters.numerics_pathintegral_squared_threshold = 1E4;
            Log::L3( "- [PathIntegralCorrelation] Adjusted Cutoff Threshold to {} (x{})\n", std::sqrt( s.parameters.numerics_pathintegral_squared_threshold ), ratio );
        }
        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
        }
    }
    Log::L3( "- [PathIntegralCorrelation] Correlation Done!\n" );
    return true;
}