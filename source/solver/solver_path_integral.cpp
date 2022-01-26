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
    ret[0][0] = calculate_propagator_single( s, tensor_dim, t0, t_step, 0, 0, output, one ); // pathint_propagator[-1][0][0] );
// Calculate the remaining propagators
#pragma omp parallel for num_threads( s.parameters.numerics_phonons_maximum_threads )
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( i == 0 && j == 0 ) continue;
            ret[i][j] = calculate_propagator_single( s, tensor_dim, t0, t_step, i, j, output, one ); // pathint_propagator[-1][i][j] );
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

    int adm_multithreading_cores = s.parameters.numerics_phonons_maximum_threads; // std::min<int>( s.parameters.numerics_phonons_maximum_threads, different_dimensions.size() );
    Log::L2( "[PathIntegral] Using #{} cpu cores for the ADM propagation.\n", adm_multithreading_cores );
    Log::L2( "[PathIntegral] Anticipated Tensor size will be {} MB (elements) and {} MB (indices).\n", std::pow( tensor_dim, 2.0 ) * std::pow( different_dimensions.size(), 2 * s.parameters.p_phonon_nc - 2 ) * 16 / 1024. / 1024., std::pow( tensor_dim, 2.0 ) * std::pow( different_dimensions.size(), 2 * s.parameters.p_phonon_nc - 2 ) * 4 * pathint_tensor_dimensions.size() / 1024. / 1024. );

    Tensor<Scalar> adm_tensor( pathint_tensor_dimensions );
    Tensor<int> adm_skip( pathint_tensor_dimensions, 0 );
    // for ( auto &[sparse_index_x, map] : adm_skip.getNextValues() )
    //     for ( auto &[sparse_index_y, val] : map )
    //         Log::L3( "{} {} -> {}\n", sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), val );

    std::vector<std::tuple<iVector, iVector>> total_indices;
    for ( auto sparse_index_x : adm_tensor.getIndices() )
        for ( auto sparse_index_y : adm_tensor.getIndices() ) {
            total_indices.emplace_back( std::make_tuple( sparse_index_x, sparse_index_y ) );
        }

    // First step is just rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );

    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( QDLC::Math::abs2( rho.coeff( i, j ) ) == 0 )
                continue;
            iVector ii = iVector::Zero( pathint_tensor_dimensions.size() );
            iVector jj = iVector::Zero( pathint_tensor_dimensions.size() );
            ii( 0 ) = i;
            jj( 0 ) = j;
            adm_tensor.setTriplet( ii, jj, rho.coeff( i, j ) );
        }
    }
    adm_tensor.swap();

    Log::L2( "[PathIntegral] Filled Initial Tensor with elements: [" );
    for ( auto &[sparse_index_x, map] : adm_tensor.getCurrentValues() )
        for ( auto &[sparse_index_y, value] : map )
            if ( QDLC::Math::abs2( value ) != 0.0 )
                Log::L2( "{} ([{}],[{}]), ", value, sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ) );

    Log::L2( "] ({} total elements)\n", adm_tensor.nonZeros() );

    // Smart Tensor cutoff index. Either directly after nc steps. if the Tensor fillrate stays the same for this amount of iterations, assume zero elements will stay zero
    int tensor_smart_cutoff_index = s.parameters.p_phonon_nc + 1;
    double tensor_smart_cutoff_significant = 0.0; // number of significant digits.
    double tensor_smart_cutoff_last_fillrate = 0.0;
    int tensor_smart_cutoff_current = 0;
    Log::L2( "[PathIntegral] Smart Tensor Cutoff Index is {}.\n", tensor_smart_cutoff_index );

    // Iterate Path integral for further time steps
    double t_start_new = t_start; // + s.parameters.t_step_pathint * ( s.parameters.p_phonon_nc - 1 );
    for ( double t_t = t_start_new; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Calculate Correlation functions:
        if ( g12_settings.size() > 0 and g12_counter % s.parameters.iterations_t_skip == 0 ) {
            // auto cached_adms_dimensions = adm_tensor.dimensions;
            for ( auto &[purpose, matrices] : g12_settings ) {
                Log::L3( "[PathIntegral] Calculating sub-rk for {}\n", purpose );
                // adm_tensor.dimensions = cached_adms_dimensions;
                // adm_tensor.rescale_dimensions();
                int order = matrices.size() == 4 ? 2 : 1;
                std::vector<QDLC::SaveState> temp;
                auto &gmat = cache[purpose];
                auto &timemat = cache[purpose + "_time"];
                calculate_path_integral_correlation( adm_tensor, rho, t_t, t_end, t_step_initial, rkTimer, progressbar, total_progressbar_iterations, purpose, s, temp, do_output, matrices, adm_multithreading_cores );
                for ( int32_t j = 0; j < std::min<int32_t>( temp.size(), gmat.rows() * s.parameters.iterations_t_skip ); j += s.parameters.iterations_t_skip ) {
                    double t_tau = temp.at( j ).t;
                    // Log::L3( "Filling gmat for t = {}, tau = {}\n", t_t, t_tau );
                    if ( order == 2 )
                        // gmat( i / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = ( ( matrices[OP2] * matrices[OP3] ).cwiseProduct( temp.at( j ).mat ) ).sum();
                        gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[2] * matrices[1], t_tau );
                    else
                        // gmat( i / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = ( matrices[1].cwiseProduct( temp.at( j ).mat ) ).sum();
                        // gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[1], t_tau );
                        gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = std::conj( s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[1], t_tau ) ); // conj ammenakoyum?????
                    timemat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = Scalar( t_t, t_tau );
                }
            }
            // adm_tensor.dimensions = cached_adms_dimensions;
            // adm_tensor.rescale_dimensions();
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

        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );
        int nonzero = 0;
        // for ( auto sparse_index_x : adm_tensor.getIndices() )collapse( 2 )
        std::set<int> indices_to_delete;
#pragma omp parallel for num_threads( adm_multithreading_cores ) schedule( guided ) shared( nonzero )
        for ( int _i = 0; _i < total_indices.size(); _i++ ) {
            auto [sparse_index_x, sparse_index_y] = total_indices[_i];
            // if ( max_index >= sparse_index_x.size() and adm_skip.getNextTriplet( sparse_index_x, sparse_index_y ) > s.parameters.p_phonon_nc ) {
            //     // Log::L3( "Skipping element {} - {}, counter is {}\n", sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), adm_skip.getNextTriplet( sparse_index_x, sparse_index_y ) );
            //     continue;
            // }
            Scalar new_value = 0.0;
            // Set all indices to zero that are for timevalues larger than max_index
            for ( int i = max_index; i < sparse_index_x.size(); i++ ) {
                sparse_index_x( i ) = 0;
                sparse_index_y( i ) = 0;
            }

            // Indices
            int i_n = sparse_index_x( 0 );
            int j_n = sparse_index_y( 0 );

            // Groups:
            int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[i_n] : i_n;
            int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[j_n] : j_n;
            int gi_n_m1 = sparse_index_x( 1 );
            int gj_n_m1 = sparse_index_y( 1 );
            auto temp = adm_tensor.getTriplet( sparse_index_x, sparse_index_y );

            // Log::L3( "[PathIntegral] (T{}) calculating new [{}] - [{}] --> old value is {}\n", omp_get_thread_num(), sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), temp );

            //  Create "Old" sparse indices (vector shifted, last value will be replaced by summation) (n,n-1,n-2,...,n-m) -> (n-1,n-2,n-3,...,n-m,placeholder)
            auto sparse_index_x_old = sparse_index_x;
            auto sparse_index_y_old = sparse_index_y;
            for ( int i = 0; i < sparse_index_x_old.size() - 1; i++ ) {
                sparse_index_x_old( i ) = sparse_index_x_old( i + 1 );
                sparse_index_y_old( i ) = sparse_index_y_old( i + 1 );
            }
            // Sum over all States in group (sum_(k_n-1,kd_n-1))
            for ( int i_n_m1 : s.operatorMatrices.phononGroupToIndices[gi_n_m1] ) {
                for ( int j_n_m1 : s.operatorMatrices.phononGroupToIndices[gj_n_m1] ) {
                    // Switch first entry of "old" sparse index back to the actual state index
                    sparse_index_x_old( 0 ) = i_n_m1;
                    sparse_index_y_old( 0 ) = j_n_m1;
                    // Propagator value
                    Scalar propagator_value = propagator[i_n_m1][j_n_m1].coeff( i_n, j_n );
                    // If this propagation is not allowed, continue
                    if ( QDLC::Math::abs2( propagator_value ) == 0.0 )
                        continue;
                    // Sum over Groups lambda_n-n_m:
                    for ( int lambda_i = 0; lambda_i < different_dimensions.size(); lambda_i++ ) {
                        for ( int lambda_j = 0; lambda_j < different_dimensions.size(); lambda_j++ ) {
                            // Old index vector:
                            sparse_index_x_old( max_index - 1 ) = lambda_i;
                            sparse_index_y_old( max_index - 1 ) = lambda_j;
                            // Calculate S:
                            Scalar phonon_s = s.dgl_phonon_S_function( 0, gi_n, gj_n, gi_n, gj_n );
                            for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                                int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_x_old( 0 )] : sparse_index_x_old( 0 ) ) : sparse_index_x_old( tau ) );
                                int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_y_old( 0 )] : sparse_index_y_old( 0 ) ) : sparse_index_y_old( tau ) );
                                phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                            }
                            new_value += propagator_value * adm_tensor.getTriplet( sparse_index_x_old, sparse_index_y_old ) * std::exp( phonon_s );
                        }
                    }
                }
            }
            if ( QDLC::Math::abs2( new_value ) != 0.0 ) {
                nonzero++;
            } else {
                if ( max_index >= sparse_index_x.size() ) {
                    adm_skip.addToTriplet( sparse_index_x, sparse_index_y, 1 );
                    int next = adm_skip.getNextTriplet( sparse_index_x, sparse_index_y );
                    if ( next > tensor_smart_cutoff_index - sparse_index_x.size() ) {
#pragma omp critical
                        indices_to_delete.emplace( _i );
                    }
                }
            }
            // Log::L3( "[PathIntegral] (T{}) final value for new [{}] - [{}] {}\n", omp_get_thread_num(), sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), new_value );

            adm_tensor.setTriplet( sparse_index_x, sparse_index_y, new_value );
        }

        t1 = ( omp_get_wtime() - t1 );
        adm_tensor.swap();

        Log::L3( "[PathIntegral] Reducing ADM\n" );
        // Calculate the reduced density matrix by tracing over all past times for each entry
        auto t2 = omp_get_wtime();
        Dense newrho = Dense::Zero( tensor_dim, tensor_dim );

        for ( auto &[sparse_index_x, map] : adm_tensor.getCurrentValues() )
            for ( auto &[sparse_index_y, value] : map ) {
                int i_n = sparse_index_x( 0 );
                int j_n = sparse_index_y( 0 );
                newrho( i_n, j_n ) += value;
            }

        // rho = newrho.sparseView() / newrho.trace();
        rho = newrho.sparseView(); // / newrho.trace();
        t2 = ( omp_get_wtime() - t2 );
        double current_tensor_fillrate = 100.0 * nonzero / ( std::pow( tensor_dim, 2 ) * std::pow( different_dimensions.size(), 2 * s.parameters.p_phonon_nc - 2 ) );
        Log::L3( "[PathIntegral] Iteration: {}, time taken: [ Propagator: {:.4f}s, ADM Advancing: {:.4f}s (Partial append time: {:.4f}\%), ADM Setting (Parallel): {:.4f}s, ADM Reduction: {:.4f}s ], Trace: {}, Elements: {} ({} pct Fillrate)\n", t_t, t0, t1, 100.0 * total_append_time / total_time, ts, t2, s.getTrace<Scalar>( rho ), nonzero, current_tensor_fillrate );

        // Dynamic Cutoff
        if ( s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max > 0 ) {
            double ratio = (double)adm_tensor.nonZeros() / s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max;
            ratio = std::exp( ratio ) * std::pow( ratio, 5 );              // Scaled squared
            s.parameters.numerics_pathintegral_squared_threshold *= ratio; // Adjust Cutoff by squared ratio
            if ( s.parameters.numerics_pathintegral_squared_threshold < 1E-30 ) s.parameters.numerics_pathintegral_squared_threshold = 1E-30;
            if ( s.parameters.numerics_pathintegral_squared_threshold > 1E4 ) s.parameters.numerics_pathintegral_squared_threshold = 1E4;
            Log::L3( "[PathIntegral] Adjusted Cutoff Threshold to {} (x{})\n", std::sqrt( s.parameters.numerics_pathintegral_squared_threshold ), ratio );
        }

        // Smart Cutoff
        if ( indices_to_delete.size() > 0 ) {
            current_tensor_fillrate = std::floor( current_tensor_fillrate * std::pow( 10, tensor_smart_cutoff_significant ) ) / tensor_smart_cutoff_significant;
            if ( current_tensor_fillrate == tensor_smart_cutoff_last_fillrate )
                tensor_smart_cutoff_current++;
            else
                tensor_smart_cutoff_current = 0;
            if ( tensor_smart_cutoff_current >= tensor_smart_cutoff_index ) {
                if ( indices_to_delete.size() > 0 ) {
                    Log::L3( "Deleting {} indices\n", indices_to_delete.size() );
                    std::vector<std::tuple<iVector, iVector>> total_indices_new;
                    for ( int _i = 0; _i < total_indices.size(); _i++ ) {
                        if ( !indices_to_delete.contains( _i ) )
                            total_indices_new.emplace_back( total_indices[_i] );
                    }
                    total_indices = total_indices_new;
                    tensor_smart_cutoff_current = 0;
                }
            }
            tensor_smart_cutoff_last_fillrate = current_tensor_fillrate;
        }

        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );
        // Progress and time output
        rkTimer.iterate();
        g12_counter++;
        if ( do_output ) {
            Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
        }
    }
    Log::L3( "Done!\n" );
    return true;
}

bool QDLC::Numerics::ODESolver::calculate_path_integral_correlation( Tensor<Scalar> adm_correlation, Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, size_t total_progressbar_iterations, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output, const std::vector<Sparse> &matrices, int adm_multithreading_cores ) {
    int order = matrices.size() == 4 ? 2 : 1;
    Log::L3( "- [PathIntegralCorrelation] Setting up Path-Integral Solver...\n" );
    int tensor_dim = rho0.rows();
    output.reserve( s.parameters.iterations_t_max + 1 );
    // saveState( rho0, t_start, output );
    bool modified = false;
    // saveState(rho0,t_start,output);

    std::set<int> different_dimensions;
    for ( int i = 0; i < s.operatorMatrices.phononCouplingIndex.size(); i++ ) {
        different_dimensions.insert( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[i] : i );
    }

    Tensor<int> adm_skip( pathint_tensor_dimensions, 0 );
    std::vector<std::tuple<iVector, iVector>> total_indices;
    for ( auto sparse_index_x : adm_correlation.getIndices() )
        for ( auto sparse_index_y : adm_correlation.getIndices() ) {
            total_indices.emplace_back( std::make_tuple( sparse_index_x, sparse_index_y ) );
        }
    // Smart Tensor cutoff index. Either directly after nc steps. if the Tensor fillrate stays the same for this amount of iterations, assume zero elements will stay zero
    int tensor_smart_cutoff_index = s.parameters.p_phonon_nc + 1;
    double tensor_smart_cutoff_significant = 0.0; // number of significant digits.
    double tensor_smart_cutoff_last_fillrate = 0.0;
    int tensor_smart_cutoff_current = 0;

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

        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );
        int nonzero = 0;
        std::set<int> indices_to_delete;
#pragma omp parallel for num_threads( adm_multithreading_cores ) schedule( guided ) shared( nonzero )
        for ( int _i = 0; _i < total_indices.size(); _i++ ) {
            auto [sparse_index_x, sparse_index_y] = total_indices[_i];
            Scalar new_value = 0.0;
            // Set all indices to zero that are for timevalues larger than max_index
            for ( int i = max_index; i < sparse_index_x.size(); i++ ) {
                sparse_index_x( i ) = 0;
                sparse_index_y( i ) = 0;
            }

            // Indices
            int i_n = sparse_index_x( 0 );
            int j_n = sparse_index_y( 0 );

            // Groups:
            int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[i_n] : i_n;
            int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[j_n] : j_n;
            int gi_n_m1 = sparse_index_x( 1 );
            int gj_n_m1 = sparse_index_y( 1 );
            auto temp = adm_correlation.getTriplet( sparse_index_x, sparse_index_y );

            // Log::L3( "[PathIntegral] (T{}) calculating new [{}] - [{}] --> old value is {}\n", omp_get_thread_num(), sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), temp );

            //  Create "Old" sparse indices (vector shifted, last value will be replaced by summation) (n,n-1,n-2,...,n-m) -> (n-1,n-2,n-3,...,n-m,placeholder)
            auto sparse_index_x_old = sparse_index_x;
            auto sparse_index_y_old = sparse_index_y;
            for ( int i = 0; i < sparse_index_x_old.size() - 1; i++ ) {
                sparse_index_x_old( i ) = sparse_index_x_old( i + 1 );
                sparse_index_y_old( i ) = sparse_index_y_old( i + 1 );
            }
            // Sum over all States in group (sum_(k_n-1,kd_n-1))
            for ( int i_n_m1 : s.operatorMatrices.phononGroupToIndices[gi_n_m1] ) {
                for ( int j_n_m1 : s.operatorMatrices.phononGroupToIndices[gj_n_m1] ) {
                    // Switch first entry of "old" sparse index back to the actual state index
                    sparse_index_x_old( 0 ) = i_n_m1;
                    sparse_index_y_old( 0 ) = j_n_m1;
                    // Propagator value
                    Scalar propagator_value = propagator[i_n_m1][j_n_m1].coeff( i_n, j_n );
                    // If this propagation is not allowed, continue
                    if ( QDLC::Math::abs2( propagator_value ) == 0.0 )
                        continue;
                    // Sum over Groups lambda_n-n_m:
                    for ( int lambda_i = 0; lambda_i < different_dimensions.size(); lambda_i++ ) {
                        for ( int lambda_j = 0; lambda_j < different_dimensions.size(); lambda_j++ ) {
                            // Old index vector:
                            sparse_index_x_old( max_index - 1 ) = lambda_i;
                            sparse_index_y_old( max_index - 1 ) = lambda_j;
                            // Calculate S:
                            Scalar phonon_s = s.dgl_phonon_S_function( 0, gi_n, gj_n, gi_n, gj_n );
                            for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                                int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_x_old( 0 )] : sparse_index_x_old( 0 ) ) : sparse_index_x_old( tau ) );
                                int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingIndex[sparse_index_y_old( 0 )] : sparse_index_y_old( 0 ) ) : sparse_index_y_old( tau ) );
                                phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                            }
                            new_value += propagator_value * adm_correlation.getTriplet( sparse_index_x_old, sparse_index_y_old ) * std::exp( phonon_s );
                        }
                    }
                }
            }
            if ( QDLC::Math::abs2( new_value ) != 0.0 ) {
                nonzero++;
            } else {
                if ( max_index >= sparse_index_x.size() ) {
                    adm_skip.addToTriplet( sparse_index_x, sparse_index_y, 1 );
                    int next = adm_skip.getNextTriplet( sparse_index_x, sparse_index_y );
                    if ( next > tensor_smart_cutoff_index - sparse_index_x.size() ) {
#pragma omp critical
                        indices_to_delete.emplace( _i );
                    }
                }
            }
            // Log::L3( "[PathIntegral] (T{}) final value for new [{}] - [{}] {}\n", omp_get_thread_num(), sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), new_value );

            adm_correlation.setTriplet( sparse_index_x, sparse_index_y, new_value );
        }

        t1 = ( omp_get_wtime() - t1 );

        auto ts = omp_get_wtime();
        adm_correlation.swap();
        ts = omp_get_wtime() - ts;

        Log::L3( "- [PathIntegralCorrelation] Reducing Correlation ADM\n" );
        // Calculate the reduced density matrix by tracing over all past times for each entry
        auto t2 = omp_get_wtime();
        Dense cache = Dense::Zero( tensor_dim, tensor_dim );
        for ( auto &[sparse_index_x, map] : adm_correlation.getCurrentValues() )
            for ( auto &[sparse_index_y, value] : map ) {
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

        double current_tensor_fillrate = 100.0 * nonzero / ( std::pow( tensor_dim, 2 ) * std::pow( different_dimensions.size(), 2 * s.parameters.p_phonon_nc - 2 ) );
        // Smart Cutoff
        if ( indices_to_delete.size() > 0 ) {
            current_tensor_fillrate = std::floor( current_tensor_fillrate * std::pow( 10, tensor_smart_cutoff_significant ) ) / tensor_smart_cutoff_significant;
            if ( current_tensor_fillrate == tensor_smart_cutoff_last_fillrate )
                tensor_smart_cutoff_current++;
            else
                tensor_smart_cutoff_current = 0;
            if ( tensor_smart_cutoff_current >= tensor_smart_cutoff_index ) {
                if ( indices_to_delete.size() > 0 ) {
                    Log::L3( "Deleting {} indices\n", indices_to_delete.size() );
                    std::vector<std::tuple<iVector, iVector>> total_indices_new;
                    for ( int _i = 0; _i < total_indices.size(); _i++ ) {
                        if ( !indices_to_delete.contains( _i ) )
                            total_indices_new.emplace_back( total_indices[_i] );
                    }
                    total_indices = total_indices_new;
                    tensor_smart_cutoff_current = 0;
                }
            }
            tensor_smart_cutoff_last_fillrate = current_tensor_fillrate;
        }

        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
        }
    }
    Log::L3( "- [PathIntegralCorrelation] Correlation Done!\n" );
    return true;
}