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
    if ( pathint_propagator.contains( t0 ) ) {
        return pathint_propagator[t0];
    }
    Sparse one = Dense::Identity( tensor_dim, tensor_dim ).sparseView();
    std::vector<std::vector<Sparse>> ret( tensor_dim, { tensor_dim, Sparse( tensor_dim, tensor_dim ) } );
    // Calculate first by hand to ensure the Hamilton gets calculated correclty
    ret[0][0] = calculate_propagator_single( s, tensor_dim, t0, t_step, 0, 0, output, one ); // pathint_propagator[-1][0][0] );
// Calculate the remaining propagators
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_secondary_threads )
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

// TODO: statt sparse -> dense am besten dense -> nonzeros zu indices appenden, ab punkt dann einfach nur noch über nonzero indices iterieren. dann kann man von anfang an multithreaded machen, irgendwann gucken was is ungleich null,
// dann nur noch über diese indices summieren.

bool QDLC::Numerics::ODESolver::calculate_path_integral( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output ) {
    // Generate list of needed G1 and G2 functions. To save on the excessive RAM usage by caching the ADM for every timestep, we calculate the tau-direction for every t step here.
    std::map<std::string, std::vector<Sparse>> g12_settings;
    int g12_counter = 0;
    {
        // Order of Matrices is op1,op3,op2,op4 for g2, 1,1,op1,op2 for g1
        Sparse ident = Dense::Identity( rho0.rows(), rho0.cols() ).sparseView();
        auto &spectrum_s = s.parameters.input_correlation["Spectrum"];
        for ( int i = 0; i < spectrum_s.string_v["Modes"].size(); i++ ) {
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, spectrum_s.string_v["Modes"][i] );
            std::string g1 = get_operators_purpose( { s_creator, s_annihilator }, 1 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings.count( g1 ) == 0 )
                g12_settings[g1] = { ident, ident, creator, annihilator };
        }
        // Calculate Indist
        auto &indist_s = s.parameters.input_correlation["Indist"];
        for ( int i = 0; i < indist_s.string_v["Modes"].size(); i++ ) {
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, indist_s.string_v["Modes"][i] );
            std::string g1 = get_operators_purpose( { s_creator, s_annihilator }, 1 );
            std::string g2 = get_operators_purpose( { s_creator, s_annihilator, s_creator, s_annihilator }, 2 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings.count( g1 ) == 0 )
                g12_settings[g1] = { ident, ident, creator, annihilator };
            if ( g12_settings.count( g2 ) == 0 )
                g12_settings[g2] = { creator, annihilator, creator, annihilator };
        }
        // Calculate Conc
        auto &conc_s = s.parameters.input_correlation["Conc"];
        for ( auto &modes : conc_s.string_v["Modes"] ) {
            const auto mode = QDLC::String::splitline( modes, '-' );
            const auto &[s_creator_1, s_annihilator_1] = get_operator_strings( s, mode[0] );
            const auto &[s_creator_2, s_annihilator_2] = get_operator_strings( s, mode[1] );

            std::string g2_1111 = get_operators_purpose( { s_creator_1, s_annihilator_1, s_creator_1, s_annihilator_1 }, 2 );
            std::string g2_1122 = get_operators_purpose( { s_creator_1, s_annihilator_2, s_creator_1, s_annihilator_2 }, 2 );
            std::string g2_2121 = get_operators_purpose( { s_creator_2, s_annihilator_2, s_creator_1, s_annihilator_1 }, 2 );
            std::string g2_1221 = get_operators_purpose( { s_creator_1, s_annihilator_2, s_creator_2, s_annihilator_1 }, 2 );
            std::string g2_2112 = get_operators_purpose( { s_creator_2, s_annihilator_1, s_creator_1, s_annihilator_2 }, 2 );
            std::string g2_2222 = get_operators_purpose( { s_creator_2, s_annihilator_2, s_creator_2, s_annihilator_2 }, 2 );
            // std::string g2_1212 = get_operators_purpose( { s_creator_1, s_annihilator_1, s_creator_2, s_annihilator_2 }, 2 );
            // std::string g2_2211 = get_operators_purpose( { s_creator_2, s_annihilator_1, s_creator_2, s_annihilator_1 }, 2 );

            auto [creator_1, annihilator_1] = get_operators_matrices( s, s_creator_1, s_annihilator_1 );
            auto [creator_2, annihilator_2] = get_operators_matrices( s, s_creator_2, s_annihilator_2 );
            if ( g12_settings.count( g2_1111 ) == 0 )
                g12_settings[g2_1111] = { creator_1, annihilator_1, creator_1, annihilator_1 };
            if ( g12_settings.count( g2_2121 ) == 0 )
                g12_settings[g2_2121] = { creator_2, annihilator_2, creator_1, annihilator_1 };
            if ( g12_settings.count( g2_1221 ) == 0 )
                g12_settings[g2_1221] = { creator_1, annihilator_2, creator_2, annihilator_1 };
            if ( g12_settings.count( g2_2112 ) == 0 )
                g12_settings[g2_2112] = { creator_2, annihilator_1, creator_1, annihilator_2 };
            if ( g12_settings.count( g2_1122 ) == 0 )
                g12_settings[g2_1122] = { creator_1, annihilator_2, creator_1, annihilator_2 };
            if ( g12_settings.count( g2_2222 ) == 0 )
                g12_settings[g2_2222] = { creator_2, annihilator_2, creator_2, annihilator_2 };
        }
        // Calculate G1/G2 functions
        auto &gs_s = s.parameters.input_correlation["GFunc"];
        for ( int i = 0; i < gs_s.string_v["Modes"].size(); i++ ) {
            int order = std::abs( gs_s.numerical_v["Order"][i] );
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, gs_s.string_v["Modes"][i] );
            std::string g = order == 1 ? get_operators_purpose( { s_creator, s_annihilator }, 1 ) : get_operators_purpose( { s_creator, s_annihilator, s_creator, s_annihilator }, 2 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings.count( g ) == 0 )
                if ( order == 1 )
                    g12_settings[g] = { ident, ident, creator, annihilator };
                else
                    g12_settings[g] = { creator, annihilator, creator, annihilator };
        }
        int matdim = std::min( int( std::floor( ( t_end - t_start ) / s.parameters.t_step_pathint ) / s.parameters.iterations_t_skip ) + 1, s.parameters.grid_resolution ) + 1;
        // int matdim = std::min( int( std::floor( ( t_end - t_start ) / s.parameters.t_step ) / s.parameters.iterations_t_skip ) + 1, s.parameters.grid_resolution ) + 1;
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
    for ( int i = 0; i < s.operatorMatrices.phonon_hilbert_index_to_group_index.size(); i++ ) {
        different_dimensions.insert( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i] : i );
    }

    pathint_tensor_dimensions = { tensor_dim };
    for ( int i = 1; i < s.parameters.p_phonon_nc; i++ )
        pathint_tensor_dimensions.emplace_back( different_dimensions.size() );

    std::stringstream adm_tensor_dimensions;
    std::copy( std::begin( pathint_tensor_dimensions ), std::end( pathint_tensor_dimensions ), std::ostream_iterator<int>( adm_tensor_dimensions, ", " ) );
    Log::L2( "[PathIntegral] ADM Tensor dimensions will reach: [{}]\n", adm_tensor_dimensions.str() );

    Log::L2( "[PathIntegral] Tensor Type is {}\n", s.parameters.numerics_pathintegral_force_dense ? "Dense" : "Sparse" );
    Log::L2( "[PathIntegral] Using a maximum of #{} cpu cores for the ADM propagation.\n", s.parameters.numerics_maximum_secondary_threads );
    Log::L2( "[PathIntegral] Anticipated Tensor size will be {} MB (elements) and {} MB (indices).\n", std::pow( tensor_dim, 2.0 ) * std::pow( different_dimensions.size(), 2 * s.parameters.p_phonon_nc - 2 ) * 16 / 1024. / 1024., std::pow( tensor_dim, 2.0 ) * std::pow( different_dimensions.size(), 2 * s.parameters.p_phonon_nc - 2 ) * 4 * pathint_tensor_dimensions.size() / 1024. / 1024. );

    // Tensor Class
    bool numerics_pathintegral_force_sparse = false;
    bool numerics_pathintegral_cleanse_insteadof_switch = true;
    Tensor<Scalar> adm_tensor( pathint_tensor_dimensions, s.parameters.numerics_pathintegral_force_dense ? Tensor<Scalar>::TYPE_DENSE : Tensor<Scalar>::TYPE_SPARSE );

    // First step is just rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );

    // Fill initial Tensor
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( QDLC::Math::abs2( rho.coeff( i, j ) ) == 0 )
                continue;
            iVector ii;
            iVector jj;
            if ( adm_tensor.isDense() ) {
                ii = iVector::Zero( pathint_tensor_dimensions.size() );
                jj = iVector::Zero( pathint_tensor_dimensions.size() );
            } else {
                ii = iVector::Zero( 1 );
                jj = iVector::Zero( 1 );
            }
            ii( 0 ) = i;
            jj( 0 ) = j;
            adm_tensor.setTriplet( ii, jj, rho.coeff( i, j ) );
        }
    }
    adm_tensor.swap();

    std::stringstream els;
    for ( auto &[sparse_index_x, map] : adm_tensor.getCurrentValues() )
        for ( auto &[sparse_index_y, value] : map )
            if ( QDLC::Math::abs2( value ) != 0.0 )
                els << fmt::format( "{} ([{}],[{}]), ", value, sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ) );
    Log::L2( "[PathIntegral] Filled Initial Tensor with elements: [{}] ({} total elements)\n", els.str(), adm_tensor.nonZeros() );

    bool filled_time = false;
    int nonzero = 0;
    int new_nonzero = 0;
    int numerics_dynamic_densitychange_counter = 0;

    // #################################################
    // ############ Temp: Simple Profiling #############
    // #################################################
    std::map<std::string, std::vector<double>> profiler_time_per_thread;
    profiler_time_per_thread["PI_propagator"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_main_iteration"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_sparse_iteration"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_dense_index_gathering"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_dense_group_summation_and_correlation_function"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_dense_value_setting"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_sparse_index_gathering"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_sparse_group_summation"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_sparse_correlation_function"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_adm_reduction"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_tensor_swap"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_tensor_prune"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_total"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );

    // Iterate Path integral for further time steps
    for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        double profiler_total = omp_get_wtime();
        // Calculate Correlation functions:
        if ( g12_settings.size() > 0 and g12_counter % s.parameters.iterations_t_skip == 0 ) {
            // auto cached_adms_dimensions = adm_tensor.dimensions;
            for ( auto &[purpose, matrices] : g12_settings ) {
                Log::L3( "[PathIntegral] Calculating sub-rk for {}\n", purpose );
                // adm_tensor.dimensions = cached_adms_dimensions;
                // adm_tensor.rescale_dimensions();
                std::vector<QDLC::SaveState> temp;
                auto &gmat = cache[purpose];
                auto &timemat = cache[purpose + "_time"];
                Log::L3( "[PathIntegral] Calculating Path Integral Correlation from t_0 = {} to t_1 = {}, initial tensor is {}...\n", t_t, t_end, adm_tensor.isSparse() ? "Sparse" : "Dense" );
                calculate_path_integral_correlation( adm_tensor, rho, t_t, t_end, t_step_initial, rkTimer, progressbar, total_progressbar_iterations, purpose, s, temp, do_output, matrices, s.parameters.numerics_maximum_secondary_threads );
                // Log::L3( "[PathIntegral] Interpolating Results...\n" );
                // if ( temp.size() > 1 )
                //     temp = Numerics::interpolate_curve( temp, t_t, s.parameters.t_end, s.parameters.grid_values, s.parameters.grid_steps, s.parameters.grid_value_indices, false, s.parameters.numerics_interpolate_method_tau );
                Log::L3( "[PathIntegral] Writing {} values to G matrix...\n", temp.size() );
                for ( int32_t j = 0; j < std::min<int32_t>( temp.size(), gmat.rows() * s.parameters.iterations_t_skip ); j += s.parameters.iterations_t_skip ) {
                    double t_tau = temp.at( j ).t;
                    // Log::L3( "Filling gmat for t = {}, tau = {}\n", t_t, t_tau );
                    // if ( matrices.size() == 4 ) {
                    // gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = ( ( matrices[2] * matrices[1] ).cwiseProduct( temp.at( j ).mat ) ).sum();
                    gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[2] * matrices[1], t_tau );
                    //} else {
                    // gmat( i / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = ( matrices[1].cwiseProduct( temp.at( j ).mat ) ).sum();
                    // gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = 5.0; // s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[1], t_tau );
                    // gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = std::conj( s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[1], t_tau ) ); // conj ammenakoyum?????
                    //}
                }
                if ( not filled_time )
                    for ( int32_t i = 0; i < gmat.rows() * s.parameters.iterations_t_skip; i += s.parameters.iterations_t_skip ) {
                        double t1 = i < temp.size() ? temp[i].t : temp.back().t + ( 1. + i - temp.size() ) * s.parameters.t_step_pathint;
                        for ( int32_t j = 0; j < gmat.rows() * s.parameters.iterations_t_skip; j += s.parameters.iterations_t_skip ) {
                            double t2 = j < temp.size() ? temp[j].t : temp.back().t + ( 1. + j - temp.size() ) * s.parameters.t_step_pathint;
                            timemat( i / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = Scalar( t1, t1 + t2 );
                        }
                    }
            }
            filled_time = true;
        }

        // Actual Path-Integral iteration

        // Calculate Propagators for current time
        /* PROFILER */ double profiler_time = omp_get_wtime();
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );
        /* PROFILER */ profiler_time_per_thread["PI_propagator"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        double cur_min = 1;
        new_nonzero = 0;

        // Main Iteration loop
        Log::L3( "[PathIntegral] Current Tensor Size: {} elements at a total of {} MB\n", adm_tensor.nonZeros(), adm_tensor.size() / 1024. / 1024. );

        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );
        // Iterate the tensor
        /* PROFILER */ profiler_time = omp_get_wtime();
        // Dense Tensor Iteration
        if ( adm_tensor.isDense() ) {
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_secondary_threads ) schedule( guided ) shared( nonzero )
            for ( auto &index : adm_tensor.getIndices() ) {
                /* PROFILER */ double profiler_d = omp_get_wtime();
                auto [sparse_index_x, sparse_index_y] = index;
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
                int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i_n] : i_n;
                int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[j_n] : j_n;
                int gi_n_m1 = sparse_index_x( 1 );
                int gj_n_m1 = sparse_index_y( 1 );

                //  Create "Old" sparse indices (vector shifted, last value will be replaced by summation) (n,n-1,n-2,...,n-m) -> (n-1,n-2,n-3,...,n-m,placeholder)
                auto sparse_index_x_old = sparse_index_x;
                auto sparse_index_y_old = sparse_index_y;
                for ( int i = 0; i < sparse_index_x_old.size() - 1; i++ ) {
                    sparse_index_x_old( i ) = sparse_index_x_old( i + 1 );
                    sparse_index_y_old( i ) = sparse_index_y_old( i + 1 );
                }
                /* PROFILER */ profiler_time_per_thread["PI_dense_index_gathering"][omp_get_thread_num()] += omp_get_wtime() - profiler_d;

                /* PROFILER */ profiler_d = omp_get_wtime();
                // Sum over all States in group (sum_(k_n-1,kd_n-1))
                for ( int i_n_m1 : s.operatorMatrices.phonon_group_index_to_hilbert_indices[gi_n_m1] ) {
                    for ( int j_n_m1 : s.operatorMatrices.phonon_group_index_to_hilbert_indices[gj_n_m1] ) {
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
                                    int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_x_old( 0 )] : sparse_index_x_old( 0 ) ) : sparse_index_x_old( tau ) );
                                    int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_y_old( 0 )] : sparse_index_y_old( 0 ) ) : sparse_index_y_old( tau ) );
                                    phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                                }
                                new_value += propagator_value * adm_tensor.getTriplet( sparse_index_x_old, sparse_index_y_old ) * std::exp( phonon_s );
                            }
                        }
                    }
                }
                /* PROFILER */ profiler_time_per_thread["PI_dense_group_summation_and_correlation_function"][omp_get_thread_num()] += omp_get_wtime() - profiler_d;
                if ( QDLC::Math::abs2( new_value ) != 0.0 ) {
                    new_nonzero++;
                }

                /* PROFILER */ profiler_d = omp_get_wtime();
                // Log::L3( "[PathIntegral] (T{}) final value for new [{}] - [{}] {}\n", omp_get_thread_num(), sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), new_value );
                adm_tensor.setTriplet( sparse_index_x, sparse_index_y, new_value );
                /* PROFILER */ profiler_time_per_thread["PI_dense_value_setting"][omp_get_thread_num()] += omp_get_wtime() - profiler_d;
            }
        } else {
            bool addeddimension = std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) < ( s.parameters.p_phonon_nc - 1 );
            for ( auto &[sparse_index_x, map] : adm_tensor.getCurrentValues() )
                for ( auto &[sparse_index_y, value] : map ) {
                    /* PROFILER */ double profiler_id = omp_get_wtime();
                    int i_n_m1 = sparse_index_x( 0 );
                    int j_n_m1 = sparse_index_y( 0 );
                    for ( int l = 0; l < propagator[i_n_m1][j_n_m1].outerSize(); ++l )
                        for ( Sparse::InnerIterator M( propagator[i_n_m1][j_n_m1], l ); M; ++M ) {
                            /* PROFILER */ double profiler_d = omp_get_wtime();
                            int i_n = M.row();
                            int j_n = M.col();
                            int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i_n] : i_n;
                            int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[j_n] : j_n;
                            if ( QDLC::Math::abs2( value ) == 0 ) continue;
                            // Log::L3( "[PathIntegral] (T{}) handling ({} > {}),({} > {}) --> {}\n", omp_get_thread_num(), gi_n, sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), gj_n, sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), value );
                            // for ( int l = 0; l < propagator[sparse_index_x( 0 )][sparse_index_y( 0 )].outerSize(); ++l )
                            //     for ( Sparse::InnerIterator M( propagator[sparse_index_x( 0 )][sparse_index_y( 0 )], l ); M; ++M ) {

                            /* PROFILER */ profiler_time_per_thread["PI_sparse_index_gathering"][omp_get_thread_num()] += omp_get_wtime() - profiler_d;

                            /* PROFILER */ profiler_d = omp_get_wtime();
                            // Log::L3( "[PathIntegral] --- Correlation Indices for tau = 0: i = {}, i' = {}, j = {}, j' = {}\n", gi_n, gi_n, gj_n, gj_n );
                            Scalar phonon_s = s.dgl_phonon_S_function( 0, gi_n, gj_n, gi_n, gj_n );
                            for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                                int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_x( 0 )] : sparse_index_x( 0 ) ) : sparse_index_x( tau ) );
                                int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_y( 0 )] : sparse_index_y( 0 ) ) : sparse_index_y( tau ) );
                                phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd ); // TODO: das hier in map cachen
                                // Log::L3( "[PathIntegral] --- Correlation Indices for tau = {}: i = {}, i' = {}, j = {}, j' = {}\n", tau + 1, gi_n, gi_nd, gj_n, gj_nd );
                            }
                            /* PROFILER */ profiler_time_per_thread["PI_sparse_correlation_function"][omp_get_thread_num()] += omp_get_wtime() - profiler_d;

                            /* PROFILER */ profiler_d = omp_get_wtime();
                            Scalar val = M.value() * value * std::exp( phonon_s );
                            double abs = QDLC::Math::abs2( val );
                            // Add Element to Triplet list if they are diagonal elements or if they surpass the given threshold.
                            if ( i_n == j_n || abs >= s.parameters.numerics_pathintegral_squared_threshold ) {
                                // Add new indices to vectors:
                                if ( addeddimension ) {
                                    iVector new_sparse_index_x = iVector::Zero( sparse_index_x.size() + 1 );
                                    iVector new_sparse_index_y = iVector::Zero( sparse_index_y.size() + 1 );
                                    for ( int i = 0; i < sparse_index_x.size(); i++ ) {
                                        new_sparse_index_x( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_x( i )] : sparse_index_x( i ) );
                                        new_sparse_index_y( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_y( i )] : sparse_index_y( i ) );
                                    }
                                    new_sparse_index_x( 0 ) = i_n;
                                    new_sparse_index_y( 0 ) = j_n;
                                    adm_tensor.addToTriplet( new_sparse_index_x, new_sparse_index_y, val );
                                } else {
                                    cur_min = cur_min != 0.0 && cur_min < abs ? cur_min : abs;
                                    if ( s.parameters.numerics_pathint_partially_summed )
                                        adm_tensor.addToTriplet( sparse_index_x, sparse_index_y, val, i_n, j_n, s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_x( 0 )], s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_y( 0 )] );
                                    else
                                        adm_tensor.addToTriplet( sparse_index_x, sparse_index_y, val, i_n, j_n );
                                }
                            }
                            /* PROFILER */ profiler_time_per_thread["PI_sparse_group_summation"][omp_get_thread_num()] += omp_get_wtime() - profiler_d;
                        }
                    /* PROFILER */ profiler_time_per_thread["PI_sparse_iteration"][omp_get_thread_num()] += omp_get_wtime() - profiler_id;
                    new_nonzero = adm_tensor.nonZeros();
                }
        }
        // Dynamically switching Tensor to Dense or Prune the Dense Tensor to enable fast multithreading
        Log::L3( "[PathIntegral] Current Fillchange Factor is {}\n", std::abs<double>( 1. - 1.0 * new_nonzero / nonzero ) );
        /* PROFILER */ profiler_time = omp_get_wtime();
        if ( ( ( adm_tensor.isSparse() and not numerics_pathintegral_force_sparse ) or ( adm_tensor.isDense() and numerics_pathintegral_cleanse_insteadof_switch ) ) and ( std::abs<double>( 1. - 1.0 * new_nonzero / nonzero ) <= s.parameters.numerics_pathintegral_sparse_to_dense_threshold or new_nonzero < nonzero ) ) {
            numerics_dynamic_densitychange_counter++;
            if ( numerics_dynamic_densitychange_counter >= s.parameters.numerics_dynamic_densitychange_limit and t_t / s.parameters.t_step_pathint > s.parameters.p_phonon_nc ) {
                if ( numerics_pathintegral_cleanse_insteadof_switch ) {
                    numerics_pathintegral_cleanse_insteadof_switch = false;
                    Log::L2( "[PathIntegral] Pruning Dense Tensor at t = {} with {} elements\n", t_t, adm_tensor.nonZeros() );
                    numerics_dynamic_densitychange_counter = 0;
                    adm_tensor.prune();
                    Log::L2( "[PathIntegral] Dense Tensor remaining at {} elements\n", adm_tensor.nonZeros() );
                } else {
                    Log::L2( "[PathIntegral] Switching to Dense Tensor at t = {} with {} elements\n", t_t, adm_tensor.nonZeros() );
                    adm_tensor.convertToDense();
                }
            }
        } else {
            numerics_dynamic_densitychange_counter = 0;
        }
        /* PROFILER */ profiler_time_per_thread["PI_tensor_prune"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;
        nonzero = new_nonzero;

        /* PROFILER */ profiler_time_per_thread["PI_main_iteration"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        /* PROFILER */ profiler_time = omp_get_wtime();
        adm_tensor.swap();
        /* PROFILER */ profiler_time_per_thread["PI_tensor_swap"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        Log::L3( "[PathIntegral] Reducing ADM\n" );
        // Calculate the reduced density matrix by tracing over all past times for each entry
        /* PROFILER */ profiler_time = omp_get_wtime();
        Dense newrho = Dense::Zero( tensor_dim, tensor_dim );
        for ( auto &[sparse_index_x, map] : adm_tensor.getCurrentValues() )
            for ( auto &[sparse_index_y, value] : map ) {
                // for ( auto &index : adm_tensor.getIndices() ) {
                //     auto &[sparse_index_x, sparse_index_y] = index;
                //     auto &value = adm_tensor.getCurrentValues()[sparse_index_x][sparse_index_y];
                int i_n = sparse_index_x( 0 );
                int j_n = sparse_index_y( 0 );
                newrho( i_n, j_n ) += value;
            }
        rho = newrho.sparseView(); // / newrho.trace();
        /* PROFILER */ profiler_time_per_thread["PI_adm_reduction"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        double current_tensor_fillrate = 100.0 * nonzero / ( std::pow( tensor_dim, 2 ) * std::pow( different_dimensions.size(), 2 * s.parameters.p_phonon_nc - 2 ) );
        // Log::L3( "[PathIntegral] Iteration: {}, time taken: [ Propagator: {:.4f}s, ADM Advancing: {:.4f}s (Partial append time: {:.4f}\%), ADM Setting (Parallel): {:.4f}s, ADM Reduction: {:.4f}s ], Trace: {}, Elements: {} ({} pct Fillrate)\n", t_t, t0, t1, 100.0 * total_append_time / total_time, ts, t2, s.get_trace<Scalar>( rho ), nonzero, current_tensor_fillrate );

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
        // saveState( s.dgl_timetrafo( rho, t_t + s.parameters.t_step_pathint ), t_t + s.parameters.t_step_pathint, output );
        saveState( rho, t_t + s.parameters.t_step_pathint, output );

        // Output Profiler Information
        /* PROFILER */ profiler_time_per_thread["PI_total"][omp_get_thread_num()] += omp_get_wtime() - profiler_total;
        Log::L3( "[PathIntegral-Profiler] Time Profile for t = {}, #Nonzeros: {}, TensorFillrate: {}\%\n", t_t, nonzero, current_tensor_fillrate );
        double total_time = std::accumulate( profiler_time_per_thread["PI_total"].begin(), profiler_time_per_thread["PI_total"].end(), 0.0 );
        for ( auto &[name, vec] : profiler_time_per_thread ) {
            double cur_time = 0;
            std::for_each( vec.begin(), vec.end(), [&]( double &num ) { cur_time+=num; num = 0; } );
            Log::L3( "[PathIntegral-Profiler]     Name: {} - Time taken: {}s - {}\% of total time\n", name, cur_time, 100.0 * cur_time / total_time );
        }

        // Progress and time output
        rkTimer.iterate();
        g12_counter++;
        if ( do_output ) {
            Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
        }
    }
    return true;
}

// TODO: tensor nicht kopieren, sondern bei erster iteration cache von neuem und values von altem tensor nehmen.
bool QDLC::Numerics::ODESolver::calculate_path_integral_correlation( Tensor<Scalar> adm_correlation, Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, size_t total_progressbar_iterations, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output, const std::vector<Sparse> &matrices, int adm_multithreading_cores ) {
    Log::L3( "- [PathIntegralCorrelation] Setting up Path-Integral Solver...\n" );
    int tensor_dim = rho0.rows();
    output.reserve( s.parameters.iterations_t_max + 1 );
    // saveState( rho0, t_start, output );
    bool modified = false;
    int nonzero = 0;
    // saveState(rho0,t_start,output);

    if ( not s.parameters.numerics_pathintegral_force_dense )
        adm_correlation.convertToSparse();

    std::set<int> different_dimensions;
    for ( int i = 0; i < s.operatorMatrices.phonon_hilbert_index_to_group_index.size(); i++ ) {
        different_dimensions.insert( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i] : i );
    }
    int numerics_dynamic_densitychange_counter = 0;
    for ( double t_t = t_start; t_t <= t_end; t_t += s.parameters.t_step_pathint ) {
        // Path-Integral iteration
        // Calculate Propagators for current time
        auto t0 = omp_get_wtime();
        auto propagator = calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );
        t0 = ( omp_get_wtime() - t0 );
        double cur_min = 1;

        // DAS hier funktioniert mit propagator als ref nicht, weil der dann bearbeitet wird!!
        if ( not modified ) {
            modified = true;
            for ( int i = 0; i < propagator.size(); i++ ) {
                for ( int j = 0; j < propagator[i].size(); j++ ) {
                    // if ( order == 2 )
                    propagator[i][j] = s.dgl_timetrafo( matrices[3] * propagator[i][j] * matrices[0], t_start );
                    // propagator[i][j] = matrices[3] * ( propagator[i][j] * matrices[0] );
                    // Sparse newPropagator = Sparse( tensor_dim, tensor_dim );
                    // for ( int n = 0; n < tensor_dim; n++ ) {
                    //    for ( int m = 0; m < tensor_dim; m++ ) {
                    //        for ( int nd = 0; nd < tensor_dim; nd++ ) {
                    //            for ( int md = 0; md < tensor_dim; md++ ) {
                    //                newPropagator.coeffRef( n, m ) += matrices[3].coeff( n, nd ) * propagator[i][j].coeffRef( nd, md ) * matrices[0].coeff( md, m );
                    //            }
                    //        }
                    //    }
                    //}
                    // propagator[i][j] = newPropagator;
                    // else
                    //     propagator[i][j] = propagator[i][j] * s.dgl_timetrafo( matrices[0], t_start );
                }
            }
            // Dense cache = Dense::Zero( tensor_dim, tensor_dim );
            // auto ds1 = Dense( s.dgl_timetrafo( matrices[3], t_start ) );
            // auto ds2 = Dense( s.dgl_timetrafo( matrices[0], t_start ) );
            // for ( auto &[sparse_index_x, map] : adm_correlation.getCurrentValues() )
            //     for ( auto &[sparse_index_y, value] : map ) {
            //         int i_n = sparse_index_x( 0 );
            //         int j_n = sparse_index_y( 0 );
            //         Scalar propagator_sum = 0;
            //         for ( int di = 0; di < ds1.rows(); di++ )
            //             for ( int dj = 0; dj < ds1.rows(); dj++ )
            //                 propagator_sum += ds1( i_n, di ) * ds2( dj, j_n );
            //         cache( i_n, j_n ) += value * propagator_sum;
            //     }
            // auto rho = cache.sparseView();
            //  Save Rho
            //  saveState( rho, t_t, output );
            //  if ( adm_correlation.isSparseTensor() ) {
            //      use_dense_tensor = true;
            //      adm_correlation.convertToDense();
            //      adm_multithreading_cores = s.parameters.numerics_maximum_secondary_threads;
            //  }
        }

        // Main Iteration loop
        Log::L3( "- [PathIntegralCorrelation] Current Correlation Tensor Size: {} elements at a total of {} MB, number of indices: {}\r", adm_correlation.nonZeros(), adm_correlation.size() / 1024. / 1024., adm_correlation.getIndices().size() );
        auto t1 = omp_get_wtime();
        double total_append_time = 0;
        double total_time = 0;

        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );
        // Iterate the tensor
        if ( adm_correlation.isDense() ) {
            nonzero = 0;
#pragma omp parallel for num_threads( adm_multithreading_cores ) schedule( guided ) shared( nonzero )
            for ( auto &index : adm_correlation.getIndices() ) {
                auto [sparse_index_x, sparse_index_y] = index;
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
                int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i_n] : i_n;
                int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[j_n] : j_n;
                int gi_n_m1 = sparse_index_x( 1 );
                int gj_n_m1 = sparse_index_y( 1 );

                //  Create "Old" sparse indices (vector shifted, last value will be replaced by summation) (n,n-1,n-2,...,n-m) -> (n-1,n-2,n-3,...,n-m,placeholder)
                auto sparse_index_x_old = sparse_index_x;
                auto sparse_index_y_old = sparse_index_y;
                for ( int i = 0; i < sparse_index_x_old.size() - 1; i++ ) {
                    sparse_index_x_old( i ) = sparse_index_x_old( i + 1 );
                    sparse_index_y_old( i ) = sparse_index_y_old( i + 1 );
                }

                // Sum over all States in group (sum_(k_n-1,kd_n-1))
                for ( int i_n_m1 : s.operatorMatrices.phonon_group_index_to_hilbert_indices[gi_n_m1] ) {
                    for ( int j_n_m1 : s.operatorMatrices.phonon_group_index_to_hilbert_indices[gj_n_m1] ) {
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
                                    int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_x_old( 0 )] : sparse_index_x_old( 0 ) ) : sparse_index_x_old( tau ) );
                                    int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_y_old( 0 )] : sparse_index_y_old( 0 ) ) : sparse_index_y_old( tau ) );
                                    phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                                }
                                new_value += propagator_value * adm_correlation.getTriplet( sparse_index_x_old, sparse_index_y_old ) * std::exp( phonon_s );
                            }
                        }
                    }
                }
                if ( QDLC::Math::abs2( new_value ) != 0.0 ) {
                    nonzero++;
                }

                // Log::L3( "[PathIntegralCorrelation] (T{}) final value for new [{}] - [{}] {}\n", omp_get_thread_num(), sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), new_value );
                adm_correlation.setTriplet( sparse_index_x, sparse_index_y, new_value );
            }
        } else {
            bool addeddimension = std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) < ( s.parameters.p_phonon_nc - 1 );
            for ( auto &[sparse_index_x, map] : adm_correlation.getCurrentValues() )
                for ( auto &[sparse_index_y, value] : map ) {
                    int i_n_m1 = sparse_index_x( 0 );
                    int j_n_m1 = sparse_index_y( 0 );
                    for ( int l = 0; l < propagator[i_n_m1][j_n_m1].outerSize(); ++l )
                        for ( Sparse::InnerIterator M( propagator[i_n_m1][j_n_m1], l ); M; ++M ) {
                            int i_n = M.row();
                            int j_n = M.col();
                            int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i_n] : i_n;
                            int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[j_n] : j_n;
                            if ( QDLC::Math::abs2( value ) == 0 ) continue;
                            // Log::L3( "[PathIntegralCorrelation] (T{}) handling ({} > {}),({} > {}) --> {}\n", omp_get_thread_num(), gi_n, sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), gj_n, sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), value );
                            // for ( int l = 0; l < propagator[sparse_index_x( 0 )][sparse_index_y( 0 )].outerSize(); ++l )
                            //     for ( Sparse::InnerIterator M( propagator[sparse_index_x( 0 )][sparse_index_y( 0 )], l ); M; ++M ) {

                            // Log::L3( "[PathIntegralCorrelation] --- Correlation Indices for tau = 0: i = {}, i' = {}, j = {}, j' = {}\n", gi_n, gi_n, gj_n, gj_n );
                            Scalar phonon_s = s.dgl_phonon_S_function( 0, gi_n, gj_n, gi_n, gj_n );
                            for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                                int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_x( 0 )] : sparse_index_x( 0 ) ) : sparse_index_x( tau ) );
                                int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_y( 0 )] : sparse_index_y( 0 ) ) : sparse_index_y( tau ) );
                                phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd ); // TODO: das hier in map cachen
                                // Log::L3( "[PathIntegralCorrelation] --- Correlation Indices for tau = {}: i = {}, i' = {}, j = {}, j' = {}\n", tau + 1, gi_n, gi_nd, gj_n, gj_nd );
                            }

                            Scalar val = M.value() * value * std::exp( phonon_s );
                            double abs = QDLC::Math::abs2( val );
                            // Add Element to Triplet list if they are diagonal elements or if they surpass the given threshold.
                            if ( i_n == j_n || abs >= s.parameters.numerics_pathintegral_squared_threshold ) {
                                // Add new indices to vectors:
                                if ( addeddimension ) {
                                    iVector new_sparse_index_x = iVector::Zero( sparse_index_x.size() + 1 );
                                    iVector new_sparse_index_y = iVector::Zero( sparse_index_y.size() + 1 );
                                    for ( int i = 0; i < sparse_index_x.size(); i++ ) {
                                        new_sparse_index_x( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_x( i )] : sparse_index_x( i ) );
                                        new_sparse_index_y( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_y( i )] : sparse_index_y( i ) );
                                    }
                                    new_sparse_index_x( 0 ) = i_n;
                                    new_sparse_index_y( 0 ) = j_n;
                                    adm_correlation.addToTriplet( new_sparse_index_x, new_sparse_index_y, val );
                                } else {
                                    cur_min = cur_min != 0.0 && cur_min < abs ? cur_min : abs;
                                    if ( s.parameters.numerics_pathint_partially_summed )
                                        adm_correlation.addToTriplet( sparse_index_x, sparse_index_y, val, i_n, j_n, s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_x( 0 )], s.operatorMatrices.phonon_hilbert_index_to_group_index[sparse_index_y( 0 )] );
                                    else
                                        adm_correlation.addToTriplet( sparse_index_x, sparse_index_y, val, i_n, j_n );
                                }
                            }
                        }
                }
            auto new_nonzero = adm_correlation.nonZeros();
            // Dynamically switching Tensor to Dense to enable fast multithreading
            Log::L3( "[PathIntegralCorrelation] Current Fillchange Factor is {} ({}, {})\n", std::abs<double>( 1. - 1.0 * new_nonzero / nonzero ), new_nonzero, nonzero );
            if ( adm_correlation.isSparse() and std::abs<double>( 1. - 1.0 * new_nonzero / nonzero ) <= s.parameters.numerics_pathintegral_sparse_to_dense_threshold ) {
                numerics_dynamic_densitychange_counter++;
                if ( numerics_dynamic_densitychange_counter >= s.parameters.numerics_dynamic_densitychange_limit and t_t / s.parameters.t_step_pathint > s.parameters.p_phonon_nc ) {
                    Log::L2( "[PathIntegralCorrelation] Switching to Dense Tensor with {} elements\n", adm_correlation.nonZeros() );
                    adm_correlation.convertToDense();
                    adm_multithreading_cores = s.parameters.numerics_maximum_secondary_threads;
                }
            } else {
                numerics_dynamic_densitychange_counter = 0;
            }
            nonzero = new_nonzero;
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
        Sparse rho = cache.sparseView();
        t2 = ( omp_get_wtime() - t2 );
        Log::L3( "- [PathIntegralCorrelation] RhoNonZ = {}, TensorNonZ = {}, Correlation Iteration: {}, time taken: [ Propagator: {:.4f}s, ADM Advancing: {:.4f}s (Partial append time: {:.4f}\%), ADM Setting: {:.4f}s, ADM Reduction: {:.4f}s ]\n", rho.nonZeros(), nonzero, t_t, t0, t1, 100.0 * total_append_time / total_time, ts, t2 );

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
        saveState( rho, t_t, output );
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
        }
    }
    return true;
}