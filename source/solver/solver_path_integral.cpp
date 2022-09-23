#include "solver/solver_ode.h"

Sparse QDLC::Numerics::ODESolver::calculate_propagator_single( System &s, size_t tensor_dim, double t0, double t_step, int i, int j, std::vector<QDLC::SaveState> &output, const Sparse &one ) {
    Log::L3( "[PathIntegral] Calculating Single Propagator at t = {} for i = {}, j = {}\n", t0, i, j );
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

// TODO: maybe we can only evaluate the oupper triangle matrix, then write a getter function that checks if j>i, return mat[i,j].dagger(). Test with samples if (i,j) = (j,i).dagger()!!
std::vector<std::vector<Sparse>> &QDLC::Numerics::ODESolver::calculate_propagator_vector( System &s, size_t tensor_dim, double t0, double t_step, std::vector<QDLC::SaveState> &output ) {
    if ( pathint_propagator.contains( t0 ) ) {
        return pathint_propagator[t0];
    }
    Sparse one = Dense::Identity( tensor_dim, tensor_dim ).sparseView();
    std::vector<std::vector<Sparse>> ret( tensor_dim, { tensor_dim, Sparse( tensor_dim, tensor_dim ) } );
    // Calculate first by hand to ensure the Hamilton gets calculated correclty
    ret[0][0] = calculate_propagator_single( s, tensor_dim, t0, t_step, 0, 0, output, one ); // pathint_propagator[-1][0][0] );
// Calculate the remaining propagators
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_secondary_threads )
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) { // <= i
            if ( i == 0 && j == 0 ) continue;
            ret[i][j] = calculate_propagator_single( s, tensor_dim, t0, t_step, i, j, output, one ); // pathint_propagator[-1][i][j] );
        }
    }
    // Only save propagator vector if correlation functions are calculated.
    if ( not cache.empty() ) {
        Log::L3( "[PathIntegral] Caching propagator vector for t = {}\n", t0 );
#pragma omp critical
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
    std::map<std::string, std::vector<Sparse>> g12_settings_map;
    int g12_counter = 0;
    {
        // Order of Matrices is op1,op3,op2,op4 for g2, 1,1,op1,op2 for g1
        Sparse ident = Dense::Identity( rho0.rows(), rho0.cols() ).sparseView();
        auto &spectrum_s = s.parameters.input_correlation["Spectrum"];
        for ( int i = 0; i < spectrum_s.string_v["Modes"].size(); i++ ) {
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, spectrum_s.string_v["Modes"][i] );
            std::string g1 = get_operators_purpose( { s_creator, s_annihilator }, 1 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings_map.count( g1 ) == 0 )
                g12_settings_map[g1] = { ident, ident, creator, annihilator };
        }
        // Calculate Indist
        auto &indist_s = s.parameters.input_correlation["Indist"];
        for ( int i = 0; i < indist_s.string_v["Modes"].size(); i++ ) {
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, indist_s.string_v["Modes"][i] );
            std::string g1 = get_operators_purpose( { s_creator, s_annihilator }, 1 );
            std::string g2 = get_operators_purpose( { s_creator, s_annihilator, s_creator, s_annihilator }, 2 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings_map.count( g1 ) == 0 )
                g12_settings_map[g1] = { ident, ident, creator, annihilator };
            if ( g12_settings_map.count( g2 ) == 0 )
                g12_settings_map[g2] = { creator, annihilator, creator, annihilator };
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
            if ( g12_settings_map.count( g2_1111 ) == 0 )
                g12_settings_map[g2_1111] = { creator_1, annihilator_1, creator_1, annihilator_1 };
            if ( g12_settings_map.count( g2_2121 ) == 0 )
                g12_settings_map[g2_2121] = { creator_2, annihilator_2, creator_1, annihilator_1 };
            if ( g12_settings_map.count( g2_1221 ) == 0 )
                g12_settings_map[g2_1221] = { creator_1, annihilator_2, creator_2, annihilator_1 };
            if ( g12_settings_map.count( g2_2112 ) == 0 )
                g12_settings_map[g2_2112] = { creator_2, annihilator_1, creator_1, annihilator_2 };
            if ( g12_settings_map.count( g2_1122 ) == 0 )
                g12_settings_map[g2_1122] = { creator_1, annihilator_2, creator_1, annihilator_2 };
            if ( g12_settings_map.count( g2_2222 ) == 0 )
                g12_settings_map[g2_2222] = { creator_2, annihilator_2, creator_2, annihilator_2 };
        }
        // Calculate G1/G2 functions
        auto &gs_s = s.parameters.input_correlation["GFunc"];
        for ( int i = 0; i < gs_s.string_v["Modes"].size(); i++ ) {
            int order = std::abs( gs_s.numerical_v["Order"][i] );
            const auto &[s_creator, s_annihilator] = get_operator_strings( s, gs_s.string_v["Modes"][i] );
            std::string g = order == 1 ? get_operators_purpose( { s_creator, s_annihilator }, 1 ) : get_operators_purpose( { s_creator, s_annihilator, s_creator, s_annihilator }, 2 );
            auto [creator, annihilator] = get_operators_matrices( s, s_creator, s_annihilator );
            if ( g12_settings_map.count( g ) == 0 )
                if ( order == 1 )
                    g12_settings_map[g] = { ident, ident, creator, annihilator };
                else
                    g12_settings_map[g] = { creator, annihilator, creator, annihilator };
        }
        // int matdim = std::min( int( std::floor( ( t_end - t_start ) / s.parameters.t_step_pathint ) / s.parameters.iterations_t_skip ) + 1, s.parameters.grid_resolution ) + 1;
        int matdim = std::floor( ( t_end - t_start ) / s.parameters.t_step_pathint ) / s.parameters.iterations_t_skip + 1;
        // int matdim = std::min( int( std::floor( ( t_end - t_start ) / s.parameters.t_step ) / s.parameters.iterations_t_skip ) + 1, s.parameters.grid_resolution ) + 1;
        for ( auto &[purpose, matrices] : g12_settings_map ) {
            Log::L2( "[PathIntegral] Calculating G-Function with purpose {} in place with path integral.\n", purpose );
            cache[purpose] = Dense::Zero( matdim, matdim );
            cache[purpose + "_time"] = Dense::Zero( matdim, matdim );
        }
    }
    std::vector<std::pair<std::string, std::vector<Sparse>>> g12_settings( g12_settings_map.begin(), g12_settings_map.end() );

    Log::L2( "[PathIntegral] Setting up Path-Integral Solver...\n" );
    output.reserve( s.parameters.iterations_t_max + 1 );
    const int total_progressbar_iterations = std::floor( ( s.parameters.t_end - s.parameters.t_start ) / s.parameters.t_step_pathint * ( 1 + 0.5 * g12_settings.size() * ( s.parameters.t_end - s.parameters.t_start ) / s.parameters.t_step_pathint ) );

    // Gather unique Coupling indices, where the number of unique indices equals N(M)x for x >= 1
    std::set<int> set_different_dimensions;
    for ( int i = 0; i < s.operatorMatrices.phonon_hilbert_index_to_group_index.size(); i++ ) {
        set_different_dimensions.insert( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i] : i );
    }
    const int different_dimensions = set_different_dimensions.size();

    // Alias tensor Dimension for readability
    const Tensor::Index tensor_dim = rho0.rows();

    // Build Tensor Dimensions. Structure is N0,M0,N1,M1,...,NN,MM
    const Tensor::Index tensor_subdim = different_dimensions;
    Tensor::IndexVector pathint_tensor_dimensions = { tensor_dim, tensor_dim }; // TODO Remove shadow?
    for ( int i = 1; i < s.parameters.p_phonon_nc; i++ ) {
        pathint_tensor_dimensions.emplace_back( tensor_subdim );
        pathint_tensor_dimensions.emplace_back( tensor_subdim );
    }
    const auto index_vector_length = pathint_tensor_dimensions.size();

    // Output Tensor Dimensions
    std::stringstream adm_tensor_dimensions;
    std::ranges::copy( std::begin( pathint_tensor_dimensions ), std::end( pathint_tensor_dimensions ), std::ostream_iterator<int>( adm_tensor_dimensions, ", " ) );
    Log::L2( "[PathIntegral] ADM Tensor dimensions: [{}]\n", adm_tensor_dimensions.str() );
    Log::L2( "[PathIntegral] Using a maximum of #{} cpu cores for the ADM propagation.\n", s.parameters.numerics_maximum_secondary_threads );

    // Tensor Object
    Tensor adm_tensor( pathint_tensor_dimensions );
    Log::L2( "[PathIntegral] Fixed Tensor size will be {} MB.\n", adm_tensor.nonZeros(), adm_tensor.size() );

    // First step is just rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );

    // Fill initial Tensor
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( QDLC::Math::abs2( rho.coeff( i, j ) ) == 0 )
                continue;
            Tensor::IndexVector index( index_vector_length, 0 );
            index[0] = i;
            index[1] = j;
            adm_tensor( index ) += rho.coeff( i, j );
            Log::L2( "[PathIntegral] Added element ({},{}) = {} to adm tensor\n", i, j, rho.coeff( i, j ) );
        }
    }

    bool filled_time = false;
    int nonzero = 0;
    int new_nonzero = 0;

    std::map<std::string, std::vector<double>> profiler_time_per_thread;
    profiler_time_per_thread["PI_total"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_propagator"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_iteration"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_reduction"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_new_tensor"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );

    // Precalculate and cache all propagators if G-functions are evaluated to allow for easy multithreading
    if ( not g12_settings.empty() ) {
        Log::L2( "[PathIntegral] Precalculating all Propagator Vectors\n" );
        //#pragma omp parallel for schedule( static )
        for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
            std::cout << "Current: " << 100. * t_t / t_end << "\r";
            calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );
        }
    }

    // Iterate Path integral for further time steps
    for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Calculate Correlation functions:
        if ( g12_settings.size() > 0 and g12_counter % s.parameters.iterations_t_skip == 0 ) {
            // Allow nesting of parallel sections
            omp_set_nested( 1 );
            omp_set_max_active_levels( 3 );
            // Only parallize sections if all Hamiltons have been cached: NOT NEEDED
            int cores = t_t == t_start ? 1 : s.parameters.numerics_maximum_secondary_threads;
            // Parallel evaluation of G1/2 functions
#pragma omp parallel for schedule( static ) // num_threads( cores )
            for ( const auto &[purpose, matrices] : g12_settings ) {
                Log::L3( "[PathIntegral] Calculating sub-rk for {} with {} ({}) nested calls\n", purpose, omp_get_nested(), cores );
                std::vector<QDLC::SaveState> temp;
                auto &gmat = cache[purpose];
                auto &timemat = cache[purpose + "_time"];
                calculate_path_integral_correlation( adm_tensor, rho, t_t, t_end, t_step_initial, rkTimer, progressbar, total_progressbar_iterations, purpose, s, temp, do_output, matrices, s.parameters.numerics_maximum_secondary_threads, different_dimensions );
                Log::L3( "[PathIntegral] Writing {} values to G matrix...\n", temp.size() );
                for ( int32_t j = 0; j < std::min<int32_t>( temp.size(), gmat.rows() * s.parameters.iterations_t_skip ); j += s.parameters.iterations_t_skip ) {
                    double t_tau = temp.at( j ).t;
                    gmat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = s.dgl_expectationvalue<Sparse, Scalar>( temp.at( j ).mat, matrices[2] * matrices[1], t_tau );
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

        double profiler_total = omp_get_wtime();

        // Calculate Propagators for current time
        /* PROFILER */ double profiler_time = omp_get_wtime();
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );
        /* PROFILER */ profiler_time_per_thread["PI_propagator"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        new_nonzero = 0;

        // Main Iteration loop
        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );
        if ( max_index < s.parameters.p_phonon_nc )
            Log::L3( "[PathIntegral] max_index = {}, nc = {}\n", max_index, s.parameters.p_phonon_nc );
        // Iterate the tensor

        // Dense Tensor Iteration
        /* PROFILER */ profiler_time = omp_get_wtime();
        Tensor adm_tensor_next( adm_tensor.nonZeros() );
        /* PROFILER */ profiler_time_per_thread["PI_new_tensor"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;
        /* PROFILER */ profiler_time = omp_get_wtime();
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_secondary_threads ) schedule( static )
        for ( const Tensor::IndexVector &index : adm_tensor.get_indices() ) {
            Scalar new_value = 0.0;
            // Extract N0,M0 indices
            int i_n = index[0];
            int j_n = index[1];

            // Extract Group indices and extract N1,M1 indices:
            int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i_n] : i_n;
            int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[j_n] : j_n;
            int gi_n_m1 = index[2];
            int gj_n_m1 = index[3];

            //  Create "Old" sparse indices (vector shifted, last value will be replaced by summation) (n,n-1,n-2,...,n-m) -> (n-1,n-2,n-3,...,n-m,placeholder)
            auto index_old = index;
            // Set all indices to zero that are for timevalues larger than max_index
            for ( int i = max_index; i < s.parameters.p_phonon_nc; i++ ) {
                index_old[2 * i] = 0;
                index_old[2 * i + 1] = 0;
            }
            // Shift all indices once to the right (twice because both index vectors are chained)
            for ( int i = 0; i < index_old.size() - 2; i++ ) {
                index_old[i] = index_old[i + 2];
            }

            // Sum over all States in group (sum_(k_n-1,kd_n-1))
            for ( int i_n_m1 : s.operatorMatrices.phonon_group_index_to_hilbert_indices[gi_n_m1] ) {
                for ( int j_n_m1 : s.operatorMatrices.phonon_group_index_to_hilbert_indices[gj_n_m1] ) {
                    // Switch first entry of "old" sparse index back to the actual state index
                    index_old[0] = i_n_m1;
                    index_old[1] = j_n_m1;
                    // Propagator value
                    Scalar propagator_value;
                    // if ( j_n_m1 > i_n_m1 )
                    //     propagator_value = std::conj( propagator[j_n_m1][i_n_m1].coeff( i_n, j_n ) );
                    // else
                    propagator_value = propagator[i_n_m1][j_n_m1].coeff( i_n, j_n );
                    // If this propagation is not allowed, continue
                    if ( QDLC::Math::abs2( propagator_value ) == 0.0 )
                        continue;
                    // Sum over Groups lambda_n-n_m:
                    for ( int lambda_i = 0; lambda_i < different_dimensions; lambda_i++ ) {
                        for ( int lambda_j = 0; lambda_j < different_dimensions; lambda_j++ ) {
                            // Old index vector:
                            index_old[index_old.size() - 2] = lambda_i;
                            index_old[index_old.size() - 1] = lambda_j;
                            // Calculate S:
                            Scalar phonon_s = s.dgl_phonon_memory_function( 0, gi_n, gj_n, gi_n, gj_n );
                            for ( int tau = 0; tau < max_index; tau++ ) {
                                int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[index_old[0]] : index_old[0] ) : index_old[2 * tau] );
                                int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[index_old[1]] : index_old[1] ) : index_old[2 * tau + 1] );
                                phonon_s += s.dgl_phonon_memory_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                            }
                            new_value += propagator_value * adm_tensor( index_old ) * std::exp( phonon_s );
                        }
                    }
                }
            }
            if ( QDLC::Math::abs2( new_value ) != 0.0 ) {
                new_nonzero++;
            }
            // Because the indices are read by-reference, we need to set-zero the trailing indices for smaller times.
            if ( max_index == s.parameters.p_phonon_nc )
                adm_tensor_next( index ) = new_value;
            else {
                auto index_new = index;
                // Set all indices to zero that are for timevalues larger than max_index
                for ( int i = 2 * max_index; i < index_new.size(); i++ ) {
                    index_new[i] = 0;
                }
                adm_tensor_next( index_new ) = new_value;
            }
        }
        Log::L3( "[PathIntegral] Non-Zeros: {}\n", new_nonzero );
        nonzero = new_nonzero;
        /* PROFILER */ profiler_time_per_thread["PI_iteration"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        Log::L3( "[PathIntegral] Reducing ADM\n" );
        /* PROFILER */ profiler_time = omp_get_wtime();
        // Calculate the reduced density matrix by tracing over all past times for each entry
        Dense newrho = Dense::Zero( tensor_dim, tensor_dim );
        for ( const Tensor::IndexVector &index : adm_tensor.get_indices() ) {
            int i_n = index[0];
            int j_n = index[1];
            newrho( i_n, j_n ) += adm_tensor_next( index );
        }
        rho = newrho.sparseView(); // / newrho.trace();
        /* PROFILER */ profiler_time_per_thread["PI_reduction"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );

        adm_tensor = adm_tensor_next;

        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
        }

        // Output Profiler Information
        Log::L3( "[PathIntegral-Profiler] Time Profile for t = {}, #Nonzeros: {}\n", t_t, nonzero );
        /* PROFILER */ profiler_time_per_thread["PI_total"][omp_get_thread_num()] += omp_get_wtime() - profiler_total;
        double total_time = std::accumulate( profiler_time_per_thread["PI_total"].begin(), profiler_time_per_thread["PI_total"].end(), 0.0 );
        for ( auto &[name, vec] : profiler_time_per_thread ) {
            double cur_time = 0;
            std::ranges::for_each( vec, [&]( double &num ) { cur_time+=num; num = 0; } );
            Log::L3( "[PathIntegral-Profiler]     Name: {} - Time taken: {}s - {}\% of total time\n", name, cur_time, 100.0 * cur_time / total_time );
        }

        g12_counter++;
    }
    return true;
}

bool QDLC::Numerics::ODESolver::calculate_path_integral_correlation( Tensor adm_correlation, Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, size_t total_progressbar_iterations, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output, const std::vector<Sparse> &matrices, int adm_multithreading_cores, int different_dimensions ) {
    Log::L3( "- [PathIntegralCorrelation] Setting up Path-Integral Solver...\n" );
    const auto tensor_dim = rho0.rows();
    output.reserve( s.parameters.iterations_t_max + 1 );
    // saveState( rho0, t_start, output );
    bool modified = false;
    int nonzero = 0;

    // Profiler
    std::map<std::string, std::vector<double>> profiler_time_per_thread;
    profiler_time_per_thread["PI_total"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_propagator"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_iteration"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_reduction"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );
    profiler_time_per_thread["PI_new_tensor"] = std::vector<double>( s.parameters.numerics_maximum_secondary_threads, 0 );

    std::vector<std::vector<Sparse>> initial_propagator = calculate_propagator_vector( s, tensor_dim, t_start, s.parameters.numerics_subiterator_stepsize, output );
    for ( auto &vec : initial_propagator )
        for ( auto &el : vec )
            el = s.dgl_timetrafo( matrices[3] * el * matrices[0], t_start );

    // Main Iteration loop
    // Iterate Path integral for further time steps
    for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        double profiler_total = omp_get_wtime();

        // Calculate Propagators for current time
        /* PROFILER */ double profiler_time = omp_get_wtime();
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );
        /* PROFILER */ profiler_time_per_thread["PI_propagator"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        int new_nonzero = 0;

        // Main Iteration loop
        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );
        if ( max_index < s.parameters.p_phonon_nc )
            Log::L3( "[PathIntegral] max_index = {}, nc = {}\n", max_index, s.parameters.p_phonon_nc );
        // Initial Propagator

        // Dense Tensor Iteration
        /* PROFILER */ profiler_time = omp_get_wtime();
        Tensor adm_correlation_next( adm_correlation.nonZeros() );
        /* PROFILER */ profiler_time_per_thread["PI_new_tensor"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;
        /* PROFILER */ profiler_time = omp_get_wtime();
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_secondary_threads ) schedule( static )
        for ( const Tensor::IndexVector &index : adm_correlation.get_indices() ) {
            Scalar new_value = 0.0;
            // Extract N0,M0 indices
            int i_n = index[0];
            int j_n = index[1];

            // Extract Group indices and extract N1,M1 indices:
            int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i_n] : i_n;
            int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[j_n] : j_n;
            int gi_n_m1 = index[2];
            int gj_n_m1 = index[3];

            //  Create "Old" sparse indices (vector shifted, last value will be replaced by summation) (n,n-1,n-2,...,n-m) -> (n-1,n-2,n-3,...,n-m,placeholder)
            auto index_old = index;
            // Set all indices to zero that are for timevalues larger than max_index
            for ( int i = max_index; i < s.parameters.p_phonon_nc; i++ ) {
                index_old[2 * i] = 0;
                index_old[2 * i + 1] = 0;
            }
            // Shift all indices once to the right (twice because both index vectors are chained)
            for ( int i = 0; i < index_old.size() - 2; i++ ) {
                index_old[i] = index_old[i + 2];
            }

            // Sum over all States in group (sum_(k_n-1,kd_n-1))
            for ( int i_n_m1 : s.operatorMatrices.phonon_group_index_to_hilbert_indices[gi_n_m1] ) {
                for ( int j_n_m1 : s.operatorMatrices.phonon_group_index_to_hilbert_indices[gj_n_m1] ) {
                    // Switch first entry of "old" sparse index back to the actual state index
                    index_old[0] = i_n_m1;
                    index_old[1] = j_n_m1;
                    // Propagator value
                    Scalar propagator_value;
                    if ( not modified ) {
                        if ( j_n_m1 > i_n_m1 )
                            propagator_value = std::conj( initial_propagator[j_n_m1][i_n_m1].coeff( i_n, j_n ) );
                        else
                            propagator_value = initial_propagator[i_n_m1][j_n_m1].coeff( i_n, j_n );
                    } else {
                        if ( j_n_m1 > i_n_m1 )
                            propagator_value = std::conj( propagator[j_n_m1][i_n_m1].coeff( i_n, j_n ) );
                        else
                            propagator_value = propagator[i_n_m1][j_n_m1].coeff( i_n, j_n );
                    }
                    // If this propagation is not allowed, continue
                    if ( QDLC::Math::abs2( propagator_value ) == 0.0 )
                        continue;
                    // Sum over Groups lambda_n-n_m:
                    for ( int lambda_i = 0; lambda_i < different_dimensions; lambda_i++ ) {
                        for ( int lambda_j = 0; lambda_j < different_dimensions; lambda_j++ ) {
                            // Old index vector:
                            index_old[index_old.size() - 2] = lambda_i;
                            index_old[index_old.size() - 1] = lambda_j;
                            // Calculate S:
                            Scalar phonon_s = s.dgl_phonon_memory_function( 0, gi_n, gj_n, gi_n, gj_n );
                            for ( int tau = 0; tau < max_index; tau++ ) {
                                int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[index_old[0]] : index_old[0] ) : index_old[2 * tau] );
                                int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[index_old[1]] : index_old[1] ) : index_old[2 * tau + 1] );
                                phonon_s += s.dgl_phonon_memory_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                            }
                            new_value += propagator_value * adm_correlation( index_old ) * std::exp( phonon_s );
                        }
                    }
                }
            }
            if ( QDLC::Math::abs2( new_value ) != 0.0 ) {
                new_nonzero++;
            }
            // Because the indices are read by-reference, we need to set-zero the trailing indices for smaller times.
            if ( max_index == s.parameters.p_phonon_nc )
                adm_correlation_next( index ) = new_value;
            else {
                auto index_new = index;
                // Set all indices to zero that are for timevalues larger than max_index
                for ( int i = 2 * max_index; i < index_new.size(); i++ ) {
                    index_new[i] = 0;
                }
                adm_correlation_next( index_new ) = new_value;
            }
        }
        modified = true;
        Log::L3( "[PathIntegral] Non-Zeros: {}\n", new_nonzero );
        nonzero = new_nonzero;
        /* PROFILER */ profiler_time_per_thread["PI_iteration"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        Log::L3( "[PathIntegral] Reducing ADM\n" );
        /* PROFILER */ profiler_time = omp_get_wtime();
        // Calculate the reduced density matrix by tracing over all past times for each entry
        Dense newrho = Dense::Zero( tensor_dim, tensor_dim );
        for ( const Tensor::IndexVector &index : adm_correlation_next.get_indices() ) {
            int i_n = index[0];
            int j_n = index[1];
            newrho( i_n, j_n ) += adm_correlation( index );
        }
        Sparse rho = newrho.sparseView(); // / newrho.trace();
        /* PROFILER */ profiler_time_per_thread["PI_reduction"][omp_get_thread_num()] += omp_get_wtime() - profiler_time;

        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );

        adm_correlation = adm_correlation_next;

        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
        }
        // Output Profiler Information
        Log::L3( "[PathIntegralG12-Profiler] Time Profile for G1/2 t = {}, #Nonzeros: {}\n", t_t, nonzero );
        /* PROFILER */ profiler_time_per_thread["PI_total"][omp_get_thread_num()] += omp_get_wtime() - profiler_total;
        double total_time = std::accumulate( profiler_time_per_thread["PI_total"].begin(), profiler_time_per_thread["PI_total"].end(), 0.0 );
        for ( auto &[name, vec] : profiler_time_per_thread ) {
            double cur_time = 0;
            std::ranges::for_each( vec, [&]( double &num ) { cur_time+=num; num = 0; } );
            Log::L3( "[PathIntegralG12-Profiler]     Name: {} - Time taken: {}s - {}\% of total time\n", name, cur_time, 100.0 * cur_time / total_time );
        }
    }
    return true;
}