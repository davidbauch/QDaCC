#include "solver/solver_ode.h"

Scalar ODESolver::path_integral_recursive( const Sparse &rho0, System &s, std::vector<std::vector<std::vector<Sparse>>> &iterates, std::vector<SaveState> &output, FixedSizeSparseMap<Scalar> &adm, bool fillADM, int max_deph, int i_n, int j_n, const std::vector<int> &indicesX, const std::vector<int> &indicesY, Scalar adm_value, int current_deph ) {
    int tensor_dim = rho0.rows();
    Scalar result = 0;
    if ( current_deph == max_deph - 1 ) {
        for ( int k = 0; k < rho0.outerSize(); ++k ) {
            for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
                int i_n_m = it.row();
                int j_n_m = it.col();
                auto new_indicesX = indicesX;
                auto new_indicesY = indicesY;
                new_indicesX.emplace_back( i_n_m );
                new_indicesY.emplace_back( j_n_m );
                Scalar phonon_s = 0;
                for ( int lt = 0; lt < new_indicesX.size(); lt++ ) {
                    int ii_n = new_indicesX.at( lt );
                    int ij_n = new_indicesY.at( lt );
                    for ( int ld = 0; ld <= lt; ld++ ) {
                        int ii_nd = new_indicesX.at( ld );
                        int ij_nd = new_indicesY.at( ld );
                        if ( s.parameters.numerics_pathint_partially_summed )
                            phonon_s += s.dgl_phonon_S_function( lt - ld, s.operatorMatrices.phononCouplingFactor( ii_n, ii_n ), s.operatorMatrices.phononCouplingFactor( ij_n, ij_n ), s.operatorMatrices.phononCouplingFactor( ii_nd, ii_nd ), s.operatorMatrices.phononCouplingFactor( ij_nd, ij_nd ) );
                        else
                            phonon_s += s.dgl_phonon_S_function( lt - ld, ii_n, ij_n, ii_nd, ij_nd );
                    }
                }
                Scalar val = rho0.coeff( i_n_m, j_n_m ) * iterates.at( max_deph - 1 - current_deph ).at( i_n_m ).at( j_n_m ).coeff( i_n, j_n ) * std::exp( phonon_s );
                if ( fillADM && ( new_indicesX[0] == new_indicesY[0] || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
                    if ( s.parameters.numerics_pathint_partially_summed ) {
                        for ( int i = 1; i < new_indicesX.size(); i++ ) {
                            new_indicesX[i] = s.operatorMatrices.phononCouplingFactor( new_indicesX[i], new_indicesX[i] );
                            new_indicesY[i] = s.operatorMatrices.phononCouplingFactor( new_indicesY[i], new_indicesY[i] );
                        }
                    }

                    adm.addTriplet( Eigen::Map<Eigen::VectorXi>( new_indicesX.data(), new_indicesX.size() ), Eigen::Map<Eigen::VectorXi>( new_indicesY.data(), new_indicesY.size() ), adm_value * val, omp_get_thread_num() );
                }
                result += val;
            }
        }
    } else {
        for ( int i_n_m = 0; i_n_m < tensor_dim; i_n_m++ ) {
            for ( int j_n_m = 0; j_n_m < tensor_dim; j_n_m++ ) {
                Scalar val = iterates[max_deph - 1 - current_deph][i_n_m][j_n_m].coeff( i_n, j_n ); //TODO: in der backwards methode: über hier i_n_m iterieren statt über i_n, dann über nonzeros von i_n iterieren, vector umgekehrt aufbauen.
                if ( abs2( val ) == 0 ) continue;
                Scalar phonon_s = 0;
                auto new_indicesX = indicesX;
                auto new_indicesY = indicesY;
                new_indicesX.emplace_back( i_n_m );
                new_indicesY.emplace_back( j_n_m );
                result += val * path_integral_recursive( rho0, s, iterates, output, adm, fillADM, max_deph, i_n_m, j_n_m, new_indicesX, new_indicesY, adm_value * val, current_deph + 1 );
            }
        }
    }
    return result;
}

Scalar ODESolver::path_integral_recursive_backwards( const Sparse &rho0, System &s, std::vector<std::vector<std::vector<Sparse>>> &iterates, std::vector<SaveState> &output, FixedSizeSparseMap<Scalar> &adm, bool fillADM, int max_deph, int i_n, int j_n, const std::vector<int> &indicesX, const std::vector<int> &indicesY, Scalar adm_value, int current_deph ) {
    if ( current_deph == -1 )
        current_deph = max_deph - 1;
    int tensor_dim = rho0.rows();
    Scalar result = 0;
    if ( current_deph == max_deph - 1 ) {
        for ( int k = 0; k < rho0.outerSize(); ++k ) {
            for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
                int i_n_m = it.row();
                int j_n_m = it.col();
                auto new_indicesX = indicesX;
                auto new_indicesY = indicesY;
                new_indicesX.emplace_back( i_n_m ); // [i_0, , i_nc]
                new_indicesY.emplace_back( j_n_m ); // [j_0, , j_nc]
                Scalar val = rho0.coeff( i_n_m, j_n_m ) * iterates.at( max_deph - 1 - current_deph ).at( i_n_m ).at( j_n_m ).coeff( i_n, j_n );
                if ( current_deph == 0 ) {
                    result += val;
                } else {
                    result += val * path_integral_recursive_backwards( rho0, s, iterates, output, adm, fillADM, max_deph, i_n_m, j_n_m, new_indicesX, new_indicesY, adm_value * val, current_deph - 1 );
                }
            }
        }
    } else {
        for ( int i_n_m = 0; i_n_m < tensor_dim; i_n_m++ ) {
            for ( int j_n_m = 0; j_n_m < tensor_dim; j_n_m++ ) {
                Scalar phonon_s = 0;
                auto new_indicesX = indicesX;
                auto new_indicesY = indicesY;
                new_indicesX.insert( new_indicesX.begin() + 1, i_n_m ); // [i_0, i_n_m, , i_nc]
                new_indicesY.insert( new_indicesY.begin() + 1, j_n_m ); // [j_0, j_n_m, , j_nc]
                Scalar val = iterates.at( max_deph - 1 - current_deph ).at( i_n_m ).at( j_n_m ).coeff( i_n, j_n );
                if ( abs2( val ) == 0 ) continue;
                if ( current_deph == 0 ) {
                    Scalar phonon_s = 0;
                    for ( int lt = 0; lt < new_indicesX.size(); lt++ ) {
                        int ii_n = new_indicesX.at( lt );
                        int ij_n = new_indicesY.at( lt );
                        for ( int ld = 0; ld <= lt; ld++ ) {
                            int ii_nd = new_indicesX.at( ld );
                            int ij_nd = new_indicesY.at( ld );
                            phonon_s += s.dgl_phonon_S_function( lt - ld, ii_n, ij_n, ii_nd, ij_nd );
                        }
                    }
                    val *= std::exp( phonon_s );
                    result += val;
                    if ( fillADM && ( ( new_indicesX[0] == new_indicesY[0] && abs2( val ) != 0.0 ) || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
                        adm.addTriplet( Eigen::Map<Eigen::VectorXi>( new_indicesX.data(), new_indicesX.size() ), Eigen::Map<Eigen::VectorXi>( new_indicesY.data(), new_indicesY.size() ), adm_value * val, omp_get_thread_num() );
                    }
                } else {
                    result += val * path_integral_recursive_backwards( rho0, s, iterates, output, adm, fillADM, max_deph, i_n_m, j_n_m, new_indicesX, new_indicesY, adm_value * val, current_deph - 1 );
                }
            }
        }
    }
    return result;
}

Sparse ODESolver::path_integral( const Sparse &rho0, System &s, std::vector<std::vector<std::vector<Sparse>>> &iterates, std::vector<SaveState> &output, FixedSizeSparseMap<Scalar> &adm, bool fillADM, int max_deph ) {
    int tensor_dim = rho0.rows();
    Dense rho_ret = Dense::Zero( tensor_dim, tensor_dim );
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_phonons_maximum_threads )
    for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
        for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
            rho_ret( i_n, j_n ) = path_integral_recursive( rho0, s, iterates, output, adm, fillADM, max_deph, i_n, j_n, { i_n }, { j_n } );
        }
    }
    if ( fillADM ) {
        adm.reduceDublicates();
    }
    return rho_ret.sparseView();
}

Sparse ODESolver::calculate_propagator_single( System &s, size_t tensor_dim, double t0, double t_step, int i, int j, std::vector<SaveState> &output, const Sparse &one ) {
    Sparse projector = Sparse( tensor_dim, tensor_dim );
    projector.coeffRef( i, j ) = 1;
    //auto H = getHamilton( s, t0 );
    //Sparse M = one + t_step * s.dgl_rungeFunction( projector, H, t0, output );
    Sparse M = iterate( projector, s, t0, t_step, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );

    //Sparse M2 = iterate( projector, s, t0, 2.0 * t_step, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
    //Sparse M11 = iterate( M, s, t0 + t_step, t_step, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
    //if ( i == 1 && j == 1 )
    //    Log::L1( "2dt:\n{}\ndt*dt:\n{}\n", Dense( M2 ), Dense( M11 ) );

    Sparse map;
    if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
        map = project_matrix_sparse( M );
    }
    for ( double tau = t_step; tau < s.parameters.t_step_pathint; tau += t_step ) {
        //auto H = getHamilton( s, tau );
        //Sparse MM = ( one + t_step * s.dgl_rungeFunction( projector, H, tau, output ) ) * M;
        //M = MM;
        M = iterate( M, s, t0 + tau, t_step, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
        //Sparse MM = iterate( projector, s, tau, t_step, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
        //Sparse DM = MM.cwiseProduct( M );
        //M = DM;
    }
    //M = Dense( M ).exp().sparseView();
    if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
        return M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
    }
    //if ( i == 1 && j == 1 )
    //    Log::L1( "KEK:\n{}\nM:\n{}\n", Dense( KEK ), Dense( M ) );
    //return s.dgl_timetrafo( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ), t0 + t_step ); //MAYBE??
    M.makeCompressed();
    return M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
}

std::vector<std::vector<Sparse>> &ODESolver::calculate_propagator_vector( System &s, size_t tensor_dim, double t0, double t_step, std::vector<SaveState> &output ) {
    if ( pathint_propagator.count( t0 ) > 0 ) {
        return pathint_propagator[t0];
    }
    Sparse one = Dense::Identity( tensor_dim, tensor_dim ).sparseView();
    std::vector<std::vector<Sparse>> ret( tensor_dim, { tensor_dim, Sparse( tensor_dim, tensor_dim ) } );
    // Calculate first by hand to ensure the Hamilton gets calculated correclty
    ret[0][0] = calculate_propagator_single( s, tensor_dim, t0, t_step, 0, 0, output, one );
// Calculate the remaining propagators
#pragma omp parallel for num_threads( s.parameters.numerics_phonons_maximum_threads )
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( i == 0 && j == 0 ) continue;
            ret[i][j] = calculate_propagator_single( s, tensor_dim, t0, t_step, i, j, output, one );
        }
    }
    // Only save propagator vector if correlation functions are calculated.
    if ( cache.size() > 0 ) {
        Log::L3( "[PATHINT] Caching propagator vector for t = {}\n", t0 );
        pathint_propagator[t0] = ret;
        return pathint_propagator[t0];
    } else {
        pathint_propagator[-1] = ret;
    }
    return pathint_propagator[-1];
}

// TODO: cache all propagators in giant map. recalculation for correlation functions scales with number of functions otherwise.
bool ODESolver::calculate_path_integral( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output ) {
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
            const auto mode = splitline( modes, '-' );
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
            Log::L3( "[PATHINT] Calculating G-Function with purpose {} in place with path integral.\n", purpose );
            cache[purpose] = Dense::Zero( matdim, matdim );
        }
        if ( cache.size() > 0 )
            cache["Time"] = Dense::Zero( matdim, matdim );
    }

    Log::L3( "[PATHINT] Setting up Path-Integral Solver...\n" );
    output.reserve( s.parameters.iterations_t_max + 1 );
    int tensor_dim = rho0.rows();

    std::set<int> different_dimensions;
    for ( int i = 0; i < s.operatorMatrices.phononCouplingFactor.rows(); i++ ) {
        different_dimensions.insert( s.operatorMatrices.phononCouplingFactor( i, i ) );
    }
    pathint_tensor_dimensions = { tensor_dim };
    for ( int i = 1; i < s.parameters.p_phonon_nc; i++ )
        pathint_tensor_dimensions.emplace_back( different_dimensions.size() );

    FixedSizeSparseMap<Scalar> adms = FixedSizeSparseMap<Scalar>( { tensor_dim }, s.parameters.numerics_phonons_maximum_threads );

    // Calculate first n_c timesteps directly
    //Log::L3( "Calculating first n_c = {} timesteps of the pathintegral directly...\n", s.parameters.p_phonon_nc );

    // First step is just rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );
    // Next s.parameters.p_phonon_nc Steps and ADM initialization

    // Precalculate all iterates per index and timestep once. Saveorder is j+j_max*i
    //Log::L3( "Pre-Calculating propagator matrices...\n" );
    //auto t_prop = omp_get_wtime();
    //std::vector<std::vector<std::vector<Sparse>>> propagators;
    //// Calculate all propagators
    //for ( int t = 0; t < s.parameters.p_phonon_nc - 1; t++ ) {
    //    propagators.emplace_back( calculate_propagator_vector( s, tensor_dim, s.parameters.t_step_pathint * t, t_step_initial, output ) );
    //}
    //Log::L3( "Done! Took {}s\n", omp_get_wtime() - t_prop );

    // Calculate the first n_c steps by recursion:
    //for ( int n = 1; n < s.parameters.p_phonon_nc; n++ ) {
    //    Log::L3( "Calculating rho(t{})...\n", n );
    //    rho = path_integral( rho0, s, propagators, output, adms, n == s.parameters.p_phonon_nc - 1, n );
    //    saveState( rho, s.parameters.t_step_pathint * n, output );
    //    rkTimer.iterate();
    //    if ( do_output ) {
    //        Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
    //    }
    //}

    //Log::L3( "Size of ADMs Tensor: {} bytes / {}\% filled\n", adms.getSizeOfValues(), adms.getFillRatio() * 100.0 );
    //Log::L3( "Calculating Path-Integral Loop for the remaining {} timesteps...\n", std::floor( ( t_end - t_start ) / s.parameters.t_step_pathint ) - s.parameters.p_phonon_nc );

    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( abs2( rho.coeff( i, j ) ) == 0 )
                continue;
            Eigen::VectorXi ii = Eigen::VectorXi( 1 );
            Eigen::VectorXi jj = Eigen::VectorXi( 1 );
            ii( 0 ) = i;
            jj( 0 ) = j;
            adms.addTriplet( ii, jj, rho.coeff( i, j ), omp_get_thread_num() );
        }
    }
    adms.reduceDublicates();
    Log::L3( "[PATHINT] Filled Initial Tensor with elements: [" );
    for ( auto &[sparse_index_x, outer_map] : adms.get() ) {
        for ( auto &[sparse_index_y, value] : outer_map ) {
            Log::L3( "{} ({},{}), ", value, sparse_index_x, sparse_index_y );
        }
    }
    Log::L3( "] ({} elements)\n", adms.nonZeros() );

    // Iterate Path integral for further time steps
    double t_start_new = t_start; // + s.parameters.t_step_pathint * ( s.parameters.p_phonon_nc - 1 );
    for ( double t_t = t_start_new; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Calculate Correlation functions:
        if ( g12_settings.size() > 0 and g12_counter % s.parameters.iterations_t_skip == 0 ) {
            //auto cached_adms_dimensions = adms.dimensions;
            bool filled_time = false;
            auto &timemat = cache["Time"];
            for ( auto &[purpose, matrices] : g12_settings ) {
                Log::L3( "[PATHINT] Calculating sub-rk for {}\n", purpose );
                //adms.dimensions = cached_adms_dimensions;
                //adms.rescale_dimensions();
                int order = matrices.size() == 4 ? 2 : 1;
                std::vector<SaveState> temp;
                auto &gmat = cache[purpose];
                calculate_path_integral_correlation( adms, rho, t_t, t_end, t_step_initial, rkTimer, progressbar, purpose, s, temp, do_output, matrices );
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
                    if ( not filled_time ) {
                        timemat( g12_counter / s.parameters.iterations_t_skip, j / s.parameters.iterations_t_skip ) = Scalar( t_t, t_tau );
                    }
                }
                filled_time = true;
            }
            //adms.dimensions = cached_adms_dimensions;
            //adms.rescale_dimensions();
        }
        // Path-Integral iteration

        // Calculate Propagators for current time
        auto t0 = omp_get_wtime();
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, t_step_initial, output );
        t0 = ( omp_get_wtime() - t0 );

        double cur_min = 1;

        // Main Iteration loop
        Log::L3( "[PATHINT] Current Tensor Size: {} elements at a total of {} MB\n", adms.nonZeros(), adms.size() / 1024. / 1024. );
        auto t1 = omp_get_wtime();
        double total_append_time = 0;
        double total_time = 0;

        bool addeddimension = false;
        if ( adms.dimensions_scaled.size() < s.parameters.p_phonon_nc ) {
            auto toadd = pathint_tensor_dimensions.at( t_t / s.parameters.t_step_pathint + 1 );
            Log::L3( "[PATHINT] Added tensor dimension {} to the path integral, total dimensions: {}\n", toadd, adms.dimensions_scaled.size() );
            adms.add_dimension( toadd );
            addeddimension = true;
        }

        ////for ( auto &threadvector : adms.get() )
        //auto tensor = adms.getIndices();
        //// Find starting and endpoint for threads:
        //std::vector<std::pair<size_t, size_t>> indices;
        //int index_size = 0; //tensor.front().first.size() - 1;
        //size_t current = tensor.front().first( index_size );
        //for ( size_t i = 1; i < tensor.size(); i++ ) {
        //    if ( tensor[i].first( index_size ) > tensor[i - 1].first( index_size ) ) {
        //        indices.emplace_back( std::make_pair( current, i ) );
        //        current = i;
        //    }
        //}
        //if ( current < tensor.size() )
        //    indices.emplace_back( std::make_pair( current, tensor.size() ) );

        //for ( auto &[a, b] : indices ) {
        //    Log::L3( "[{} {}], ", a, b );
        //}
        //    Log::L3( "\n" );

        //#pragma omp parallel for num_threads( s.parameters.numerics_phonons_maximum_threads )
        //for ( auto &[a, b] : indices )
        //for ( int index = a; index < b; index++ ) {
        //auto &[sparse_index_x, outer_map] = tensor[index];
        for ( auto &[sparse_index_x, outer_map] : adms.get() ) {
            //#pragma omp parallel for num_threads( s.parameters.numerics_phonons_maximum_threads )
            //            for ( auto i = 0; i < adms.get().size(); i++ ) {
            //                auto iterator = adms.get().begin();
            //                std::advance( iterator, i );
            //                auto &sparse_index_x = iterator->first;
            //                auto &outer_map = iterator->second;
            for ( auto &[sparse_index_y, value] : outer_map ) {
                if ( abs2( value ) == 0 ) continue;
                auto tt = omp_get_wtime();
                //Log::L3( "[PATHINT] handling ({}),({}) -> {}\n", sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), value );
                for ( int l = 0; l < propagator[sparse_index_x( 0 )][sparse_index_y( 0 )].outerSize(); ++l ) {
                    for ( Sparse::InnerIterator M( propagator[sparse_index_x( 0 )][sparse_index_y( 0 )], l ); M; ++M ) {
                        int i_n = M.row();
                        int j_n = M.col();
                        int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( i_n, i_n ) : i_n; //TODO: nicht couplingFactor, sondern couplingIndex = set(couplingFactor), damit auch non-int couplings gehen!
                        int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( j_n, j_n ) : j_n;

                        Scalar phonon_s = s.dgl_phonon_S_function( 0, gi_n, gj_n, gi_n, gj_n );
                        for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                            int gi_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( sparse_index_x( 0 ), sparse_index_x( 0 ) ) : sparse_index_x( 0 ) ) : sparse_index_x( tau ) );
                            int gj_nd = ( tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( sparse_index_y( 0 ), sparse_index_y( 0 ) ) : sparse_index_y( 0 ) ) : sparse_index_y( tau ) );
                            phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                        }
                        Scalar val = M.value() * value * std::exp( phonon_s );
                        double abs = abs2( val );
                        auto appendtime = omp_get_wtime();
                        // Add Element to Triplet list if they are diagonal elements or if they surpass the given threshold.
                        if ( i_n == j_n || abs >= s.parameters.numerics_pathintegral_squared_threshold ) {
                            // Add new indices to vectors:
                            if ( addeddimension ) {
                                Eigen::VectorXi new_sparse_index_x = Eigen::VectorXi::Zero( sparse_index_x.size() + 1 );
                                Eigen::VectorXi new_sparse_index_y = Eigen::VectorXi::Zero( sparse_index_y.size() + 1 );
                                for ( int i = 0; i < sparse_index_x.size(); i++ ) {
                                    new_sparse_index_x( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( sparse_index_x( i ), sparse_index_x( i ) ) : sparse_index_x( i ) );
                                    new_sparse_index_y( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( sparse_index_y( i ), sparse_index_y( i ) ) : sparse_index_y( i ) );
                                }
                                new_sparse_index_x( 0 ) = i_n;
                                new_sparse_index_y( 0 ) = j_n;
                                adms.addTriplet( new_sparse_index_x, new_sparse_index_y, val, omp_get_thread_num() );
                            } else {
                                cur_min = cur_min != 0.0 && cur_min < abs ? cur_min : abs;
                                if ( s.parameters.numerics_pathint_partially_summed )
                                    adms.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n, s.operatorMatrices.phononCouplingFactor( sparse_index_x( 0 ), sparse_index_x( 0 ) ), s.operatorMatrices.phononCouplingFactor( sparse_index_y( 0 ), sparse_index_y( 0 ) ) );
                                else
                                    adms.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n );
                            }
                        }
                        total_append_time += omp_get_wtime() - appendtime;
                    }
                }
                total_time += omp_get_wtime() - tt;
            }
        }

        t1 = ( omp_get_wtime() - t1 );

        auto ts = omp_get_wtime();
        adms.reduceDublicates();
        ts = omp_get_wtime() - ts;

        Log::L3( "[PATHINT] Reducing ADM\n" );
        // Calculate the reduced density matrix by tracing over all past times for each entry
        auto t2 = omp_get_wtime();
        Dense newrho = Dense::Zero( tensor_dim, tensor_dim );
        for ( auto &[sparse_index_x, outer_map] : adms.get() ) {
            for ( auto &[sparse_index_y, value] : outer_map ) {
                int i_n = sparse_index_x( 0 );
                int j_n = sparse_index_y( 0 );
                newrho( i_n, j_n ) += value;
            }
        }
        //rho = newrho.sparseView() / newrho.trace();
        rho = newrho.sparseView() / newrho.trace();
        t2 = ( omp_get_wtime() - t2 );
        Log::L3( "[PATHINT] Iteration: {}, time taken: [ Propagator: {:.4f}s, ADM Advancing: {:.4f}s (Partial append time: {:.4f}\%), ADM Setting: {:.4f}s, ADM Reduction: {:.4f}s ], Trace: {}, Elements: {}\n", t_t, t0, t1, 100.0 * total_append_time / total_time, ts, t2, s.getTrace<Scalar>( rho ), adms.nonZeros() );

        // Dynamic Cutoff
        if ( s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max > 0 ) {
            double ratio = (double)adms.nonZeros() / s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max;
            ratio = std::exp( ratio ) * std::pow( ratio, 5 );              // Scaled squared
            s.parameters.numerics_pathintegral_squared_threshold *= ratio; // Adjust Cutoff by squared ratio
            if ( s.parameters.numerics_pathintegral_squared_threshold < 1E-30 ) s.parameters.numerics_pathintegral_squared_threshold = 1E-30;
            if ( s.parameters.numerics_pathintegral_squared_threshold > 1E4 ) s.parameters.numerics_pathintegral_squared_threshold = 1E4;
            Log::L3( "[PATHINT] Adjusted Cutoff Threshold to {} (x{})\n", std::sqrt( s.parameters.numerics_pathintegral_squared_threshold ), ratio );
        }

        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );
        // Progress and time output
        rkTimer.iterate();
        g12_counter++;
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
        }
    }
    //Log::L3( "Final Size of ADMs Tensor: {} bytes / {}\% filled\n", adms.getSizeOfValues(), adms.getFillRatio() * 100.0 );
    Log::L3( "Done!\n" );
    return true;
}

bool ODESolver::calculate_path_integral_correlation( FixedSizeSparseMap<Scalar> adm_correlation, Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output, const std::vector<Sparse> &matrices ) {
    int order = matrices.size() == 4 ? 2 : 1;
    Log::L3( "| [PATHINT SUB] Setting up Path-Integral Solver...\n" );
    int tensor_dim = rho0.rows();
    output.reserve( s.parameters.iterations_t_max + 1 );
    //saveState( rho0, t_start, output );
    bool modified = false;

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

        bool addeddimension = false;
        if ( adm_correlation.dimensions_scaled.size() < s.parameters.p_phonon_nc ) {
            auto toadd = pathint_tensor_dimensions.at( t_t / s.parameters.t_step_pathint + 1 );
            Log::L3( "| [PATHINT SUB] Added tensor dimension {} to the path integral, total dimensions: {}\n", toadd, adm_correlation.dimensions_scaled.size() );
            adm_correlation.add_dimension( toadd );
            addeddimension = true;
        }

        // Main Iteration loop
        Log::L3( "| [PATHINT SUB] Current Correlation Tensor Size: {} elements at a total of {} MB\r", adm_correlation.nonZeros(), adm_correlation.size() / 1024. / 1024. );
        auto t1 = omp_get_wtime();
        double total_append_time = 0;
        double total_time = 0;

        for ( auto &[sparse_index_x, outer_map] : adm_correlation.get() ) {
            for ( auto &[sparse_index_y, value] : outer_map ) {
                if ( abs2( value ) == 0 ) continue;
                auto tt = omp_get_wtime();
                for ( int l = 0; l < propagator[sparse_index_x( 0 )][sparse_index_y( 0 )].outerSize(); ++l ) {
                    for ( Sparse::InnerIterator M( propagator[sparse_index_x( 0 )][sparse_index_y( 0 )], l ); M; ++M ) {
                        int i_n = M.row();
                        int j_n = M.col();
                        int gi_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( i_n, i_n ) : i_n;
                        int gj_n = s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( j_n, j_n ) : j_n;
                        Scalar phonon_s = s.dgl_phonon_S_function( 0, gi_n, gj_n, gi_n, gj_n );
                        for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                            int gi_nd = tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( sparse_index_x( 0 ), sparse_index_x( 0 ) ) : sparse_index_x( 0 ) ) : sparse_index_x( tau );
                            int gj_nd = tau == 0 ? ( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( sparse_index_y( 0 ), sparse_index_y( 0 ) ) : sparse_index_y( 0 ) ) : sparse_index_y( tau );
                            phonon_s += s.dgl_phonon_S_function( tau + 1, gi_n, gj_n, gi_nd, gj_nd );
                        }
                        Scalar val = M.value() * value * std::exp( phonon_s );
                        double abs = abs2( val );
                        auto appendtime = omp_get_wtime();
                        // Add Element to Triplet list if they are diagonal elements or if they surpass the given threshold.
                        if ( i_n == j_n || abs >= s.parameters.numerics_pathintegral_squared_threshold ) {
                            // Add new indices to vectors:
                            if ( addeddimension ) {
                                Eigen::VectorXi new_sparse_index_x = Eigen::VectorXi::Zero( sparse_index_x.size() + 1 );
                                Eigen::VectorXi new_sparse_index_y = Eigen::VectorXi::Zero( sparse_index_y.size() + 1 );
                                for ( int i = 0; i < sparse_index_x.size(); i++ ) {
                                    new_sparse_index_x( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( sparse_index_x( i ), sparse_index_x( i ) ) : sparse_index_x( i ) );
                                    new_sparse_index_y( i + 1 ) = ( i == 0 and s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phononCouplingFactor( sparse_index_y( i ), sparse_index_y( i ) ) : sparse_index_y( i ) );
                                }
                                new_sparse_index_x( 0 ) = i_n;
                                new_sparse_index_y( 0 ) = j_n;
                                //Log::L3( "[PATHINT SUB] saving ({}),({}) -> {}\n", new_sparse_index_x.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), new_sparse_index_y.format( Eigen::IOFormat( 0, 0, ", ", " ", "", "" ) ), val );
                                adm_correlation.addTriplet( new_sparse_index_x, new_sparse_index_y, val, omp_get_thread_num() );
                            } else {
                                cur_min = cur_min != 0.0 && cur_min < abs ? cur_min : abs;
                                if ( s.parameters.numerics_pathint_partially_summed )
                                    adm_correlation.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n, s.operatorMatrices.phononCouplingFactor( sparse_index_x( 0 ), sparse_index_x( 0 ) ), s.operatorMatrices.phononCouplingFactor( sparse_index_y( 0 ), sparse_index_y( 0 ) ) );
                                else
                                    adm_correlation.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n );
                            }
                        }
                        total_append_time += omp_get_wtime() - appendtime;
                        //}
                    }
                }
                total_time += omp_get_wtime() - tt;
            }
        }
        //modified = true;

        t1 = ( omp_get_wtime() - t1 );

        auto ts = omp_get_wtime();
        adm_correlation.reduceDublicates();
        ts = omp_get_wtime() - ts;

        Log::L3( "| [PATHINT SUB] Reducing Correlation ADM\n" );
        // Calculate the reduced density matrix by tracing over all past times for each entry
        auto t2 = omp_get_wtime();
        Dense cache = Dense::Zero( tensor_dim, tensor_dim );
        for ( auto &[sparse_index_x, outer_map] : adm_correlation.get() ) {
            for ( auto &[sparse_index_y, value] : outer_map ) {
                int i_n = sparse_index_x( 0 );
                int j_n = sparse_index_y( 0 );
                cache( i_n, j_n ) += value;
            }
        }
        auto rho = cache.sparseView();
        t2 = ( omp_get_wtime() - t2 );
        Log::L3( "| [PATHINT SUB] Correlation Iteration: {}, time taken: [ Propagator: {:.4f}s, ADM Advancing: {:.4f}s (Partial append time: {:.4f}\%), ADM Setting: {:.4f}s, ADM Reduction: {:.4f}s ]\n", t_t, t0, t1, 100.0 * total_append_time / total_time, ts, t2 );

        // Dynamic Cutoff
        if ( s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max > 0 ) {
            double ratio = (double)adm_correlation.nonZeros() / s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max;
            ratio = std::exp( ratio ) * std::pow( ratio, 5 );              // Scaled squared
            s.parameters.numerics_pathintegral_squared_threshold *= ratio; // Adjust Cutoff by squared ratio
            if ( s.parameters.numerics_pathintegral_squared_threshold < 1E-30 ) s.parameters.numerics_pathintegral_squared_threshold = 1E-30;
            if ( s.parameters.numerics_pathintegral_squared_threshold > 1E4 ) s.parameters.numerics_pathintegral_squared_threshold = 1E4;
            Log::L3( "| [PATHINT SUB] Adjusted Cutoff Threshold to {} (x{})\n", std::sqrt( s.parameters.numerics_pathintegral_squared_threshold ), ratio );
        }
        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
        }
    }
    Log::L3( "| [PATHINT SUB] Correlation Done!\n" );
    return true;
}

//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_phonons_maximum_threads )
//for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
//    for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
//        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
//            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
//                if ( std::abs( iterates.at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) ) == 0.0 ) continue;
//                for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
//                    for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
//                        for ( int i_n_m3 = 0; i_n_m3 < tensor_dim; i_n_m3++ ) {
//                            for ( int j_n_m3 = 0; j_n_m3 < tensor_dim; j_n_m3++ ) {
//                                Scalar sum = 0;
//                                for ( int i_n_m4 = 0; i_n_m4 < tensor_dim; i_n_m4++ ) {
//                                    for ( int j_n_m4 = 0; j_n_m4 < tensor_dim; j_n_m4++ ) {
//                                        sum += adms.get( i_n_m1, i_n_m2, i_n_m3, i_n_m4, j_n_m1, j_n_m2, j_n_m3, j_n_m4 );
//                                    }
//                                }
//                                sum *= iterates.at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n );
//                                if ( std::abs( sum ) > 0 ) { //1E-15 ) {
//                                    #pragma omp critical
//                                    adms.addTriplet( i_n, i_n_m1, i_n_m2, i_n_m3, j_n, j_n_m1, j_n_m2, j_n_m3, sum );
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//}

//for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
//    for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
//        Scalar sum = 0;
//        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
//            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
//                for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
//                    for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
//                        for ( int i_n_m3 = 0; i_n_m3 < tensor_dim; i_n_m3++ ) {
//                            for ( int j_n_m3 = 0; j_n_m3 < tensor_dim; j_n_m3++ ) {
//                                sum += adms.get( i_n, i_n_m1, i_n_m2, i_n_m3, j_n, j_n_m1, j_n_m2, j_n_m3 );
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        triplets.emplace_back( i_n, j_n, sum );
//    }
//}
//rho.setFromTriplets( triplets.begin(), triplets.end() );

// Temp, brute force function
//Sparse ODESolver::path_integral_bruteforce( const Sparse &rho0, System &s, std::vector<std::vector<Sparse>> &iterates, std::vector<SaveState> &output, FixedSizeSparseMap<Scalar> &adm, bool fillADM, int deph ) {
//    int tensor_dim = rho0.rows();
//    Dense rho_ret = Dense::Zero( tensor_dim, tensor_dim );
//
//    if ( deph == 1 ) {
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_phonons_maximum_threads )
//        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
//            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
//                Scalar result = 0;
//                for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
//                    for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
//                        Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 );
//                        phonon_s = 0;
//                        result += rho0.coeff( i_n_m1, j_n_m1 ) * 1.0 * iterates.at( 0 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * std::exp( phonon_s );
//                    }
//                }
//                rho_ret( i_n, j_n ) = result;
//            }
//        }
//    } else if ( deph == 2 ) {
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_phonons_maximum_threads )
//        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
//            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
//                Scalar result = 0;
//                for ( int k = 0; k < rho0.outerSize(); ++k ) {
//                    for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
//                        auto i_n_m2 = (int)it.row();
//                        auto j_n_m2 = (int)it.col();
//                        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
//                            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
//                                if ( abs2( iterates.at( 0 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
//                                Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 0, i_n_m2, j_n_m2, i_n_m2, j_n_m2 );
//                                phonon_s += s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n_m1, j_n_m1, i_n_m2, j_n_m2 );
//                                phonon_s += s.dgl_phonon_S_function( 2, i_n, j_n, i_n_m2, j_n_m2 );
//                                phonon_s = 0;
//                                Scalar val = it.value() * 1.0 * iterates.at( 1 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates.at( 0 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * std::exp( phonon_s );
//                                result += val;
//                                if ( fillADM && ( i_n == j_n || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
//                                    //#pragma omp critical
//                                    adm.addTriplet( { i_n, i_n_m1, i_n_m2 }, { j_n, j_n_m1, j_n_m2 }, val, omp_get_thread_num() );
//                                }
//                            }
//                        }
//                    }
//                }
//                rho_ret( i_n, j_n ) = result;
//            }
//        }
//    } else if ( deph == 3 ) {
//        int a = 0;
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_phonons_maximum_threads )
//        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
//            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
//                //fmt::print( "N_c = 3 -- i = {}, j = {} -> {:3.0f}\%\n", i_n, j_n, 100.0 * a / 36 / 36 );
//                //#pragma omp critical
//                a++;
//                Scalar result = 0;
//                for ( int k = 0; k < rho0.outerSize(); ++k ) {
//                    for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
//                        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
//                            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
//                                if ( abs2( iterates.at( 2 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) ) == 0.0 ) continue;
//                                for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
//                                    for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
//                                        if ( abs2( iterates.at( 1 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
//                                        auto i_n_m3 = (int)it.row();
//                                        auto j_n_m3 = (int)it.col();
//                                        Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 0, i_n_m2, j_n_m2, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 0, i_n_m3, j_n_m3, i_n_m3, j_n_m3 );
//                                        phonon_s += s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n_m1, j_n_m1, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 1, i_n_m2, j_n_m2, i_n_m3, j_n_m3 );
//                                        phonon_s += s.dgl_phonon_S_function( 2, i_n, j_n, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 2, i_n_m1, j_n_m1, i_n_m3, j_n_m3 );
//                                        phonon_s += s.dgl_phonon_S_function( 3, i_n, j_n, i_n_m3, j_n_m3 );
//                                        Scalar val = it.value() * 1.0 * iterates.at( 2 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates.at( 1 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * iterates.at( 0 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) * std::exp( phonon_s );
//                                        phonon_s = 0;
//                                        result += val;
//                                        // Indexing is time-upwards, meaning t_n, t_n-1, t_n-2, t_n-3
//                                        if ( fillADM && ( i_n == j_n || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
//                                            //#pragma omp critical
//                                            adm.addTriplet( { i_n, i_n_m1, i_n_m2, i_n_m3 }, { j_n, j_n_m1, j_n_m2, j_n_m3 }, val, omp_get_thread_num() );
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                rho_ret( i_n, j_n ) = result;
//            }
//        }
//    } else if ( deph == 4 ) {
//        int a = 0;
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_phonons_maximum_threads )
//        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
//            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
//                //fmt::print( "N_c = 4 -- i = {}, j = {} -> {:3.0f}\%\n", i_n, j_n, 100.0 * a / 36 / 36 );
//                //#pragma omp critical
//                a++;
//                Scalar result = 0;
//                for ( int k = 0; k < rho0.outerSize(); ++k ) {
//                    for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
//                        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
//                            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
//                                if ( abs2( iterates.at( 3 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) ) == 0.0 ) continue;
//                                for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
//                                    for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
//                                        if ( abs2( iterates.at( 2 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
//                                        for ( int i_n_m3 = 0; i_n_m3 < tensor_dim; i_n_m3++ ) {
//                                            for ( int j_n_m3 = 0; j_n_m3 < tensor_dim; j_n_m3++ ) {
//                                                if ( abs2( iterates.at( 1 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) ) == 0.0 ) continue;
//                                                auto i_n_m4 = (int)it.row();
//                                                auto j_n_m4 = (int)it.col();
//                                                Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 0, i_n_m2, j_n_m2, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 0, i_n_m3, j_n_m3, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 0, i_n_m4, j_n_m4, i_n_m4, j_n_m4 );
//                                                phonon_s += s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n_m1, j_n_m1, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 1, i_n_m2, j_n_m2, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 1, i_n_m3, j_n_m3, i_n_m4, j_n_m4 );
//                                                phonon_s += s.dgl_phonon_S_function( 2, i_n, j_n, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 2, i_n_m1, j_n_m1, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 2, i_n_m2, j_n_m2, i_n_m4, j_n_m4 );
//                                                phonon_s += s.dgl_phonon_S_function( 3, i_n, j_n, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 3, i_n_m1, j_n_m1, i_n_m4, j_n_m4 );
//                                                phonon_s += s.dgl_phonon_S_function( 4, i_n, j_n, i_n_m4, j_n_m4 );
//                                                Scalar val = it.value() * 1.0 * iterates.at( 3 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates.at( 2 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * iterates.at( 1 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) * iterates.at( 0 ).at( j_n_m4 + tensor_dim * i_n_m4 ).coeff( i_n_m3, j_n_m3 ) * std::exp( phonon_s );
//                                                phonon_s = 0;
//                                                result += val;
//                                                // Indexing is time-upwards, meaning t_n, t_n-1, t_n-2, t_n-3
//                                                if ( fillADM && ( i_n == j_n || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
//                                                    //#pragma omp critical
//                                                    adm.addTriplet( { i_n, i_n_m1, i_n_m2, i_n_m3, i_n_m4 }, { j_n, j_n_m1, j_n_m2, j_n_m3, j_n_m4 }, val, omp_get_thread_num() );
//                                                }
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                rho_ret( i_n, j_n ) = result;
//            }
//        }
//    } else if ( deph == 5 ) {
//        int a = 0;
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_phonons_maximum_threads )
//        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
//            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
//                //fmt::print( "N_c = 5 --  i = {}, j = {} -> {:3.0f}\%\n", i_n, j_n, 100.0 * a / 36 / 36 );
//                a++;
//                Scalar result = 0;
//                for ( int k = 0; k < rho0.outerSize(); ++k ) {
//                    for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
//                        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
//                            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
//                                if ( abs2( iterates.at( 4 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) ) == 0.0 ) continue;
//                                for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
//                                    for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
//                                        if ( abs2( iterates.at( 3 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
//                                        for ( int i_n_m3 = 0; i_n_m3 < tensor_dim; i_n_m3++ ) {
//                                            for ( int j_n_m3 = 0; j_n_m3 < tensor_dim; j_n_m3++ ) {
//                                                if ( abs2( iterates.at( 2 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) ) == 0.0 ) continue;
//                                                for ( int i_n_m4 = 0; i_n_m4 < tensor_dim; i_n_m4++ ) {
//                                                    for ( int j_n_m4 = 0; j_n_m4 < tensor_dim; j_n_m4++ ) {
//                                                        if ( std::abs( iterates.at( 1 ).at( j_n_m4 + tensor_dim * i_n_m4 ).coeff( i_n_m3, j_n_m3 ) ) == 0.0 ) continue;
//                                                        auto i_n_m5 = (int)it.row();
//                                                        auto j_n_m5 = (int)it.col();
//                                                        Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 0, i_n_m2, j_n_m2, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 0, i_n_m3, j_n_m3, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 0, i_n_m4, j_n_m4, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 0, i_n_m5, j_n_m5, i_n_m5, j_n_m5 );
//                                                        phonon_s += s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n_m1, j_n_m1, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 1, i_n_m2, j_n_m2, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 1, i_n_m3, j_n_m3, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 1, i_n_m4, j_n_m4, i_n_m5, j_n_m5 );
//                                                        phonon_s += s.dgl_phonon_S_function( 2, i_n, j_n, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 2, i_n_m1, j_n_m1, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 2, i_n_m2, j_n_m2, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 2, i_n_m3, j_n_m3, i_n_m5, j_n_m5 );
//                                                        phonon_s += s.dgl_phonon_S_function( 3, i_n, j_n, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 3, i_n_m1, j_n_m1, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 3, i_n_m2, j_n_m2, i_n_m5, j_n_m5 );
//                                                        phonon_s += s.dgl_phonon_S_function( 4, i_n, j_n, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 4, i_n_m1, j_n_m1, i_n_m5, j_n_m5 );
//                                                        phonon_s += s.dgl_phonon_S_function( 5, i_n, j_n, i_n_m5, j_n_m5 );
//                                                        phonon_s = 0;
//                                                        Scalar val = it.value() * 1.0 * iterates.at( 4 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates.at( 3 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * iterates.at( 2 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) * iterates.at( 1 ).at( j_n_m4 + tensor_dim * i_n_m4 ).coeff( i_n_m3, j_n_m3 ) * iterates.at( 0 ).at( j_n_m5 + tensor_dim * i_n_m5 ).coeff( i_n_m4, j_n_m4 ) * std::exp( phonon_s );
//                                                        result += val;
//                                                        // Indexing is time-upwards, meaning t_n, t_n-1, t_n-2, t_n-3
//                                                        if ( fillADM && ( i_n == j_n || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
//                                                            //#pragma omp critical
//                                                            adm.addTriplet( { i_n, i_n_m1, i_n_m2, i_n_m3, i_n_m4, i_n_m5 }, { j_n, j_n_m1, j_n_m2, j_n_m3, j_n_m4, j_n_m5 }, val, omp_get_thread_num() );
//                                                        }
//                                                    }
//                                                }
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                rho_ret( i_n, j_n ) = result;
//            }
//        }
//    }
//
//    if ( fillADM ) {
//        adm.reduceDublicates();
//    }
//    return rho_ret.sparseView();
//}