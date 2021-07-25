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
                        phonon_s += s.dgl_phonon_S_function( lt - ld, ii_n, ij_n, ii_nd, ij_nd );
                    }
                }
                Scalar val = rho0.coeff( i_n_m, j_n_m ) * iterates.at( max_deph - 1 - current_deph ).at( i_n_m ).at( j_n_m ).coeff( i_n, j_n ) * std::exp( phonon_s );
                if ( fillADM && ( new_indicesX[0] == new_indicesY[0] || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
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
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
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
    return M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
}

std::vector<std::vector<Sparse>> ODESolver::calculate_propagator_vector( System &s, size_t tensor_dim, double t0, double t_step, std::vector<SaveState> &output ) {
    Sparse one = Dense::Identity( tensor_dim, tensor_dim ).sparseView();
    std::vector<std::vector<Sparse>> ret( tensor_dim, { tensor_dim, Sparse( tensor_dim, tensor_dim ) } );
    // Calculate first by hand to ensure the Hamilton gets calculated correclty
    ret[0][0] = calculate_propagator_single( s, tensor_dim, t0, t_step, 0, 0, output, one );
// Calculate the remaining propagators
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( i == 0 && j == 0 ) continue;
            ret[i][j] = calculate_propagator_single( s, tensor_dim, t0, t_step, i, j, output, one );
        }
    }
    return ret;
}

bool ODESolver::calculate_path_integral( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output ) {
    Log::L3( "Setting up Path-Integral Solver...\n" );
    output.reserve( s.parameters.iterations_t_max + 1 );
    size_t tensor_dim = rho0.rows();

    FixedSizeSparseMap<Scalar> adms = FixedSizeSparseMap<Scalar>( std::vector<int>( s.parameters.p_phonon_nc, tensor_dim ), s.parameters.numerics_maximum_threads );

    // Calculate first n_c timesteps directly
    Log::L3( "Calculating first n_c = {} timesteps of the pathintegral directly...\n", s.parameters.p_phonon_nc );

    // First step is just rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );
    // Next s.parameters.p_phonon_nc Steps and ADM initialization

    // Precalculate all iterates per index and timestep once. Saveorder is j+j_max*i
    Log::L3( "Pre-Calculating propagator matrices...\n" );
    auto t_prop = omp_get_wtime();
    std::vector<std::vector<std::vector<Sparse>>> propagators;
    // Calculate all propagators
    for ( int t = 0; t < s.parameters.p_phonon_nc - 1; t++ ) {
        propagators.emplace_back( calculate_propagator_vector( s, tensor_dim, s.parameters.t_step_pathint * t, t_step_initial, output ) );
    }
    Log::L3( "Done! Took {}s\n", omp_get_wtime() - t_prop );

    // Calculate the first n_c steps by recursion:
    for ( int n = 1; n < s.parameters.p_phonon_nc; n++ ) {
        Log::L3( "Calculating rho(t{})...\n", n );
        rho = path_integral( rho0, s, propagators, output, adms, n == s.parameters.p_phonon_nc - 1, n );
        saveState( rho, s.parameters.t_step_pathint * n, output );
    }

    //Log::L3( "Size of ADMs Tensor: {} bytes / {}\% filled\n", adms.getSizeOfValues(), adms.getFillRatio() * 100.0 );
    Log::L3( "Calculating Path-Integral Loop for the remaining {} timesteps...\n", std::floor( ( t_end - t_start ) / s.parameters.t_step_pathint ) - s.parameters.p_phonon_nc );
    // Iterate Path integral for further time steps

    double t_start_new = t_start + s.parameters.t_step_pathint * ( s.parameters.p_phonon_nc - 1 );
    for ( double t_t = t_start_new; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Path-Integral iteration
        if ( true ) {
            // Calculate Propagators for current time
            auto t0 = omp_get_wtime();
            auto propagator = calculate_propagator_vector( s, tensor_dim, t_t, t_step_initial, output );
            t0 = ( omp_get_wtime() - t0 );

            double cur_min = 1;

            // Main Iteration loop
            Log::L3( "Current Tensor Size: {} elements at a total of {} MB\r", adms.nonZeros(), adms.size() / 1024. / 1024. );
            auto t1 = omp_get_wtime();
            double total_append_time = 0;
            double total_time = 0;

            //for ( auto &threadvector : adms.get() )
            //auto tensor = adms.getIndices();
            //#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
            for ( auto &[sparse_index_x, outer_map] : adms.get() ) {
                for ( auto &[sparse_index_y, value] : outer_map ) {
                    if ( abs2( value ) == 0 ) continue;
                    auto tt = omp_get_wtime();
                    for ( int l = 0; l < propagator[sparse_index_x( 0 )][sparse_index_y( 0 )].outerSize(); ++l ) {
                        for ( Sparse::InnerIterator M( propagator[sparse_index_x( 0 )][sparse_index_y( 0 )], l ); M; ++M ) {
                            int i_n = M.row();
                            int j_n = M.col();

                            Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n );
                            for ( int tau = 0; tau < sparse_index_x.size(); tau++ ) {
                                int i_nd = sparse_index_x( tau );
                                int j_nd = sparse_index_y( tau );
                                phonon_s += s.dgl_phonon_S_function( tau + 1, i_n, j_n, i_nd, j_nd );
                            }

                            Scalar val = M.value() * value * std::exp( phonon_s );
                            double abs = abs2( val );
                            auto appendtime = omp_get_wtime();
                            // Add Element to Triplet list if they are diagonal elements or if they surpass the given threshold.
                            if ( i_n == j_n || abs >= s.parameters.numerics_pathintegral_squared_threshold ) {
                                cur_min = cur_min != 0.0 && cur_min < abs ? cur_min : abs;
                                adms.addTriplet( sparse_index_x, sparse_index_y, val, omp_get_thread_num(), i_n, j_n );
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

            Log::L3( "Reducing ADM\n" );
            // Calculate the reduced density matrix by tracing over all past times for each entry
            auto t2 = omp_get_wtime();
            // Calculate rho by tracing over the n_c last steps
            Dense cache = Dense::Zero( tensor_dim, tensor_dim );
            //#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
            //for ( auto &threadvector : adms.get() )
            for ( auto &[sparse_index_x, outer_map] : adms.get() ) {
                for ( auto &[sparse_index_y, value] : outer_map ) {
                    int i_n = sparse_index_x( 0 );
                    int j_n = sparse_index_y( 0 );
                    cache( i_n, j_n ) += value;
                }
            }
            rho = cache.sparseView();
            t2 = ( omp_get_wtime() - t2 );
            Log::L3( "Iteration: {}, time taken: [ Propagator: {:.4f}s, ADM Advancing: {:.4f}s (Partial append time: {:.4f}\%), ADM Setting: {:.4f}s, ADM Reduction: {:.4f}s ]\n", t_t, t0, t1, 100.0 * total_append_time / total_time, ts, t2 );

            // Dynamic Cutoff
            if ( s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max > 0 ) {
                double ratio = (double)adms.nonZeros() / s.parameters.numerics_pathintegral_dynamiccutoff_iterations_max;
                ratio = std::exp( ratio ) * std::pow( ratio, 5 );              // Scaled squared
                s.parameters.numerics_pathintegral_squared_threshold *= ratio; // Adjust Cutoff by squared ratio
                if ( s.parameters.numerics_pathintegral_squared_threshold < 1E-30 ) s.parameters.numerics_pathintegral_squared_threshold = 1E-30;
                if ( s.parameters.numerics_pathintegral_squared_threshold > 1E4 ) s.parameters.numerics_pathintegral_squared_threshold = 1E4;
                Log::L3( "Adjusted Cutoff Threshold to {} (x{})\n", std::sqrt( s.parameters.numerics_pathintegral_squared_threshold ), ratio );
            }
        } else {
            //propagators.emplace_back( empty_outer );
            ////propagators.erase( propagators.begin() );
            //
            //// The first iteration of this is not threadsafe due to the hamilton function caching itself
            //Sparse projector = Sparse( tensor_dim, tensor_dim );
            //projector.coeffRef( 0, 0 ) = 1;
            //Sparse M = iterate( projector, s, t_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
            //Sparse map;
            //if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
            //    map = project_matrix_sparse( M );
            //}
            //for ( double sub_t = t_t + t_step_initial * substepsize; sub_t < t_t + t_step_initial - t_step_initial * substepsize * 0.5; sub_t += t_step_initial * substepsize ) {
            //    M = iterate( M, s, sub_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
            //}
            //if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
            //    propagators.back().at( 0 ).at( 0 ) = ( M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
            //} else {
            //    propagators.back().at( 0 ).at( 0 ) = ( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
            //}
            //#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
            //for ( int t = 1; t < s.parameters.p_phonon_nc; t++ ) {
            //    for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
            //        for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
            //            if ( i_n == 0 && i_n == j_n ) continue;
            //            Sparse projector = Sparse( tensor_dim, tensor_dim );
            //            projector.coeffRef( i_n, j_n ) = 1;
            //            Sparse M = iterate( projector, s, t_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
            //            Sparse map;
            //            if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
            //                map = project_matrix_sparse( M );
            //            }
            //            for ( double sub_t = t_t + t_step_initial * substepsize; sub_t < t_t + t_step_initial - t_step_initial * substepsize * 0.5; sub_t += t_step_initial * substepsize ) {
            //                M = iterate( M, s, sub_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
            //            }
            //            if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
            //                propagators.back().at( i_n ).at( j_n ) = ( M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
            //            } else {
            //                propagators.back().at( i_n ).at( j_n ) = ( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
            //            }
            //        }
            //    }
            //}
            ////rho = path_integral( rho, s, propagators, output, adms, s.parameters.numerics_pathintegral_squared_threshold, false, s.parameters.p_phonon_nc - 1 );
            //rho = path_integral( rho0, s, propagators, output, adms, false, (int)std::floor( t_t / s.parameters.t_step ) );
        }
        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );
        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
        }
    }
    //Log::L3( "Final Size of ADMs Tensor: {} bytes / {}\% filled\n", adms.getSizeOfValues(), adms.getFillRatio() * 100.0 );
    Log::L3( "Done!\n" );
    return true;
}

//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
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
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
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
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
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
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
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
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
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
//#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
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