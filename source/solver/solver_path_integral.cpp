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
                    adm.addTriplet( new_indicesX, new_indicesY, adm_value * val, omp_get_thread_num() );
                }
                result += val;
            }
        }
    } else {
        for ( int i_n_m = 0; i_n_m < tensor_dim; i_n_m++ ) {
            for ( int j_n_m = 0; j_n_m < tensor_dim; j_n_m++ ) {
                Scalar phonon_s = 0;
                auto new_indicesX = indicesX;
                auto new_indicesY = indicesY;
                new_indicesX.emplace_back( i_n_m );
                new_indicesY.emplace_back( j_n_m );
                Scalar val = iterates.at( max_deph - 1 - current_deph ).at( i_n_m ).at( j_n_m ).coeff( i_n, j_n );
                if ( abs2( val ) == 0 ) continue;
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
                    if ( fillADM && ( new_indicesX[0] == new_indicesY[0] || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
                        adm.addTriplet( new_indicesX, new_indicesY, adm_value * val, omp_get_thread_num() );
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

bool ODESolver::calculate_path_integral( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output ) {
    Log::L3( "Setting up Path-Integral Solver...\n" );
    output.reserve( s.parameters.iterations_t_max + 1 );
    int tensor_dim = rho0.rows();
    // Configurables
    double substepsize = s.parameters.numerics_pathintegral_stepsize_iterator / t_step_initial;
    Log::L2( "Path Integral Substepsize is 100fs -> {} delta\n", substepsize );

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
    std::vector<Sparse> empty_inner( tensor_dim, Sparse( tensor_dim, tensor_dim ) );
    std::vector<std::vector<Sparse>> empty_outer( tensor_dim, empty_inner );
    for ( int t = 1; t < s.parameters.p_phonon_nc; t++ ) {
        propagators.emplace_back( empty_outer );
    }
    for ( int t = 1; t < s.parameters.p_phonon_nc; t++ ) {
        // The first iteration of this is not threadsafe due to the hamilton function caching itself
        Sparse projector = Sparse( tensor_dim, tensor_dim );
        projector.coeffRef( 0, 0 ) = 1;
        Sparse M = iterate( projector, s, t_step_initial * t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
        Sparse map;
        if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
            map = project_matrix_sparse( M );
        }
        for ( double sub_t = t_step_initial * t + t_step_initial * substepsize; sub_t < t_step_initial * ( t + 1 ) - t_step_initial * substepsize * 0.5; sub_t += t_step_initial * substepsize ) {
            M = iterate( M, s, sub_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
        }
        if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
            propagators.at( t - 1 ).at( 0 ).at( 0 ) = ( M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
        } else {
            propagators.at( t - 1 ).at( 0 ).at( 0 ) = ( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
        }
    }
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
    for ( int t = 1; t < s.parameters.p_phonon_nc; t++ ) {
        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                if ( i_n == 0 && i_n == j_n ) continue;
                Sparse projector = Sparse( tensor_dim, tensor_dim );
                projector.coeffRef( i_n, j_n ) = 1;
                Sparse M = iterate( projector, s, t_step_initial * t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
                Sparse map;
                if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                    map = project_matrix_sparse( M );
                }
                for ( double sub_t = t_step_initial * t + t_step_initial * substepsize; sub_t < t_step_initial * ( t + 1 ) - t_step_initial * substepsize * 0.5; sub_t += t_step_initial * substepsize ) {
                    M = iterate( M, s, sub_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
                }
                if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                    propagators.at( t - 1 ).at( i_n ).at( j_n ) = ( M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
                } else {
                    propagators.at( t - 1 ).at( i_n ).at( j_n ) = ( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
                }
            }
        }
    }
    Log::L3( "Done! Took {}s\n", omp_get_wtime() - t_prop );

    // Calculate the first n_c steps by recursion:
    for ( int n = 1; n < s.parameters.p_phonon_nc; n++ ) {
        Log::L3( "Calculating rho(t{})...\n", n );
        //rho = path_integral_bruteforce( rho0, s, iterates, output, adms, s.parameters.numerics_pathintegral_squared_threshold, n == s.parameters.p_phonon_nc - 1, n );
        rho = path_integral( rho0, s, propagators, output, adms, n == s.parameters.p_phonon_nc - 1, n );
        saveState( rho, t_step_initial * n, output );
    }

    //Log::L3( "Size of ADMs Tensor: {} bytes / {}\% filled\n", adms.getSizeOfValues(), adms.getFillRatio() * 100.0 );
    Log::L3( "Calculating Path-Integral Loop for the remaining {} timesteps...\n", std::floor( ( t_end - t_start ) / t_step_initial ) - s.parameters.p_phonon_nc );
    // Iterate Path integral for further time steps

    std::ofstream test( s.parameters.subfolder + "test.txt", std::ofstream::out );
    test << "Time\tADM_Size(#)\tADM_Size(Bytes)\tMin\tCachevectorsizes\n";
    double t_start_new = t_start + t_step_initial * s.parameters.p_phonon_nc;
    for ( double t_t = t_start_new; t_t <= t_end; t_t += t_step_initial ) {
        // Path-Integral iteration
        if ( true ) {
            std::vector<std::vector<Sparse>> propagator;
            propagator.reserve( tensor_dim );
            for ( int i = 0; i < tensor_dim; i++ ) {
                propagator.emplace_back( std::vector<Sparse>() );
            }
            // The first iteration of this is not threadsafe due to the hamilton function caching itself
            Sparse projector = Sparse( tensor_dim, tensor_dim );
            projector.coeffRef( 0, 0 ) = 1;
            Sparse M = iterate( projector, s, t_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
            Sparse map;
            if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                map = project_matrix_sparse( M );
            }
            std::string nonzerros = fmt::format( "[{:d} ", (int)M.nonZeros() );
            for ( double sub_t = t_t + t_step_initial * substepsize; sub_t < t_t + t_step_initial - t_step_initial * substepsize / 2.0; sub_t += t_step_initial * substepsize ) {
                M = iterate( M, s, sub_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ); // The first iteration of this is not threadsafe due to the hamilton function caching itself
                nonzerros = fmt::format( "{}{:d} ", nonzerros, (int)M.nonZeros() );
            }
            if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                propagator.at( 0 ).emplace_back( M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
            } else {
                propagator.at( 0 ).emplace_back( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
            }
            nonzerros = fmt::format( "{} -> {:d}]", nonzerros, (int)propagator[0].back().nonZeros() );
            Log::L3( "Non-Zeros evolution for 0,0 projector propagation: {}\n", nonzerros );

            //Log::L3( "Curren time t_t = {}, subtimes are '{}'\n", t_t, times );
            auto t0 = omp_get_wtime();
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
            for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
                for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                    if ( !( i_n == j_n && i_n == 0 ) ) {
                        Sparse projector = Sparse( tensor_dim, tensor_dim );
                        projector.coeffRef( i_n, j_n ) = 1;
                        Sparse M = iterate( projector, s, t_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
                        Sparse map;
                        if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                            map = project_matrix_sparse( M );
                        }
                        for ( double sub_t = t_t + t_step_initial * substepsize; sub_t < t_t + t_step_initial - t_step_initial * substepsize / 2.0; sub_t += t_step_initial * substepsize ) {
                            M = iterate( M, s, sub_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ); // The first iteration of this is not threadsafe due to the hamilton function caching itself
                        }
                        if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                            propagator.at( i_n ).emplace_back( M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
                        } else {
                            propagator.at( i_n ).emplace_back( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
                        }
                    }
                }
            }
            t0 = ( omp_get_wtime() - t0 );

            double cur_min = 1;
            test << t_t << "\t" << adms.nonZeros() << "\t" << adms.getSizeOfTensor() << "\t";

            // Main Iteration loop
            Log::L3( "current non zeros in adm: {}   \r", adms.nonZeros() );
            auto t1 = omp_get_wtime();
            double total_append_time = 0;
            double total_time = 0;

            //for ( size_t b = 0; b < adms.get().bucket_count(); b++ )
            //for ( auto bi = adms.get().begin( b ); bi != adms.get().end( b ); bi++ ) {
            //auto &triplet = bi->second;
            //#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
            for ( auto &outer : adms.get() )
                for ( auto &inner : outer ) {
                    auto &triplet = inner.second;
                    auto tt = omp_get_wtime();
                    //TODO: TEST propagator hier nochmal ausrechnen, damits auch definitiv der richtige ist.
                    // Dafür: code oben auslagern in calculate_propagator_single (die dann auch hier testweise benutzen), und calculate_propagators (die dann in vektor)
                    for ( int l = 0; l < propagator[triplet.indicesX[0]][triplet.indicesY[0]].outerSize(); ++l ) {
                        for ( Sparse::InnerIterator it( propagator[triplet.indicesX[0]][triplet.indicesY[0]], l ); it; ++it ) {
                            int i_n = it.row();
                            int j_n = it.col();

                            Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n );
                            for ( int tau = 0; tau < triplet.indicesX.size(); tau++ ) {
                                int i_nd = triplet.indicesX[tau];
                                int j_nd = triplet.indicesY[tau];
                                phonon_s += s.dgl_phonon_S_function( tau + 1, i_n, j_n, i_nd, j_nd );
                            }

                            Scalar val = it.value() * triplet.value * std::exp( phonon_s );
                            double abs = abs2( val );
                            auto appendtime = omp_get_wtime();
                            // Add Element to Triplet list if they are diagonal elements or if they surpass the given threshold.
                            if ( i_n == j_n || abs >= s.parameters.numerics_pathintegral_squared_threshold ) {
                                cur_min = cur_min != 0.0 && cur_min < abs ? cur_min : abs;
                                adms.addTriplet( triplet.indicesX, triplet.indicesY, val, omp_get_thread_num(), i_n, j_n );
                            }
                            total_append_time += omp_get_wtime() - appendtime;
                        }
                    }
                    total_time += omp_get_wtime() - tt;
                }
            t1 = ( omp_get_wtime() - t1 );

            test << cur_min << "\t" << adms.getSizesOfCache() << std::endl;

            auto ts = omp_get_wtime();
            adms.reduceDublicates();
            ts = omp_get_wtime() - ts;

            Log::L3( "Reducing ADM\n" );
            // Calculate the reduced density matrix by tracing over all past times for each entry
            auto t2 = omp_get_wtime();
            // Calculate rho by tracing over the n_c last steps
            Dense cache = Dense::Zero( tensor_dim, tensor_dim );
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
            for ( auto &outer : adms.get() )
                for ( auto &inner : outer ) {
                    auto &triplet = inner.second;
                    int i_n = triplet.indicesX[0];
                    int j_n = triplet.indicesY[0];
#pragma omp critical
                    cache( i_n, j_n ) += triplet.value;
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
            propagators.emplace_back( empty_outer );
            //propagators.erase( propagators.begin() );

            // The first iteration of this is not threadsafe due to the hamilton function caching itself
            Sparse projector = Sparse( tensor_dim, tensor_dim );
            projector.coeffRef( 0, 0 ) = 1;
            Sparse M = iterate( projector, s, t_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
            Sparse map;
            if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                map = project_matrix_sparse( M );
            }
            for ( double sub_t = t_t + t_step_initial * substepsize; sub_t < t_t + t_step_initial - t_step_initial * substepsize * 0.5; sub_t += t_step_initial * substepsize ) {
                M = iterate( M, s, sub_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
            }
            if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                propagators.back().at( 0 ).at( 0 ) = ( M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
            } else {
                propagators.back().at( 0 ).at( 0 ) = ( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
            }
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
            for ( int t = 1; t < s.parameters.p_phonon_nc; t++ ) {
                for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
                    for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                        if ( i_n == 0 && i_n == j_n ) continue;
                        Sparse projector = Sparse( tensor_dim, tensor_dim );
                        projector.coeffRef( i_n, j_n ) = 1;
                        Sparse M = iterate( projector, s, t_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
                        Sparse map;
                        if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                            map = project_matrix_sparse( M );
                        }
                        for ( double sub_t = t_t + t_step_initial * substepsize; sub_t < t_t + t_step_initial - t_step_initial * substepsize * 0.5; sub_t += t_step_initial * substepsize ) {
                            M = iterate( M, s, sub_t, t_step_initial * substepsize, output ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold );
                        }
                        if ( s.parameters.numerics_pathintegral_docutoff_propagator ) {
                            propagators.back().at( i_n ).at( j_n ) = ( M.cwiseProduct( map ).pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
                        } else {
                            propagators.back().at( i_n ).at( j_n ) = ( M.pruned( s.parameters.numerics_pathintegral_sparse_prune_threshold ) );
                        }
                    }
                }
            }
            //rho = path_integral( rho, s, propagators, output, adms, s.parameters.numerics_pathintegral_squared_threshold, false, s.parameters.p_phonon_nc - 1 );
            rho = path_integral( rho0, s, propagators, output, adms, false, (int)std::floor( t_t / s.parameters.t_step ) );
        }
        // Save Rho
        saveState( rho, t_t, output );
        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
        }
    }
    test.close();
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
Sparse ODESolver::path_integral_bruteforce( const Sparse &rho0, System &s, std::vector<std::vector<Sparse>> &iterates, std::vector<SaveState> &output, FixedSizeSparseMap<Scalar> &adm, bool fillADM, int deph ) {
    int tensor_dim = rho0.rows();
    Dense rho_ret = Dense::Zero( tensor_dim, tensor_dim );

    if ( deph == 1 ) {
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                Scalar result = 0;
                for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
                    for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
                        Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 );
                        phonon_s = 0;
                        result += rho0.coeff( i_n_m1, j_n_m1 ) * 1.0 * iterates.at( 0 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * std::exp( phonon_s );
                    }
                }
                rho_ret( i_n, j_n ) = result;
            }
        }
    } else if ( deph == 2 ) {
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                Scalar result = 0;
                for ( int k = 0; k < rho0.outerSize(); ++k ) {
                    for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
                        auto i_n_m2 = (int)it.row();
                        auto j_n_m2 = (int)it.col();
                        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
                            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
                                if ( abs2( iterates.at( 0 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
                                Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 0, i_n_m2, j_n_m2, i_n_m2, j_n_m2 );
                                phonon_s += s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n_m1, j_n_m1, i_n_m2, j_n_m2 );
                                phonon_s += s.dgl_phonon_S_function( 2, i_n, j_n, i_n_m2, j_n_m2 );
                                phonon_s = 0;
                                Scalar val = it.value() * 1.0 * iterates.at( 1 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates.at( 0 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * std::exp( phonon_s );
                                result += val;
                                if ( fillADM && ( i_n == j_n || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
                                    //#pragma omp critical
                                    adm.addTriplet( { i_n, i_n_m1, i_n_m2 }, { j_n, j_n_m1, j_n_m2 }, val, omp_get_thread_num() );
                                }
                            }
                        }
                    }
                }
                rho_ret( i_n, j_n ) = result;
            }
        }
    } else if ( deph == 3 ) {
        int a = 0;
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                //fmt::print( "N_c = 3 -- i = {}, j = {} -> {:3.0f}\%\n", i_n, j_n, 100.0 * a / 36 / 36 );
                //#pragma omp critical
                a++;
                Scalar result = 0;
                for ( int k = 0; k < rho0.outerSize(); ++k ) {
                    for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
                        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
                            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
                                if ( abs2( iterates.at( 2 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) ) == 0.0 ) continue;
                                for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
                                    for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
                                        if ( abs2( iterates.at( 1 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
                                        auto i_n_m3 = (int)it.row();
                                        auto j_n_m3 = (int)it.col();
                                        Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 0, i_n_m2, j_n_m2, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 0, i_n_m3, j_n_m3, i_n_m3, j_n_m3 );
                                        phonon_s += s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n_m1, j_n_m1, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 1, i_n_m2, j_n_m2, i_n_m3, j_n_m3 );
                                        phonon_s += s.dgl_phonon_S_function( 2, i_n, j_n, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 2, i_n_m1, j_n_m1, i_n_m3, j_n_m3 );
                                        phonon_s += s.dgl_phonon_S_function( 3, i_n, j_n, i_n_m3, j_n_m3 );
                                        Scalar val = it.value() * 1.0 * iterates.at( 2 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates.at( 1 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * iterates.at( 0 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) * std::exp( phonon_s );
                                        phonon_s = 0;
                                        result += val;
                                        // Indexing is time-upwards, meaning t_n, t_n-1, t_n-2, t_n-3
                                        if ( fillADM && ( i_n == j_n || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
                                            //#pragma omp critical
                                            adm.addTriplet( { i_n, i_n_m1, i_n_m2, i_n_m3 }, { j_n, j_n_m1, j_n_m2, j_n_m3 }, val, omp_get_thread_num() );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                rho_ret( i_n, j_n ) = result;
            }
        }
    } else if ( deph == 4 ) {
        int a = 0;
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                //fmt::print( "N_c = 4 -- i = {}, j = {} -> {:3.0f}\%\n", i_n, j_n, 100.0 * a / 36 / 36 );
                //#pragma omp critical
                a++;
                Scalar result = 0;
                for ( int k = 0; k < rho0.outerSize(); ++k ) {
                    for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
                        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
                            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
                                if ( abs2( iterates.at( 3 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) ) == 0.0 ) continue;
                                for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
                                    for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
                                        if ( abs2( iterates.at( 2 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
                                        for ( int i_n_m3 = 0; i_n_m3 < tensor_dim; i_n_m3++ ) {
                                            for ( int j_n_m3 = 0; j_n_m3 < tensor_dim; j_n_m3++ ) {
                                                if ( abs2( iterates.at( 1 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) ) == 0.0 ) continue;
                                                auto i_n_m4 = (int)it.row();
                                                auto j_n_m4 = (int)it.col();
                                                Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 0, i_n_m2, j_n_m2, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 0, i_n_m3, j_n_m3, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 0, i_n_m4, j_n_m4, i_n_m4, j_n_m4 );
                                                phonon_s += s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n_m1, j_n_m1, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 1, i_n_m2, j_n_m2, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 1, i_n_m3, j_n_m3, i_n_m4, j_n_m4 );
                                                phonon_s += s.dgl_phonon_S_function( 2, i_n, j_n, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 2, i_n_m1, j_n_m1, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 2, i_n_m2, j_n_m2, i_n_m4, j_n_m4 );
                                                phonon_s += s.dgl_phonon_S_function( 3, i_n, j_n, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 3, i_n_m1, j_n_m1, i_n_m4, j_n_m4 );
                                                phonon_s += s.dgl_phonon_S_function( 4, i_n, j_n, i_n_m4, j_n_m4 );
                                                Scalar val = it.value() * 1.0 * iterates.at( 3 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates.at( 2 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * iterates.at( 1 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) * iterates.at( 0 ).at( j_n_m4 + tensor_dim * i_n_m4 ).coeff( i_n_m3, j_n_m3 ) * std::exp( phonon_s );
                                                phonon_s = 0;
                                                result += val;
                                                // Indexing is time-upwards, meaning t_n, t_n-1, t_n-2, t_n-3
                                                if ( fillADM && ( i_n == j_n || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
                                                    //#pragma omp critical
                                                    adm.addTriplet( { i_n, i_n_m1, i_n_m2, i_n_m3, i_n_m4 }, { j_n, j_n_m1, j_n_m2, j_n_m3, j_n_m4 }, val, omp_get_thread_num() );
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                rho_ret( i_n, j_n ) = result;
            }
        }
    } else if ( deph == 5 ) {
        int a = 0;
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                //fmt::print( "N_c = 5 --  i = {}, j = {} -> {:3.0f}\%\n", i_n, j_n, 100.0 * a / 36 / 36 );
                a++;
                Scalar result = 0;
                for ( int k = 0; k < rho0.outerSize(); ++k ) {
                    for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
                        for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
                            for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
                                if ( abs2( iterates.at( 4 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) ) == 0.0 ) continue;
                                for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
                                    for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
                                        if ( abs2( iterates.at( 3 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
                                        for ( int i_n_m3 = 0; i_n_m3 < tensor_dim; i_n_m3++ ) {
                                            for ( int j_n_m3 = 0; j_n_m3 < tensor_dim; j_n_m3++ ) {
                                                if ( abs2( iterates.at( 2 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) ) == 0.0 ) continue;
                                                for ( int i_n_m4 = 0; i_n_m4 < tensor_dim; i_n_m4++ ) {
                                                    for ( int j_n_m4 = 0; j_n_m4 < tensor_dim; j_n_m4++ ) {
                                                        if ( std::abs( iterates.at( 1 ).at( j_n_m4 + tensor_dim * i_n_m4 ).coeff( i_n_m3, j_n_m3 ) ) == 0.0 ) continue;
                                                        auto i_n_m5 = (int)it.row();
                                                        auto j_n_m5 = (int)it.col();
                                                        Scalar phonon_s = s.dgl_phonon_S_function( 0, i_n, j_n, i_n, j_n ) + s.dgl_phonon_S_function( 0, i_n_m1, j_n_m1, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 0, i_n_m2, j_n_m2, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 0, i_n_m3, j_n_m3, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 0, i_n_m4, j_n_m4, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 0, i_n_m5, j_n_m5, i_n_m5, j_n_m5 );
                                                        phonon_s += s.dgl_phonon_S_function( 1, i_n, j_n, i_n_m1, j_n_m1 ) + s.dgl_phonon_S_function( 1, i_n_m1, j_n_m1, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 1, i_n_m2, j_n_m2, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 1, i_n_m3, j_n_m3, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 1, i_n_m4, j_n_m4, i_n_m5, j_n_m5 );
                                                        phonon_s += s.dgl_phonon_S_function( 2, i_n, j_n, i_n_m2, j_n_m2 ) + s.dgl_phonon_S_function( 2, i_n_m1, j_n_m1, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 2, i_n_m2, j_n_m2, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 2, i_n_m3, j_n_m3, i_n_m5, j_n_m5 );
                                                        phonon_s += s.dgl_phonon_S_function( 3, i_n, j_n, i_n_m3, j_n_m3 ) + s.dgl_phonon_S_function( 3, i_n_m1, j_n_m1, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 3, i_n_m2, j_n_m2, i_n_m5, j_n_m5 );
                                                        phonon_s += s.dgl_phonon_S_function( 4, i_n, j_n, i_n_m4, j_n_m4 ) + s.dgl_phonon_S_function( 4, i_n_m1, j_n_m1, i_n_m5, j_n_m5 );
                                                        phonon_s += s.dgl_phonon_S_function( 5, i_n, j_n, i_n_m5, j_n_m5 );
                                                        phonon_s = 0;
                                                        Scalar val = it.value() * 1.0 * iterates.at( 4 ).at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates.at( 3 ).at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * iterates.at( 2 ).at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 ) * iterates.at( 1 ).at( j_n_m4 + tensor_dim * i_n_m4 ).coeff( i_n_m3, j_n_m3 ) * iterates.at( 0 ).at( j_n_m5 + tensor_dim * i_n_m5 ).coeff( i_n_m4, j_n_m4 ) * std::exp( phonon_s );
                                                        result += val;
                                                        // Indexing is time-upwards, meaning t_n, t_n-1, t_n-2, t_n-3
                                                        if ( fillADM && ( i_n == j_n || abs2( val ) > s.parameters.numerics_pathintegral_squared_threshold ) ) {
                                                            //#pragma omp critical
                                                            adm.addTriplet( { i_n, i_n_m1, i_n_m2, i_n_m3, i_n_m4, i_n_m5 }, { j_n, j_n_m1, j_n_m2, j_n_m3, j_n_m4, j_n_m5 }, val, omp_get_thread_num() );
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                rho_ret( i_n, j_n ) = result;
            }
        }
    }

    if ( fillADM ) {
        adm.reduceDublicates();
    }
    return rho_ret.sparseView();
}