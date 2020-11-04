#include "solver/solver_ode.h"

Sparse ODESolver::iterateRungeKutta4( const Sparse &rho, System &s, const double t, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-4 ausrechnen
    Sparse H_calc_k1 = getHamilton( s, t );
    Sparse H_calc_k23 = getHamilton( s, t + s.parameters.t_step * 0.5 );
    Sparse H_calc_k4 = getHamilton( s, t + s.parameters.t_step );
    // k1-4 ausrechnen
    Sparse rk1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    Sparse rk2 = s.dgl_rungeFunction( rho + s.parameters.t_step * 0.5 * rk1, H_calc_k23, t + s.parameters.t_step * 0.5, savedStates );
    Sparse rk3 = s.dgl_rungeFunction( rho + s.parameters.t_step * 0.5 * rk2, H_calc_k23, t + s.parameters.t_step * 0.5, savedStates );
    Sparse rk4 = s.dgl_rungeFunction( rho + s.parameters.t_step * rk3, H_calc_k4, t + s.parameters.t_step, savedStates );
    // Dichtematrix
    return rho + s.parameters.t_step / 6.0 * ( rk1 + 2. * rk2 + 2. * rk3 + rk4 );
}

Sparse ODESolver::iterateRungeKutta5( const Sparse &rho, System &s, const double t, std::vector<SaveState> &savedStates ) {
    // Verschiedene H's fuer k1-6 ausrechnen
    Sparse H_calc_k1 = getHamilton( s, t );
    Sparse H_calc_k2 = getHamilton( s, t + a2 * s.parameters.t_step );
    Sparse H_calc_k3 = getHamilton( s, t + a3 * s.parameters.t_step );
    Sparse H_calc_k4 = getHamilton( s, t + a4 * s.parameters.t_step );
    Sparse H_calc_k5 = getHamilton( s, t + a5 * s.parameters.t_step );
    Sparse H_calc_k6 = getHamilton( s, t + a6 * s.parameters.t_step );
    // k1-6 ausrechnen
    Sparse k1 = s.dgl_rungeFunction( rho, H_calc_k1, t, savedStates );
    Sparse k2 = s.dgl_rungeFunction( rho + s.parameters.t_step * b11 * k1, H_calc_k2, t + a2 * s.parameters.t_step, savedStates );
    Sparse k3 = s.dgl_rungeFunction( rho + s.parameters.t_step * ( b21 * k1 + b22 * k2 ), H_calc_k3, t + a3 * s.parameters.t_step, savedStates );
    Sparse k4 = s.dgl_rungeFunction( rho + s.parameters.t_step * ( b31 * k1 + b32 * k2 + b33 * k3 ), H_calc_k4, t + a4 * s.parameters.t_step, savedStates );
    Sparse k5 = s.dgl_rungeFunction( rho + s.parameters.t_step * ( b41 * k1 + b42 * k2 + b43 * k3 + b44 * k4 ), H_calc_k5, t + a5 * s.parameters.t_step, savedStates );
    Sparse k6 = s.dgl_rungeFunction( rho + s.parameters.t_step * ( b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4 + b55 * k5 ), H_calc_k6, t + a6 * s.parameters.t_step, savedStates );
    // Dichtematrix
    return rho + s.parameters.t_step * ( b61 * k1 + b63 * k3 + b64 * k4 + b65 * k5 + b66 * k6 );
}

Sparse ODESolver::iterate( const Sparse &rho, System &s, const double t, std::vector<SaveState> &savedStates, const int dir ) {
    int order = dir == DIR_T ? s.parameters.numerics_order_t : s.parameters.numerics_order_tau;
    if ( order == 4 ) {
        return iterateRungeKutta4( rho, s, t, savedStates );
    }
    return iterateRungeKutta5( rho, s, t, savedStates );
}

bool ODESolver::calculate_runge_kutta( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output ) {
    Log::L3( "Setting up Runge-Kutta Solver...\n" );
    // Reserve Output Vector
    output.reserve( s.parameters.iterations_t_max + 1 );
    // Save initial value
    saveState( rho0, t_start, output );
    // Calculate Remaining
    Sparse rho = rho0;
    Log::L3( "Calculating Runge-Kutta Loop...\n" );
    for ( double t_t = t_start + t_step_initial; t_t <= t_end; t_t += t_step_initial ) {
        // Runge-Kutta iteration
        rho = iterate( rho, s, t_t, output );
        // Save Rho
        saveState( rho, t_t, output );
        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
        }
    }
    Log::L3( "Done!\n" );
    return true;
}

//Scalar ODESolver::path_integral_recursive( const Sparse &rho0, System &s, std::vector<SaveState> &output, int max_deph, int outer_i, int outer_j, int current_deph ) {
//    Scalar result = 0;
//    int tensor_dim = rho0.rows();
//
//    for ( int i = 0; i < tensor_dim; i++ ) {
//        for ( int j = 0; j < tensor_dim; j++ ) {
//            Dense projector = Dense::Zero( tensor_dim, tensor_dim );
//            projector( i, j ) = 1;
//            Dense M = iterate( projector.sparseView(), s, s.parameters.t_step * ( max_deph - current_deph ), output ).toDense();
//            if ( std::abs( M( outer_i, outer_j ) ) <= 1E-15 ) continue;
//            if ( current_deph == max_deph - 1 ) {
//                result += rho0.coeff( i, j ) * M( outer_i, outer_j );
//            } else {
//                result += M( outer_i, outer_j ) * path_integral_recursive( rho0, s, output, max_deph, i, j, current_deph + 1 );
//            }
//        }
//    }
//
//    return result;
//}
//Dense rho_ret = Dense::Zero( tensor_dim, tensor_dim );
//int ncmax = 3;
//int max_deph = std::min( ncmax, t_index );
//Dense r0 = output.at( output.size() - max_deph ).mat;
//#pragma omp parallel for collapse( 2 )
//for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
//    for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
//        Scalar result = path_integral_recursive( rho0, s, output, max_deph, i_n, j_n );
//        rho_ret( i_n, j_n ) = result;
//    }
//}
//rho = rho_ret.sparseView();

bool ODESolver::calculate_path_integral( Sparse &rho0, double t_start, double t_end, double t_step_initial, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<SaveState> &output, bool do_output ) {
    // Erstmal alle M, S sachen hier lokal machen, dann später outsourcen wie bei der RK PMEQ auch (S in system precalculaten, M als weiteren saveState vector)
    // FIXME: Assume zeros photons for now
    Log::L3( "Setting up Path-Integral Solver...\n" );
    output.reserve( s.parameters.iterations_t_max + 1 );
    int n_c_max = 4;
    int tensor_dim = rho0.rows();

    FixedSizeSparseMap<Scalar> adms = FixedSizeSparseMap<Scalar>( {tensor_dim, tensor_dim, tensor_dim, tensor_dim} );

    Log::L3( "Size of ADMs Tensor: {} bytes\n", sizeof( adms ) );
    // Calculate first n_c timesteps directly
    Log::L3( "Calculating first n_c = {} timesteps of the pathintegral directly...\n", n_c_max );

    // First step is just rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );
    // Next 3 Steps and ADM initialization

    // Precalculate all iterates per index and timestep once. Saveorder is j+j_max*i
    Log::L3( "Pre-Calculating propagator matrices...\n" );
    std::vector<Sparse> iterates_t0, iterates_t1, iterates_t2;
    for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
        for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
            Sparse projector = Sparse( tensor_dim, tensor_dim );
            projector.coeffRef( i_n, j_n ) = 1;
            Sparse M_t0 = iterate( projector, s, 0.333 * t_step_initial, output );
            Sparse M_t1 = iterate( projector, s, 0.666 * t_step_initial, output );
            Sparse M_t2 = iterate( projector, s, 1.0 * t_step_initial, output );
            iterates_t0.emplace_back( M_t0 );
            iterates_t1.emplace_back( M_t1 );
            iterates_t2.emplace_back( M_t2 );
        }
    }

    // Calculate the first 3 steps by hand:
    Log::L3( "Calculating rho(t1)...\n" );
    Dense rho_ret = Dense::Zero( tensor_dim, tensor_dim );
    Dense r0 = rho0.toDense();
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
        for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
            Scalar result = 0;
            for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
                for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
                    result += r0( i_n_m1, j_n_m1 ) * 1.0 * iterates_t0.at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n );
                }
            }
            rho_ret( i_n, j_n ) = result;
        }
    }
    rho = rho_ret.sparseView();
    saveState( rho, t_step_initial, output );

    Log::L3( "Calculating rho(t2)...\n" );
    rho_ret = Dense::Zero( tensor_dim, tensor_dim );
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
        for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
            Scalar result = 0;
            for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
                for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
                    for ( int k = 0; k < rho0.outerSize(); ++k ) {
                        for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
                            auto i_n_m2 = it.row();
                            auto j_n_m2 = it.col();
                            result += it.value() * 1.0 * iterates_t1.at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates_t0.at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 );
                        }
                    }
                }
            }
            rho_ret( i_n, j_n ) = result;
        }
    }
    rho = rho_ret.sparseView();
    saveState( rho, 2.0 * t_step_initial, output );

    Log::L3( "Calculating rho(t3) while building ADM Tensor...\n" );
    rho_ret = Dense::Zero( tensor_dim, tensor_dim );

    std::vector<std::vector<Eigen::Triplet<Scalar>>> cpu_caches;
    for ( int i = 0; i < s.parameters.numerics_maximum_threads; i++ ) {
        cpu_caches.emplace_back( std::vector<Eigen::Triplet<Scalar>>() );
    }
#pragma omp parallel for collapse( 2 ) num_threads( s.parameters.numerics_maximum_threads )
    for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
        for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
            Scalar result = 0;
            for ( int i_n_m1 = 0; i_n_m1 < tensor_dim; i_n_m1++ ) {
                for ( int j_n_m1 = 0; j_n_m1 < tensor_dim; j_n_m1++ ) {
                    if ( std::abs( iterates_t2.at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) ) == 0.0 ) continue;
                    for ( int i_n_m2 = 0; i_n_m2 < tensor_dim; i_n_m2++ ) {
                        for ( int j_n_m2 = 0; j_n_m2 < tensor_dim; j_n_m2++ ) {
                            if ( std::abs( iterates_t1.at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) ) == 0.0 ) continue;
                            for ( int k = 0; k < rho0.outerSize(); ++k ) {
                                for ( Sparse::InnerIterator it( rho0, k ); it; ++it ) {
                                    auto i_n_m3 = it.row();
                                    auto j_n_m3 = it.col();
                                    Scalar val = it.value() * 1.0 * iterates_t2.at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * iterates_t1.at( j_n_m2 + tensor_dim * i_n_m2 ).coeff( i_n_m1, j_n_m1 ) * iterates_t0.at( j_n_m3 + tensor_dim * i_n_m3 ).coeff( i_n_m2, j_n_m2 );
                                    result += val;
                                    // Indexing is time-upwards, meaning t_n, t_n-1, t_n-2, t_n-3
                                    if ( std::abs( val ) != 0 ) {
                                        adms.addTripletTo( i_n, i_n_m1, i_n_m2, i_n_m3, j_n, j_n_m1, j_n_m2, j_n_m3, val, cpu_caches.at( omp_get_thread_num() ) );
                                    }
                                }
                            }
                        }
                    }
                }
            }
            //if ( omp_get_thread_num() == 0 )
            std::cout << "Current progress: " << 100.0 * ( i_n * tensor_dim + j_n ) / tensor_dim / tensor_dim << "\%\r";
            rho_ret( i_n, j_n ) = result;
        }
    }
    for ( int i = 0; i < s.parameters.numerics_maximum_threads; i++ ) {
        adms.getTriplets().insert( adms.getTriplets().begin(), cpu_caches.at( i ).begin(), cpu_caches.at( i ).end() );
    }

    adms.setFromTripletList();
    rho = rho_ret.sparseView();
    saveState( rho, 3.0 * t_step_initial, output );

    double t_start_new = t_start + t_step_initial * n_c_max;
    rho = output.back().mat;
    Log::L3( "Size of ADMs Tensor: {} bytes / {}\% filled\n", adms.getSizeOfCache(), adms.getFillRatio() * 100.0 );
    Log::L3( "Calculating Path-Integral Loop for the remaining {} timesteps...\n", std::floor( ( t_end - t_start ) / t_step_initial ) - n_c_max );
    // Iterate Path integral for further time steps

    std::ofstream test( s.parameters.subfolder + "test.txt", std::ofstream::out );
    for ( double t_t = t_start_new; t_t <= t_end; t_t += t_step_initial ) {
        int n_c = (int)( std::min( n_c_max, (int)( t_t / s.parameters.t_step ) ) );
        int t_index = (int)( t_t / s.parameters.t_step );

        double cur_min = 1;
        auto c1 = clock();
        // Path-Integral iteration
        // Again, precalculate the iterated M for every timestep, because we only want N^2 instead of N^2*N^2 iteration calls.
        std::vector<Sparse> iterates;
        iterates.reserve( tensor_dim * tensor_dim );
        for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
            for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                Sparse projector = Sparse( tensor_dim, tensor_dim );
                projector.coeffRef( i_n, j_n ) = 1;
                Sparse M = iterate( projector, s, t_t, output );
                iterates.emplace_back( M );
            }
        }

        // Main Iteration loop
        Log::L3( "current non zeros in adm: {:.5f}\r", adms.getFillRatio() * 100.0 );
        test << t_t << "\t" << adms.getFillRatio() * 100.0 << "\t" << adms.getSizeOfCache() << "\t";

        for ( int k = 0; k < adms.mat().outerSize(); ++k ) {
            for ( Sparse::InnerIterator it( adms.mat(), k ); it; ++it ) {
                auto indices_i = adms.indexToIndices( it.row() );
                auto indices_j = adms.indexToIndices( it.col() );
                int i_n_m1 = indices_i.at( 0 );
                int i_n_m2 = indices_i.at( 1 );
                int i_n_m3 = indices_i.at( 2 );
                int i_n_m4 = indices_i.at( 3 );
                int j_n_m1 = indices_j.at( 0 );
                int j_n_m2 = indices_j.at( 1 );
                int j_n_m3 = indices_j.at( 2 );
                int j_n_m4 = indices_j.at( 3 );
                for ( int i_n = 0; i_n < tensor_dim; i_n++ ) {
                    for ( int j_n = 0; j_n < tensor_dim; j_n++ ) {
                        Scalar phonon_s = 0.0;
                        for ( int l = 0; l <= n_c_max; l++ ) {
                            int i_nd = ( l == 0 ? i_n : indices_i.at( l - 1 ) );
                            int j_nd = ( l == 0 ? j_n : indices_j.at( l - 1 ) );
                            phonon_s += s.dgl_phonon_S_function( l, i_n, j_n, i_nd, j_nd );
                        }
                        //if ( std::real( phonon_s ) != 0 ) std::cout << "t = " << t_t << " -> ps = " << phonon_s << " -> " << std::exp( phonon_s ) << std::endl;
                        Scalar val = iterates.at( j_n_m1 + tensor_dim * i_n_m1 ).coeff( i_n, j_n ) * it.value() * std::exp( phonon_s );
                        double abs = std::abs( val );
                        if ( abs <= 1E-10 ) continue;
                        cur_min = cur_min < abs ? cur_min : abs;
                        adms.addTriplet( i_n, i_n_m1, i_n_m2, i_n_m3, j_n, j_n_m1, j_n_m2, j_n_m3, val );
                    }
                }
            }
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

        Log::L3( "Iteration t = {}, time taken for propagation: {}s\n", t_t, ( clock() - c1 ) / CLOCKS_PER_SEC );
        test << cur_min << std::endl;
        // Copy new ADM into old variable for consistency
        adms.setFromTripletList();

        // Calculate the reduced density matrix by tracing over all past times for each entry
        std::vector<Eigen::Triplet<Scalar>> triplets;
        // Calculate rho by tracing over the n_c last steps
        c1 = clock();

        Dense cache = Dense::Zero( tensor_dim, tensor_dim );
        //#pragma omp parallel for num_threads( s.parameters.numerics_maximum_threads )
        for ( int k = 0; k < adms.mat().outerSize(); ++k ) {
            for ( Sparse::InnerIterator it( adms.mat(), k ); it; ++it ) {
                auto indices_i = adms.indexToIndices( it.row() );
                auto indices_j = adms.indexToIndices( it.col() );
                int i_n = indices_i.at( 0 );
                int j_n = indices_j.at( 0 );
                cache( i_n, j_n ) += it.value();
            }
        }
        rho = cache.sparseView();

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

        Log::L3( "Iteration t = {}, time taken for adm trace: {}s\n", t_t, ( clock() - c1 ) / CLOCKS_PER_SEC );

        // Save Rho
        saveState( rho, t_t, output );
        // Progress and time output
        rkTimer.iterate();
        if ( do_output ) {
            Timers::outputProgress( s.parameters.output_handlerstrings, rkTimer, progressbar, s.parameters.iterations_t_max, progressbar_name );
        }
    }
    test.close();
    Log::L3( "Final Size of ADMs Tensor: {} bytes / {}\% filled\n", adms.getSizeOfCache(), adms.getFillRatio() * 100.0 );
    Log::L3( "Done!\n" );
    return true;
}
