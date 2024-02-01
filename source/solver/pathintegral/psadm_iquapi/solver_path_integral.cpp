#include "solver/solver_ode.h"

// TODO: statt sparse -> dense am besten dense -> nonzeros zu indices appenden, ab punkt dann einfach nur noch über nonzero indices iterieren. dann kann man von anfang an multithreaded machen, irgendwann gucken was is ungleich null,
// dann nur noch über diese indices summieren.

QDACC::Type::MatrixMain _reduce_adm_tensor( QDACC::Numerics::Tensor &tensor ) {
    const auto tensor_dim = tensor.primary_dimensions();
    QDACC::Type::MatrixMain ret( tensor_dim, tensor_dim );
    ret.setZero();
    for ( QDACC::Numerics::Tensor::IndexVector &index : tensor.get_indices() ) {
        const int &i_n = index[0];
        const int &j_n = index[1];
        ret.coeffRef( i_n, j_n ) += tensor( index );
    }
    return ret;
}

void _fill_tensor( QDACC::Numerics::Tensor &tensor, const QDACC::Type::MatrixMain &rho0 ) {
    const auto tensor_dim = tensor.primary_dimensions();
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( QDACC::Math::abs2( rho0.coeff( i, j ) ) == 0 )
                continue;
            QDACC::Numerics::Tensor::IndexVector index( tensor.index_size(), 0 );
            index[0] = i;
            index[1] = j;
            tensor( index ) += rho0.coeff( i, j );
        }
    }
}

// for ( int i = 0; i < tensor_dim; i++ ) {
//         for ( int j = 0; j < tensor_dim; j++ ) {
//             if ( QDACC::Math::abs2( rho0.coeff( i, j ) ) == 0 )
//                 continue;
//             QDACC::Numerics::Tensor::IndexVector index( index_vector_length, 0 );
//             index[0] = i;
//             index[1] = j;
//             adm_tensor( index ) += rho0.coeff( i, j );
//             Log::L2( "[PathIntegral] Added element ({},{}) = {} to adm tensor\n", i, j, rho0.coeff( i, j ) );
//         }
//     }

QDACC::Numerics::Tensor _generate_initial_tensor( QDACC::System &s, const QDACC::Type::MatrixMain &rho0 ) {
    // Gather unique Coupling indices, where the number of unique indices equals N(M)x for x >= 1
    std::set<int> set_different_dimensions;
    for ( int i = 0; i < s.operatorMatrices.phonon_hilbert_index_to_group_index.size(); i++ ) {
        set_different_dimensions.insert( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i] : i );
    }
    const int different_dimensions = set_different_dimensions.size();

    // Alias tensor Dimension for readability
    const QDACC::Numerics::Tensor::Index tensor_dim = rho0.rows();

    // Build Tensor Dimensions. Structure is N0,M0,N1,M1,...,NN,MM
    QDACC::Numerics::Tensor::IndexVector pathint_tensor_dimensions = { tensor_dim, tensor_dim }; // TODO Remove shadow?
    for ( int i = 1; i < s.parameters.p_phonon_nc; i++ ) {
        pathint_tensor_dimensions.emplace_back( different_dimensions );
        pathint_tensor_dimensions.emplace_back( different_dimensions );
    }
    const auto index_vector_length = pathint_tensor_dimensions.size();

    // Output Tensor Dimensions
    std::stringstream adm_tensor_dimensions;
    std::ranges::copy( std::begin( pathint_tensor_dimensions ), std::end( pathint_tensor_dimensions ), std::ostream_iterator<int>( adm_tensor_dimensions, ", " ) );
    Log::L2( "[PathIntegral] ADM Tensor dimensions: [{}]\n", adm_tensor_dimensions.str() );
    Log::L2( "[PathIntegral] Using a maximum of #{} cpu cores for the ADM propagation.\n", s.parameters.numerics_maximum_secondary_threads );

    // Tensor Object
    QDACC::Numerics::Tensor adm_tensor( pathint_tensor_dimensions );
    Log::L2( "[PathIntegral] Fixed Tensor size will be {} MB.\n", adm_tensor.nonZeros(), adm_tensor.size() );

    // Fill initial Tensor
    _fill_tensor( adm_tensor, rho0 );

    return adm_tensor;
}

QDACC::Numerics::Tensor QDACC::Numerics::ODESolver::iterate_path_integral( System &s, QDACC::Numerics::Tensor &adm_tensor, std::vector<std::vector<MatrixMain>> &propagator, const int max_index ) {
    // Create New Tensor
    Tensor adm_tensor_next( adm_tensor.nonZeros() );
    const auto different_dimensions = adm_tensor.secondary_dimensions();
    // Iterate Tensor
#pragma omp parallel for num_threads( s.parameters.numerics_maximum_secondary_threads ) schedule( static )
    for ( Tensor::IndexVector &index : adm_tensor.get_indices() ) {
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
                Scalar propagator_value = propagator[i_n_m1][j_n_m1].coeff( i_n, j_n );
                // If this propagation is not allowed, continue
                if ( QDACC::Math::abs2( propagator_value ) == 0.0 )
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

    // Return New ADM Tensor
    return adm_tensor_next;
}

bool QDACC::Numerics::ODESolver::calculate_path_integral( MatrixMain &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDACC::SaveState> &output, bool do_output ) {
    Log::L2( "[PathIntegral] Setting up Path-Integral Solver...\n" );
    output.reserve( s.parameters.iterations_t_max + 1 );
    const int total_progressbar_iterations = std::floor( ( s.parameters.t_end - s.parameters.t_start ) / s.parameters.t_step_pathint );

    auto adm_tensor = _generate_initial_tensor( s, rho0 );
    const auto tensor_dim = adm_tensor.primary_dimensions();
    adm_tensor.save_to_file( 0 );

    // Only save the tensors to files if we calculate correlations AND if we do not use the QRT
    const bool save_tensor = not s.parameters.input_correlation.empty() and not s.parameters.numerics_pathintegral_use_qrt;

    // First step is just rho0
    MatrixMain rho = rho0;
    saveState( rho0, t_start, output );

    // Iterate Path integral for further time steps
    for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Calculate the maximum Vector index. For times < cutofftime, the number of active elements is reduced
        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );

        // Calculate Propagators for current time
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );
        // Iterate the tensor
        // adm_tensor = iterate_path_integral_gpu( s, adm_tensor, propagator, max_index );
        adm_tensor = iterate_path_integral( s, adm_tensor, propagator, max_index );
        rho = _reduce_adm_tensor( adm_tensor );
        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );

        // Save ADM Tensor to file
        if ( save_tensor ) {
            const int tensor_int_time = std::floor( 1.001 * t_t / s.parameters.t_step_pathint + 1. );
            adm_tensor.save_to_file( tensor_int_time );
        }

        // Increase Tmax if required
        if ( s.parameters.numerics_calculate_till_converged and t_t + s.parameters.t_step_pathint >= t_end and t_t < s.parameters.numerics_hard_t_max and std::real( output.back().mat.coeff( s.parameters.numerics_groundstate, s.parameters.numerics_groundstate ) ) < 0.999 ) {
            t_end += 10.0 * s.parameters.t_step_pathint;
            Log::L3( "[PathIntegral] Adjusted Calculation end to {}\n", t_end );
        }

        // If prune == auto: prune tensor if non_zero elements dont change below difference of 10 elements and NC > NC_max
        // if prune = manual: Xps:prune:[threshold];Yps:extend;Zps:prune:[threshold] where prune calls prune, extend cleas pruned vector.
        // if prune is set to auto, automatically determine prunes and extends using the pulse, chirp
        // FOR NOW: disable it
        // if (t_t / s.parameters.t_step_pathint > 6. and not adm_tensor.is_pruned())
        //    adm_tensor.make_indices_sparse();

        // Progress and time output
        rkTimer.iterate();
        if ( do_output )
            Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
    }

    if ( s.parameters.numerics_calculate_till_converged ) {
        s.parameters.numerics_calculate_till_converged = false;
        s.parameters.t_end = t_end;
        Log::L1( "[PathIntegral] Adjusted t_end to {}.\n", s.parameters.t_end );
    }

    return true;
}

bool QDACC::Numerics::ODESolver::calculate_path_integral_correlation( MatrixMain &rho0, double t_start, double t_end, Timer &rkTimer, System &s, std::vector<QDACC::SaveState> &output, const MatrixMain &op_l, const MatrixMain &op_i, int adm_multithreading_cores ) {
    output.reserve( s.parameters.iterations_t_max + 1 );
    const int tensor_int_time = std::floor( 1.001 * t_start / s.parameters.t_step_pathint ); // This assumes t_start for the main direction is zero!
    auto adm_correlation = Tensor( Tensor::max_size() );
    if ( not s.parameters.numerics_pathintegral_use_qrt )
        adm_correlation.load_from_file( tensor_int_time );
    else
        _fill_tensor( adm_correlation, rho0 );
    const auto tensor_dim = adm_correlation.primary_dimensions();

    saveState( rho0, t_start, output );

    // Modify Initial Propagator
    if ( not s.parameters.numerics_pathintegral_use_qrt ) {
        std::vector<std::vector<MatrixMain>> initial_propagator = calculate_propagator_vector( s, tensor_dim, t_start, s.parameters.numerics_subiterator_stepsize, output );
        for ( auto &vec : initial_propagator )
            for ( auto &el : vec )
                el = s.dgl_timetrafo( op_l, t_start ) * el * s.dgl_timetrafo( op_i, t_start );
        // First step by hand
        adm_correlation = iterate_path_integral( s, adm_correlation, initial_propagator, 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_start / s.parameters.t_step_pathint ) ) );
        saveState( _reduce_adm_tensor( adm_correlation ), t_start + s.parameters.t_step_pathint, output );
        t_start += s.parameters.t_step_pathint;
    }

    // Iterate Path integral for further time steps
    for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Calculate the maximum Vector index. For times < cutofftime, the number of active elements is reduced
        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );

        // Calculate Propagators for current time
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );

        // Iterate the tensor
        adm_correlation = iterate_path_integral( s, adm_correlation, propagator, max_index );

        const auto rho = _reduce_adm_tensor( adm_correlation );

        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );

        rkTimer.iterate();
    }

    return true;
}