#include "solver/solver_ode.h"

// TODO: statt sparse -> dense am besten dense -> nonzeros zu indices appenden, ab punkt dann einfach nur noch über nonzero indices iterieren. dann kann man von anfang an multithreaded machen, irgendwann gucken was is ungleich null,
// dann nur noch über diese indices summieren.

// TODO
// Nur eine Funktion für PI bzw weitestgehend aufteilen.
// Unterfunktion für G1 und G2 die die selbe syntax wie die normalen funktionen haben.
// PathInt wie RKXY nutzen. Dafür Tensor in T-Richtung auf Disk oder im Ram cachen wenn G1 oder G2 gerechnet werden.
// Dann rechnen andere funktionen G1 und G2 aus. Dafür dann dort statt nur "calculate_runge_kutta" fallunterscheidung mit RK und PI.

QDLC::Type::Sparse _reduce_adm_tensor( QDLC::Numerics::Tensor &tensor ) {
    const auto tensor_dim = tensor.primary_dimensions();
    QDLC::Type::Dense ret = QDLC::Type::Dense::Zero( tensor_dim, tensor_dim );
    for ( QDLC::Numerics::Tensor::IndexVector &index : tensor.get_indices() ) {
        const int &i_n = index[0];
        const int &j_n = index[1];
        ret( i_n, j_n ) += tensor( index );
    }
    return ret.sparseView();
}

QDLC::Numerics::Tensor _generate_initial_tensor( QDLC::System &s, const QDLC::Type::Sparse &rho0 ) {
    // Gather unique Coupling indices, where the number of unique indices equals N(M)x for x >= 1
    std::set<int> set_different_dimensions;
    for ( int i = 0; i < s.operatorMatrices.phonon_hilbert_index_to_group_index.size(); i++ ) {
        set_different_dimensions.insert( s.parameters.numerics_pathint_partially_summed ? s.operatorMatrices.phonon_hilbert_index_to_group_index[i] : i );
    }
    const int different_dimensions = set_different_dimensions.size();

    // Alias tensor Dimension for readability
    const QDLC::Numerics::Tensor::Index tensor_dim = rho0.rows();

    // Build Tensor Dimensions. Structure is N0,M0,N1,M1,...,NN,MM
    QDLC::Numerics::Tensor::IndexVector pathint_tensor_dimensions = { tensor_dim, tensor_dim }; // TODO Remove shadow?
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
    QDLC::Numerics::Tensor adm_tensor( pathint_tensor_dimensions );
    Log::L2( "[PathIntegral] Fixed Tensor size will be {} MB.\n", adm_tensor.nonZeros(), adm_tensor.size() );

    // Fill initial Tensor
    for ( int i = 0; i < tensor_dim; i++ ) {
        for ( int j = 0; j < tensor_dim; j++ ) {
            if ( QDLC::Math::abs2( rho0.coeff( i, j ) ) == 0 )
                continue;
            QDLC::Numerics::Tensor::IndexVector index( index_vector_length, 0 );
            index[0] = i;
            index[1] = j;
            adm_tensor( index ) += rho0.coeff( i, j );
            Log::L2( "[PathIntegral] Added element ({},{}) = {} to adm tensor\n", i, j, rho0.coeff( i, j ) );
        }
    }

    return adm_tensor;
}

QDLC::Numerics::Tensor QDLC::Numerics::ODESolver::iterate_path_integral( System &s, QDLC::Numerics::Tensor &adm_tensor, std::vector<std::vector<Sparse>> &propagator, const int max_index ) {
    // Create New Tensor
    Tensor adm_tensor_next( adm_tensor.nonZeros() );
    const auto different_dimensions = adm_tensor.secondary_dimensions();

    // Iterate Tensor
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

bool QDLC::Numerics::ODESolver::calculate_path_integral( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, ProgressBar &progressbar, std::string progressbar_name, System &s, std::vector<QDLC::SaveState> &output, bool do_output ) {
    Log::L2( "[PathIntegral] Setting up Path-Integral Solver...\n" );
    output.reserve( s.parameters.iterations_t_max + 1 );
    const int total_progressbar_iterations = std::floor( ( s.parameters.t_end - s.parameters.t_start ) / s.parameters.t_step_pathint );

    auto adm_tensor = _generate_initial_tensor( s, rho0 );
    const auto tensor_dim = adm_tensor.primary_dimensions();

    // First step is just rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );

    // Iterate Path integral for further time steps
    for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Calculate the maximum Vector index. For times < cutofftime, the number of active elements is reduced
        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );
        if ( max_index < s.parameters.p_phonon_nc ) {
            Log::L3( "[PathIntegral] max_index = {}, nc = {}\n", max_index, s.parameters.p_phonon_nc );
        }

        // Calculate Propagators for current time
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );
        // Iterate the tensor
        adm_tensor = iterate_path_integral( s, adm_tensor, propagator, max_index );
        rho = _reduce_adm_tensor( adm_tensor );
        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );

        // Save ADM Tensor to file
        const int tensor_int_time = t_t / s.parameters.t_step_pathint;
        adm_tensor.save_to_file( tensor_int_time );

        // Progress and time output
        rkTimer.iterate();
        if ( do_output )
            Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
    }
    return true;
}

bool QDLC::Numerics::ODESolver::calculate_path_integral_correlation( Sparse &rho0, double t_start, double t_end, Timer &rkTimer, System &s, std::vector<QDLC::SaveState> &output, const Sparse &op_l, const Sparse &op_i, int adm_multithreading_cores ) {
    output.reserve( s.parameters.iterations_t_max + 1 );
    
    const int tensor_int_time = t_start / s.parameters.t_step_pathint; // This assumes t_start for the main direction is zero!
    auto adm_correlation = Tensor( Tensor::max_size() );
    adm_correlation.load_from_file( tensor_int_time );
    const auto tensor_dim = adm_correlation.primary_dimensions();

    // Modify Initial Propagator
    std::vector<std::vector<Sparse>> initial_propagator = calculate_propagator_vector( s, tensor_dim, t_start, s.parameters.numerics_subiterator_stepsize, output );
    for ( auto &vec : initial_propagator )
        for ( auto &el : vec )
            el = s.dgl_timetrafo( op_l * el * op_i, t_start );

    // First step is just rho0. TODO: iterate once by hand to include new prop, then save that as rho0
    Sparse rho = rho0;
    saveState( rho0, t_start, output );

    // Iterate Path integral for further time steps
    for ( double t_t = t_start; t_t < t_end; t_t += s.parameters.t_step_pathint ) {
        // Calculate the maximum Vector index. For times < cutofftime, the number of active elements is reduced
        int max_index = 2 + std::min<int>( s.parameters.p_phonon_nc - 2, std::floor( 1.001 * t_t / s.parameters.t_step_pathint ) );
        if ( max_index < s.parameters.p_phonon_nc ) {
            Log::L3( "[PathIntegral] max_index = {}, nc = {}\n", max_index, s.parameters.p_phonon_nc );
        }

        // Calculate Propagators for current time
        auto &propagator = calculate_propagator_vector( s, tensor_dim, t_t, s.parameters.numerics_subiterator_stepsize, output );
        // Iterate the tensor
        if ( t_t == t_start )
            adm_correlation = iterate_path_integral( s, adm_correlation, initial_propagator, max_index );
        else
            adm_correlation = iterate_path_integral( s, adm_correlation, propagator, max_index );
        rho = _reduce_adm_tensor( adm_correlation );
        // Save Rho
        saveState( rho, t_t + s.parameters.t_step_pathint, output );

        // Progress and time output
        rkTimer.iterate();
        //Timers::outputProgress( rkTimer, progressbar, rkTimer.getTotalIterationNumber(), total_progressbar_iterations, progressbar_name );
    }
    return true;
}