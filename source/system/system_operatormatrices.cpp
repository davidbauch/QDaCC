#include "system/operatormatrices.h"
#include "misc/helperfunctions_matrix.h"
//#include "system/operatormatrices_text.h"

OperatorMatrices::OperatorMatrices( Parameters &p ) {
    Timer &timer_operatormatrices = Timers::create( "Operator Matrices", true, false );
    timer_operatormatrices.start();
    Log::L2( "[System-OperatorMatrices] Generating operator matrices...\n" );
    if ( !generate_operators( p ) ) {
        Log::L2( "[System-OperatorMatrices] Generating operator matrices failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    timer_operatormatrices.end();
    Log::L2( "[System-OperatorMatrices] Generating operator matrices was successful. Elapsed time is {}ms\n", timer_operatormatrices.getWallTime( Timers::MILLISECONDS ) );
}

bool OperatorMatrices::generate_operators( Parameters &p ) {
    output_format = Eigen::IOFormat( 4, 0, ", ", "\n", "[", "]" );
    // Zeroing Hamiltons (redundant at this point)
    Log::L2( "[System-OperatorMatrices] Creating operator matrices, dimension = {}, creating base matrices...\n", p.maxStates );

    // Generate Electronic and Photonic Base States and Self-Hilbert Matrices
    // Electronic
    Log::L2( "[System-OperatorMatrices] Creating Self-Hilbert Electronic states...\n" );
    std::vector<std::string> base_electronic;
    base_selfhilbert.emplace_back( Dense::Identity( p.input_electronic.size(), p.input_electronic.size() ) );
    int index = 0;
    int maximum_electronic_levels = p.input_electronic.size();
    for ( auto &mat : p.input_electronic ) {
        base_electronic.emplace_back( mat.first );
        auto &state = el_states[mat.first];
        state.ket = Dense::Zero( maximum_electronic_levels, 1 );
        state.ket( index++ ) = 1;
        state.bra = state.ket.transpose();
        state.self_hilbert = state.ket * state.bra;
        state.base = 0;
        state.energy = mat.second.numerical["Energy"];
    }

    // Photonic
    Log::L2( "[System-OperatorMatrices] Creating Self-Hilbert Photonic states...\n" );
    std::vector<std::vector<std::string>> base_photonic( p.input_photonic.size() );
    int curcav = 0;
    for ( auto &cav : p.input_photonic ) {
        int max_photons = cav.second.numerical["MaxPhotons"];
        base_selfhilbert.emplace_back( Dense::Identity( max_photons + 1, max_photons + 1 ) );
        for ( int n = 0; n <= max_photons; n++ )
            base_photonic[curcav].emplace_back( std::to_string( n ) + cav.first );
        auto &state = ph_states[cav.first];
        state.self_hilbert = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_STATE, max_photons );
        state.base = ++curcav;
        state.energy = cav.second.numerical["Energy"];
        // std::cout << "Matrix " << cav.first << "\n"
        //           << state.self_hilbert << std::endl;
    }

    // Tensor all matrices into the total Hilbert Space.
    // String Base Tensors
    auto subb = base_photonic;
    subb.insert( subb.begin(), base_electronic );
    base = QDLC::Matrix::tensor( subb );

    // Create Universal Identity Matrix
    identity = Sparse( base.size(), base.size() );

    Log::L2( "[System-OperatorMatrices] Evaluating Kronecker Delta for the electronic states...\n" );
    // Sparse total Hilbert tensors
    for ( auto &bel : el_states ) {
        auto current = base_selfhilbert;
        current.front() = bel.second.self_hilbert; // Electronic states have basis 0 for now. Maybe change when TODO: add tensoring of multiple electronic bases
        bel.second.hilbert = QDLC::Matrix::tensor( current ).sparseView();
        bel.second.projector = QDLC::Matrix::sparse_projector( bel.second.hilbert );
        // std::cout << "Mat = " << bel.first << "  " << bel.second.hilbert.rows() << "\n"
        //           << Dense( bel.second.hilbert ) << std::endl;
    }
    Log::L2( "[System-OperatorMatrices] Evaluating Kronecker Delta for the photonic states...\n" );
    for ( auto &phot : ph_states ) {
        auto current = base_selfhilbert;
        current[phot.second.base] = phot.second.self_hilbert;
        phot.second.hilbert = QDLC::Matrix::tensor( current ).sparseView();
        phot.second.projector = QDLC::Matrix::sparse_projector( phot.second.hilbert );
        // std::cout << "Mat = " << phot.first << "  " << phot.second.hilbert.rows() << "\n"
        //           << Dense( phot.second.hilbert ) << std::endl;
    }
    Log::L2( "[System-OperatorMatrices] Creating Electronic Transition Matrices...\n" );
    // Generate Transition Matrices
    for ( auto &bel : el_states ) {
        for ( auto &trans_to : p.input_electronic[bel.first].string_v["CoupledTo"] ) {
            if ( trans_to.front() == '-' ) continue;
            std::string transition = bel.first + trans_to;
            el_transitions[transition].self_hilbert = bel.second.ket * el_states[trans_to].bra;
            el_transitions[transition].direction = -1; // DOWN
            el_transitions[transition].energy = std::abs( p.input_electronic[trans_to].numerical["Energy"] - p.input_electronic[bel.first].numerical["Energy"] );
            std::string transition_transposed = trans_to + bel.first;
            el_transitions[transition_transposed].self_hilbert = el_states[trans_to].ket * bel.second.bra;
            el_transitions[transition_transposed].direction = 1; // UP
            el_transitions[transition_transposed].energy = std::abs( p.input_electronic[trans_to].numerical["Energy"] - p.input_electronic[bel.first].numerical["Energy"] );
        }
    }
    Log::L2( "[System-OperatorMatrices] Creating Photonic Transition Matrices...\n" );
    curcav = 1;
    for ( auto &cav : p.input_photonic ) {
        int max_photons = cav.second.numerical["MaxPhotons"];
        auto &state_b = ph_transitions[cav.first + "b"];
        state_b.self_hilbert = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_ANNIHILATE, max_photons );
        state_b.base = curcav;
        state_b.direction = -1;
        state_b.energy = cav.second.numerical["Energy"];
        auto &state_bd = ph_transitions[cav.first + "bd"];
        state_bd.self_hilbert = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_CREATE, max_photons );
        state_bd.base = curcav++;
        state_bd.direction = 1;
        state_bd.energy = cav.second.numerical["Energy"];
    }
    Log::L2( "[System-OperatorMatrices] Creating Electronic Projector Matrices...\n" );
    // Tensor remaining transition
    for ( auto &bel : el_transitions ) {
        auto current = base_selfhilbert;
        current.front() = bel.second.self_hilbert;
        bel.second.hilbert = QDLC::Matrix::tensor( current ).sparseView();
        Log::L3( "[System-OperatorMatrices] Electronic Transition Matrix for Transition {} in self-Hilbert space:\n{}\n", bel.first, Dense( bel.second.self_hilbert ).format( output_format ) );
        Log::L3( "[System-OperatorMatrices] Electronic Transition Matrix for Transition {} in total-Hilbert space:\n{}\n", bel.first, Dense( bel.second.hilbert ).format( output_format ) );
        bel.second.projector = QDLC::Matrix::sparse_projector( bel.second.hilbert );
    }
    Log::L2( "[System-OperatorMatrices] Creating Photonic Projector Matrices...\n" );
    for ( auto &phot : ph_transitions ) {
        auto current = base_selfhilbert;
        current[phot.second.base] = phot.second.self_hilbert;
        phot.second.hilbert = QDLC::Matrix::tensor( current ).sparseView();
        Log::L3( "[System-OperatorMatrices] Cavity Transition Matrix for Transition {} in self-Hilbert space:\n{}\n", phot.first, Dense( phot.second.self_hilbert ).format( output_format ) );
        Log::L3( "[System-OperatorMatrices] Cavity Transition Matrix for Transition {} in total-Hilbert space:\n{}\n", phot.first, Dense( phot.second.hilbert ).format( output_format ) );
        phot.second.projector = QDLC::Matrix::sparse_projector( phot.second.hilbert );
    }

    Log::L2( "[System-OperatorMatrices] Creating Hilbert Space Index Matrices...\n" );
    // Create total Hilbert space indices. This routine assumes the individual bases are added in ascending order, e.g. 0,1,2,3...:
    for ( int i = 0; i < base_selfhilbert.size(); i++ ) {
        auto current = base_selfhilbert;
        int index = 1;
        for ( int k = 0; k < current[i].rows(); k++ )
            for ( int j = 0; j < current[i].cols(); j++ )
                current[i]( k, j ) = ( k + 1 ) + 1.0i * ( j + 1 );

        base_hilbert_index.emplace_back( QDLC::Matrix::tensor( current ) );
        // Log::L2( "Tensor for base i = {}:\nCurrent:\n{}\n\n{}\n\n", i, current[i], base_hilbert_index.back() );
    }

    Log::L2( "[System-OperatorMatrices] Creating Pulse Cache Matrices...\n" );
    // Generate prechaced pulse and chirp matrices. 2 matrices are generated per pulse for Omega and Omega^*
    // Generate transition matrices from the initial state matrices to allow for addidional transitions and the cavity modes to be included.
    // Create a cache Vector for the cavity pulses to allow for the phonon scaling to still apply, meaning the cavity pulse matrices will be added after the phonon scaling was applied:
    std::vector<QDLC::Type::Sparse> pulse_mat_cavity_cache;
    for ( auto &pulse : p.input_pulse ) {
        Sparse pulsemat = Sparse( base.size(), base.size() );
        Sparse pulsemat_star = Sparse( base.size(), base.size() );
        pulse_mat_cavity_cache.emplace_back( pulsemat );
        pulse_mat_cavity_cache.emplace_back( pulsemat_star );
        for ( int i = 0; i < pulse.second.string_v["CoupledTo"].size(); i++ ) {
            std::string transition = pulse.second.string_v["CoupledTo"][i];
            std::string transition_transposed = transition;
            std::reverse( transition_transposed.begin(), transition_transposed.end() );
            if ( el_transitions.count( transition ) > 0 ) {
                pulsemat += el_transitions[transition].hilbert;
                pulsemat_star += el_transitions[transition_transposed].hilbert;
            } else if ( std::isupper( transition.front() ) ) {
                Log::L2( "[System-OperatorMatrices] Electronic Pulse transition {} is not in the list of allowed electronic transitions, recreating transition matrices...\n", transition );
                auto ket1 = el_states[transition.substr( 0, 1 )].ket;
                auto bra1 = el_states[transition.substr( 0, 1 )].bra;
                auto ket2 = el_states[transition.substr( 1, 1 )].ket;
                auto bra2 = el_states[transition.substr( 1, 1 )].bra;
                auto current = base_selfhilbert;
                current.front() = ket1 * bra2;
                Sparse transition_hilbert = QDLC::Matrix::tensor( current ).sparseView();
                {
                    auto &state = extra_transitions[transition];
                    state.ket = ket1;
                    state.bra = state.ket.transpose();
                    state.self_hilbert = state.ket * state.bra;
                    state.base = 0;
                    state.hilbert = transition_hilbert;
                }

                current.front() = ket2 * bra1;
                Sparse transition_transposed_hilbert = transition_hilbert.adjoint();
                {
                    auto &state = extra_transitions[transition_transposed];
                    state.ket = ket2;
                    state.bra = state.ket.transpose();
                    state.self_hilbert = state.ket * state.bra;
                    state.base = 0;
                    state.hilbert = transition_transposed_hilbert;
                }

                pulsemat += transition_hilbert;
                pulsemat_star += transition_transposed_hilbert;
            } else {
                Log::L2( "[System-OperatorMatrices] Pulse transition {} is cavity...\n", transition );
                // pulsemat += ph_transitions[transition + "b"].hilbert;
                // pulsemat_star += ph_transitions[transition + "bd"].hilbert;
                pulse_mat_cavity_cache[pulse_mat_cavity_cache.size() - 2] += ph_transitions[transition + "b"].hilbert;
                pulse_mat_cavity_cache[pulse_mat_cavity_cache.size() - 1] += ph_transitions[transition + "bd"].hilbert;
            }
        }
        pulse_mat.emplace_back( pulsemat );
        pulse_mat.emplace_back( pulsemat_star );
        Log::L3( "[System-OperatorMatrices] Added Pulse Matrix for Pulse {} in total-Hilbert space (normal+transposed):\n{}\n", pulse.first, Dense( pulsemat + pulse_mat_cavity_cache[pulse_mat_cavity_cache.size() - 2] + pulsemat_star + pulse_mat_cavity_cache[pulse_mat_cavity_cache.size() - 1] ).format( output_format ) );
    }

    // TODO: chirp cavity
    Log::L2( "[System-OperatorMatrices] Creating Chirp Cache Matrices...\n" );
    // 1 matrix is generated per chirp
    for ( auto &chirp : p.input_chirp ) {
        Sparse chirpmat = Sparse( base.size(), base.size() );
        Sparse chirpmat_star = Sparse( base.size(), base.size() );
        for ( int i = 0; i < chirp.second.string_v["CoupledTo"].size(); i++ ) {
            std::string state = chirp.second.string_v["CoupledTo"][i];
            chirpmat += chirp.second.numerical_v["AmpFactor"][i] * el_states[state].hilbert; // TODO: chirp cavity similar to pulse
        }
        chirp_mat.emplace_back( chirpmat );
    }

    Log::L2( "[System-OperatorMatrices] Creating Base Index Maps...\n" );
    // Create Index Map
    index = 0;
    std::stringstream ss;
    for ( auto &b : base ) {
        b = "|" + b + ">";
        base_index_map[b] = index++;
        ss << b << " ";
    }
    Log::L1( "System Base ({}): {}\n", base.size(), ss.str() );
    p.maxStates = base.size();

    H_0 = Sparse( p.maxStates, p.maxStates );
    H_I = Sparse( p.maxStates, p.maxStates );
    H_used = Sparse( p.maxStates, p.maxStates );
    rho = Sparse( p.maxStates, p.maxStates );

    // Analytical Time Trafo Matrix
    // TODO: Puls mit in trafo? Trafohamilton dann H_el + H_phot + H_pulse
    Log::L2( "[System-OperatorMatrices] Creating Analytical Time Transformation Cache Matrix...\n" );
    timetrafo_cachematrix = Dense::Zero( base.size(), base.size() );
    for ( auto &iss : base_index_map ) {
        std::string is = iss.first;
        is.pop_back();
        is.erase( is.begin() );
        auto i = QDLC::String::splitline( is, '|' );
        for ( auto &jss : base_index_map ) {
            Scalar val = 0;
            std::string js = jss.first;
            js.pop_back();
            js.erase( js.begin() );
            auto j = QDLC::String::splitline( js, '|' );
            val += p.input_electronic[i.front()].numerical["Energy"] - p.input_electronic[j.front()].numerical["Energy"];
            for ( int a = 1; a < i.size(); a++ ) {
                std::string ai = i[a];
                std::string aj = j[a];
                std::string mode = "";
                mode.append( 1, i[a].back() );
                ai.pop_back();
                aj.pop_back();
                double ni = std::stod( ai.c_str() );
                double nj = std::stod( aj.c_str() );
                // std::cout << iss.first << ", " << jss.first << " --> " << ni << " " << nj << ", modes = " << mode << std::endl;
                val += p.input_photonic[mode].numerical["Energy"] * ( ni - nj );
            }
            timetrafo_cachematrix( iss.second, jss.second ) = 1.0i * val;
        }
    }

    Log::L2( "[System-OperatorMatrices] Creating Polaron Cache Matrices...\n" );
    // Precalculate Polaron Matrices.
    polaron_factors.emplace_back( Sparse( base.size(), base.size() ) );
    // Transition is always |0><1| (annihilator), hence reversed is |1><0| (creator)
    for ( auto &[mode, param] : p.input_photonic ) {
        int i = 0;
        for ( auto transition : param.string_v["CoupledTo"] ) {
            Log::L2( "[System-PME] Adding Polaron Cavity Transition |{}><{}|b_{}\n", transition.front(), transition.back(), mode );
            std::reverse( transition.begin(), transition.end() );
            polaron_factors[0] += el_transitions[transition].hilbert * p.p_omega_coupling * param.numerical_v["CouplingScaling"][i] * ph_transitions[mode + "b"].hilbert;
            i++;
        }
    }
    for ( auto &[mode, param] : p.input_pulse ) {
        Sparse temp = Sparse( base.size(), base.size() );
        for ( auto transition : param.string_v["CoupledTo"] ) {
            polaron_pulse_factors_explicit_time.emplace_back( Sparse( base.size(), base.size() ) );
            std::string transition_transposed = transition;
            std::reverse( transition_transposed.begin(), transition_transposed.end() );
            Log::L2( "[System-PME] Adding Polaron Pulse Transition |{}><{}|Omega_{}\n", transition.front(), transition.back(), mode );
            if ( el_transitions.count( transition_transposed ) > 0 ) {
                temp += el_transitions[transition_transposed].hilbert;
                polaron_pulse_factors_explicit_time.back() += el_transitions[transition_transposed].projector;
            } else if ( std::isupper( transition.front() ) ) {
                temp += extra_transitions[transition_transposed].hilbert;
                polaron_pulse_factors_explicit_time.back() += QDLC::Matrix::sparse_projector( extra_transitions[transition_transposed].hilbert );
            } else {
                Log::L2( "[System-OperatorMatrices] Pulse transition {} is cavity...\n", transition );
                temp += ph_transitions[transition + "bd"].hilbert;
                polaron_pulse_factors_explicit_time.back() += ph_transitions[transition + "bd"].projector;
            }
        }
        polaron_factors.emplace_back( temp );
    }

    Log::L2( "[System-OperatorMatrices] Creating Hamiltonoperator...\n" );
    // Generate Self Action Hamilton H_0:
    H_0 = Sparse( base.size(), base.size() );
    for ( auto &state : el_states ) {
        double energy = p.input_electronic[state.first].numerical["Energy"];
        H_0 += energy * state.second.hilbert;
    }
    for ( auto &state : ph_states ) {
        double energy = p.input_photonic[state.first].numerical["Energy"];
        H_0 += energy * state.second.hilbert;
    }
    // Generate Interaction Hamilton H_I:
    H_I = Sparse( base.size(), base.size() );
    Sparse H_I_a = H_I; // RWA Part
    Sparse H_I_b = H_I; // Non RWA Part
    for ( auto &cav : p.input_photonic ) {
        for ( int i = 0; i < cav.second.string_v["CoupledTo"].size(); i++ ) {
            std::string transition = cav.second.string_v["CoupledTo"][i];
            std::string transition_transposed = transition;
            std::reverse( transition_transposed.begin(), transition_transposed.end() );
            H_I_a += cav.second.numerical_v["CouplingScaling"][i] * el_transitions[transition_transposed].hilbert * ph_transitions[cav.first + "b"].hilbert;
            H_I_a += cav.second.numerical_v["CouplingScaling"][i] * el_transitions[transition].hilbert * ph_transitions[cav.first + "bd"].hilbert;
            H_I_b += cav.second.numerical_v["CouplingScaling"][i] * el_transitions[transition_transposed].hilbert * ph_transitions[cav.first + "bd"].hilbert;
            H_I_b += cav.second.numerical_v["CouplingScaling"][i] * el_transitions[transition].hilbert * ph_transitions[cav.first + "b"].hilbert;
        }
    }

    Log::L2( "[System-OperatorMatrices] Creating Path Integral Sorting Vectors...\n" );
    // The Path Integral can be partially summed by either: the electronic state index; or the electronic state coupling factor.
    // Because multiple states can have the same coupling factor, the latter will usually be faster, but at least as fast as the index summation.
    // As a default value, "factor" should be used.
    std::string sorting = "factor";
    if ( sorting == "index" ) {
        std::map<std::string, int> temp_base_indices;
        int new_index = 0;
        for ( int i = 0; i < base.size(); i++ ) {
            std::string state1 = base.at( i ).substr( 1, 1 );
            auto factor = (double)p.input_electronic[state1].numerical["PhononCoupling"].get();
            auto index = state1;
            if ( !temp_base_indices.count( index ) > 0 ) {
                temp_base_indices[index] = new_index++;
                phonon_group_index_to_coupling_value.emplace_back( factor );
            }
            phonon_hilbert_index_to_group_index.emplace_back( temp_base_indices[index] );
        }
        phonon_group_index_to_hilbert_indices = std::vector<std::vector<int>>( temp_base_indices.size() );
        for ( int i = 0; i < base.size(); i++ ) {
            phonon_group_index_to_hilbert_indices[phonon_hilbert_index_to_group_index[i]].emplace_back( i );
        }
    } else if ( sorting == "factor" ) {
        std::map<double, int> temp_base_indices;
        int new_index = 0;
        for ( int i = 0; i < base.size(); i++ ) {
            std::string state1 = base.at( i ).substr( 1, 1 );
            double factor = (double)p.input_electronic[state1].numerical["PhononCoupling"].get();
            auto index = factor; // state1;
            if ( !temp_base_indices.count( index ) > 0 ) {
                temp_base_indices[index] = new_index++;
                phonon_group_index_to_coupling_value.emplace_back( factor );
            }
            phonon_hilbert_index_to_group_index.emplace_back( temp_base_indices[index] );
        }
        phonon_group_index_to_hilbert_indices = std::vector<std::vector<int>>( temp_base_indices.size() );
        for ( int i = 0; i < base.size(); i++ ) {
            phonon_group_index_to_hilbert_indices[phonon_hilbert_index_to_group_index[i]].emplace_back( i );
        }
    } else {
        Log::L2( "[System-OperatorMatrices] Sorting Parameter for Partially Summed ADM mismatch!\n" );
    }
    Log::L2( "[System-OperatorMatrices] Phonon Coupling Index Vector: {}\n", phonon_hilbert_index_to_group_index );
    Log::L2( "[System-OperatorMatrices] Phonon Coupling Value Vector: {}\n", phonon_group_index_to_coupling_value );

    // Creating Phonon Coupling Matrix
    polaron_phonon_coupling_matrix = Sparse( base.size(), base.size() );
    for ( auto i = 0; i < polaron_phonon_coupling_matrix.rows(); i++ )
        for ( auto j = 0; j < polaron_phonon_coupling_matrix.cols(); j++ ) {
            Scalar coupling = phonon_group_index_to_coupling_value[phonon_hilbert_index_to_group_index[i]] * phonon_group_index_to_coupling_value[phonon_hilbert_index_to_group_index[j]];
            if ( QDLC::Math::abs2( coupling ) == 0.0 )
                coupling = std::max( phonon_group_index_to_coupling_value[phonon_hilbert_index_to_group_index[i]], phonon_group_index_to_coupling_value[phonon_hilbert_index_to_group_index[j]] );
            if ( QDLC::Math::abs2( coupling ) != 0.0 )
                polaron_phonon_coupling_matrix.coeffRef( i, j ) = coupling;
        }
    // Map this matrix to every nonzero entries in any of the polaron transition matrices, all other entries can be discarded
    Sparse projector = QDLC::Matrix::sparse_projector( std::accumulate( polaron_factors.begin(), polaron_factors.end(), Sparse( base.size(), base.size() ) ) );
    Sparse projector_adjoint = projector.adjoint();
    // Remove to-be-zero elements from the Coupling Matrix
    // polaron_phonon_coupling_matrix = polaron_phonon_coupling_matrix.cwiseProduct( projector + projector_adjoint );
    // Scale Interaction Hamilton Operator and Pulse matrices with <B> if PME is used. Note that because K = exp(-0.5*integral J), the scaling has to be exponential in B, e.g. local_b^factor
    // The Scaling for the Polaron Frame is squared, because the coupling factor n_Level is applied twice in the PI and otherwise only once in the PME.
    if ( p.numerics_phonon_approximation_order != PHONON_PATH_INTEGRAL ) {
        Sparse b_matrix = polaron_phonon_coupling_matrix.unaryExpr( [&]( Scalar val ) { return std::pow( p.p_phonon_b.get(), 1.0 * val ); } );
        Log::L2( "[System-OperatorMatrices] Scaling H_I,Cavity and H_I,Pulse with <B> = {}\n", p.p_phonon_b );
        Log::L2( "[System-OperatorMatrices] <B>-Matrix:\n{}\n", Dense( b_matrix ).format( output_format ) );
        H_I_a = H_I_a.cwiseProduct( b_matrix );
        H_I_b = H_I_b.cwiseProduct( b_matrix );
        // Scale Pulse Matrices for Pulse evaluations for the Hamilton Operators. Note: If the pulse matrix consists of cavity transitions, the <B>-scaling must not be applied there!
        std::for_each( pulse_mat.begin(), pulse_mat.end(), [&]( auto &mat ) { mat = mat.cwiseProduct( b_matrix ); } );
        // Scale Polaron Cache Matrices for Chi evaluations
        std::for_each( polaron_factors.begin(), polaron_factors.end(), [&]( auto &mat ) { mat = mat.cwiseProduct( b_matrix ); } );
    }
    // Add Cavity Pulse Transitions after scaling of the electronic transitions. Note that this way, the pure cavity driving is not influenced by the phonon coupling. IDK if this is correct.
    std::for_each( pulse_mat.begin(), pulse_mat.end(), [&, indx = 0]( auto &mat ) mutable { mat += pulse_mat_cavity_cache[indx]; indx++; } );

    for ( auto &a : phonon_group_index_to_hilbert_indices )
        Log::L2( "[System-OperatorMatrices] Phonon Group Index Vector: {}\n", a );

    Log::L2( "[System-OperatorMatrices] Creating Initial State Vector...\n" );
    // Create Initial State
    // Split starting state into superposition. States can be passed as "|...>+|...>" with amplitudes
    initial_state_vector_ket = Dense::Zero( base.size(), 1 );
    for ( auto state : QDLC::String::splitline( p.p_initial_state_s, '+' ) ) {
        // For a coherent and normal input the syntax is |base>. replace number in base with 'alpha#' to indicate a coherent state. replace number in base with 'xiy&' to indicate a squeezed state. no spaces.
        // Scalar amp = state[0] == '|' ? 1.0 : ( state.back() == 'i' ? 1.0i * stod( QDLC::String::splitline( state.substr( 0, state.size() - 1 ), '|' ).front() ) : stod( QDLC::String::splitline( state, '|' ).front() ) );
        Scalar amp = state[0] == '|' ? 1.0 : stod( QDLC::String::splitline( state, '|' ).front() );
        std::string pure_state = state.substr( state.find( '|' ) );
        // Coherent State
        if ( pure_state.find( "#" ) != std::string::npos ) {
            auto modev = QDLC::String::splitline( pure_state, '#' );
            auto front = QDLC::String::splitline( modev.front(), '|' );
            std::string front_base = "|";
            for ( int i = 0; i < front.size() - 1; i++ ) {
                front_base += front[i] + "|";
            }
            double alpha = std::stod( front.back() );
            std::string mode = modev.back().substr( 0, 1 );
            Log::L2( "[System-OperatorMatrices] Creating superpositioned coherent state {} for mode {} with alpha = {} and scaled amplitude {}\n", pure_state, mode, alpha, amp );
            for ( int n = 0; n < p.input_photonic[mode].numerical["MaxPhotons"]; n++ ) {
                std::string current = front_base + std::to_string( n ) + modev.back();
                Log::L2( "[System-OperatorMatrices] Creating coherent substate {} ({}) with amplitude {}\n", current, base_index_map[current], QDLC::Math::getCoherent( alpha, n ) * amp );
                // Add initial state with amplitudes
                initial_state_vector_ket( base_index_map[current] ) += amp * QDLC::Math::getCoherent( alpha, n );
            }
        }
        // Squeezed State
        else if ( pure_state.find( "&" ) != std::string::npos ) {
            auto modev = QDLC::String::splitline( pure_state, '&' );
            auto front = QDLC::String::splitline( modev.front(), '|' );
            std::string front_base = "|";
            for ( int i = 0; i < front.size() - 1; i++ ) {
                front_base += front[i] + "|";
            }
            auto cmplx = QDLC::String::splitline( front.back(), '_' );
            double r = std::stod( cmplx.front() );
            double phi = std::stod( cmplx.back() );
            std::string mode = modev.back().substr( 0, 1 );
            Log::L2( "[System-OperatorMatrices] Creating superpositioned squeezed state {} for mode {} with r = {}, phi = {} and scaled amplitude {}\n", pure_state, mode, r, phi, amp );
            for ( int n = 0; n < p.input_photonic[mode].numerical["MaxPhotons"]; n++ ) {
                if ( n % 2 != 0 )
                    continue;
                std::string current = front_base + std::to_string( n ) + modev.back();
                Log::L2( "[System-OperatorMatrices] Creating squeezed substate {} ({}) with amplitude {}\n", current, base_index_map[current], QDLC::Math::getSqueezed( r, phi, n ) * amp );
                // Add initial state with amplitudes
                initial_state_vector_ket( base_index_map[current] ) += amp * QDLC::Math::getSqueezed( r, phi, n );
            }
        }
        // Thermal State
        else if ( pure_state.find( "@" ) != std::string::npos ) {
            auto modev = QDLC::String::splitline( pure_state, '@' );
            auto front = QDLC::String::splitline( modev.front(), '|' );
            std::string front_base = "|";
            for ( int i = 0; i < front.size() - 1; i++ ) {
                front_base += front[i] + "|";
            }
            double alpha = std::stod( front.back() );
            std::string mode = modev.back().substr( 0, 1 );
            Log::L2( "[System-OperatorMatrices] Creating superpositioned Thermal state {} for mode {} with alpha = {} and scaled amplitude {}\n", pure_state, mode, alpha, amp );
            for ( int n = 0; n < p.input_photonic[mode].numerical["MaxPhotons"]; n++ ) {
                std::string current = front_base + std::to_string( n ) + modev.back();
                Log::L2( "[System-OperatorMatrices] Creating thermal substate {} ({}) with amplitude {}\n", current, base_index_map[current], QDLC::Math::getThermal( alpha, n ) * amp );
                // Add initial state with amplitudes
                initial_state_vector_ket( base_index_map[current] ) += amp * QDLC::Math::getThermal( alpha, n );
            }
        } else {
            Log::L2( "[System-OperatorMatrices] Creating superpositioned state {} ({}) with amplitude {}\n", pure_state, base_index_map[pure_state], amp );
            // Add initial state with amplitudes
            initial_state_vector_ket( base_index_map[pure_state] ) += amp;
        }
    }

    initial_state_vector_ket.normalize();
    Log::L2( "[System-OperatorMatrices] Initial State Vector: [{}]\n", initial_state_vector_ket.format( Eigen::IOFormat( 4, 0, ", ", " ", "", "" ) ) );
    rho = ( initial_state_vector_ket * initial_state_vector_ket.transpose() ).sparseView();

    // Choose Final Hamilton. The Coupling scalings are incorporateed in H_I_a and H_I_b. The PME scaling <B> is incorporated in H_I_a/b and the Pulse Matrices.
    Log::L2( "[System-OperatorMatrices] Choosing final Hamilton Operator...\n" );
    if ( p.numerics_use_rwa )
        H_I = p.p_omega_coupling * H_I_a;
    else
        H_I = p.p_omega_coupling * ( H_I_a + H_I_b );
    if ( p.numerics_use_interactionpicture )
        H_used = H_I;
    else
        H_used = H_0 + H_I;

    Log::L2( "[System-OperatorMatrices] Hamilton Eigenvalues:\n" );
    Log::L2( "[System-OperatorMatrices] H_0: [{}]\n", ( Dense( H_0 ).eigenvalues() * QDLC::Math::ev_conversion ).format( Eigen::IOFormat( -1, 0, ", ", ", ", "", "" ) ) );
    Log::L2( "[System-OperatorMatrices] H_I: [{}]\n", ( Dense( H_I ).eigenvalues() * QDLC::Math::ev_conversion ).format( Eigen::IOFormat( -1, 0, ", ", ", ", "", "" ) ) );
    Log::L2( "[System-OperatorMatrices] H: [{}]\n", ( Dense( H_0 + H_I ).eigenvalues() * QDLC::Math::ev_conversion ).format( Eigen::IOFormat( -1, 0, ", ", ", ", "", "" ) ) );
    return true;
}

void OperatorMatrices::output_operators( Parameters &p ) {
    if ( p.output_operators > 0 ) {
        std::ostringstream out;
        Eigen::IOFormat CleanFmt( 4, 0, ", ", "\n", "[", "]" );
        out << "[System-OperatorMatrices] Hamilton and Rho:\n";
        out << "H_0\n"
            << Dense( H_0 ).format( CleanFmt ) << std::endl;
        out << "H_I\n"
            << Dense( H_I ).format( CleanFmt ) << std::endl;
        out << "H_used\n"
            << Dense( H_used ).format( CleanFmt ) << std::endl;
        out << "rho\n"
            << Dense( rho ).format( CleanFmt ) << std::endl;
        // out << "test1\n" << test1.format(CleanFmt)<< "\ntest2\n" << test2.format(CleanFmt) << std::endl;
        Log::L2( out.str() );
        // Log::L2( "[System-OperatorMatrices] Outputting String Matrices...\n" );
        //  OperatorMatricesText test = OperatorMatricesText();
        //  test.generate_operators( p );
        if ( p.output_operators == 3 )
            exit( 0 );
    }
}