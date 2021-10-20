#include "system/operatormatrices.h"
//#include "system/operatormatrices_text.h"

OperatorMatrices::OperatorMatrices( Parameters &p ) {
    Timer &timer_operatormatrices = Timers::create( "Operator Matrices", true, false );
    timer_operatormatrices.start();
    Log::L2( "Generating operator matrices... " );
    if ( !generateOperators( p ) ) {
        Log::L2( "Generating operator matrices failed! Exitting program...\n" );
        Log::close();
        exit( EXIT_FAILURE );
    }
    timer_operatormatrices.end();
    Log::L2( "successful. Elapsed time is {}ms\n", timer_operatormatrices.getWallTime( TIMER_MILLISECONDS ) );
}

bool OperatorMatrices::generateOperators( Parameters &p ) {
    Eigen::IOFormat CleanFmt( 4, 0, ", ", "\n", "[", "]" );
    // Zeroing Hamiltons (redundant at this point)
    Log::L2( "Creating operator matrices, dimension = {}\nCreating base matrices...\n", p.maxStates );

    // Generate Electronic and Photonic Base States and Self-Hilbert Matrices
    // Electronic
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
        //std::cout << "Matrix " << cav.first << "\n"
        //          << state.self_hilbert << std::endl;
    }

    // Tensor all matrices into the total Hilbert Space.
    // String Base Tensors
    auto subb = base_photonic;
    subb.insert( subb.begin(), base_electronic );
    base = tensor( subb );

    // Sparse total Hilbert tensors
    for ( auto &bel : el_states ) {
        auto current = base_selfhilbert;
        current.front() = bel.second.self_hilbert; // Electronic states have basis 0 for now. Maybe change when TODO: add tensoring of multiple electronic bases
        bel.second.hilbert = tensor( current ).sparseView();
        bel.second.projector = project_matrix_sparse( bel.second.hilbert );
        //std::cout << "Mat = " << bel.first << "  " << bel.second.hilbert.rows() << "\n"
        //          << Dense( bel.second.hilbert ) << std::endl;
    }
    for ( auto &phot : ph_states ) {
        auto current = base_selfhilbert;
        current[phot.second.base] = phot.second.self_hilbert;
        phot.second.hilbert = tensor( current ).sparseView();
        phot.second.projector = project_matrix_sparse( phot.second.hilbert );
        //std::cout << "Mat = " << phot.first << "  " << phot.second.hilbert.rows() << "\n"
        //          << Dense( phot.second.hilbert ) << std::endl;
    }
    // Generate Transition Matrices
    for ( auto &bel : el_states ) {
        for ( auto &trans_to : p.input_electronic[bel.first].string_v["CoupledTo"] ) {
            if ( trans_to.front() == '-' ) continue;
            std::string transition = bel.first + trans_to;
            el_transitions[transition].self_hilbert = bel.second.ket * el_states[trans_to].bra;
            el_transitions[transition].direction = -1; //DOWN
            el_transitions[transition].energy = std::abs( p.input_electronic[trans_to].numerical["Energy"] - p.input_electronic[bel.first].numerical["Energy"] );
            std::string transition_transposed = trans_to + bel.first;
            el_transitions[transition_transposed].self_hilbert = el_states[trans_to].ket * bel.second.bra;
            el_transitions[transition_transposed].direction = 1; //UP
            el_transitions[transition_transposed].energy = std::abs( p.input_electronic[trans_to].numerical["Energy"] - p.input_electronic[bel.first].numerical["Energy"] );
        }
    }
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
    // Tensor remaining transition
    for ( auto &bel : el_transitions ) {
        auto current = base_selfhilbert;
        current.front() = bel.second.self_hilbert;
        bel.second.hilbert = tensor( current ).sparseView();
        bel.second.projector = project_matrix_sparse( bel.second.hilbert );
    }
    for ( auto &phot : ph_transitions ) {
        auto current = base_selfhilbert;
        current[phot.second.base] = phot.second.self_hilbert;
        phot.second.hilbert = tensor( current ).sparseView();
        phot.second.projector = project_matrix_sparse( phot.second.hilbert );
    }

    // Create total Hilbert space indices. This routine assumes the individual bases are added in ascending order, e.g. 0,1,2,3...:
    for ( int i = 0; i < base_selfhilbert.size(); i++ ) {
        auto current = base_selfhilbert;
        int index = 1;
        for ( int k = 0; k < current[i].rows(); k++ )
            for ( int j = 0; j < current[i].cols(); j++ )
                current[i]( k, j ) = ( k + 1 ) + 1i * ( j + 1 );

        base_hilbert_index.emplace_back( tensor( current ) );
        //Log::L2( "Tensor for base i = {}:\nCurrent:\n{}\n\n{}\n\n", i, current[i], base_hilbert_index.back() );
    }

    // Generate prechaced pulse and chirp matrices. 2 matrices are generated per pulse for Omega and Omega^*
    // Generate transition matrices from the initial state matrices to allow for addidional transitions and the cavity modes to be included.
    for ( auto &pulse : p.input_pulse ) {
        Sparse pulsemat = Sparse( base.size(), base.size() );
        Sparse pulsemat_star = Sparse( base.size(), base.size() );
        for ( int i = 0; i < pulse.second.string_v["CoupledTo"].size(); i++ ) {
            std::string transition = pulse.second.string_v["CoupledTo"][i];
            std::string transition_transposed = transition;
            std::reverse( transition_transposed.begin(), transition_transposed.end() );
            if ( el_transitions.count( transition ) > 0 ) {
                pulsemat += el_transitions[transition].hilbert;
                pulsemat_star += el_transitions[transition_transposed].hilbert;
            } else if ( std::isupper( transition.front() ) ) {
                Log::L2( "Electronic Pulse transition {} is not in the list of allowed electronic transitions, recreating transition matrices...\n", transition );
                auto ket1 = el_states[transition.substr( 0, 1 )].ket;
                auto bra1 = el_states[transition.substr( 0, 1 )].bra;
                auto ket2 = el_states[transition.substr( 1, 1 )].ket;
                auto bra2 = el_states[transition.substr( 1, 1 )].bra;
                auto current = base_selfhilbert;
                current.front() = ket2 * bra1;
                Sparse transition_hilbert = tensor( current ).sparseView();
                current.front() = ket1 * bra2;
                Sparse transition_transposed_hilbert = transition_hilbert.adjoint();//tensor( current ).sparseView();
                pulsemat += transition_hilbert;
                pulsemat_star += transition_transposed_hilbert;
            } else {
                Log::L2( "Pulse transition {} is cavity...\n", transition );
                pulsemat += ph_transitions[transition + "b"].hilbert;
                pulsemat_star += ph_transitions[transition + "bd"].hilbert;
            }
        }
        pulse_mat.emplace_back( pulsemat );
        pulse_mat.emplace_back( pulsemat_star );
    }
    // 1 matrix is generated per chirp
    for ( auto &chirp : p.input_chirp ) {
        Sparse chirpmat = Sparse( base.size(), base.size() );
        Sparse chirpmat_star = Sparse( base.size(), base.size() );
        for ( int i = 0; i < chirp.second.string_v["CoupledTo"].size(); i++ ) {
            std::string state = chirp.second.string_v["CoupledTo"][i];
            chirpmat += chirp.second.numerical_v["AmpFactor"][i] * el_states[state].hilbert; //TODO: chirp cavity
        }
        chirp_mat.emplace_back( chirpmat );
    }

    // Create Index Map
    index = 0;
    Log::L1( "System Base ({}): ", base.size() );
    for ( auto &b : base ) {
        b = "|" + b + ">";
        base_index_map[b] = index++;
        Log::L1( "{} ", b );
        //std::cout << b << " ";
    }
    Log::L1( "\n" );

    //FIXME:
    p.maxStates = base.size();

    H = Sparse( p.maxStates, p.maxStates );
    H_0 = Sparse( p.maxStates, p.maxStates );
    H_I = Sparse( p.maxStates, p.maxStates );
    H_used = Sparse( p.maxStates, p.maxStates );
    rho = Sparse( p.maxStates, p.maxStates );

    // Analytical Time Trafo Matrix
    timetrafo_cachematrix = Dense::Zero( base.size(), base.size() );
    for ( auto &iss : base_index_map ) {
        std::string is = iss.first;
        is.pop_back();
        is.erase( is.begin() );
        auto i = splitline( is, '|' );
        for ( auto &jss : base_index_map ) {
            Scalar val = 0;
            std::string js = jss.first;
            js.pop_back();
            js.erase( js.begin() );
            auto j = splitline( js, '|' );
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
                //std::cout << iss.first << ", " << jss.first << " --> " << ni << " " << nj << ", modes = " << mode << std::endl;
                val += p.input_photonic[mode].numerical["Energy"] * ( ni - nj );
            }
            timetrafo_cachematrix( iss.second, jss.second ) = 1i * val;
        }
    }
    // Precalculate Lindbladians
    lindblad_factors = std::vector<Sparse>( 6, Sparse( base.size(), base.size() ) );
    // Cavity Terms //TODO: Testen das die hier richtig sind!
    //for ( auto &cav : ph_transitions ) {
    //    if ( cav.second.direction == 1 )
    //        continue;
    //    std::string mode = cav.first.substr( 0, 1 );
    //    auto &params = p.input_photonic[mode];
    //    lindblad_factors[0] += 2.0 * 0.5 * p.p_omega_cavity_loss * params.numerical["DecayScaling"] * ph_transitions[mode + "b"].hilbert;
    //    lindblad_factors[1] += ph_transitions[mode + "bd"].hilbert;
    //    lindblad_factors[2] += 0.5 * p.p_omega_cavity_loss * params.numerical["DecayScaling"] * ph_transitions[mode + "bd"].hilbert * ph_transitions[mode + "b"].hilbert;
    //    lindblad_factors[3] += 0.5 * p.p_omega_cavity_loss * params.numerical["DecayScaling"] * ph_transitions[mode + "bd"].hilbert * ph_transitions[mode + "b"].hilbert;
    //}
    //// Radiative decay //TODO: Testen das die hier richtig sind!
    //for ( auto &trans : el_transitions ) {
    //    if ( trans.second.direction == 1 )
    //        continue;
    //    std::string transition = trans.first;
    //    std::string state = transition.substr( transition.size() - 1, 1 );
    //    auto &params = p.input_electronic[state];
    //    //std::cout << "Transition " << transition << " with state info from " << state << " --> " << params.numerical["DecayScaling"] << std::endl;
    //    std::string trans_transposed = transition;
    //    std::reverse( trans_transposed.begin(), trans_transposed.end() );
    //    lindblad_factors[0] += 0.5 * p.p_omega_decay * params.numerical["DecayScaling"] * el_transitions[transition].hilbert;
    //    lindblad_factors[1] += el_transitions[trans_transposed].hilbert;
    //    lindblad_factors[2] += 0.5 * p.p_omega_decay * params.numerical["DecayScaling"] * el_transitions[trans_transposed].hilbert * el_transitions[transition].hilbert;
    //    lindblad_factors[3] += 0.5 * p.p_omega_decay * params.numerical["DecayScaling"] * el_transitions[transition].hilbert * el_transitions[trans_transposed].hilbert;
    //}
    //// Electronic Dephasing
    //for ( auto &state_a : el_states ) {
    //    for ( auto &state_b : el_states ) {
    //        if ( state_a.first.compare( state_b.first ) == 0 )
    //            continue;
    //        lindblad_factors[4] += 0.5 * p.p_omega_pure_dephasing * p.input_electronic[state_a.first].numerical["DephasingScaling"] * el_states[state_a.first].hilbert;
    //        lindblad_factors[5] += p.input_electronic[state_b.first].numerical["DephasingScaling"] * el_states[state_b.first].hilbert;
    //    }
    //}

    // Precalculate Polaron Matrices
    polaron_factors.emplace_back( Sparse( base.size(), base.size() ) );
    for ( auto &[mode, param] : p.input_photonic ) {
        int i = 0;
        for ( auto transition : param.string_v["CoupledTo"] ) {
            std::reverse( transition.begin(), transition.end() );
            polaron_factors.back() += el_transitions[transition].hilbert * p.p_omega_coupling * param.numerical_v["CouplingScaling"][i] * ph_transitions[mode + "b"].hilbert; // TODO: was ist mit back coupling?
            i++;
        }
    }
    for ( auto &[mode, param] : p.input_pulse ) {
        Sparse temp = Sparse( base.size(), base.size() );
        for ( auto transition : param.string_v["CoupledTo"] ) {
            std::reverse( transition.begin(), transition.end() );
            if ( el_transitions.count( transition ) > 0 ) {
                temp += el_transitions[transition].hilbert;
            } else if ( std::isupper( transition.front() ) ) {
                Log::L2( "Electronic Pulse transition {} is not in the list of allowed electronic transitions, recreating transition matrices...\n", transition );
                auto ket1 = el_states[transition.substr( 0, 1 )].ket;
                auto bra1 = el_states[transition.substr( 0, 1 )].bra;
                auto ket2 = el_states[transition.substr( 1, 1 )].ket;
                auto bra2 = el_states[transition.substr( 1, 1 )].bra;
                auto current = base_selfhilbert;
                current.front() = ket2 * bra1;
                Sparse transition_transposed_hilbert = tensor( current ).sparseView();
                temp += transition_transposed_hilbert;
            } else {
                Log::L2( "Pulse transition {} is cavity...\n", transition );
                temp += ph_transitions[transition + "bd"].hilbert;
            }
        }
        polaron_factors.emplace_back( temp );
    }

    Log::L2( "Done! Creating Hamiltonoperator...\n" );
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
    if ( p.numerics_use_rwa )
        H_I = p.p_omega_coupling * H_I_a;
    else
        H_I = p.p_omega_coupling * ( H_I_a + H_I_b );
    H = H_0 + H_I;
    if ( p.numerics_use_interactionpicture )
        H_used = H_I;
    else
        H_used = H;

    phononCouplingFactor = dDense::Zero( base.size(), base.size() );
    std::map<double, int> temp_base_indices;
    int new_index = 0;
    for ( int i = 0; i < base.size(); i++ ) {
        for ( int j = 0; j < base.size(); j++ ) { // base is |el|...>
            std::string state1 = base.at( i ).substr( 1, 1 );
            std::string state2 = base.at( j ).substr( 1, 1 );
            if ( i == j ) {
                phononCouplingFactor( i, j ) = (double)std::min( p.input_electronic[state1].numerical["PhononCoupling"].get() * p.input_electronic[state2].numerical["PhononCoupling"].get(), std::max( p.input_electronic[state1].numerical["PhononCoupling"].get(), p.input_electronic[state2].numerical["PhononCoupling"].get() ) );
                if ( !temp_base_indices.contains( phononCouplingFactor( i, j ) ) ) {
                    temp_base_indices[phononCouplingFactor( i, j )] = new_index++;
                    phononCouplingIndexValue.emplace_back( phononCouplingFactor( i, j ) );
                }
                phononCouplingIndex.emplace_back( temp_base_indices[phononCouplingFactor( i, j )] );
            }
        }
    }
    //Log::L2( "Phonon Coupling Matrix:\n{}\n", phononCouplingFactor.format( CleanFmt ) );
    Log::L2( "Phonon Coupling Index Vector: {}\n", phononCouplingIndex );
    Log::L2( "Phonon Coupling Value Vector: {}\n", phononCouplingIndexValue );

    Log::L2( "Hamiltonoperator done!\n" );

    // Create Initial State
    // Split starting state into superposition. States can be passed as "|...>+|...>" with amplitudes
    initial_state_vector_ket = Dense::Zero( base.size(), 1 );
    for ( auto state : splitline( p.p_initial_state_s, '+' ) ) {
        // For a coherent and normal input the syntax is |base>. replace number in base with 'alpha#' to indicate a coherent state. replace number in base with 'xiy&' to indicate a squeezed state. no spaces.
        double amp = state[0] == '|' ? 1.0 : stod( splitline( state, '|' ).front() );
        std::string pure_state = state.substr( state.find( '|' ) );
        if ( state.find( "#" ) != std::string::npos ) {
            auto modev = splitline( state, '#' );
            auto front = splitline( modev.front(), '|' );
            std::string front_base = "|";
            for ( int i = 0; i < front.size() - 1; i++ ) {
                front_base += front[i] + "|";
            } 
            double alpha = std::stod( front.back() );
            std::string mode = modev.back().substr( 0, 1 );
            Log::L2( "Creating superpositioned coherent state {} for mode {} with alpha = {} and scaled amplitude {}\n", pure_state, mode, alpha, amp );
            for ( int n = 0; n < p.input_photonic[mode].numerical["MaxPhotons"]; n++ ) {
                std::string current = front_base + std::to_string( n ) + modev.back();
                Log::L2( "Creating coherent substate {} ({}) with amplitude {}\n", current, base_index_map[current], getCoherent( alpha, n ) * amp );
                // Add initial state with amplitudes
                initial_state_vector_ket( base_index_map[current] ) += amp * getCoherent( alpha, n );
            }
        } else if ( state.find( "&" ) != std::string::npos ) {
            auto modev = splitline( state, '&' );
            auto front = splitline( modev.front(), '|' );
            std::string front_base = "|";
            for ( int i = 0; i < front.size() - 1; i++ ) {
                front_base += front[i] + "|";
            }
            auto cmplx = splitline( front.back(), 'i' );
            double r = std::stod( cmplx.front() );
            double phi = std::stod( cmplx.back() );
            std::string mode = modev.back().substr( 0, 1 );
            Log::L2( "Creating superpositioned coherent state {} for mode {} with r = {}, phi = {} and scaled amplitude {}\n", pure_state, mode, r,phi, amp );
            for ( int n = 0; n < p.input_photonic[mode].numerical["MaxPhotons"]; n++ ) {
                if ( n % 2 != 0 )
                    continue;
                std::string current = front_base + std::to_string( n ) + modev.back();
                Log::L2( "Creating coherent substate {} ({}) with amplitude {}\n", current, base_index_map[current], getSqueezed( r, phi, n ) * amp );
                // Add initial state with amplitudes
                initial_state_vector_ket( base_index_map[current] ) += amp * getSqueezed( r, phi, n );
            }
        } else {
            Log::L2( "Creating superpositioned state {} ({}) with amplitude {}\n", pure_state, base_index_map[pure_state], amp );
            // Add initial state with amplitudes
            initial_state_vector_ket( base_index_map[pure_state] ) += amp;
        }
    }

    initial_state_vector_ket.normalize();
    Log::L2( "Initial State Vector: {}\n", initial_state_vector_ket );
    rho = ( initial_state_vector_ket * initial_state_vector_ket.transpose() ).sparseView();

    // if ( !p.startCoherent || true ) {
    //p.p_initial_state = base_index_map[p.p_initial_state_s];
    //Log::L2( "Setting initial rho as pure state with rho_0 = {}... \n", p.p_initial_state );
    //rho.coeffRef( p.p_initial_state, p.p_initial_state ) = 1;
    //std::cout << "initial state = " << p.p_initial_state_s << " -> " << p.p_initial_state << std::endl;
    //} else {
    //Log::L2( "Setting initial rho as pure coherent state with alpha_h = {}, alpha_v = ... \n", p.p_initial_state_photon_h, p.p_initial_state_photon_v );
    //double trace_rest_h = 1.0;
    //double trace_rest_v = 1.0;
    //for ( int i = 0; i < p.p_max_photon_number; i++ ) {
    //    auto coherent_h = getCoherent( std::sqrt( p.p_initial_state_photon_h ), i );
    //    auto coherent_v = 0.0; //getCoherent( std::sqrt(p.p_initial_state_photon_v), i );
    //    int state_h = 4 * ( p.p_max_photon_number + 1 ) * i;
    //    int state_v = 4 * i;
    //    rho.coeffRef( state_h, state_h ) = coherent_h; // Remember, <n> = alpha^2 -> alpha = sqrt(n) !!
    //    //rho.coeffRef( state_v, state_v ) = coherent_v; // Remember, <n> = alpha^2 -> alpha = sqrt(n) !!
    //    trace_rest_h -= coherent_h;
    //    trace_rest_v -= coherent_v;
    //    Log::L3( "Coherent state at N = {} with coefficient H = {}, V = {}\n", i, coherent_h, coherent_v );
    //}
    //int final_h = 4 * ( p.p_max_photon_number + 1 ) * p.p_max_photon_number;
    //int final_v = 4 * p.p_max_photon_number;
    //rho.coeffRef( final_h, final_h ) = trace_rest_h;
    ////rho.coeffRef( final_v, final_v ) = trace_rest_v;
    //Log::L3( "Coherent state at N = {} with coefficient H = {}, V = {}\n", p.p_max_photon_number, trace_rest_h, trace_rest_v );
    //}
    Log::L2( "Done!\n" );
    return true;
}

void OperatorMatrices::outputOperators( Parameters &p ) {
    if ( p.output_operators > 0 ) {
        std::ostringstream out;
        Eigen::IOFormat CleanFmt( 4, 0, ", ", "\n", "[", "]" );
        out << "Hamilton and Rho:\nH=H_0+H_I (no RWA)\n"
            << Dense( H ).format( CleanFmt ) << std::endl;
        out << "H_0\n"
            << Dense( H_0 ).format( CleanFmt ) << std::endl;
        out << "H_I\n"
            << Dense( H_I ).format( CleanFmt ) << std::endl;
        out << "H_used\n"
            << Dense( H_used ).format( CleanFmt ) << std::endl;
        out << "rho\n"
            << Dense( rho ).format( CleanFmt ) << std::endl;
        //out << "test1\n" << test1.format(CleanFmt)<< "\ntest2\n" << test2.format(CleanFmt) << std::endl;
        Log::L2( out.str() );
        Log::L2( "Outputting String Matrices...\n" );
        //OperatorMatricesText test = OperatorMatricesText();
        //test.generateOperators( p );
        if ( p.output_operators == 3 )
            exit( 0 );
    }
}