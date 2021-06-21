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
    Dense m_base_electronic = Dense::Identity( p.input_electronic.size(), p.input_electronic.size() );
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
    std::vector<Dense> m_base_photonic;
    int curcav = 0;
    for ( auto &cav : p.input_photonic ) {
        int max_photons = cav.second.numerical["MaxPhotons"];
        m_base_photonic.emplace_back( Dense::Identity( max_photons + 1, max_photons + 1 ) );
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
    std::vector<Dense> all_bases = m_base_photonic;
    all_bases.insert( all_bases.begin(), m_base_electronic );
    for ( auto &bel : el_states ) {
        auto current = all_bases;
        current.front() = bel.second.self_hilbert;
        bel.second.hilbert = tensor( current ).sparseView();
        bel.second.projector = project_matrix_sparse( bel.second.hilbert );
        //std::cout << "Mat = " << bel.first << "  " << bel.second.hilbert.rows() << "\n"
        //          << Dense( bel.second.hilbert ) << std::endl;
    }
    for ( auto &phot : ph_states ) {
        auto current = all_bases;
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
            el_transitions[transition].direction = -1;
            el_transitions[transition].energy = std::abs( p.input_electronic[trans_to].numerical["Energy"] - p.input_electronic[bel.first].numerical["Energy"] );
            std::string transition_transposed = trans_to + bel.first;
            el_transitions[transition_transposed].self_hilbert = el_states[trans_to].ket * bel.second.bra;
            el_transitions[transition_transposed].direction = 1;
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
        auto current = all_bases;
        current.front() = bel.second.self_hilbert;
        bel.second.hilbert = tensor( current ).sparseView();
        bel.second.projector = project_matrix_sparse( bel.second.hilbert );
    }
    for ( auto &phot : ph_transitions ) {
        auto current = all_bases;
        current[phot.second.base] = phot.second.self_hilbert;
        phot.second.hilbert = tensor( current ).sparseView();
        phot.second.projector = project_matrix_sparse( phot.second.hilbert );
    }

    // Generate prechaced pulse and chirp matrices. 2 matrices are generated per pulse for Omega and Omega^*
    for ( auto &pulse : p.input_pulse ) {
        Sparse pulsemat = Sparse( base.size(), base.size() );
        Sparse pulsemat_star = Sparse( base.size(), base.size() );
        for ( int i = 0; i < pulse.second.string_v["CoupledTo"].size(); i++ ) {
            std::string transition = pulse.second.string_v["CoupledTo"][i];
            std::string transition_transposed = transition;
            std::reverse( transition_transposed.begin(), transition_transposed.end() );
            pulsemat += el_transitions[transition].hilbert;
            pulsemat_star += el_transitions[transition_transposed].hilbert;
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
            chirpmat += chirp.second.numerical_v["AmpFactor"][i] * el_states[state].hilbert;
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

    //std::cout << std::endl;

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
            polaron_factors.back() += el_transitions[transition].hilbert * p.p_omega_coupling * param.numerical_v["CouplingScaling"][i] * ph_transitions[mode + "b"].hilbert;
            i++;
        }
    }
    for ( auto &[mode, param] : p.input_pulse ) {
        Sparse temp = Sparse( base.size(), base.size() );
        for ( auto transition : param.string_v["CoupledTo"] ) {
            std::reverse( transition.begin(), transition.end() );
            temp += el_transitions[transition].hilbert;
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
            double scaling = cav.second.numerical_v["CouplingScaling"][i];
            std::string transition_transposed = transition;
            std::reverse( transition_transposed.begin(), transition_transposed.end() );
            H_I_a += scaling * el_transitions[transition_transposed].hilbert * ph_transitions[cav.first + "b"].hilbert;
            H_I_a += scaling * el_transitions[transition].hilbert * ph_transitions[cav.first + "bd"].hilbert;
            H_I_b += scaling * el_transitions[transition_transposed].hilbert * ph_transitions[cav.first + "bd"].hilbert;
            H_I_b += scaling * el_transitions[transition].hilbert * ph_transitions[cav.first + "b"].hilbert;
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

    // Create Initial State

    // REMOVE: OLD CODE
    // Create Base Matrices and Base Vector:
    //Dense m_base1 = Dense::Identity( p.p_max_photon_number + 1, p.p_max_photon_number + 1 );
    //Dense m_base2 = Dense::Identity( p.p_max_photon_number + 1, p.p_max_photon_number + 1 );
    //Dense m_base3 = Dense::Identity( 4, 4 );
    //std::vector<std::string> base1;
    //std::vector<std::string> base2;
    //std::vector<std::string> base3 = { "G", "X_H", "X_V", "B" };
    //for ( int i = 0; i <= p.p_max_photon_number; i++ ) {
    //    base1.emplace_back( std::to_string( i ) + "_H" );
    //    base2.emplace_back( std::to_string( i ) + "_V" );
    //}
    //auto base_old = tensor( tensor( base1, base2 ), base3 );

    Log::L3( "Operator Base: (size {})\n", base.size() );
    phononCouplingFactor = Dense::Zero( base.size(), base.size() );
    for ( int i = 0; i < base.size(); i++ ) {
        for ( int j = 0; j < base.size(); j++ ) {
            //if ( i == j ) {
            if ( base.at( i ).back() == 'G' or base.at( j ).back() == 'G' )
                continue;
            if ( base.at( i ).back() == 'B' or base.at( j ).back() == 'B' ) {
                //phononCouplingFactor( i, j ) = 2.0;
            } else { //if ( base.at( i ).back() == 'H' or base.at( i ).back() == 'V' or base.at( j ).back() == 'H' or base.at( j ).back() == 'V' ) {
                //phononCouplingFactor( i, j ) = 1.0;
            }

            //}
        }
    }

    //Log::L3( "Coupling Matrix:\n{}\n", phononCouplingFactor );
    //std::string bb = "";
    //for ( auto b : base_old ) {
    //    bb = fmt::format( "{}|{}> ", bb, b );
    //}
    //Log::L3( "Base: {}\n", bb );
    //Dense ket = Dense::Zero( 4, 1 );
    //ket( 0 ) = 1;
    //Dense bra = ket.transpose();
    //std::cout << "ket: " << ket.rows() << " " << ket.cols() << std::endl
    //          << ket << std::endl;
    //std::cout << "bra: " << bra.rows() << " " << bra.cols() << std::endl
    //          << bra << std::endl
    //          << std::endl;
    //;
    //Dense comb = bra * ket;
    //Dense comb2 = ket * bra;
    //std::cout << "bra*ket" << std::endl
    //          << comb << std::endl;
    //std::cout << "ket*bra" << std::endl
    //          << comb2 << std::endl;

    // Initializing bare matrices:
    Log::L2( "Initializing base matrices...\n" );
    // Atomic state operators
    //bare_atom_state_ground = Dense::Zero( 4, 4 );
    //bare_atom_state_ground << 1, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0;
    //bare_atom_state_biexciton = Dense::Zero( 4, 4 );
    //bare_atom_state_biexciton << 0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 1;
    //bare_atom_state_H = Dense::Zero( 4, 4 );
    //bare_atom_state_H << 0, 0, 0, 0,
    //    0, 1, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0;
    //bare_atom_state_V = Dense::Zero( 4, 4 );
    //bare_atom_state_V << 0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 1, 0,
    //    0, 0, 0, 0;
    //bare_atom_sigmaplus_G_H = Dense::Zero( 4, 4 );
    //bare_atom_sigmaplus_G_H << 0, 0, 0, 0,
    //    1, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0;
    //bare_atom_sigmaminus_G_H = Dense::Zero( 4, 4 );
    //bare_atom_sigmaminus_G_H << 0, 1, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0;
    //bare_atom_sigmaplus_G_V = Dense::Zero( 4, 4 );
    //bare_atom_sigmaplus_G_V << 0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    1, 0, 0, 0,
    //    0, 0, 0, 0;
    //bare_atom_sigmaminus_G_V = Dense::Zero( 4, 4 );
    //bare_atom_sigmaminus_G_V << 0, 0, 1, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0;
    //bare_atom_sigmaplus_H_B = Dense::Zero( 4, 4 );
    //bare_atom_sigmaplus_H_B << 0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 1, 0, 0;
    //bare_atom_sigmaminus_H_B = Dense::Zero( 4, 4 );
    //bare_atom_sigmaminus_H_B << 0, 0, 0, 0,
    //    0, 0, 0, 1,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0;
    //bare_atom_sigmaplus_V_B = Dense::Zero( 4, 4 );
    //bare_atom_sigmaplus_V_B << 0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 1, 0;
    //bare_atom_sigmaminus_V_B = Dense::Zero( 4, 4 );
    //bare_atom_sigmaminus_V_B << 0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 1,
    //    0, 0, 0, 0;
    //bare_atom_inversion_G_H = Dense::Zero( 4, 4 );
    //bare_atom_inversion_G_H << -1, 0, 0, 0,
    //    0, 1, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0;
    //bare_atom_inversion_G_V = Dense::Zero( 4, 4 );
    //bare_atom_inversion_G_V << -1, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 1, 0,
    //    0, 0, 0, 0;
    //bare_atom_inversion_H_B = Dense::Zero( 4, 4 );
    //bare_atom_inversion_H_B << 0, 0, 0, 0,
    //    0, -1, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 1;
    //bare_atom_inversion_V_B = Dense::Zero( 4, 4 );
    //bare_atom_inversion_V_B << 0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, -1, 0,
    //    0, 0, 0, 1;
    //bare_atom_inversion_G_B = Dense::Zero( 4, 4 );
    //bare_atom_inversion_G_B << -1, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 1;
    //bare_atom_sigmaminus_G_B = Dense::Zero( 4, 4 );
    //bare_atom_sigmaminus_G_B << 0, 0, 0, 1,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0,
    //    0, 0, 0, 0;
    //
    //// Photon operators
    //bare_photon_create_H = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number );
    //bare_photon_annihilate_H = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number );
    //bare_photon_n_H = bare_photon_create_H * bare_photon_annihilate_H;
    //bare_photon_create_V = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number );
    //bare_photon_annihilate_V = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number );
    //bare_photon_n_V = bare_photon_create_V * bare_photon_annihilate_V;
    //
    //// Expanding both states
    //Log::L2( "Expanding single state matrices...\n" );
    //atom_state_ground = tensor( m_base1, m_base2, bare_atom_state_ground ).sparseView();
    //atom_state_biexciton = tensor( m_base1, m_base2, bare_atom_state_biexciton ).sparseView();
    //atom_state_H = tensor( m_base1, m_base2, bare_atom_state_H ).sparseView();
    //atom_state_V = tensor( m_base1, m_base2, bare_atom_state_V ).sparseView();
    //atom_sigmaplus_G_H = tensor( m_base1, m_base2, bare_atom_sigmaplus_G_H ).sparseView();
    //atom_sigmaminus_G_H = tensor( m_base1, m_base2, bare_atom_sigmaminus_G_H ).sparseView();
    //atom_sigmaplus_H_B = tensor( m_base1, m_base2, bare_atom_sigmaplus_H_B ).sparseView();
    //atom_sigmaminus_H_B = tensor( m_base1, m_base2, bare_atom_sigmaminus_H_B ).sparseView();
    //atom_sigmaplus_G_V = tensor( m_base1, m_base2, bare_atom_sigmaplus_G_V ).sparseView();
    //atom_sigmaminus_G_V = tensor( m_base1, m_base2, bare_atom_sigmaminus_G_V ).sparseView();
    //atom_sigmaplus_V_B = tensor( m_base1, m_base2, bare_atom_sigmaplus_V_B ).sparseView();
    //atom_sigmaminus_V_B = tensor( m_base1, m_base2, bare_atom_sigmaminus_V_B ).sparseView();
    //atom_inversion_G_H = tensor( m_base1, m_base2, bare_atom_inversion_G_H ).sparseView();
    //atom_inversion_G_V = tensor( m_base1, m_base2, bare_atom_inversion_G_V ).sparseView();
    //atom_inversion_H_B = tensor( m_base1, m_base2, bare_atom_inversion_H_B ).sparseView();
    //atom_inversion_V_B = tensor( m_base1, m_base2, bare_atom_inversion_V_B ).sparseView();
    //atom_inversion_G_B = tensor( m_base1, m_base2, bare_atom_inversion_G_B ).sparseView();
    //atom_sigmaminus_G_B = tensor( m_base1, m_base2, bare_atom_sigmaminus_G_B ).sparseView();
    //
    //photon_create_H = tensor( bare_photon_create_H, m_base2, m_base3 ).sparseView();
    //photon_create_V = tensor( m_base1, bare_photon_create_V, m_base3 ).sparseView();
    //photon_annihilate_H = tensor( bare_photon_annihilate_H, m_base2, m_base3 ).sparseView();
    //photon_annihilate_V = tensor( m_base1, bare_photon_annihilate_V, m_base3 ).sparseView();
    //photon_n_H = tensor( bare_photon_n_H, m_base2, m_base3 ).sparseView();
    //photon_n_V = tensor( m_base1, bare_photon_n_V, m_base3 ).sparseView();
    //
    //// Compressing
    //atom_state_ground.makeCompressed();
    //atom_state_biexciton.makeCompressed();
    //atom_state_H.makeCompressed();
    //atom_state_V.makeCompressed();
    //atom_sigmaplus_G_H.makeCompressed();
    //atom_sigmaminus_G_H.makeCompressed();
    //atom_sigmaplus_H_B.makeCompressed();
    //atom_sigmaminus_H_B.makeCompressed();
    //atom_sigmaplus_G_V.makeCompressed();
    //atom_sigmaminus_G_V.makeCompressed();
    //atom_sigmaplus_V_B.makeCompressed();
    //atom_sigmaminus_V_B.makeCompressed();
    //atom_inversion_G_H.makeCompressed();
    //atom_inversion_G_V.makeCompressed();
    //atom_inversion_H_B.makeCompressed();
    //atom_inversion_V_B.makeCompressed();
    //atom_inversion_G_B.makeCompressed();
    //atom_sigmaminus_G_B.makeCompressed();
    //
    //photon_create_H.makeCompressed();
    //photon_create_V.makeCompressed();
    //photon_annihilate_H.makeCompressed();
    //photon_annihilate_V.makeCompressed();
    //photon_n_H.makeCompressed();
    //photon_n_V.makeCompressed();
    //
    //// Projector Matrices
    //projector_atom_sigmaplus_G_H = project_matrix_sparse( atom_sigmaplus_G_H );
    //projector_atom_sigmaminus_G_H = project_matrix_sparse( atom_sigmaminus_G_H );
    //projector_atom_sigmaplus_H_B = project_matrix_sparse( atom_sigmaplus_H_B );
    //projector_atom_sigmaminus_H_B = project_matrix_sparse( atom_sigmaminus_H_B );
    //projector_atom_sigmaplus_G_V = project_matrix_sparse( atom_sigmaplus_G_V );
    //projector_atom_sigmaminus_G_V = project_matrix_sparse( atom_sigmaminus_G_V );
    //projector_atom_sigmaplus_V_B = project_matrix_sparse( atom_sigmaplus_V_B );
    //projector_atom_sigmaminus_V_B = project_matrix_sparse( atom_sigmaminus_V_B );
    //projector_photon_create_H = project_matrix_sparse( photon_create_H );
    //projector_photon_annihilate_H = project_matrix_sparse( photon_annihilate_H );
    //projector_photon_create_V = project_matrix_sparse( photon_create_V );
    //projector_photon_annihilate_V = project_matrix_sparse( photon_annihilate_V );

    Sparse H_0_old, H_I_old, H_old, H_used_old;
    // All possible Hamiltonions
    // H_0
    //H_0_old = p.p_omega_atomic_G_H * atom_state_H + p.p_omega_atomic_G_V * atom_state_V + p.p_omega_atomic_B * atom_state_biexciton + p.p_omega_cavity_H * photon_n_H + p.p_omega_cavity_V * photon_n_V;
    //H_0 = p.p_omega_atomic_G_H / 2.0 * atom_inversion_G_H + p.p_omega_atomic_G_V / 2.0 * atom_inversion_G_V + p.p_omega_atomic_H_B / 2.0 * atom_inversion_H_B + p.p_omega_atomic_V_B / 2.0 * atom_inversion_V_B + p.p_omega_cavity_H * photon_n_H + p.p_omega_cavity_V * photon_n_V;
    // H_I
    // RWA
    //if ( p.numerics_use_rwa ) {
    //    Log::L2( "using RWA... " );
    //    H_I_old = p.p_omega_coupling * ( atom_sigmaplus_G_H * photon_annihilate_H + atom_sigmaminus_G_H * photon_create_H + atom_sigmaplus_G_V * photon_annihilate_V + atom_sigmaminus_G_V * photon_create_V + atom_sigmaplus_H_B * photon_annihilate_H + atom_sigmaminus_H_B * photon_create_H + atom_sigmaplus_V_B * photon_annihilate_V + atom_sigmaminus_V_B * photon_create_V );
    //}
    //// non RWA
    //if ( !p.numerics_use_rwa ) {
    //    Log::L2( "NOT using RWA... " );
    //    H_I_old = p.p_omega_coupling * ( atom_sigmaplus_G_H * photon_create_H + atom_sigmaplus_G_H * photon_annihilate_H + atom_sigmaminus_G_H * photon_create_H + atom_sigmaminus_G_H * photon_annihilate_H + atom_sigmaplus_G_V * photon_create_V + atom_sigmaplus_G_V * photon_annihilate_V + atom_sigmaminus_G_V * photon_create_V + atom_sigmaminus_G_V * photon_annihilate_V + atom_sigmaplus_H_B * photon_create_H + atom_sigmaplus_H_B * photon_annihilate_H + atom_sigmaminus_H_B * photon_create_H + atom_sigmaminus_H_B * photon_annihilate_H + atom_sigmaplus_V_B * photon_create_V + atom_sigmaplus_V_B * photon_annihilate_V + atom_sigmaminus_V_B * photon_create_V + atom_sigmaminus_V_B * photon_annihilate_V );
    //}
    //// H
    //H_old = H_0_old + H_I_old;
    //// Interaction picture
    //if ( p.numerics_use_interactionpicture ) {
    //    Log::L2( "using interaction picture...\n" );
    //    H_used_old = H_I_old;
    //}
    //if ( !p.numerics_use_interactionpicture ) {
    //    Log::L2( "NOT using interaction picture...\n" );
    //    H_used_old = H_old;
    //}
    Log::L2( "Hamiltonoperator done!\n" );
    // rho, experimental: start with coherent state. in this case, always start in ground state.
    if ( !p.startCoherent || true ) {
        p.p_initial_state = base_index_map[p.p_initial_state_s];
        Log::L2( "Setting initial rho as pure state with rho_0 = {}... \n", p.p_initial_state );
        rho.coeffRef( p.p_initial_state, p.p_initial_state ) = 1;
        std::cout << "initial state = " << p.p_initial_state_s << " -> " << p.p_initial_state << std::endl;
    } else {
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
    }
    Log::L2( "Done!\n" );
    return true;
}

void OperatorMatrices::outputOperators( Parameters &p ) {
    if ( p.output_operators > 0 ) {
        std::ostringstream out;
        Eigen::IOFormat CleanFmt( 4, 0, ", ", "\n", "[", "]" );
        if ( p.output_operators > 1 ) {
            out << "General Operators:\natom_state_ground\n"
                << Dense( atom_state_ground ).format( CleanFmt ) << std::endl;
            out << "atom_state_biexciton\n"
                << Dense( atom_state_biexciton ).format( CleanFmt ) << std::endl;
            out << "atom_state_H\n"
                << Dense( atom_state_H ).format( CleanFmt ) << std::endl;
            out << "atom_state_V\n"
                << Dense( atom_state_V ).format( CleanFmt ) << std::endl;
            out << "atom_sigmaplus_G_H\n"
                << Dense( atom_sigmaplus_G_H ).format( CleanFmt ) << std::endl;
            out << "atom_sigmaminus_G_H\n"
                << Dense( atom_sigmaminus_G_H ).format( CleanFmt ) << std::endl;
            out << "atom_sigmaplus_H_B\n"
                << Dense( atom_sigmaplus_H_B ).format( CleanFmt ) << std::endl;
            out << "atom_sigmaminus_H_B\n"
                << Dense( atom_sigmaminus_H_B ).format( CleanFmt ) << std::endl;
            out << "atom_sigmaplus_G_V\n"
                << Dense( atom_sigmaplus_G_V ).format( CleanFmt ) << std::endl;
            out << "atom_sigmaminus_G_V\n"
                << Dense( atom_sigmaminus_G_V ).format( CleanFmt ) << std::endl;
            out << "atom_sigmaplus_V_B\n"
                << Dense( atom_sigmaplus_V_B ).format( CleanFmt ) << std::endl;
            out << "atom_sigmaminus_V_B\n"
                << Dense( atom_sigmaminus_V_B ).format( CleanFmt ) << std::endl;
            out << "atom_inversion_G_H\n"
                << Dense( atom_inversion_G_H ).format( CleanFmt ) << std::endl;
            out << "atom_inversion_G_V\n"
                << Dense( atom_inversion_G_V ).format( CleanFmt ) << std::endl;
            out << "atom_inversion_H_B\n"
                << Dense( atom_inversion_H_B ).format( CleanFmt ) << std::endl;
            out << "atom_inversion_V_B\n"
                << Dense( atom_inversion_V_B ).format( CleanFmt ) << std::endl;
            out << "atom_inversion_G_B\n"
                << Dense( atom_inversion_G_B ).format( CleanFmt ) << std::endl;
            out << "photon_create_H\n"
                << Dense( photon_create_H ).format( CleanFmt ) << std::endl;
            out << "photon_create_V\n"
                << Dense( photon_create_V ).format( CleanFmt ) << std::endl;
            out << "photon_annihilate_H\n"
                << Dense( photon_annihilate_H ).format( CleanFmt ) << std::endl;
            out << "photon_annihilate_V\n"
                << Dense( photon_annihilate_V ).format( CleanFmt ) << std::endl;
            out << "photon_n_H\n"
                << Dense( photon_n_H ).format( CleanFmt ) << std::endl;
            out << "photon_n_V\n"
                << Dense( photon_n_V ).format( CleanFmt ) << std::endl;
        }
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