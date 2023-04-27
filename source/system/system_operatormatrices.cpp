#include "system/operatormatrices.h"
#include "misc/helperfunctions_matrix.h"
// #include "system/operatormatrices_text.h"

using namespace QDLC;
/**
 * TODO: function add_state, add_resontor die dann ket =, bra = , bname = , name^T = usw übernehmen!
 */

// TODO for this: change initial state to initial state matrix element --> |X> to |X><X|, or with new syntax: :X; to :X;X:
// This way, we can also initialize the off-diagonals, like :X;G:.
// Function to parse states like "|><|" aka :X;X:
// add scalings -> a:X;X: or i:X;X: or ai:X;X:

// Separates amp:state; or amp:state;state: strings and returns the isolated Scalar amplitude and the string state
std::tuple<Scalar, std::string> separate_state(const std::string& input) {
    // Edge Case; No Amplitude
    if (input.front() == ':')
        return {1,input};
    // Split Amp and State at first ":"
    auto [amp, state] = QDLC::String::split_pair(input, ":");
    Scalar res = 1;
    if (amp.back() == 'i') {
        amp.pop_back();
        res = 1i;
    }
    res *= std::stod(amp);
    return {res, state};

}

OperatorMatrices::OperatorMatrices( Parameters &p ) {
    Timer &timer_operatormatrices = Timers::create( "Operator Matrices", true, false ).start();
    Log::L2( "[System-OperatorMatrices] Generating operator matrices...\n" );
    if ( !generate_operators( p ) ) {
        Log::L2( "[System-OperatorMatrices] Generating operator matrices failed! Exitting program...\n" );
        Log::Logger::close();
        exit( EXIT_FAILURE );
    }
    timer_operatormatrices.end();
    Log::L2( "[System-OperatorMatrices] Generating operator matrices was successful. Elapsed time is {}ms\n", timer_operatormatrices.getWallTime( Timers::MILLISECONDS ) );
}

bool OperatorMatrices::generate_operators( Parameters &p ) {
    output_format = Eigen::IOFormat( 4, 0, ", ", "\n", "[", "]" );
    const auto output_operators = p.output_dict.contains( "operators" );
    // Zeroing Hamiltons (redundant at this point)
    Log::L2( "[System-OperatorMatrices] Creating operator matrices, dimension = {}, creating base matrices...\n", p.maxStates );
    // Generate Electronic and Photonic Base States and Self-Hilbert Matrices
    // Electronic
    Log::L2( "[System-OperatorMatrices] Creating Self-Hilbert Electronic states...\n" );
    std::vector<std::string> base_electronic;
    base_selfhilbert.emplace_back( Dense::Identity( p.input_electronic.size(), p.input_electronic.size() ) );
    int index = 0;
    auto maximum_electronic_levels = p.input_electronic.size();
    for ( auto &[name, data] : p.input_electronic ) {
        base_electronic.emplace_back( name );
        auto &state = el_states[name];
        state.ket = Dense::Zero( maximum_electronic_levels, 1 );
        state.ket( index ) = 1;
        state.bra = state.ket.transpose();
        state.self_hilbert = state.ket * state.bra;
        state.base = 0;
        state.energy = data.get( "Energy" );
        state.name = name;
        // Increase Electronic Index
        index++;
        if ( output_operators )
            Log::L2( "[System-OperatorMatrices] Added electronic state {}\n", name );
    }

    // Photonic
    Log::L2( "[System-OperatorMatrices] Creating Self-Hilbert Photonic states...\n" );
    std::vector<std::vector<std::string>> base_photonic( p.input_photonic.size() );
    int curcav = 0;
    for ( auto &[name, data] : p.input_photonic ) {
        auto max_photons = int( data.get( "MaxPhotons" ) );
        base_selfhilbert.emplace_back( Dense::Identity( max_photons + 1, max_photons + 1 ) );
        for ( int n = 0; n <= max_photons; n++ )
            base_photonic[curcav].emplace_back( std::to_string( n ) + name );
        auto &state = ph_states[name];
        state.self_hilbert = create_photonic_operator<Dense>( QDLC::PhotonicOperator::State, max_photons );
        state.energy = data.get( "Energy" );
        state.name = name;
        // Increase Cavity Index
        curcav++;
        state.base = curcav;
        if ( output_operators )
            Log::L2( "[System-OperatorMatrices] Added photonic resonator {}\n", name );
    }

    // Tensor all matrices into the total Hilbert Space.
    // String Base Tensors
    auto subb = base_photonic;
    subb.insert( subb.begin(), base_electronic );
    base = QDLC::Matrix::tensor( subb );

    // Create Universal Identity Matrix
    identity = Dense::Identity( base.size(), base.size() ).sparseView();

    Log::L2( "[System-OperatorMatrices] Evaluating Kronecker Delta for the electronic states...\n" );
    // Sparse total Hilbert tensors
    for ( auto &[name, data] : el_states ) {
        auto current = base_selfhilbert;
        current.front() = data.self_hilbert; // Electronic states have basis 0 for now. Maybe change when TODO: add tensoring of multiple electronic bases
        data.hilbert = QDLC::Matrix::tensor( current ).sparseView();
        data.projector = QDLC::Matrix::sparse_projector( data.hilbert );
    }
    Log::L2( "[System-OperatorMatrices] Evaluating Kronecker Delta for the photonic states...\n" );
    for ( auto &[name, data] : ph_states ) {
        auto current = base_selfhilbert;
        current[data.base] = data.self_hilbert;
        data.hilbert = QDLC::Matrix::tensor( current ).sparseView();
        data.projector = QDLC::Matrix::sparse_projector( data.hilbert );
    }
    // TODO: ggf nur mit key "up" speichern, dann "up", "down" namen und matritzen dazu.
    Log::L2( "[System-OperatorMatrices] Creating Electronic Transition Matrices...\n" );
    // Generate Transition Matrices
    for ( const auto &[name, data] : el_states ) {
        for ( const auto &trans_to : p.input_electronic[name].string_v["CoupledTo"] ) {
            if ( trans_to == "-" ) continue;
            std::string transition = name + p.transition_delimiter + trans_to;
            el_transitions[transition].self_hilbert = data.ket * el_states[trans_to].bra;
            el_transitions[transition].direction = -1; // DOWN / sigma_-
            el_transitions[transition].energy = std::abs( p.input_electronic[trans_to].get( "Energy" ) - p.input_electronic[name].get( "Energy" ) );
            std::string transition_transposed = trans_to + p.transition_delimiter + name;
            el_transitions[transition_transposed].self_hilbert = el_states[trans_to].ket * data.bra;
            el_transitions[transition_transposed].direction = 1; // UP / sigma_+
            el_transitions[transition_transposed].energy = std::abs( p.input_electronic[trans_to].get( "Energy" ) - p.input_electronic[name].get( "Energy" ) );
            // Names. Should be equal to the key.
            el_transitions[transition].name = transition;
            el_transitions[transition_transposed].name = transition_transposed;
            // Transposed names
            el_transitions[transition].name_transposed = transition_transposed;
            el_transitions[transition_transposed].name_transposed = transition;
            // From-To Notes
            el_transitions[transition].from = name;
            el_transitions[transition].to = trans_to;
            el_transitions[transition_transposed].to = name;
            el_transitions[transition_transposed].from = trans_to;
            if ( output_operators ) {
                Log::L2( "[System-OperatorMatrices] Added electronic transition {}\n", transition );
                Log::L2( "[System-OperatorMatrices] Added electronic transition {}\n", transition_transposed );
            }
        }
    }
    Log::L2( "[System-OperatorMatrices] Creating Photonic Transition Matrices...\n" );
    curcav = 1;
    for ( auto &[name, data] : p.input_photonic ) {
        auto max_photons = int( data.get( "MaxPhotons" ) );
        auto &state_b = ph_transitions[name + "b"];
        state_b.self_hilbert = create_photonic_operator<Dense>( QDLC::PhotonicOperator::Annihilate, max_photons );
        state_b.base = curcav;
        state_b.direction = -1;
        state_b.energy = data.get( "Energy" );
        auto &state_bd = ph_transitions[name + "bd"];
        state_bd.self_hilbert = create_photonic_operator<Dense>( QDLC::PhotonicOperator::Create, max_photons );
        state_bd.base = curcav++;
        state_bd.direction = 1;
        state_bd.energy = data.get( "Energy" );
        // Names. Should be equal to the key
        state_b.name = name + "b";
        state_bd.name = name + "bd";
        state_b.name_transposed = name + "bd";
        state_bd.name_transposed = name + "b";
        // From-To. For the photonic states, from and to just point to the corresponding cavity
        state_b.from = name;
        state_b.to = name;
        state_bd.from = name;
        state_bd.to = name;
        if ( output_operators ) {
            Log::L2( "[System-OperatorMatrices] Added photonic transition {}\n", name + "b" );
            Log::L2( "[System-OperatorMatrices] Added photonic transition {}\n", name + "bd" );
        }
    }
    Log::L2( "[System-OperatorMatrices] Creating Electronic Projector Matrices...\n" );
    // Tensor remaining transition
    for ( auto &[name, data] : el_transitions ) {
        auto current = base_selfhilbert;
        current.front() = data.self_hilbert;
        data.hilbert = QDLC::Matrix::tensor( current ).sparseView();
        if ( output_operators ) {
            Log::L2( "[System-OperatorMatrices] Electronic Transition Matrix for Transition {} in self-Hilbert space:\n{}\n", name, Dense( data.self_hilbert ).format( output_format ) );
            Log::L2( "[System-OperatorMatrices] Electronic Transition Matrix for Transition {} in total-Hilbert space:\n{}\n", name, Dense( data.hilbert ).format( output_format ) );
        }
        data.projector = QDLC::Matrix::sparse_projector( data.hilbert );
    }
    Log::L2( "[System-OperatorMatrices] Creating Photonic Projector Matrices...\n" );
    for ( auto &[name, data] : ph_transitions ) {
        auto current = base_selfhilbert;
        current[data.base] = data.self_hilbert;
        data.hilbert = QDLC::Matrix::tensor( current ).sparseView();
        if ( output_operators ) {
            Log::L2( "[System-OperatorMatrices] Cavity Transition Matrix for Transition {} in self-Hilbert space:\n{}\n", name, Dense( data.self_hilbert ).format( output_format ) );
            Log::L2( "[System-OperatorMatrices] Cavity Transition Matrix for Transition {} in total-Hilbert space:\n{}\n", name, Dense( data.hilbert ).format( output_format ) );
        }
        data.projector = QDLC::Matrix::sparse_projector( data.hilbert );
    }

    Log::L2( "[System-OperatorMatrices] Creating Hilbert Space Index Matrices...\n" );
    // Create total Hilbert space indices. This routine assumes the individual bases are added in ascending order, e.g. 0,1,2,3...:
    for ( int i = 0; i < base_selfhilbert.size(); i++ ) {
        auto current = base_selfhilbert;
        for ( int k = 0; k < current[i].rows(); k++ )
            for ( int j = 0; j < current[i].cols(); j++ )
                current[i]( k, j ) = ( k + 1 ) + 1.0i * ( j + 1 );

        base_hilbert_index.emplace_back( QDLC::Matrix::tensor( current ) );
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
            std::string transition = pulse.second.string_v["CoupledTo"][i]; // sigma_-
            if ( el_transitions.contains( transition ) ) {
                Log::L2( "[System-OperatorMatrices] Electronic Pulse transition {} added to pulse {}...\n", transition, i );
                pulsemat += el_transitions[transition].hilbert;
                std::string transition_transposed = el_transitions[transition].name_transposed; // sigma_+
                pulsemat_star += el_transitions[transition_transposed].hilbert;
            } else if ( std::isupper( transition.front() ) ) { // This will require the electronic states to be capitalized and not start with a numerical.
                Log::L2( "[System-OperatorMatrices] Electronic Pulse transition {} is not in the list of allowed electronic transitions, recreating transition matrices and adding to pulse {}\n", transition, i );
                auto [from, to] = QDLC::String::split_pair( transition, p.transition_delimiter );
                auto ket1 = el_states[from].ket;
                auto bra1 = el_states[from].bra;
                auto ket2 = el_states[to].ket;
                auto bra2 = el_states[to].bra;
                auto current = base_selfhilbert;
                current.front() = ket1 * bra2;
                Sparse transition_hilbert = QDLC::Matrix::tensor( current ).sparseView();
                std::string transition_transposed = QDLC::String::split_and_reverse( transition, p.transition_delimiter ); // sigma_+
                {
                    auto &state = extra_transitions[transition];
                    state.ket = ket1;
                    state.bra = state.ket.transpose();
                    state.self_hilbert = state.ket * state.bra;
                    state.base = 0;
                    state.hilbert = transition_hilbert;
                    state.name = transition;
                    state.name_transposed = transition_transposed;
                    state.from = from;
                    state.to = to;
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
                    state.name = transition_transposed;
                    state.name_transposed = transition;
                    state.from = to;
                    state.to = from;
                }

                pulsemat += transition_hilbert;
                pulsemat_star += transition_transposed_hilbert;
            } else {
                Log::L2( "[System-OperatorMatrices] Pulse transition {} is cavity...\n", transition );
                // pulsemat += ph_transitions[transition + "b"].hilbert;
                // pulsemat_star += ph_transitions[transition + "bd"].hilbert;
                const auto pindx = pulse_mat_cavity_cache.size();
                pulse_mat_cavity_cache[pindx - 2] += ph_transitions[transition + "b"].hilbert;
                pulse_mat_cavity_cache[pindx - 1] += ph_transitions[transition + "bd"].hilbert;
            }
        }
        pulse_mat.emplace_back( pulsemat );      // sigma_-
        pulse_mat.emplace_back( pulsemat_star ); // sigma_+
        const auto pindx = pulse_mat_cavity_cache.size();

        if ( output_operators )
            Log::L2( "[System-OperatorMatrices] Added Pulse Matrix for Pulse {} in total-Hilbert space (normal+transposed):\n{}\n", pulse.first, Dense( pulsemat + pulse_mat_cavity_cache[pindx - 2] + pulsemat_star + pulse_mat_cavity_cache[pindx - 1] ).format( output_format ) );
    }

    // TODO: chirp cavity
    Log::L2( "[System-OperatorMatrices] Creating Chirp Cache Matrices...\n" );
    // 1 matrix is generated per chirp
    for ( auto &chirp : p.input_chirp ) {
        auto chirpmat = Sparse( base.size(), base.size() );
        auto chirpmat_star = Sparse( base.size(), base.size() );
        for ( int i = 0; i < chirp.second.string_v["CoupledTo"].size(); i++ ) {
            std::string state = chirp.second.string_v["CoupledTo"][i];
            chirpmat += chirp.second.get( "AmpFactor", i ) * el_states[state].hilbert; // TODO: chirp cavity similar to pulse
        }
        chirp_mat.emplace_back( chirpmat );
    }

    Log::L2( "[System-OperatorMatrices] Creating Base Index Maps...\n" );
    // Create Index Map
    index = 0;
    std::stringstream ss;
    for ( auto &b : base ) {
        b = ":" + b + ";";
        base_index_map[b] = index++;
        ss << b << ", ";
    }
    Log::L2( "System Base ({}): {}\n", base.size(), ss.str() );
    p.maxStates = int( base.size() );

    if ( p.numerics_groundstate_string.front() == ':' ) {
        if ( not base_index_map.contains( p.numerics_groundstate_string ) ) {
            Log::L1( "[ERROR] Provided Groundstate {} is unknown! Using Groundstate index 0!\n", p.numerics_groundstate_string );
            p.numerics_groundstate = 0;
        } else {
            p.numerics_groundstate = base_index_map[p.numerics_groundstate_string];
            Log::L2( "[System-OperatorMatrices] Setting Groundstate index from groundstate string {} to {}\n", p.numerics_groundstate_string, p.numerics_groundstate );
        }
    } else {
        p.numerics_groundstate = std::stoi( p.numerics_groundstate_string.c_str() );
        Log::L2( "[System-OperatorMatrices] Setting Groundstate index to {}\n", p.numerics_groundstate );
    }

    H_0 = Sparse( p.maxStates, p.maxStates );
    H_I = Sparse( p.maxStates, p.maxStates );
    H_used = Sparse( p.maxStates, p.maxStates );
    rho = Sparse( p.maxStates, p.maxStates );

    // Analytical Time Trafo Matrix
    // TODO: Puls mit in trafo? Trafohamilton dann H_el + H_phot + H_pulse
    Log::L2( "[System-OperatorMatrices] Creating Analytical Time Transformation Cache Matrix...\n" );
    timetrafo_cachematrix = Dense::Zero( base.size(), base.size() );
    for ( const auto &[name_i, index_i] : base_index_map ) {
        std::string is = name_i;
        is.pop_back();
        is.erase( is.begin() );
        auto i = QDLC::String::splitline( is, ':' );
        for ( const auto &[name_j, index_j] : base_index_map ) {
            Scalar val = 0;
            std::string js = name_j;
            js.pop_back();
            js.erase( js.begin() );
            auto j = QDLC::String::splitline( js, ':' );
            val += p.input_electronic[i.front()].get( "Energy" ) - p.input_electronic[j.front()].get( "Energy" );
            for ( int a = 1; a < i.size(); a++ ) {
                std::string ai = i[a];
                std::string aj = j[a];
                std::string sni = "";
                std::string snj = "";
                while ( std::isdigit( ai.at( sni.size() ) ) ) {
                    sni += ai.substr( 0, 1 );
                }
                while ( std::isdigit( aj.at( snj.size() ) ) ) {
                    snj += aj.substr( 0, 1 );
                }
                ai = ai.substr( sni.size() );
                aj = aj.substr( snj.size() );
                std::string mode = ai;
                double ni = std::stod( sni.c_str() );
                double nj = std::stod( snj.c_str() );
                // std::cout << iss.first << ", " << jss.first << " --> " << ni << " " << nj << ", modes = " << mode << std::endl;
                val += p.input_photonic[mode].get( "Energy" ) * ( ni - nj );
            }
            timetrafo_cachematrix( index_i, index_j ) = 1.0i * val;
        }
    }
    if ( output_operators )
        Log::L2( "[System-OperatorMatrices] Timetrafo Cachematrix:\n{}\n", timetrafo_cachematrix.format( output_format ) );

    Log::L2( "[System-OperatorMatrices] Creating Polaron Cache Matrices...\n" );
    // Precalculate Polaron Matrices.
    polaron_factors.emplace_back( base.size(), base.size() );
    // Transition is always |0><1| (annihilator), hence reversed is |1><0| (creator)
    for ( auto &[mode, param] : p.input_photonic ) { // TODO: das hier auf input matritzen ändern kekw.
        int i = 0;
        for ( const auto &transition : param.string_v["CoupledTo"] ) {
            Log::L2( "[System-PME] Adding Polaron Cavity Transition |{}><{}|b_{}\n", el_transitions[transition].to, el_transitions[transition].from, mode );
            auto transition_transposed = el_transitions[transition].name_transposed;
            polaron_factors[0] += el_transitions[transition_transposed].hilbert * p.p_omega_coupling * param.get( "CouplingScaling", i ) * ph_transitions[mode + "b"].hilbert;
            i++;
        }
    }
    for ( int current = 0; auto const &[mode, param] : p.input_pulse ) {
        // auto temp = Sparse( base.size(), base.size() );
        // for ( const auto &transition : param.string_v["CoupledTo"] ) {
        //     polaron_pulse_factors_explicit_time.emplace_back( base.size(), base.size() );
        //     if ( el_transitions.contains( transition ) ) {
        //         Log::L2( "[System-PME] Adding Polaron Pulse Transition |{}><{}|Omega_{}\n", el_transitions[transition].to, el_transitions[transition].from, mode );
        //         std::string transition_transposed = el_transitions[transition].name_transposed; // QDLC::String::split_and_reverse( transition, p.transition_delimiter );
        //         temp += el_transitions[transition_transposed].hilbert;
        //         polaron_pulse_factors_explicit_time.back() += el_transitions[transition_transposed].projector;
        //         Log::L2( "[System-PME] Transition Matrix:\n{}\n", Dense( el_transitions[transition_transposed].hilbert ).format( output_format ) );
        //     } else if ( extra_transitions.contains( transition ) ) {
        //         Log::L2( "[System-PME] Adding Polaron Pulse Transition |{}><{}|Omega_{} which is a not-allowed electronic transition.\n", extra_transitions[transition].to, extra_transitions[transition].from, mode );
        //         std::string transition_transposed = extra_transitions[transition].name_transposed; // QDLC::String::split_and_reverse( transition, p.transition_delimiter );
        //         temp += extra_transitions[transition_transposed].hilbert;
        //         polaron_pulse_factors_explicit_time.back() += QDLC::Matrix::sparse_projector( extra_transitions[transition_transposed].hilbert );
        //         Log::L2( "[System-PME] Transition Matrix:\n{}\n", Dense( extra_transitions[transition_transposed].hilbert ).format( output_format ) );
        //     } else {
        //         Log::L2( "[System-OperatorMatrices] Pulse transition {} is cavity...\n", transition );
        //         temp += ph_transitions[transition + "bd"].hilbert;
        //         polaron_pulse_factors_explicit_time.back() += ph_transitions[transition + "bd"].projector;
        //         Log::L2( "[System-PME] Transition Matrix:\n{}\n", Dense( ph_transitions[transition + "bd"].hilbert ).format( output_format ) );
        //     }
        // }
        // Log::L2( "[System-PME] Total Polaron Factor:\n{}\n", Dense( temp ).format( output_format ) );
        polaron_factors.emplace_back( pulse_mat[2 * current + 1] );
        polaron_pulse_factors_explicit_time.emplace_back( QDLC::Matrix::sparse_projector( pulse_mat[2 * current + 1] ) );
        current++;
    }

    Log::L2( "[System-OperatorMatrices] Creating Hamiltonoperator...\n" );
    // Generate Self Action Hamilton H_0:
    H_0 = Sparse( base.size(), base.size() );
    for ( const auto &[name, data] : el_states ) {
        double energy = p.input_electronic[name].get( "Energy" );
        H_0 += energy * data.hilbert;
    }
    for ( const auto &[name, data] : ph_states ) {
        double energy = p.input_photonic[name].get( "Energy" );
        H_0 += energy * data.hilbert;
    }
    // Generate Interaction Hamilton H_I:
    H_I = Sparse( base.size(), base.size() );
    Sparse H_I_a = H_I; // RWA Part
    Sparse H_I_b = H_I; // Non RWA Part
    for ( auto &[name, data] : p.input_photonic ) {
        for ( int i = 0; i < data.string_v["CoupledTo"].size(); i++ ) {
            std::string transition = data.string_v["CoupledTo"][i];
            std::string transition_transposed = QDLC::String::split_and_reverse( transition, p.transition_delimiter );
            const auto coupling_scaling = data.get( "CouplingScaling", i );
            H_I_a += coupling_scaling * el_transitions[transition_transposed].hilbert * ph_transitions[name + "b"].hilbert;
            H_I_a += coupling_scaling * el_transitions[transition].hilbert * ph_transitions[name + "bd"].hilbert;
            H_I_b += coupling_scaling * el_transitions[transition_transposed].hilbert * ph_transitions[name + "bd"].hilbert;
            H_I_b += coupling_scaling * el_transitions[transition].hilbert * ph_transitions[name + "b"].hilbert;
        }
    }
    if ( output_operators )
        Log::L2( "[System-OperatorMatrices] Used Hamilton:\n{}\n", Dense( H_I_a ) );

    Log::L2( "[System-OperatorMatrices] Creating Path Integral Sorting Vectors...\n" );
    // The Path Integral can be partially summed by either: the electronic state index; or the electronic state coupling factor.
    // Because multiple states can have the same coupling factor, the latter will usually be faster, but at least as fast as the index summation.
    // As a default value, "factor" should be used.
    std::string sorting = "factor";
    if ( sorting == "index" ) {
        std::map<std::string, int> temp_base_indices;
        int new_index = 0;
        for ( int i = 0; i < base.size(); i++ ) {
            auto state1 = QDLC::String::split( base.at( i ).substr( 1 ), ":" ).front();
            double factor = p.input_electronic[state1].get( "PhononCoupling" ).get();
            if ( not temp_base_indices.contains( state1 ) ) {
                temp_base_indices[state1] = new_index++;
                phonon_group_index_to_coupling_value.emplace_back( factor );
            }
            phonon_hilbert_index_to_group_index.emplace_back( temp_base_indices[state1] );
        }
        phonon_group_index_to_hilbert_indices = std::vector<std::vector<int>>( temp_base_indices.size() );
        for ( int i = 0; i < base.size(); i++ ) {
            phonon_group_index_to_hilbert_indices[phonon_hilbert_index_to_group_index[i]].emplace_back( i );
        }
    } else if ( sorting == "factor" ) {
        std::map<double, int> temp_base_indices;
        int new_index = 0;
        for ( const auto &current : base ) {
            auto state1 = QDLC::String::split( current.substr( 1 ), ":" ).front();
            if ( state1.back() == ';' )
                state1.pop_back();
            auto factor = (double)p.input_electronic[state1].get( "PhononCoupling" ).get();
            auto index = factor; // state1;
            if ( !temp_base_indices.contains( index ) ) {
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
    if ( p.numerics_phonon_approximation_order != QDLC::PhononApproximation::PathIntegral ) {
        Sparse b_matrix = polaron_phonon_coupling_matrix.unaryExpr( [&p]( Scalar val ) { return std::pow( p.p_phonon_b.get(), 1.0 * val ); } );
        Log::L2( "[System-OperatorMatrices] Scaling H_I,Cavity and H_I,Pulse with <B> = {}\n", p.p_phonon_b );

        if ( output_operators )
            Log::L2( "[System-OperatorMatrices] <B>-Matrix:\n{}\n", Dense( b_matrix ).format( output_format ) );

        // H_I Scaled once by <B>
        H_I_a = H_I_a.cwiseProduct( b_matrix );
        H_I_b = H_I_b.cwiseProduct( b_matrix );
        // Scale Pulse Matrices for Pulse evaluations for the Hamilton Operators. Note: If the pulse matrix consists of cavity transitions, the <B>-scaling must not be applied there! This is the <B> from the <B>*H_I scaling for the pulse.
        Log::L2( "[System-PME] Scaled Pulse Transition Matrices to:\n" );
        std::ranges::for_each( pulse_mat, [&, indx = 0]( auto &mat ) mutable { mat = mat.cwiseProduct( b_matrix ); if (output_operators) Log::L2("[System-PME] Matrix {}:\n{}\n",indx,Dense(mat).format(output_format)); indx++; } );

        // Scale Polaron Cache Matrices for Chi evaluations. This is the <B> from the <B>^2 Term in front of the polaron green functions. The <B>^2 term in front of the green functions is the result from the commutator multiplication [<B>A,<B>B]=<B>^2[A,B]
        Log::L2( "[System-PME] Scaled Polaron Factor Matrices to:\n" );
        std::ranges::for_each( polaron_factors, [&, indx = 0]( auto &mat ) mutable { mat = mat.cwiseProduct( b_matrix ); if (output_operators) Log::L2("[System-PME] Matrix {}:\n{}\n",indx,Dense(mat).format(output_format)); indx++; } );

        // Add Cavity Pulse Transitions after scaling of the electronic transitions. Note that this way, the pure cavity driving is not influenced by the phonon coupling. IDK if this is correct.
        Log::L2( "[System-PME] Modified Pulse Transition Matrices to include cavity driving to:\n" );
        std::ranges::for_each( pulse_mat, [&, indx = 0]( auto &mat ) mutable { mat += pulse_mat_cavity_cache[indx]; if (output_operators) Log::L2("[System-PME] Matrix {}:\n{}\n",indx,Dense(mat).format(output_format)); indx++; } );
    }

    for ( const auto &a : phonon_group_index_to_hilbert_indices )
        Log::L2( "[System-OperatorMatrices] Phonon Group Index Vector: {}\n", a );

    Log::L2( "[System-OperatorMatrices] Creating Initial State Vector...\n" );
    // Create Initial State
    // Split starting state into superposition. States can be passed as "|...>+|...>" with amplitudes
    initial_state_vector_ket = Dense::Zero( base.size(), 1 );
    for ( auto amped_state : QDLC::String::splitline( p.p_initial_state_s, '+' ) ) {
        // State Syntax: |state> where state can be an actual system state or coherent(alpha)mode, squeezed(x,y)mode or thermal(alpha)mode
        //Scalar amp = state[0] == ':' ? 1.0 : std::stod( QDLC::String::splitline( state, ':' ).front() );
        auto [amp, state] = separate_state(amped_state);
        std::string pure_state = state.substr( state.find( ':' ) );
        // Coherent State
        if ( pure_state.find( "coherent" ) != std::string::npos ) {
            auto state_left = pure_state.substr( 0, pure_state.find( "coherent(" ) );
            auto state_right = pure_state.substr( pure_state.find( ")" ) + 1 );
            auto alpha = std::stod( pure_state.substr( state_left.size() + 9, pure_state.size() - state_left.size() - 9 - state_right.size() ).c_str() );
            auto mode = state_right.find( ":" ) != std::string::npos ? QDLC::String::splitline( state_right, ':' ).front() : state_right.substr( 0, state_right.size() - 1 );
            Log::L2( "[System-OperatorMatrices] Creating superpositioned coherent state {} for mode {} with alpha = {} and scaled amplitude {}\n", pure_state, mode, alpha, amp );
            for ( int n = 0; n <= p.input_photonic[mode].get( "MaxPhotons" ); n++ ) {
                std::string current = state_left + std::to_string( n ) + state_right;
                Log::L2( "[System-OperatorMatrices] Creating coherent substate {} ({}) with amplitude {}\n", current, base_index_map[current], QDLC::Math::getCoherent( alpha, n ) * amp );
                // Add initial state with amplitudes
                initial_state_vector_ket( base_index_map[current] ) += amp * QDLC::Math::getCoherent( alpha, n );
            }
        }
        // Squeezed State
        else if ( pure_state.find( "squeezed" ) != std::string::npos ) {
            auto state_left = pure_state.substr( 0, pure_state.find( "squeezed(" ) );
            auto state_right = pure_state.substr( pure_state.find( ")" ) + 1 );
            auto cmplx = QDLC::String::splitline( pure_state.substr( state_left.size() + 9, pure_state.size() - state_left.size() - 9 - state_right.size() ), ',' );
            auto mode = state_right.find( ":" ) != std::string::npos ? QDLC::String::splitline( state_right, ':' ).front() : state_right.substr( 0, state_right.size() - 1 );
            double r = std::stod( cmplx.front() );
            double phi = std::stod( cmplx.back() );
            Log::L2( "[System-OperatorMatrices] Creating superpositioned squeezed state {} for mode {} with r = {}, phi = {} and scaled amplitude {}\n", pure_state, mode, r, phi, amp );
            for ( int n = 0; n <= p.input_photonic[mode].get( "MaxPhotons" ); n++ ) {
                if ( n % 2 != 0 )
                    continue;
                std::string current = state_left + std::to_string( n ) + state_right;
                Log::L2( "[System-OperatorMatrices] Creating squeezed substate {} ({}) with amplitude {}\n", current, base_index_map[current], QDLC::Math::getSqueezed( r, phi, n ) * amp );
                // Add initial state with amplitudes
                initial_state_vector_ket( base_index_map[current] ) += amp * QDLC::Math::getSqueezed( r, phi, n );
            }
        }
        // Thermal State
        else if ( pure_state.find( "thermal" ) != std::string::npos ) {
            auto state_left = pure_state.substr( 0, pure_state.find( "thermal(" ) );
            auto state_right = pure_state.substr( pure_state.find( ")" ) + 1 );
            auto alpha = std::stod( pure_state.substr( state_left.size() + 8, pure_state.size() - state_left.size() - 8 - state_right.size() ).c_str() );
            auto mode = state_right.find( ":" ) != std::string::npos ? QDLC::String::splitline( state_right, ':' ).front() : state_right.substr( 0, state_right.size() - 1 );
            Log::L2( "[System-OperatorMatrices] Creating superpositioned Thermal state {} for mode {} with alpha = {} and scaled amplitude {}\n", pure_state, mode, alpha, amp );
            for ( int n = 0; n < p.input_photonic[mode].get( "MaxPhotons" ); n++ ) {
                std::string current = state_left + std::to_string( n ) + state_right;
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
    if ( p.output_dict.contains( "operators" ) ) {
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
        // if ( p.output_dict_operators == 3 )
        //    exit( 0 );
    }
}