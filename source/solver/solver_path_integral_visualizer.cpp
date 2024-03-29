#include "solver/solver_ode.h"

bool QDACC::Numerics::ODESolver::visualize_path( MatrixMain &rho0, System &s ) {
    Log::L2( "[Dot-Visualizer] Starting to output dot file to path.dot\n" );
    std::vector<QDACC::SaveState> dummy{ { rho0, 0.0 } };
    auto &fp_dot = FileOutput::add_file( "path", "dot" );
    fp_dot << "digraph G{{\ngraph [pad=\"0.5\", nodesep=\"0.1\", ranksep=\"3\", rankdir=\"TB\"]\n";

    // Cache Parameters
    auto omega_cavity_loss = s.parameters.p_omega_cavity_loss;
    auto omega_decay = s.parameters.p_omega_decay;
    auto omega_pure_dephasing = s.parameters.p_omega_pure_dephasing;
    auto H_used = s.operatorMatrices.H_used;
    auto sh = s.parameters.numerics_use_saved_hamiltons;
    auto dgl_pulse = s.pulse;
    auto dgl_chirp = s.chirp;
    s.parameters.numerics_use_saved_hamiltons = false;

    s.parameters.p_omega_cavity_loss = 0;
    s.parameters.p_omega_decay = 0;
    s.parameters.p_omega_pure_dephasing = 0;
    s.operatorMatrices.H_used = MatrixMain( H_used.rows(), H_used.cols() );
    s.pulse = std::vector<Pulse>();
    s.chirp = std::vector<Chirp>();

    // Calculate Propagator for Cavity Loss only
    s.parameters.p_omega_cavity_loss = omega_cavity_loss;
    auto propagators_cavity = calculate_propagator_vector( s, rho0.rows(), s.parameters.t_step_pathint, s.parameters.t_step_pathint, dummy );
    s.parameters.p_omega_cavity_loss = 0;

    // Calculate Propagator for Radiative Loss only
    s.parameters.p_omega_decay = omega_decay;
    auto propagators_rad = calculate_propagator_vector( s, rho0.rows(), s.parameters.t_step_pathint * 2., s.parameters.t_step_pathint, dummy );
    s.parameters.p_omega_decay = 0;

    // Calculate Propagator for Dephasing only
    s.parameters.p_omega_pure_dephasing = omega_pure_dephasing;
    auto propagators_dephasing = calculate_propagator_vector( s, rho0.rows(), s.parameters.t_step_pathint * 3., s.parameters.t_step_pathint, dummy );
    s.parameters.p_omega_pure_dephasing = 0;

    // Calculate Propagator for Hamilton only
    s.operatorMatrices.H_used = H_used;
    auto propagators_hamilton = calculate_propagator_vector( s, rho0.rows(), 0, s.parameters.t_step_pathint, dummy );
    s.operatorMatrices.H_used = MatrixMain( H_used.rows(), H_used.cols() );

    // Calculate Propagator for Pulse only
    s.pulse = dgl_pulse;
    auto propagators_pulse = calculate_propagator_vector( s, rho0.rows(), 0, s.parameters.t_step_pathint, dummy );
    s.pulse = std::vector<Pulse>();

    // Calculate Propagator for Chirp only
    s.chirp = dgl_chirp;
    auto propagators_chirp = calculate_propagator_vector( s, rho0.rows(), 0, s.parameters.t_step_pathint, dummy );
    s.chirp = std::vector<Chirp>();

    auto propagators = std::vector{ propagators_hamilton, propagators_pulse, propagators_chirp, propagators_cavity, propagators_rad, propagators_dephasing };
    auto colors = std::vector{ "black", "red", "green", "orange", "blue", "purple" };
    auto names = std::vector{ "Hamilton", "Pulse", "Chirp", "Cavity", "Radiative", "Dephasing" };

    // Reset parameters
    s.parameters.p_omega_cavity_loss = omega_cavity_loss;
    s.parameters.p_omega_decay = omega_decay;
    s.parameters.p_omega_pure_dephasing = omega_pure_dephasing;
    s.operatorMatrices.H_used = H_used;
    s.parameters.numerics_use_saved_hamiltons = sh;
    s.chirp = dgl_chirp;
    s.pulse = dgl_pulse;
    std::set<std::string> labels;
    std::map<std::string, double> total_weights;

    // Weights
    for ( size_t k = 0; k < propagators.size(); k++ ) {
        auto &propagator = propagators[k];
        auto &color = colors[k];
        for ( auto i = 0; i < propagator.size(); i++ ) {
            for ( auto j = 0; j < propagator.size(); j++ ) {
                //for ( int l = 0; l < propagator[i][j].outerSize(); ++l )
                //    for ( Sparse::InnerIterator M( propagator[i][j], l ); M; ++M ) {
                //        int i_n = M.row();
                //        int j_n = M.col();
                for ( int i_n = 0; i_n < propagator[i][j].rows(); i_n++ )
                    for ( int j_n = 0; j_n < propagator[i][j].cols(); j_n++ ) {
                        std::string name = std::format( "{},{}->{},{}", i, j, i_n, j_n );
                        if ( !( total_weights.contains( name ) ) ) total_weights[name] = 0.0;
                        total_weights[name] += std::abs( propagator[i][j].coeff( i_n, j_n ) );
                    }
            }
        }
    }
    // Output Propagator Mappings
    for ( size_t k = 0; k < propagators.size(); k++ ) {
        fp_dot << std::format( "\n{} [pos=\"0,{}!\", color=\"{}\"]\n", names[k], 5 - k, colors[k] );
        auto &propagator = propagators[k];
        auto &color = colors[k];
        for ( auto i = 0; i < propagator.size(); i++ ) {
            for ( auto j = 0; j < propagator.size(); j++ ) {
                if ( !labels.contains( std::format( "{}_0,{}_0", i, j ) ) ) {
                    fp_dot << std::format( "\"{1}_{0},{2}_{0}\" [label=\"|{3}><{4}|\" pos=\"{5},{6}!\"];\n", 0, i, j, s.operatorMatrices.base.at( i ).substr( 1, s.operatorMatrices.base.at( i ).size() - 2 ),
                                           s.operatorMatrices.base.at( j ).substr( 1, s.operatorMatrices.base.at( j ).size() - 2 ), i * propagator.size() * 2 + j * 2 + 3, 5 );
                    labels.emplace( std::format( "{}_0,{}_0", i, j ) );
                }
                //for ( int l = 0; l < propagator[i][j].outerSize(); ++l )
                //    for ( Sparse::InnerIterator M( propagator[i][j], l ); M; ++M ) {
                //        int i_n = M.row();
                //        int j_n = M.col();
                for ( int i_n = 0; i_n < propagator[i][j].rows(); i_n++ )
                    for ( int j_n = 0; j_n < propagator[i][j].cols(); j_n++ ) {
                        if ( !labels.contains( std::format( "{}_1,{}_1", i_n, j_n ) ) ) {
                            fp_dot << std::format( "\"{1}_{0},{2}_{0}\" [label=\"|{3}><{4}|\" pos=\"{5},{6}!\"];\n", 1, i_n, j_n, s.operatorMatrices.base.at( i_n ).substr( 1, s.operatorMatrices.base.at( i_n ).size() - 2 ),
                                                   s.operatorMatrices.base.at( j_n ).substr( 1, s.operatorMatrices.base.at( j_n ).size() - 2 ), i_n * propagator.size() * 2 + j_n * 2 + 3, 0 );
                            labels.emplace( std::format( "{}_1,{}_1", i_n, j_n ) );
                        }
                        auto element = propagator[i][j].coeff( i_n, j_n );
                        double stroke =
                            ( i == i_n and j == j_n and std::real( element ) == 1.0 ) ? 0.2 : 2.0; // std::min( 2.0, std::max( 0.2, std::abs( element ) / total_weights[std::format( "{},{}->{},{}", i, j, i_n, j_n )] ) );
                        double arrowsize = ( i == i_n and j == j_n and std::real( element ) == 1.0 ) ? 0.05 : 0.1;
                        fp_dot << std::format( "\"{1}_{0},{2}_{0}\"->\"{4}_{3},{5}_{3}\" [color=\"{8}\" penwidth=\"{9}\" arrowsize=\"{10}\" edgetooltip=\"{11} Value: ({6},{7})\" fontsize=\"5\"];\n", 0, i, j, 1, i_n, j_n,
                                               std::real( element ), std::imag( element ), color, stroke, arrowsize, names[k] );
                    }
            }
        }
    }

    // Path Integral Mappings
    try {
        if ( s.parameters.p_phonon_T < 0 or s.parameters.numerics_phonon_approximation_order != QDACC::PhononApproximation::PathIntegral )
            throw std::runtime_error( "Phonon Path Integral Kernel is only available for Path Integral Phonons." );
        fp_dot << std::format( "\n\"Phonon Kernel\" [pos=\"0,-1!\", color=\"dodgerblue2\"]\n" );
        // Phonon correlation function:
        if ( s.phi_vector_int.size() == 0 ) s.initialize_path_integral_functions();
        for ( size_t i = 0; i < rho0.rows(); i++ )
            for ( size_t j = 0; j < rho0.rows(); j++ )
                for ( size_t i_n = 0; i_n < rho0.rows(); i_n++ )
                    for ( size_t j_n = 0; j_n < rho0.rows(); j_n++ ) {
                        Scalar val = 0;
                        for ( int tau = 0; tau < s.parameters.p_phonon_nc; tau++ )
                            if ( s.parameters.numerics_pathint_partially_summed )
                                val += s.dgl_phonon_memory_function( tau, s.operatorMatrices.phonon_hilbert_index_to_group_index[i], s.operatorMatrices.phonon_hilbert_index_to_group_index[j],
                                                                     s.operatorMatrices.phonon_hilbert_index_to_group_index[i_n], s.operatorMatrices.phonon_hilbert_index_to_group_index[j_n] );
                            else
                                val += s.dgl_phonon_memory_function( tau, i, j, i_n, j_n );
                        if ( std::abs( val ) > 1E-15 )
                            fp_dot << std::format( "\"{1}_{0},{2}_{0}\"->\"{4}_{3},{5}_{3}\" [color=\"{8}\" penwidth=\"{9}\" arrowsize=\"{10}\" edgetooltip=\"{11} Value: ({6},{7})\" fontsize=\"5\"];\n", 0, i, j, 1, i_n, j_n,
                                                   std::real( val ), std::imag( val ), "dodgerblue2", 1.0, 0.1, "Phonon Kernel" );
                        // fmt::print( "i = {}, j = {}, id = {}, jd = {} --> converted i = {}, j = {}, id = {}, jd = {} --> {}\n", i, j, i_n, j_n, s.operatorMatrices.phonon_hilbert_index_to_group_index[i], s.operatorMatrices.phonon_hilbert_index_to_group_index[j], s.operatorMatrices.phonon_hilbert_index_to_group_index[i_n], s.operatorMatrices.phonon_hilbert_index_to_group_index[j_n], val );
                    }
    } catch ( const std::exception &e ) {
        Log::L2( "[Dot-Visualizer] Error when calculating and outputting Path Integral Kernel. Exception: {}\n", e.what() );
    }
    // Polaron Frame Mapping
    try {
        if ( s.parameters.p_phonon_T < 0 or s.parameters.numerics_phonon_approximation_order == QDACC::PhononApproximation::PathIntegral ) throw std::runtime_error( "Phonon PME Kernel is only available for PME Phonons." );
        s.initialize_polaron_frame_functions();
        fp_dot << std::format( "\n\"Polaron Mapping\" [pos=\"3,-1!\", color=\"firebrick4\"]\n" );
        std::vector<std::vector<MatrixMain>> polaron( (size_t)rho0.rows(), { (size_t)rho0.rows(), MatrixMain( rho0.rows(), rho0.rows() ) } );
        // Calculate the remaining propagators
        for ( int i = 0; i < rho0.rows(); i++ ) {
            for ( int j = 0; j < rho0.rows(); j++ ) {
                MatrixMain projector = MatrixMain( rho0.rows(), rho0.rows() );
                projector.coeffRef( i, j ) = 1;
                polaron[i][j] = s.dgl_phonons_pmeq( projector, s.parameters.p_phonon_tcutoff, dummy );
            }
        }
        for ( auto i = 0; i < polaron.size(); i++ )
            for ( auto j = 0; j < polaron.size(); j++ )
                //for ( int l = 0; l < polaron[i][j].outerSize(); ++l )
                //    for ( Sparse::InnerIterator M( polaron[i][j], l ); M; ++M ) {
                //        int i_n = M.row();
                //        int j_n = M.col();
                for ( int i_n = 0; i_n < polaron[i][j].rows(); i_n++ )
                    for ( int j_n = 0; j_n < polaron[i][j].cols(); j_n++ ) {
                        auto element = polaron[i][j].coeff( i_n, j_n );
                        double stroke =
                            ( i == i_n and j == j_n and std::real( element ) == 1.0 ) ? 0.2 : 1.0; // std::min( 2.0, std::max( 0.2, std::abs( element ) / total_weights[std::format( "{},{}->{},{}", i, j, i_n, j_n )] ) );
                        double arrowsize = ( i == i_n and j == j_n and std::real( element ) == 1.0 ) ? 0.05 : 0.1;
                        fp_dot << std::format( "\"{1}_{0},{2}_{0}\"->\"{4}_{3},{5}_{3}\" [color=\"{8}\" penwidth=\"{9}\" arrowsize=\"{10}\" edgetooltip=\"{11} Value: ({6},{7})\" fontsize=\"5\"];\n", 0, i, j, 1, i_n, j_n,
                                               std::real( element ), std::imag( element ), "firebrick4", stroke, arrowsize, "Polaron Mapping" );
                    }
    } catch ( const std::exception &e ) {
        Log::L2( "[Dot-Visualizer] Error when calculating and outputting PME Kernel.\n" );
    }
    // Always End the File
    fp_dot << "}}";
    // Close File to ensure Flushing
    fp_dot.close();
    // RIP. TODO: no system call; for now dont care
    Log::L2( "[Dot-Visualizer] Outputting dot file to path.dot\n" );
    Log::L2( "[Dot-Visualizer] dot -v -Kneato -Tsvg \"{0}path.dot\" -o \"{0}path.svg\"\n", s.parameters.working_directory );
    Log::L2( "[Dot-Visualizer] Return value: {}\n", std::system( std::format( "dot -v -Kneato -Tsvg \"{0}path.dot\" -o \"{0}path.svg\"", s.parameters.working_directory ).c_str() ) );
    Log::L2( "[Dot-Visualizer] dot -v -Kneato -Tpng \"{0}path.dot\" -o \"{0}path.png\"\n", s.parameters.working_directory );
    Log::L2( "[Dot-Visualizer] Return value: {}\n", std::system( std::format( "dot -v -Kneato -Tpng \"{0}path.dot\" -o \"{0}path.png\"", s.parameters.working_directory ).c_str() ) );
    Log::L2( "[Dot-Visualizer] dot -v -Tsvg \"{0}path.dot\" -o \"{0}path_unordered.svg\"\n", s.parameters.working_directory );
    Log::L2( "[Dot-Visualizer] Return value: {}\n", std::system( std::format( "dot -v -Tsvg \"{0}path.dot\" -o \"{0}path_unordered.svg\"", s.parameters.working_directory ).c_str() ) );
    return true;
}