#include "system/operatormatrices.h"
#include "system/operatormatrices_text.h"

OperatorMatrices::OperatorMatrices( const Parameters &p ) {
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

bool OperatorMatrices::generateOperators( const Parameters &p ) {
    // Zeroing Hamiltons (redundant at this point)
    Log::L2( "Creating operator matrices, dimension = {}\nCreating base matrices... ", p.maxStates );
    H = Sparse( p.maxStates, p.maxStates );
    H_0 = Sparse( p.maxStates, p.maxStates );
    H_I = Sparse( p.maxStates, p.maxStates );
    H_used = Sparse( p.maxStates, p.maxStates );
    rho = Sparse( p.maxStates, p.maxStates );

    // Create Base Matrices and Base Vector:
    Dense m_base1 = Dense::Identity( p.p_max_photon_number + 1, p.p_max_photon_number + 1 );
    Dense m_base2 = Dense::Identity( p.p_max_photon_number + 1, p.p_max_photon_number + 1 );
    Dense m_base3 = Dense::Identity( 4, 4 );
    std::vector<std::string> base1;
    std::vector<std::string> base2;
    std::vector<std::string> base3 = {"G", "X_H", "X_V", "B"};
    for ( int i = 0; i <= p.p_max_photon_number; i++ ) {
        base1.emplace_back( std::to_string( i ) + "_H" );
        base2.emplace_back( std::to_string( i ) + "_V" );
    }
    base = tensor( tensor( base1, base2 ), base3 );

    Log::L3( "Operator Base: (size {}) ", base.size() );
    for ( auto b : base )
        Log::L3( "|{}> ", b );
    Log::L3( "\n" );

    // Initializing bare matrices:
    Log::L2( "Initializing base matrices... " );
    // Atomic state operators
    bare_atom_state_ground = Dense::Zero( 4, 4 );
    bare_atom_state_ground << 1, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    bare_atom_state_biexciton = Dense::Zero( 4, 4 );
    bare_atom_state_biexciton << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 1;
    bare_atom_state_H = Dense::Zero( 4, 4 );
    bare_atom_state_H << 0, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    bare_atom_state_V = Dense::Zero( 4, 4 );
    bare_atom_state_V << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 0;
    bare_atom_sigmaplus_G_H = Dense::Zero( 4, 4 );
    bare_atom_sigmaplus_G_H << 0, 0, 0, 0,
        1, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    bare_atom_sigmaminus_G_H = Dense::Zero( 4, 4 );
    bare_atom_sigmaminus_G_H << 0, 1, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    bare_atom_sigmaplus_G_V = Dense::Zero( 4, 4 );
    bare_atom_sigmaplus_G_V << 0, 0, 0, 0,
        0, 0, 0, 0,
        1, 0, 0, 0,
        0, 0, 0, 0;
    bare_atom_sigmaminus_G_V = Dense::Zero( 4, 4 );
    bare_atom_sigmaminus_G_V << 0, 0, 1, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    bare_atom_sigmaplus_H_B = Dense::Zero( 4, 4 );
    bare_atom_sigmaplus_H_B << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 1, 0, 0;
    bare_atom_sigmaminus_H_B = Dense::Zero( 4, 4 );
    bare_atom_sigmaminus_H_B << 0, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 0, 0,
        0, 0, 0, 0;
    bare_atom_sigmaplus_V_B = Dense::Zero( 4, 4 );
    bare_atom_sigmaplus_V_B << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 1, 0;
    bare_atom_sigmaminus_V_B = Dense::Zero( 4, 4 );
    bare_atom_sigmaminus_V_B << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 0, 0;
    bare_atom_inversion_G_H = Dense::Zero( 4, 4 );
    bare_atom_inversion_G_H << -1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    bare_atom_inversion_G_V = Dense::Zero( 4, 4 );
    bare_atom_inversion_G_V << -1, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 0;
    bare_atom_inversion_H_B = Dense::Zero( 4, 4 );
    bare_atom_inversion_H_B << 0, 0, 0, 0,
        0, -1, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 1;
    bare_atom_inversion_V_B = Dense::Zero( 4, 4 );
    bare_atom_inversion_V_B << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, -1, 0,
        0, 0, 0, 1;
    bare_atom_inversion_G_B = Dense::Zero( 4, 4 );
    bare_atom_inversion_G_B << -1, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 1;
    bare_atom_sigmaminus_G_B = Dense::Zero( 4, 4 );
    bare_atom_sigmaminus_G_B << 0, 0, 0, 1,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

    // Photon operators
    bare_photon_create_H = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number );
    bare_photon_annihilate_H = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number );
    bare_photon_n_H = bare_photon_create_H * bare_photon_annihilate_H;
    bare_photon_create_V = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number );
    bare_photon_annihilate_V = create_photonic_operator<Dense>( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number );
    bare_photon_n_V = bare_photon_create_V * bare_photon_annihilate_V;

    // Expanding both states
    Log::L2( "Expanding single state matrices... " );
    atom_state_ground = tensor( m_base1, m_base2, bare_atom_state_ground ).sparseView();
    atom_state_biexciton = tensor( m_base1, m_base2, bare_atom_state_biexciton ).sparseView();
    atom_state_H = tensor( m_base1, m_base2, bare_atom_state_H ).sparseView();
    atom_state_V = tensor( m_base1, m_base2, bare_atom_state_V ).sparseView();
    atom_sigmaplus_G_H = tensor( m_base1, m_base2, bare_atom_sigmaplus_G_H ).sparseView();
    atom_sigmaminus_G_H = tensor( m_base1, m_base2, bare_atom_sigmaminus_G_H ).sparseView();
    atom_sigmaplus_H_B = tensor( m_base1, m_base2, bare_atom_sigmaplus_H_B ).sparseView();
    atom_sigmaminus_H_B = tensor( m_base1, m_base2, bare_atom_sigmaminus_H_B ).sparseView();
    atom_sigmaplus_G_V = tensor( m_base1, m_base2, bare_atom_sigmaplus_G_V ).sparseView();
    atom_sigmaminus_G_V = tensor( m_base1, m_base2, bare_atom_sigmaminus_G_V ).sparseView();
    atom_sigmaplus_V_B = tensor( m_base1, m_base2, bare_atom_sigmaplus_V_B ).sparseView();
    atom_sigmaminus_V_B = tensor( m_base1, m_base2, bare_atom_sigmaminus_V_B ).sparseView();
    atom_inversion_G_H = tensor( m_base1, m_base2, bare_atom_inversion_G_H ).sparseView();
    atom_inversion_G_V = tensor( m_base1, m_base2, bare_atom_inversion_G_V ).sparseView();
    atom_inversion_H_B = tensor( m_base1, m_base2, bare_atom_inversion_H_B ).sparseView();
    atom_inversion_V_B = tensor( m_base1, m_base2, bare_atom_inversion_V_B ).sparseView();
    atom_inversion_G_B = tensor( m_base1, m_base2, bare_atom_inversion_G_B ).sparseView();
    atom_sigmaminus_G_B = tensor( m_base1, m_base2, bare_atom_sigmaminus_G_B ).sparseView();

    photon_create_H = tensor( bare_photon_create_H, m_base2, m_base3 ).sparseView();
    photon_create_V = tensor( m_base1, bare_photon_create_V, m_base3 ).sparseView();
    photon_annihilate_H = tensor( bare_photon_annihilate_H, m_base2, m_base3 ).sparseView();
    photon_annihilate_V = tensor( m_base1, bare_photon_annihilate_V, m_base3 ).sparseView();
    photon_n_H = tensor( bare_photon_n_H, m_base2, m_base3 ).sparseView();
    photon_n_V = tensor( m_base1, bare_photon_n_V, m_base3 ).sparseView();

    // Compressing
    atom_state_ground.makeCompressed();
    atom_state_biexciton.makeCompressed();
    atom_state_H.makeCompressed();
    atom_state_V.makeCompressed();
    atom_sigmaplus_G_H.makeCompressed();
    atom_sigmaminus_G_H.makeCompressed();
    atom_sigmaplus_H_B.makeCompressed();
    atom_sigmaminus_H_B.makeCompressed();
    atom_sigmaplus_G_V.makeCompressed();
    atom_sigmaminus_G_V.makeCompressed();
    atom_sigmaplus_V_B.makeCompressed();
    atom_sigmaminus_V_B.makeCompressed();
    atom_inversion_G_H.makeCompressed();
    atom_inversion_G_V.makeCompressed();
    atom_inversion_H_B.makeCompressed();
    atom_inversion_V_B.makeCompressed();
    atom_inversion_G_B.makeCompressed();
    atom_sigmaminus_G_B.makeCompressed();

    photon_create_H.makeCompressed();
    photon_create_V.makeCompressed();
    photon_annihilate_H.makeCompressed();
    photon_annihilate_V.makeCompressed();
    photon_n_H.makeCompressed();
    photon_n_V.makeCompressed();

    // Projector Matrices
    projector_atom_sigmaplus_G_H = project_matrix_sparse( atom_sigmaplus_G_H );
    projector_atom_sigmaminus_G_H = project_matrix_sparse( atom_sigmaminus_G_H );
    projector_atom_sigmaplus_H_B = project_matrix_sparse( atom_sigmaplus_H_B );
    projector_atom_sigmaminus_H_B = project_matrix_sparse( atom_sigmaminus_H_B );
    projector_atom_sigmaplus_G_V = project_matrix_sparse( atom_sigmaplus_G_V );
    projector_atom_sigmaminus_G_V = project_matrix_sparse( atom_sigmaminus_G_V );
    projector_atom_sigmaplus_V_B = project_matrix_sparse( atom_sigmaplus_V_B );
    projector_atom_sigmaminus_V_B = project_matrix_sparse( atom_sigmaminus_V_B );
    projector_photon_create_H = project_matrix_sparse( photon_create_H );
    projector_photon_annihilate_H = project_matrix_sparse( photon_annihilate_H );
    projector_photon_create_V = project_matrix_sparse( photon_create_V );
    projector_photon_annihilate_V = project_matrix_sparse( photon_annihilate_V );

    // All possible Hamiltonions
    Log::L2( "Done! Creating Hamiltonoperator... " );
    // H_0
    H_0 = p.p_omega_atomic_G_H * atom_state_H + p.p_omega_atomic_G_V * atom_state_V + p.p_omega_atomic_B * atom_state_biexciton + p.p_omega_cavity_H * photon_n_H + p.p_omega_cavity_V * photon_n_V;
    //H_0 = p.p_omega_atomic_G_H / 2.0 * atom_inversion_G_H + p.p_omega_atomic_G_V / 2.0 * atom_inversion_G_V + p.p_omega_atomic_H_B / 2.0 * atom_inversion_H_B + p.p_omega_atomic_V_B / 2.0 * atom_inversion_V_B + p.p_omega_cavity_H * photon_n_H + p.p_omega_cavity_V * photon_n_V;
    // H_I
    // RWA
    if ( p.numerics_use_rwa ) {
        Log::L2( "using RWA... " );
        H_I = p.p_omega_coupling * ( atom_sigmaplus_G_H * photon_annihilate_H + atom_sigmaminus_G_H * photon_create_H + atom_sigmaplus_G_V * photon_annihilate_V + atom_sigmaminus_G_V * photon_create_V + atom_sigmaplus_H_B * photon_annihilate_H + atom_sigmaminus_H_B * photon_create_H + atom_sigmaplus_V_B * photon_annihilate_V + atom_sigmaminus_V_B * photon_create_V );
    }
    // non RWA
    if ( !p.numerics_use_rwa ) {
        Log::L2( "NOT using RWA... " );
        H_I = p.p_omega_coupling * ( atom_sigmaplus_G_H * photon_create_H + atom_sigmaplus_G_H * photon_annihilate_H + atom_sigmaminus_G_H * photon_create_H + atom_sigmaminus_G_H * photon_annihilate_H + atom_sigmaplus_G_V * photon_create_V + atom_sigmaplus_G_V * photon_annihilate_V + atom_sigmaminus_G_V * photon_create_V + atom_sigmaminus_G_V * photon_annihilate_V + atom_sigmaplus_H_B * photon_create_H + atom_sigmaplus_H_B * photon_annihilate_H + atom_sigmaminus_H_B * photon_create_H + atom_sigmaminus_H_B * photon_annihilate_H + atom_sigmaplus_V_B * photon_create_V + atom_sigmaplus_V_B * photon_annihilate_V + atom_sigmaminus_V_B * photon_create_V + atom_sigmaminus_V_B * photon_annihilate_V );
    }
    // H
    H = H_0 + H_I;
    // Interaction picture
    if ( p.numerics_use_interactionpicture ) {
        Log::L2( "using interaction picture... " );
        H_used = H_I;
    }
    if ( !p.numerics_use_interactionpicture ) {
        Log::L2( "NOT using interaction picture... " );
        H_used = H;
    }
    Log::L2( "Hamiltonoperator done!\n" );
    // rho, experimental: start with coherent state. in this case, always start in ground state.
    if ( !p.startCoherent ) {
        Log::L2( "Setting initial rho as pure state with rho_0 = {}... \n", p.p_initial_state );
        rho.coeffRef( p.p_initial_state, p.p_initial_state ) = 1;
    } else {
        Log::L2( "Setting initial rho as pure coherent state with alpha_h = {}, alpha_v = ... \n", p.p_initial_state_photon_h, p.p_initial_state_photon_v );
        double trace_rest_h = 1.0;
        double trace_rest_v = 1.0;
        for ( int i = 0; i < p.p_max_photon_number; i++ ) {
            auto coherent_h = getCoherent( std::sqrt( p.p_initial_state_photon_h ), i );
            auto coherent_v = 0.0; //getCoherent( std::sqrt(p.p_initial_state_photon_v), i );
            int state_h = 4 * ( p.p_max_photon_number + 1 ) * i;
            int state_v = 4 * i;
            rho.coeffRef( state_h, state_h ) = coherent_h; // Remember, <n> = alpha^2 -> alpha = sqrt(n) !!
            //rho.coeffRef( state_v, state_v ) = coherent_v; // Remember, <n> = alpha^2 -> alpha = sqrt(n) !!
            trace_rest_h -= coherent_h;
            trace_rest_v -= coherent_v;
            Log::L3( "Coherent state at N = {} with coefficient H = {}, V = {}\n", i, coherent_h, coherent_v );
        }
        int final_h = 4 * ( p.p_max_photon_number + 1 ) * p.p_max_photon_number;
        int final_v = 4 * p.p_max_photon_number;
        rho.coeffRef( final_h, final_h ) = trace_rest_h;
        //rho.coeffRef( final_v, final_v ) = trace_rest_v;
        Log::L3( "Coherent state at N = {} with coefficient H = {}, V = {}\n", p.p_max_photon_number, trace_rest_h, trace_rest_v );
    }
    Log::L2( "Done!\n" );
    return true;
}

void OperatorMatrices::outputOperators( const Parameters &p ) {
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
        OperatorMatricesText test = OperatorMatricesText();
        test.generateOperators( p );
        if ( p.output_operators == 3 )
            exit( 0 );
    }
}