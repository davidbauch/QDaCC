#include "../operatormatrices.cpp"

class OperatorMatrices : public OperatorMatrices_Parent {
   public:
    SparseMat H;
    SparseMat H_0;
    SparseMat H_I;
    SparseMat rho;
    SparseMat H_used;
    std::vector<std::string> base;

    //Operator Matrices
    SparseMat photon_create_H, photon_annihilate_H, photon_create_V, photon_annihilate_V;
    SparseMat atom_state_biexciton, atom_state_ground, atom_state_H, atom_state_V, photon_n_H, photon_n_V;
    SparseMat atom_sigmaplus_G_H, atom_sigmaminus_G_H, atom_sigmaplus_H_B, atom_sigmaminus_H_B, atom_sigmaplus_G_V, atom_sigmaminus_G_V, atom_sigmaplus_V_B, atom_sigmaminus_V_B;
    SparseMat atom_inversion_G_H, atom_inversion_G_V, atom_inversion_H_B, atom_inversion_V_B, atom_inversion_G_B;

    // Bare matrices:
    MatrixXcd bare_photon_create_H, bare_photon_annihilate_H, bare_photon_create_V, bare_photon_annihilate_V;
    MatrixXcd bare_atom_state_biexciton, bare_atom_state_ground, bare_atom_state_H, bare_atom_state_V, bare_photon_n_H, bare_photon_n_V;
    MatrixXcd bare_atom_sigmaplus_G_H, bare_atom_sigmaminus_G_H, bare_atom_sigmaplus_H_B, bare_atom_sigmaminus_H_B, bare_atom_sigmaplus_G_V, bare_atom_sigmaminus_G_V, bare_atom_sigmaplus_V_B, bare_atom_sigmaminus_V_B;
    MatrixXcd bare_atom_inversion_G_H, bare_atom_inversion_G_V, bare_atom_inversion_H_B, bare_atom_inversion_V_B, bare_atom_inversion_G_B;

    OperatorMatrices(){};
    OperatorMatrices( const Parameters &p ) {
        init( p );
    }

    bool generateOperators( const Parameters &p ) {
        // Zeroing Hamiltons (redundant at this point)
        logs.level2( "Creating operator matrices, dimension = {}\nCreating base matrices... ", p.maxStates );
        H = SparseMat( p.maxStates, p.maxStates );
        H_0 = SparseMat( p.maxStates, p.maxStates );
        H_I = SparseMat( p.maxStates, p.maxStates );
        H_used = SparseMat( p.maxStates, p.maxStates );
        rho = SparseMat( p.maxStates, p.maxStates );

        // Create Base Matrices and Base Vector:
        Eigen::MatrixXcd m_base1 = MatrixXcd::Identity( p.p_max_photon_number + 1, p.p_max_photon_number + 1 );
        Eigen::MatrixXcd m_base2 = MatrixXcd::Identity( p.p_max_photon_number + 1, p.p_max_photon_number + 1 );
        Eigen::MatrixXcd m_base3 = MatrixXcd::Identity( 4, 4 );
        std::vector<std::string> base1;
        std::vector<std::string> base2;
        std::vector<std::string> base3 = {"G", "X_H", "X_V", "B"};
        for ( int i = 0; i <= p.p_max_photon_number; i++ ) {
            base1.emplace_back( std::to_string( i ) + "_H" );
            base2.emplace_back( std::to_string( i ) + "_V" );
        }
        base = tensor( tensor( base1, base2 ), base3 );

        logs.level2( "Operator Base: (size {}) ", base.size() );
        for ( auto b : base )
            logs.level2( "|{}> ", b );
        logs.level2( "\n" );

        // Initializing bare matrices:
        logs.level2( "Initializing base matrices... " );
        // Atomic state operators
        bare_atom_state_ground = MatrixXcd::Zero( 4, 4 );
        bare_atom_state_ground << 1, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;
        bare_atom_state_biexciton = MatrixXcd::Zero( 4, 4 );
        bare_atom_state_biexciton << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 1;
        bare_atom_state_H = MatrixXcd::Zero( 4, 4 );
        bare_atom_state_H << 0, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;
        bare_atom_state_V = MatrixXcd::Zero( 4, 4 );
        bare_atom_state_V << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 0;
        bare_atom_sigmaplus_G_H = MatrixXcd::Zero( 4, 4 );
        bare_atom_sigmaplus_G_H << 0, 0, 0, 0,
            1, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;
        bare_atom_sigmaminus_G_H = MatrixXcd::Zero( 4, 4 );
        bare_atom_sigmaminus_G_H << 0, 1, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;
        bare_atom_sigmaplus_G_V = MatrixXcd::Zero( 4, 4 );
        bare_atom_sigmaplus_G_V << 0, 0, 0, 0,
            0, 0, 0, 0,
            1, 0, 0, 0,
            0, 0, 0, 0;
        bare_atom_sigmaminus_G_V = MatrixXcd::Zero( 4, 4 );
        bare_atom_sigmaminus_G_V << 0, 0, 1, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;
        bare_atom_sigmaplus_H_B = MatrixXcd::Zero( 4, 4 );
        bare_atom_sigmaplus_H_B << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 1, 0, 0;
        bare_atom_sigmaminus_H_B = MatrixXcd::Zero( 4, 4 );
        bare_atom_sigmaminus_H_B << 0, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 0, 0,
            0, 0, 0, 0;
        bare_atom_sigmaplus_V_B = MatrixXcd::Zero( 4, 4 );
        bare_atom_sigmaplus_V_B << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 1, 0;
        bare_atom_sigmaminus_V_B = MatrixXcd::Zero( 4, 4 );
        bare_atom_sigmaminus_V_B << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 0, 0;
        bare_atom_inversion_G_H = MatrixXcd::Zero( 4, 4 );
        bare_atom_inversion_G_H << -1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;
        bare_atom_inversion_G_V = MatrixXcd::Zero( 4, 4 );
        bare_atom_inversion_G_V << -1, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 0;
        bare_atom_inversion_H_B = MatrixXcd::Zero( 4, 4 );
        bare_atom_inversion_H_B << 0, 0, 0, 0,
            0, -1, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 1;
        bare_atom_inversion_V_B = MatrixXcd::Zero( 4, 4 );
        bare_atom_inversion_V_B << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, -1, 0,
            0, 0, 0, 1;
        bare_atom_inversion_G_B = MatrixXcd::Zero( 4, 4 );
        bare_atom_inversion_G_B << -1, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 1;

        // Photon operators
        bare_photon_create_H = create_photonic_operator<Eigen::MatrixXcd>( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number );
        bare_photon_annihilate_H = create_photonic_operator<Eigen::MatrixXcd>( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number );
        bare_photon_n_H = bare_photon_create_H * bare_photon_annihilate_H;
        bare_photon_create_V = create_photonic_operator<Eigen::MatrixXcd>( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number );
        bare_photon_annihilate_V = create_photonic_operator<Eigen::MatrixXcd>( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number );
        bare_photon_n_V = bare_photon_create_V * bare_photon_annihilate_V;

        // Expanding both states
        logs.level2( "Expanding single state matrices... " );
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

        photon_create_H.makeCompressed();
        photon_create_V.makeCompressed();
        photon_annihilate_H.makeCompressed();
        photon_annihilate_V.makeCompressed();
        photon_n_H.makeCompressed();
        photon_n_V.makeCompressed();


        // All possible Hamiltonions
        logs.level2( "Done! Creating Hamiltonoperator... " );
        // H_0
        H_0 = p.p_omega_atomic_G_H / 2.0 * atom_inversion_G_H + p.p_omega_atomic_G_V / 2.0 * atom_inversion_G_V + p.p_omega_atomic_H_B / 2.0 * atom_inversion_H_B + p.p_omega_atomic_V_B / 2.0 * atom_inversion_V_B + p.p_omega_cavity_H * photon_n_H + p.p_omega_cavity_V * photon_n_V;
        // H_I
        // RWA
        if ( p.numerics_use_rwa ) {
            logs.level2( "using RWA... " );
            H_I = p.p_omega_coupling * ( atom_sigmaplus_G_H * photon_annihilate_H + atom_sigmaminus_G_H * photon_create_H + atom_sigmaplus_G_V * photon_annihilate_V + atom_sigmaminus_G_V * photon_create_V + atom_sigmaplus_H_B * photon_annihilate_H + atom_sigmaminus_H_B * photon_create_H + atom_sigmaplus_V_B * photon_annihilate_V + atom_sigmaminus_V_B * photon_create_V );
        }
        // non RWA
        if ( !p.numerics_use_rwa ) {
            logs.level2( "NOT using RWA... " );
            H_I = p.p_omega_coupling * ( atom_sigmaplus_G_H * photon_create_H + atom_sigmaplus_G_H * photon_annihilate_H + atom_sigmaminus_G_H * photon_create_H + atom_sigmaminus_G_H * photon_annihilate_H + atom_sigmaplus_G_V * photon_create_V + atom_sigmaplus_G_V * photon_annihilate_V + atom_sigmaminus_G_V * photon_create_V + atom_sigmaminus_G_V * photon_annihilate_V + atom_sigmaplus_H_B * photon_create_H + atom_sigmaplus_H_B * photon_annihilate_H + atom_sigmaminus_H_B * photon_create_H + atom_sigmaminus_H_B * photon_annihilate_H + atom_sigmaplus_V_B * photon_create_V + atom_sigmaplus_V_B * photon_annihilate_V + atom_sigmaminus_V_B * photon_create_V + atom_sigmaminus_V_B * photon_annihilate_V );
        }
        // H
        H = H_0 + H_I;
        // Interaction picture
        if ( p.numerics_use_interactionpicture ) {
            logs.level2( "using interaction picture... " );
            H_used = H_I;
        }
        if ( !p.numerics_use_interactionpicture ) {
            logs.level2( "NOT using interaction picture... " );
            H_used = H;
        }
        logs.level2( "Hamiltonoperator done! Used:\n{}\nSetting initial rho as pure state with rho_0 = {}... ", H_used, p.p_initial_state );
        // rho, experimental: start with coherent state. in this case, always start in ground state.
        if ( !p.startCoherent )
            rho.coeffRef( p.p_initial_state, p.p_initial_state ) = 1;
        else {
            double trace_rest = 1.0;
            for ( int i = 0; i < p.p_max_photon_number; i++ ) {
                rho.coeffRef( i * 4, i * 4 ) = getCoherent( p.p_initial_state, i ); // Remember, <n> = alpha^2 -> alpha = sqrt(n) !!
                trace_rest -= getCoherent( p.p_initial_state, i );
                logs( "Coherent state at N = {} with coefficient {}\n", i, getCoherent( p.p_initial_state, i ) );
            }
            rho.coeffRef( p.p_max_photon_number * 2, p.p_max_photon_number * 2 ) = trace_rest;
            logs( "Coherent state at N = {} with coefficient {}\n", p.p_max_photon_number, trace_rest );
        }
        return true;
    }

    void outputOperators( const Parameters &p ) {
        if ( p.output_operators > 0 ) {
            std::ostringstream out;
            Eigen::IOFormat CleanFmt( 4, 0, ", ", "\n", "[", "]" );
            if ( p.output_operators > 1 ) {
                out << "General Operators:\natom_state_ground\n"
                    << DenseMat(atom_state_ground).format( CleanFmt ) << std::endl;
                out << "atom_state_biexciton\n"
                    << DenseMat(atom_state_biexciton).format( CleanFmt ) << std::endl;
                out << "atom_state_H\n"
                    << DenseMat(atom_state_H).format( CleanFmt ) << std::endl;
                out << "atom_state_V\n"
                    << DenseMat(atom_state_V).format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus_G_H\n"
                    << DenseMat(atom_sigmaplus_G_H).format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus_G_H\n"
                    << DenseMat(atom_sigmaminus_G_H).format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus_H_B\n"
                    << DenseMat(atom_sigmaplus_H_B).format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus_H_B\n"
                    << DenseMat(atom_sigmaminus_H_B).format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus_G_V\n"
                    << DenseMat(atom_sigmaplus_G_V).format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus_G_V\n"
                    << DenseMat(atom_sigmaminus_G_V).format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus_V_B\n"
                    << DenseMat(atom_sigmaplus_V_B).format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus_V_B\n"
                    << DenseMat(atom_sigmaminus_V_B).format( CleanFmt ) << std::endl;
                out << "atom_inversion_G_H\n"
                    << DenseMat(atom_inversion_G_H).format( CleanFmt ) << std::endl;
                out << "atom_inversion_G_V\n"
                    << DenseMat(atom_inversion_G_V).format( CleanFmt ) << std::endl;
                out << "atom_inversion_H_B\n"
                    << DenseMat(atom_inversion_H_B).format( CleanFmt ) << std::endl;
                out << "atom_inversion_V_B\n"
                    << DenseMat(atom_inversion_V_B).format( CleanFmt ) << std::endl;
                out << "atom_inversion_G_B\n"
                    << DenseMat(atom_inversion_G_B).format( CleanFmt ) << std::endl;
                out << "photon_create_H\n"
                    << DenseMat(photon_create_H).format( CleanFmt ) << std::endl;
                out << "photon_create_V\n"
                    << DenseMat(photon_create_V).format( CleanFmt ) << std::endl;
                out << "photon_annihilate_H\n"
                    << DenseMat(photon_annihilate_H).format( CleanFmt ) << std::endl;
                out << "photon_annihilate_V\n"
                    << DenseMat(photon_annihilate_V).format( CleanFmt ) << std::endl;
                out << "photon_n_H\n"
                    << DenseMat(photon_n_H).format( CleanFmt ) << std::endl;
                out << "photon_n_V\n"
                    << DenseMat(photon_n_V).format( CleanFmt ) << std::endl;
            }
            out << "Hamilton and Rho:\nH=H_0+H_I (no RWA)\n"
                << DenseMat(H).format( CleanFmt ) << std::endl;
            out << "H_0\n"
                << DenseMat(H_0).format( CleanFmt ) << std::endl;
            out << "H_I\n"
                << DenseMat(H_I).format( CleanFmt ) << std::endl;
            out << "H_used\n"
                << DenseMat(H_used).format( CleanFmt ) << std::endl;
            out << "rho\n"
                << DenseMat(rho).format( CleanFmt ) << std::endl;
            //out << "test1\n" << test1.format(CleanFmt)<< "\ntest2\n" << test2.format(CleanFmt) << std::endl;
            logs.level2( out.str() );
            if ( p.output_operators == 3 )
                exit( 0 );
        }
    }
};