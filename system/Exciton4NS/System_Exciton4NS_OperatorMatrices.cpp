#include "../operatormatrices.cpp"

class OperatorMatrices : public OperatorMatrices_Parent {
   public:
    MatrixXcd H;
    MatrixXcd H_0;
    MatrixXcd H_I;
    MatrixXcd rho;
    MatrixXcd H_used;

    //Operator Matrices
    MatrixXcd photon_create_H, photon_annihilate_H, photon_create_V, photon_annihilate_V;
    MatrixXcd atom_state_biexciton, atom_state_ground, atom_state_H, atom_state_V, photon_n_H, photon_n_V;
    MatrixXcd atom_sigmaplus_G_H, atom_sigmaminus_G_H, atom_sigmaplus_H_B, atom_sigmaminus_H_B, atom_sigmaplus_G_V, atom_sigmaminus_G_V, atom_sigmaplus_V_B, atom_sigmaminus_V_B;
    MatrixXcd atom_inversion_G_H, atom_inversion_G_V, atom_inversion_H_B, atom_inversion_V_B, atom_inversion_G_B;

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
        H = MatrixXcd::Zero( p.maxStates, p.maxStates );
        H_0 = MatrixXcd::Zero( p.maxStates, p.maxStates );
        H_I = MatrixXcd::Zero( p.maxStates, p.maxStates );
        H_used = MatrixXcd::Zero( p.maxStates, p.maxStates );
        rho = MatrixXcd::Zero( p.maxStates, p.maxStates );

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

        // Photonic mapping operators
        MatrixXcd map1 = MatrixXcd::Zero( 2, 2 );
        MatrixXcd map2 = MatrixXcd::Zero( 2, 2 );
        map1 << 1, 0, 0, 0;
        map2 << 0, 0, 0, 1;
        // Photon operators
        bare_photon_create_H = map_operator( create_photonic_operator( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number ), map1 );
        bare_photon_annihilate_H = map_operator( create_photonic_operator( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number ), map1 );
        bare_photon_n_H = bare_photon_create * bare_photon_annihilate;
        bare_photon_create_V = map_operator( create_photonic_operator( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number ), map2 );
        bare_photon_annihilate_V = map_operator( create_photonic_operator( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number ), map2 );
        bare_photon_n_V = bare_photon_create * bare_photon_annihilate;

        // Expanding both states
        logs.level2( "Expanding single state matrices... " );
        atom_state_ground = expand_atomic_operator( bare_atom_state_ground, p.p_max_photon_number );
        atom_state_biexciton = expand_atomic_operator( bare_atom_state_biexciton, p.p_max_photon_number );
        atom_state_H = expand_atomic_operator( bare_atom_state_H, p.p_max_photon_number );
        atom_state_V = expand_atomic_operator( bare_atom_state_V, p.p_max_photon_number );
        atom_sigmaplus_G_H = expand_atomic_operator( bare_atom_sigmaplus_G_H, p.p_max_photon_number );
        atom_sigmaminus_G_H = expand_atomic_operator( bare_atom_sigmaminus_G_H, p.p_max_photon_number );
        atom_sigmaplus_H_B = expand_atomic_operator( bare_atom_sigmaplus_H_B, p.p_max_photon_number );
        atom_sigmaminus_H_B = expand_atomic_operator( bare_atom_sigmaminus_H_B, p.p_max_photon_number );
        atom_sigmaplus_G_V = expand_atomic_operator( bare_atom_sigmaplus_G_V, p.p_max_photon_number );
        atom_sigmaminus_G_V = expand_atomic_operator( bare_atom_sigmaminus_G_V, p.p_max_photon_number );
        atom_sigmaplus_V_B = expand_atomic_operator( bare_atom_sigmaplus_V_B, p.p_max_photon_number );
        atom_sigmaminus_V_B = expand_atomic_operator( bare_atom_sigmaminus_V_B, p.p_max_photon_number );
        atom_inversion_G_H = expand_atomic_operator( bare_atom_inversion_G_H, p.p_max_photon_number );
        atom_inversion_G_V = expand_atomic_operator( bare_atom_inversion_G_V, p.p_max_photon_number );
        atom_inversion_H_B = expand_atomic_operator( bare_atom_inversion_H_B, p.p_max_photon_number );
        atom_inversion_V_B = expand_atomic_operator( bare_atom_inversion_V_B, p.p_max_photon_number );
        atom_inversion_G_B = expand_atomic_operator( bare_atom_inversion_G_B, p.p_max_photon_number );

        photon_create_H = expand_photonic_operator( bare_photon_create_H, 4 );
        photon_create_V = expand_photonic_operator( bare_photon_create_V, 4 );
        photon_annihilate_H = expand_photonic_operator( bare_photon_annihilate_H, 4 );
        photon_annihilate_V = expand_photonic_operator( bare_photon_annihilate_V, 4 );
        photon_n_H = expand_photonic_operator( bare_photon_n_H, 4 );
        photon_n_V = expand_photonic_operator( bare_photon_n_V, 4 );

        // All possible Hamiltonions
        logs.level2( "Done! Creating Hamiltonoperator... " );
        // H_0
        H_0 = p.p_omega_atomic_G_H / 2.0 * atom_inversion_G_H + p.p_omega_atomic_G_V / 2.0 * atom_inversion_G_V + p.p_omega_atomic_H_B / 2.0 * atom_inversion_H_B + p.p_omega_atomic_V_B / 2.0 * atom_inversion_V_B + p.p_omega_cavity_H * photon_n_H + p.p_omega_cavity_V * photon_n_V;
        // H_I
        // RWA
        if ( p.numerics_use_rwa ) {
            logs.level2( "using RWA... " );
            H_I = p.p_omega_coupling * ( atom_sigmaplus_G_H * photon_annihilate_G_H + atom_sigmaminus_G_H * photon_create_G_H + atom_sigmaplus_G_V * photon_annihilate_G_V + atom_sigmaminus_G_V * photon_create_G_V + atom_sigmaplus_H_B * photon_annihilate_H_B + atom_sigmaminus_H_B * photon_create_H_B + atom_sigmaplus_V_B * photon_annihilate_V_B + atom_sigmaminus_V_B * photon_create_V_B );
        }
        // non RWA
        if ( !p.numerics_use_rwa ) {
            logs.level2( "NOT using RWA... " );
            H_I = p.p_omega_coupling * ( atom_sigmaplus_G_H * photon_create_G_H + atom_sigmaplus_G_H * photon_annihilate_G_H + atom_sigmaminus_G_H * photon_create_G_H + atom_sigmaminus_G_H * photon_annihilate_G_H + atom_sigmaplus_G_V * photon_create_G_V + atom_sigmaplus_G_V * photon_annihilate_G_V + atom_sigmaminus_G_V * photon_create_G_V + atom_sigmaminus_G_V * photon_annihilate_G_V + atom_sigmaplus_H_B * photon_create_H_B + atom_sigmaplus_H_B * photon_annihilate_H_B + atom_sigmaminus_H_B * photon_create_H_B + atom_sigmaminus_H_B * photon_annihilate_H_B + atom_sigmaplus_V_B * photon_create_V_B + atom_sigmaplus_V_B * photon_annihilate_V_B + atom_sigmaminus_V_B * photon_create_V_B + atom_sigmaminus_V_B * photon_annihilate_V_B );
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
            rho( p.p_initial_state, p.p_initial_state ) = 1;
        else {
            double trace_rest = 1.0;
            for ( int i = 0; i < p.p_max_photon_number; i++ ) {
                rho( i * 4, i * 4 ) = getCoherent( p.p_initial_state, i ); // Remember, <n> = alpha^2 -> alpha = sqrt(n) !!
                trace_rest -= getCoherent( p.p_initial_state, i );
                logs( "Coherent state at N = {} with coefficient {}\n", i, getCoherent( p.p_initial_state, i ) );
            }
            rho( p.p_max_photon_number * 2, p.p_max_photon_number * 2 ) = trace_rest;
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
                    << atom_state_ground.format( CleanFmt ) << std::endl;
                out << "atom_state_biexciton\n"
                    << atom_state_biexciton.format( CleanFmt ) << std::endl;
                out << "atom_state_H\n"
                    << atom_state_H.format( CleanFmt ) << std::endl;
                out << "atom_state_V\n"
                    << atom_state_V.format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus_G_H\n"
                    << atom_sigmaplus_G_H.format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus_G_H\n"
                    << atom_sigmaminus_G_H.format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus_H_B\n"
                    << atom_sigmaplus_H_B.format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus_H_B\n"
                    << atom_sigmaminus_H_B.format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus_G_V\n"
                    << atom_sigmaplus_G_V.format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus_G_V\n"
                    << atom_sigmaminus_G_V.format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus_V_B\n"
                    << atom_sigmaplus_V_B.format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus_V_B\n"
                    << atom_sigmaminus_V_B.format( CleanFmt ) << std::endl;
                out << "atom_inversion_G_H\n"
                    << atom_inversion_G_H.format( CleanFmt ) << std::endl;
                out << "atom_inversion_G_V\n"
                    << atom_inversion_G_V.format( CleanFmt ) << std::endl;
                out << "atom_inversion_H_B\n"
                    << atom_inversion_H_B.format( CleanFmt ) << std::endl;
                out << "atom_inversion_V_B\n"
                    << atom_inversion_V_B.format( CleanFmt ) << std::endl;
                out << "atom_inversion_G_B\n"
                    << atom_inversion_G_B.format( CleanFmt ) << std::endl;
                out << "photon_create_H\n"
                    << photon_create_H.format( CleanFmt ) << std::endl;
                out << "photon_create_V\n"
                    << photon_create_V.format( CleanFmt ) << std::endl;
                out << "photon_annihilate_H\n"
                    << photon_annihilate_H.format( CleanFmt ) << std::endl;
                out << "photon_annihilate_V\n"
                    << photon_annihilate_V.format( CleanFmt ) << std::endl;
                out << "photon_n_H\n"
                    << photon_n_H.format( CleanFmt ) << std::endl;
                out << "photon_n_V\n"
                    << photon_n_V.format( CleanFmt ) << std::endl;
                
            }
            out << "Hamilton and Rho:\nH=H_0+H_I (no RWA)\n"
                << H.format( CleanFmt ) << std::endl;
            out << "H_0\n"
                << H_0.format( CleanFmt ) << std::endl;
            out << "H_I\n"
                << H_I.format( CleanFmt ) << std::endl;
            out << "H_used\n"
                << H_used.format( CleanFmt ) << std::endl;
            out << "rho\n"
                << rho.format( CleanFmt ) << std::endl;
            //out << "test1\n" << test1.format(CleanFmt)<< "\ntest2\n" << test2.format(CleanFmt) << std::endl;
            logs.level2( out.str() );
            if ( p.output_operators == 3 )
                exit( 0 );
        }
    }
};