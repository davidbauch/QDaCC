#include "../operatormatrices.cpp"

class OperatorMatrices : public OperatorMatrices_Parent {
   public:
    MatrixXcd H;
    MatrixXcd H_0;
    MatrixXcd H_I;
    MatrixXcd rho;
    MatrixXcd H_used;
    std::string base;

    //Operator Matrices
    MatrixXcd photon_create, photon_annihilate, atom_exited, atom_ground, photon_n, atom_sigmaplus, atom_sigmaminus, atom_inversion;

    // Bare matrices:
    MatrixXcd bare_photon_create, bare_photon_annihilate, bare_photon_n, bare_atom_exited, bare_atom_ground, bare_atom_sigmaplus, bare_atom_sigmaminus, bare_atom_inversion;

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

        std::vector<std::string> base2 = {"G", "E"};
        std::vector<std::string> base1; for (int i = 0; i <= p.p_max_photon_number; i++) {base1.emplace_back(std::to_string(i));}
        base = tensor( base1, base2 );
        
        logs.level2( "Operator Base: (size {}) ", base.size() );
        for ( auto b : base )
            logs.level2( "|{}> ", b );
        logs.level2( "\n" );

        // Initializing bare matrices:
        logs.level2( "Initializing base matrices... " );
        bare_atom_exited = MatrixXcd::Zero( 2, 2 );
        bare_atom_exited << 0, 0,
            0, 1;
        bare_atom_ground = MatrixXcd::Zero( 2, 2 );
        bare_atom_ground << 1, 0,
            0, 0;
        bare_atom_sigmaplus = MatrixXcd::Zero( 2, 2 );
        bare_atom_sigmaplus << 0, 0,
            1, 0;
        bare_atom_sigmaminus = MatrixXcd::Zero( 2, 2 );
        bare_atom_sigmaminus << 0, 1,
            0, 0;
        bare_atom_inversion = MatrixXcd::Zero( 2, 2 );
        bare_atom_inversion << -1, 0,
            0, 1;
        bare_photon_create = create_photonic_operator<Eigen::MatrixXcd>( OPERATOR_PHOTONIC_CREATE, p.p_max_photon_number );
        bare_photon_annihilate = create_photonic_operator<Eigen::MatrixXcd>( OPERATOR_PHOTONIC_ANNIHILATE, p.p_max_photon_number );
        bare_photon_n = bare_photon_create * bare_photon_annihilate;

        Eigen::MatrixXcd m_base1 = MatrixXcd::Identity(bare_photon_create.rows(), bare_photon_create.cols()); 
        Eigen::MatrixXcd m_base2 = MatrixXcd::Identity(2,2);

        // Expanding both states
        logs.level2( "Expanding single state matrices... " );
        atom_exited = tensor( m_base1, bare_atom_exited );
        atom_ground = tensor( m_base1,bare_atom_ground );
        atom_sigmaminus = tensor( m_base1,bare_atom_sigmaminus );
        atom_sigmaplus = tensor( m_base1,bare_atom_sigmaplus );
        atom_inversion = tensor( m_base1,bare_atom_inversion );

        photon_create = tensor( bare_photon_create, m_base2 );
        photon_annihilate = tensor( bare_photon_annihilate, m_base2 );
        photon_n = tensor( bare_photon_n, m_base2 );

        // All possible Hamiltonions
        logs.level2( "Done! Creating Hamiltonoperator... " );
        // H_0
        H_0 = p.p_omega_atomic / 2.0 * atom_inversion + p.p_omega_cavity * photon_n;
        // H_I
        // RWA
        if ( p.numerics_use_rwa ) {
            logs.level2( "using RWA... " );
            H_I = p.p_omega_coupling * ( atom_sigmaplus * photon_annihilate + atom_sigmaminus * photon_create );
        }
        // non RWA
        if ( !p.numerics_use_rwa ) {
            logs.level2( "NOT using RWA... " );
            H_I = p.p_omega_coupling * ( atom_sigmaplus * photon_create + atom_sigmaplus * photon_annihilate + atom_sigmaminus * photon_create + atom_sigmaminus * photon_annihilate );
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
        logs.level2( "Hamiltonoperator done! Used:\n{}\nSetting initial rho as pure state with rho_0 = {}...\n", H_used, p.p_initial_state );
        // rho, experimental: start with coherent state. in this case, always start in ground state.
        if ( !p.startCoherent )
            rho( p.p_initial_state, p.p_initial_state ) = 1;
        else {
            double trace_rest = 1.0;
            for ( int i = 0; i < p.p_max_photon_number; i++ ) {
                rho( i * 2, i * 2 ) = getCoherent( p.p_initial_state, i ); // Remember, <n> = alpha^2 -> alpha = sqrt(n) !!
                trace_rest -= getCoherent( p.p_initial_state, i );
                logs( "Coherent state at N = {} with coefficient {}\n", i, getCoherent( p.p_initial_state, i ) );
            }
            rho( p.p_max_photon_number * 2, p.p_max_photon_number * 2 ) = trace_rest;
            logs( "Coherent state at N = {} with coefficient {}\n", p.p_max_photon_number, trace_rest );
        }
        return checkMatrices2NS( p );
    }

    void outputOperators( const Parameters &p ) {
        if ( p.output_operators > 0 ) {
            std::ostringstream out;
            Eigen::IOFormat CleanFmt( 4, 0, ", ", "\n", "[", "]" );
            if ( p.output_operators > 1 ) {
                out << "General Operators:\natom_exited\n"
                    << atom_exited.format( CleanFmt ) << std::endl;
                out << "atom_ground\n"
                    << atom_ground.format( CleanFmt ) << std::endl;
                out << "atom_sigmaplus\n"
                    << atom_sigmaplus.format( CleanFmt ) << std::endl;
                out << "atom_sigmaminus\n"
                    << atom_sigmaminus.format( CleanFmt ) << std::endl;
                out << "atom_inversion\n"
                    << atom_inversion.format( CleanFmt ) << std::endl;
                out << "photon_create\n"
                    << photon_create.format( CleanFmt ) << std::endl;
                out << "photon_annihilate\n"
                    << photon_annihilate.format( CleanFmt ) << std::endl;
                out << "photon_n\n"
                    << photon_n.format( CleanFmt ) << std::endl;
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
    // Redundant
    bool checkMatrices2NS( const Parameters &p ) {
        bool valid = true;
        // Creating Matrices the old way
        MatrixXcd photon_create_check, photon_annihilate_check, atom_exited_check, atom_ground_check, photon_n_check, atom_sigmaplus_check, atom_sigmaminus_check, atom_inversion_check;
        photon_create_check = MatrixXcd::Zero( p.maxStates, p.maxStates );
        photon_annihilate_check = MatrixXcd::Zero( p.maxStates, p.maxStates );
        photon_n_check = MatrixXcd::Zero( p.maxStates, p.maxStates );
        atom_exited_check = MatrixXcd::Zero( p.maxStates, p.maxStates );
        atom_ground_check = MatrixXcd::Zero( p.maxStates, p.maxStates );
        atom_sigmaplus_check = MatrixXcd::Zero( p.maxStates, p.maxStates );
        atom_sigmaminus_check = MatrixXcd::Zero( p.maxStates, p.maxStates );
        atom_inversion_check = MatrixXcd::Zero( p.maxStates, p.maxStates );
        // Generating Matrices
        for ( int n = 0; n < p.maxStates; n++ )
            for ( int m = 0; m < p.maxStates; m++ ) {
                if ( n == m ) photon_n_check( n, m ) = (int)( ( n ) / 2 );
                photon_create_check( n, m ) = std::sqrt( (int)( m / 2 + 1 ) ) * delta( n % 2, m % 2 ) * delta( (int)( n / 2 ), (int)( m / 2 ) + 1 );
                photon_annihilate_check( n, m ) = std::sqrt( (int)( m / 2 ) ) * delta( n % 2, m % 2 ) * delta( (int)( n / 2 ), (int)( m / 2 ) - 1 );
                atom_sigmaplus_check( n, m ) = delta( (int)( n / 2 ), (int)( m / 2 ) ) * delta( (int)( n % 2 ), 1 ) * delta( (int)( m % 2 ), 0 );
                atom_sigmaminus_check( n, m ) = delta( (int)( n / 2 ), (int)( m / 2 ) ) * delta( (int)( n % 2 ), 0 ) * delta( (int)( m % 2 ), 1 );
                atom_ground_check( n, m ) = delta( (int)( n / 2 ), (int)( m / 2 ) ) * delta( (int)( n % 2 ), 0 ) * delta( (int)( m % 2 ), 0 );
                atom_exited_check( n, m ) = delta( (int)( n / 2 ), (int)( m / 2 ) ) * delta( (int)( n % 2 ), 1 ) * delta( (int)( m % 2 ), 1 );
            }
        for ( int n = 0; n < p.maxStates; n++ ) {
            atom_inversion_check( n, n ) = ( n % 2 == 0 ) ? -1 : 1;
        }
        // Checking
        if ( !atom_exited.isApprox( atom_exited_check ) ) {
            logs.level2( "Operator atom_exited is wrong!\n" );
            valid = false;
        }
        if ( !atom_ground.isApprox( atom_ground_check ) ) {
            logs.level2( "Operator atom_ground is wrong!\n" );
            valid = false;
        }
        if ( !atom_sigmaminus.isApprox( atom_sigmaminus_check ) ) {
            logs.level2( "Operator atom_sigmaminus is wrong!\n" );
            valid = false;
        }
        if ( !atom_sigmaplus.isApprox( atom_sigmaplus_check ) ) {
            logs.level2( "Operator atom_sigmaplus is wrong!\n" );
            valid = false;
        }
        if ( !atom_inversion.isApprox( atom_inversion_check ) ) {
            logs.level2( "Operator atom_inversion is wrong!\n" );
            valid = false;
        }
        if ( !photon_create.isApprox( photon_create_check ) ) {
            logs.level2( "Operator photon_create is wrong!\n" );
            valid = false;
        }
        if ( !photon_annihilate.isApprox( photon_annihilate_check ) ) {
            logs.level2( "Operator photon_annihilate is wrong!\n" );
            valid = false;
        }
        if ( !photon_n.isApprox( photon_n_check ) ) {
            logs.level2( "Operator photon_n is wrong!\n" );
            valid = false;
        }
        return valid;
    }
};