#include "system/operatormatrices.h"
#include "system/operatormatrices_text.h"
OperatorMatrix::OperatorMatrix( const OperatorMatrix &other ) {
    mat = other.mat;
}
OperatorMatrix::OperatorMatrix( int n, int m ) {
    mat.reserve( n );
    for ( int i = 0; i < n; i++ ) {
        mat.emplace_back( std::vector<std::string>() );
        mat.back().reserve( m );
        for ( int j = 0; j < m; j++ ) {
            mat.at( i ).emplace_back( "0" );
        }
    }
}
int OperatorMatrix::cols() const {
    return mat.at( 0 ).size();
}
int OperatorMatrix::rows() const {
    return mat.size();
}
// Row, Col
std::string OperatorMatrix::operator()( int i, int j ) const {
    return mat.at( i ).at( j );
}
// Row, Col
std::string OperatorMatrix::at( int i, int j ) const {
    return mat.at( i ).at( j );
}
// Row, Col
void OperatorMatrix::set( int i, int j, const std::string &val ) {
    mat.at( i ).at( j ) = val;
}
// Row, Col
std::string OperatorMatrix::getFormatted( int i, int j ) const {
    return mat.at( i ).at( j );
}
// Matrix-Matrix operations
OperatorMatrix OperatorMatrix::operator+( const OperatorMatrix &other ) {
    OperatorMatrix ret( rows(), cols() );
    for ( int i = 0; i < rows(); i++ ) {
        for ( int j = 0; j < cols(); j++ ) {
            if ( at( i, j ).front() == '0' && other.at( i, j ).front() != '0' )
                ret.set( i, j, other.at( i, j ) );
            else if ( at( i, j ).front() != '0' && other.at( i, j ).front() == '0' )
                ret.set( i, j, at( i, j ) );
            else if ( at( i, j ).front() != '0' && other.at( i, j ).front() != '0' )
                ret.set( i, j, at( i, j ) + "+" + other.at( i, j ) );
        }
    }
    return ret;
}
OperatorMatrix OperatorMatrix::operator-( const OperatorMatrix &other ) {
    OperatorMatrix ret = OperatorMatrix( rows(), cols() );
    for ( int i = 0; i < rows(); i++ ) {
        for ( int j = 0; j < cols(); j++ ) {
            if ( at( i, j ).front() == '0' && other.at( i, j ).front() != '0' )
                ret.set( i, j, "-" + other.at( i, j ) );
            else if ( at( i, j ).front() != '0' && other.at( i, j ).front() == '0' )
                ret.set( i, j, at( i, j ) );
            else if ( at( i, j ).front() != '0' && other.at( i, j ).front() != '0' )
                ret.set( i, j, at( i, j ) + "-" + other.at( i, j ) );
        }
    }
    return ret;
}
OperatorMatrix OperatorMatrix::operator*( const OperatorMatrix &other ) {
    OperatorMatrix ret = OperatorMatrix( rows(), cols() );
    for ( int i = 0; i < rows(); i++ ) {
        for ( int k = 0; k < other.cols(); k++ ) {
            std::string result = "";
            for ( int j = 0; j < cols(); j++ ) {
                if ( at( i, j ).front() != '0' && other.at( j, k ).front() != '0' ) {
                    if ( result.size() > 0 )
                        result += "+";
                    result += at( i, j ) + other.at( j, k );
                }
            }
            if ( result.size() == 0 )
                result = "0";
            ret.set( i, k, result );
        }
    }
    return ret;
}

OperatorMatrix OperatorMatrix::elementWiseMul( const OperatorMatrix &other ) {
    OperatorMatrix ret = OperatorMatrix( rows(), cols() );
    for ( int i = 0; i < rows(); i++ ) {
        for ( int j = 0; j < cols(); j++ ) {
            if ( at( i, j ).front() != '0' && other.at( i, j ).front() != '0' )
                ret.set( i, j, at( i, j ) + other.at( i, j ) );
            else
                ret.set( i, j, "0" );
        }
    }
    return ret;
}
OperatorMatrix OperatorMatrix::concat( const OperatorMatrix &other ) {
    OperatorMatrix ret = OperatorMatrix( rows(), cols() );
    for ( int i = 0; i < rows(); i++ ) {
        for ( int j = 0; j < cols(); j++ ) {
            ret.set( i, j, at( i, j ) + other.at( i, j ) );
        }
    }
    return ret;
}
OperatorMatrix OperatorMatrix::addElement( const std::string &element, int before_or_after ) {
    OperatorMatrix ret = OperatorMatrix( rows(), cols() );
    for ( int i = 0; i < rows(); i++ ) {
        for ( int j = 0; j < cols(); j++ ) {
            if ( before_or_after == -1 )
                ret.set( i, j, element + at( i, j ) );
            else
                ret.set( i, j, at( i, j ) + element );
        }
    }
    return ret;
}

OperatorMatrix operator+( const OperatorMatrix &mat, const std::string &var ) {
    OperatorMatrix ret( mat.rows(), mat.cols() );
    for ( int i = 0; i < mat.rows(); i++ ) {
        for ( int j = 0; j < mat.cols(); j++ ) {
            if ( mat.at( i, j ).front() == '0' )
                ret.set( i, j, var );
            else
                ret.set( i, j, mat.at( i, j ) + "+" + var );
        }
    }
    return ret;
}
OperatorMatrix operator-( const OperatorMatrix &mat, const std::string &var ) {
    OperatorMatrix ret( mat.rows(), mat.cols() );
    for ( int i = 0; i < mat.rows(); i++ ) {
        for ( int j = 0; j < mat.cols(); j++ ) {
            if ( mat.at( i, j ).front() == '0' )
                ret.set( i, j, "-" + var );
            else
                ret.set( i, j, mat.at( i, j ) + "-" + var );
        }
    }
    return ret;
}
OperatorMatrix operator*( const OperatorMatrix &mat, const std::string &var ) {
    OperatorMatrix ret( mat.rows(), mat.cols() );
    for ( int i = 0; i < mat.rows(); i++ ) {
        for ( int j = 0; j < mat.cols(); j++ ) {
            if ( mat.at( i, j ).front() == '0' )
                ret.set( i, j, "0" );
            else
                ret.set( i, j, mat.at( i, j ) + var );
        }
    }
    return ret;
}
OperatorMatrix operator/( const OperatorMatrix &mat, const std::string &var ) {
    OperatorMatrix ret( mat.rows(), mat.cols() );
    for ( int i = 0; i < mat.rows(); i++ ) {
        for ( int j = 0; j < mat.cols(); j++ ) {
            if ( mat.at( i, j ).front() == '0' )
                ret.set( i, j, "0" );
            else
                ret.set( i, j, mat.at( i, j ) + "/" + var );
        }
    }
    return ret;
}

OperatorMatrix operator+( const std::string &var, const OperatorMatrix &mat ) {
    OperatorMatrix ret( mat.rows(), mat.cols() );
    for ( int i = 0; i < mat.rows(); i++ ) {
        for ( int j = 0; j < mat.cols(); j++ ) {
            if ( mat.at( i, j ).front() == '0' )
                ret.set( i, j, var );
            else
                ret.set( i, j, var + "+" + mat.at( i, j ) );
        }
    }
    return ret;
}
OperatorMatrix operator-( const std::string &var, const OperatorMatrix &mat ) {
    OperatorMatrix ret( mat.rows(), mat.cols() );
    for ( int i = 0; i < mat.rows(); i++ ) {
        for ( int j = 0; j < mat.cols(); j++ ) {
            if ( mat.at( i, j ).front() == '0' )
                ret.set( i, j, var );
            else
                ret.set( i, j, var + "-" + mat.at( i, j ) );
        }
    }
    return ret;
}
OperatorMatrix operator*( const std::string &var, const OperatorMatrix &mat ) {
    OperatorMatrix ret( mat.rows(), mat.cols() );
    for ( int i = 0; i < mat.rows(); i++ ) {
        for ( int j = 0; j < mat.cols(); j++ ) {
            if ( mat.at( i, j ).front() == '0' )
                ret.set( i, j, "0" );
            else
                ret.set( i, j, var + mat.at( i, j ) );
        }
    }
    return ret;
}
OperatorMatrix operator/( const std::string &var, const OperatorMatrix &mat ) {
    OperatorMatrix ret( mat.rows(), mat.cols() );
    for ( int i = 0; i < mat.rows(); i++ ) {
        for ( int j = 0; j < mat.cols(); j++ ) {
            if ( mat.at( i, j ).front() == '0' )
                ret.set( i, j, "0" );
            else
                ret.set( i, j, var + "/" + mat.at( i, j ) );
        }
    }
    return ret;
}
OperatorMatrix OperatorMatricesText::commutator( OperatorMatrix &left, OperatorMatrix &right ) {
    auto commutator_1 = left * right;
    auto commutator_2 = right * left;
    auto result = commutator_1.concat( commutator_2.addElement( " - (", -1 ).addElement( ")", 1 ) );
    return result;
}
void OperatorMatrix::cleanup() {
    // Clean up "||"
    for ( int i = 0; i < rows(); i++ ) {
        for ( int j = 0; j < cols(); j++ ) {
            std::string current = at( i, j );
            for ( int c = current.size() - 2; c > 0; c-- ) {
                if ( current.at( c ) == current.at( c + 1 ) ) {
                    current.erase( current.begin() + c, current.begin() + c + 1 );
                }
            }
            set( i, j, current );
        }
    }
    // Clean up <i|i> = 1
    for ( int i = 0; i < rows(); i++ ) {
        for ( int j = 0; j < cols(); j++ ) {
            std::string current = at( i, j );
            for ( int c = current.size() - 4; c > 1; c-- ) {
                if ( current.at( c ) == current.at( c + 2 ) && current.at( c ) != '>' && current.at( c ) != '<' && current.at( c ) != ' ' ) {
                    current.erase( current.begin() + c - 2, current.begin() + c + 3 );
                }
            }
            set( i, j, current );
        }
    }
}
void OperatorMatrix::toTEX() {
    for ( int i = 0; i < rows(); i++ ) {
        for ( int j = 0; j < cols(); j++ ) {
            std::string newstr = "";
            std::string current = at( i, j );
            for ( int c = 0; c < current.size(); c++ ) {
                if ( current.at( c ) == '>' ) {
                    // Check for pre '|'
                    for ( int h = c; h > 0; h-- ) {
                        //if (current.at(h) == '|')
                    }
                    current.erase( current.begin() + c, current.begin() + c + 1 );
                }
            }
            set( i, j, current );
        }
    }
}

void OperatorMatricesText::generateOperators( const Parameters &p ) {
    // Base for electronic states is <G|, <H|, <V|, <B| (rows), N^2
    atom_base = OperatorMatrix( {{"|G><G|", "|H><G|", "|V><G|", "|B><G|"}, {"|G><H|", "|H><H|", "|V><H|", "|B><H|"}, {"|G><V|", "|H><V|", "|V><V|", "|B><V|"}, {"|G><B|", "|H><B|", "|V><B|", "|B><B|"}} );
    //bare_atom_base = OperatorMatrix( {{"1", "0", "0", "0"}, {"0", "1", "0", "0"}, {"0", "0", "1", "0"}, {"0", "0", "0", "1"}} );
    bare_atom_state_ground = OperatorMatrix( {{"|G><G|", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}} );
    bare_atom_state_H = OperatorMatrix( {{"0", "0", "0", "0"}, {"0", "|H><H|", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}} );
    bare_atom_state_V = OperatorMatrix( {{"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "|V><V|", "0"}, {"0", "0", "0", "0"}} );
    bare_atom_state_biexciton = OperatorMatrix( {{"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "|B><B|"}} );
    bare_atom_sigmaplus_G_H = OperatorMatrix( {{"0", "0", "0", "0"}, {"|H><G|", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}} );
    bare_atom_sigmaminus_G_H = OperatorMatrix( {{"0", "|G><H|", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}} );
    bare_atom_sigmaplus_G_V = OperatorMatrix( {{"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"|H><V|", "0", "0", "0"}, {"0", "0", "0", "0"}} );
    bare_atom_sigmaminus_G_V = OperatorMatrix( {{"0", "0", "|V><G|", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}} );
    bare_atom_sigmaplus_H_B = OperatorMatrix( {{"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "|B><H|", "0", "0"}} );
    bare_atom_sigmaminus_H_B = OperatorMatrix( {{"0", "0", "0", "0"}, {"0", "0", "0", "|H><B|"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}} );
    bare_atom_sigmaplus_V_B = OperatorMatrix( {{"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "|B><V|", "0"}} );
    bare_atom_sigmaminus_V_B = OperatorMatrix( {{"0", "0", "0", "0"}, {"0", "0", "0", "0"}, {"0", "0", "0", "|H><B|"}, {"0", "0", "0", "0"}} );
    // Base for photonic state is <0|, <1|, <2| (rows)
    // 3 Photons
    //OperatorMatrix photon_base_h( {{"|0_H><0_H|", "|0_H><1_H|", "|0_H><2_H|", "|0_H><3_H|"}, {"|0_H><1_H|", "|1_H><1_H|", "|2_H><1_H|", "|3_H><1_H|"}, {"|0_H><2_H|", "|1_H><2_H|", "|2_H><2_H|", "|3_H><2_H|"}, {"|0_H><3_H|", "|1_H><3_H|", "|2_H><3_H|", "|3_H><3_H|"}} );
    //OperatorMatrix photon_base_v( {{"|0_V><0_V|", "|0_V><1_V|", "|0_V><2_V|", "|0_V><3_V|"}, {"|0_V><1_V|", "|1_V><1_V|", "|2_V><1_V|", "|3_V><1_V|"}, {"|0_V><2_V|", "|1_V><2_V|", "|2_V><2_V|", "|3_V><2_V|"}, {"|0_V><3_V|", "|1_V><3_V|", "|2_V><3_V|", "|3_V><3_V|"}} );
    //OperatorMatrix photon_create_H( {{"0", "|1_H><0_H|", "0", "0"}, {"0", "0", "sqrt(2)*|2_H><1_H|", "0"}, {"0", "0", "0", "sqrt(3)|3_H><2_H|"}, {"0", "0", "0", "0"}} );
    //OperatorMatrix photon_annihilate_H( {{"0", "0", "0", "0"}, {"|0_H><1_H|", "0", "0", "0"}, {"0", "sqrt(2)*|1_H><2_H|", "0", "0"}, {"0", "0", "sqrt(3)|2_H><3_H|", "0"}} );
    //OperatorMatrix photon_create_V( {{"0", "|1_V><0_V|", "0", "0"}, {"0", "0", "sqrt(2)*|2_V><1_V|", "0"}, {"0", "0", "0", "sqrt(3)|3_V><2_V|"}, {"0", "0", "0", "0"}} );
    //OperatorMatrix photon_annihilate_V( {{"0", "0", "0", "0"}, {"|0_V><1_V|", "0", "0", "0"}, {"0", "sqrt(2)*|1_V><2_V|", "0", "0"}, {"0", "0", "sqrt(3)|2_V><3_V|", "0"}} );
    photon_base_h = OperatorMatrix( {{"|0_H><0_H|", "|0_H><1_H|", "|0_H><2_H|"}, {"|0_H><1_H|", "|1_H><1_H|", "|2_H><1_H|"}, {"|0_H><2_H|", "|1_H><2_H|", "|2_H><2_H|"}} );
    photon_base_v = OperatorMatrix( {{"|0_V><0_V|", "|0_V><1_V|", "|0_V><2_V|"}, {"|0_V><1_V|", "|1_V><1_V|", "|2_V><1_V|"}, {"|0_V><2_V|", "|1_V><2_V|", "|2_V><2_V|"}} );
    bare_photon_create_H = OperatorMatrix( {{"0", "|1_H><0_H|", "0"}, {"0", "0", "sqrt(2)*|2_H><1_H|"}, {"0", "0", "0"}} );
    bare_photon_annihilate_H = OperatorMatrix( {{"0", "0", "0"}, {"|0_H><1_H|", "0", "0"}, {"0", "sqrt(2)*|1_H><2_H|", "0"}} );
    bare_photon_create_V = OperatorMatrix( {{"0", "|1_V><0_V|", "0"}, {"0", "0", "sqrt(2)*|2_V><1_V|"}, {"0", "0", "0"}} );
    bare_photon_annihilate_V = OperatorMatrix( {{"0", "0", "0"}, {"|0_V><1_V|", "0", "0"}, {"0", "sqrt(2)*|1_V><2_V|", "0"}} );
    bare_photon_n_H = bare_photon_create_H * bare_photon_annihilate_V;
    bare_photon_n_V = bare_photon_create_V * bare_photon_annihilate_V;

    atom_state_ground = tensor( photon_base_h, photon_base_v, bare_atom_state_ground );
    atom_state_biexciton = tensor( photon_base_h, photon_base_v, bare_atom_state_biexciton );
    atom_state_H = tensor( photon_base_h, photon_base_v, bare_atom_state_H );
    atom_state_V = tensor( photon_base_h, photon_base_v, bare_atom_state_V );
    atom_sigmaplus_G_H = tensor( photon_base_h, photon_base_v, bare_atom_sigmaplus_G_H );
    atom_sigmaminus_G_H = tensor( photon_base_h, photon_base_v, bare_atom_sigmaminus_G_H );
    atom_sigmaplus_H_B = tensor( photon_base_h, photon_base_v, bare_atom_sigmaplus_H_B );
    atom_sigmaminus_H_B = tensor( photon_base_h, photon_base_v, bare_atom_sigmaminus_H_B );
    atom_sigmaplus_G_V = tensor( photon_base_h, photon_base_v, bare_atom_sigmaplus_G_V );
    atom_sigmaminus_G_V = tensor( photon_base_h, photon_base_v, bare_atom_sigmaminus_G_V );
    atom_sigmaplus_V_B = tensor( photon_base_h, photon_base_v, bare_atom_sigmaplus_V_B );
    atom_sigmaminus_V_B = tensor( photon_base_h, photon_base_v, bare_atom_sigmaminus_V_B );
    photon_create_H = tensor( bare_photon_create_H, photon_base_v, atom_base );
    photon_create_V = tensor( photon_base_h, bare_photon_create_V, atom_base );
    photon_annihilate_H = tensor( bare_photon_annihilate_H, photon_base_v, atom_base );

    photon_annihilate_V = tensor( photon_base_h, bare_photon_annihilate_V, atom_base );
    photon_n_H = tensor( bare_photon_n_H, photon_base_v, atom_base );
    photon_n_V = tensor( photon_base_h, bare_photon_n_V, atom_base );

    complete_base = tensor( tensor( photon_base_h, photon_base_v ), atom_base );
    auto H_0 = "E_H" * atom_state_H + "E_V" * atom_state_V + "E_B" * atom_state_biexciton + "E_C" * photon_n_H + "E_C" * photon_n_V;
    auto H_I = "g" * ( atom_sigmaplus_G_H * photon_annihilate_H + atom_sigmaminus_G_H * photon_create_H + atom_sigmaplus_G_V * photon_annihilate_V + atom_sigmaminus_G_V * photon_create_V + atom_sigmaplus_H_B * photon_annihilate_H + atom_sigmaminus_H_B * photon_create_H + atom_sigmaplus_V_B * photon_annihilate_V + atom_sigmaminus_V_B * photon_create_V );
    auto H = H_0 + H_I;
    // Test with only electronic system
    auto rho = complete_base.addElement( "rho_(", -1 );
    rho = rho.addElement( ")", 1 );
    auto result = rho.addElement( "d/dt ", -1 );
    result = result.addElement( " = -i ( ", 1 );
    auto comm = commutator( H, rho );
    result = result.concat( comm );
    result = result.addElement( " )", 1 );
    result.cleanup();
    std::cout << result.at( 0, 0 ) << std::endl;
}

void OperatorMatricesText::generateOperatorsSimple( const Parameters &p ) {
    // Base for electronic states is <G|, <E|
    auto atom_base = OperatorMatrix( {{"|G><G|", "|E><G|"}, {"|G><E|", "|E><E|"}} );
    //bare_atom_base = OperatorMatrix( {{"1", "0", "0", "0"}, {"0", "1", "0", "0"}, {"0", "0", "1", "0"}, {"0", "0", "0", "1"}} );
    auto bare_atom_state_ground = OperatorMatrix( {{"|G><G|", "0"}, {"0", "0"}} );
    auto bare_atom_state_excited = OperatorMatrix( {{"0", "0"}, {"0", "|E><E|"}} );
    auto bare_atom_sigmaplus = OperatorMatrix( {{"0", "0"}, {"|E><G|", "0"}} );
    auto bare_atom_sigmaminus = OperatorMatrix( {{"0", "|G><E|"}, {"0", "0"}} );
    // Base for photonic state is <0|, <1| (rows)
    auto photon_base = OperatorMatrix( {{"|0><0|", "|0><1|"}, {"|0><1|", "|1><1|"}} );
    auto bare_photon_create = OperatorMatrix( {{"0", "|1><0|"}, {"0", "0"}} );
    auto bare_photon_annihilate = OperatorMatrix( {{"0", "0"}, {"|0><1|", "0"}} );
    auto bare_photon_n = bare_photon_create * bare_photon_annihilate;

    auto atom_state_ground = tensor( photon_base, bare_atom_state_ground );
    auto atom_state_excited = tensor( photon_base, bare_atom_state_excited );
    auto atom_sigmaplus = tensor( photon_base, bare_atom_sigmaplus );
    auto atom_sigmaminus = tensor( photon_base, bare_atom_sigmaminus );

    auto photon_create = tensor( bare_photon_create, atom_base );
    auto photon_annihilate = tensor( bare_photon_annihilate, atom_base );
    auto photon_n = tensor( bare_photon_n, atom_base );

    auto complete_base = tensor( photon_base, atom_base );
    auto H_0 = "E" * atom_state_excited + "E_C" * photon_n;
    auto H_I = "g" * ( atom_sigmaplus * photon_annihilate + atom_sigmaminus * photon_create );
    auto H = H_0 + H_I;
    // Test with only electronic system
    auto rho = complete_base.addElement( "rho_(", -1 );
    rho = rho.addElement( ")", 1 );
    auto result = rho.addElement( "\frac{\text{d}}{\text{dt}} ", -1 );
    result = result.addElement( " = -i ( ", 1 );
    auto comm = commutator( H, rho );
    result = result.concat( comm );
    result = result.addElement( " )", 1 );
    result.cleanup();
    std::cout << result.at( 0, 0 ) << std::endl;
    //std::cout << test << std::endl;
    //std::cout << photon_n_H << std::endl;
}