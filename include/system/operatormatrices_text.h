#pragma once

#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "system/parameters.h"

class OperatorMatrix {
   public:
    std::vector<std::vector<std::string>> mat;
    // Rows, Cols
    OperatorMatrix(){};
    OperatorMatrix( const OperatorMatrix &other );
    OperatorMatrix( const std::vector<std::vector<std::string>> &mat ) : mat( mat ){};
    OperatorMatrix( int n, int m );

    int cols() const;
    int rows() const;
    // Row, Col
    std::string operator()( int i, int j ) const;
    // Row, Col
    std::string at( int i, int j ) const;
    // Row, Col
    void set( int i, int j, const std::string &val );
    // Row, Col
    std::string getFormatted( int i, int j ) const;
    void cleanup();
    void toTEX();
    // Matrix-Matrix operations
    OperatorMatrix operator+( const OperatorMatrix &other );
    OperatorMatrix operator-( const OperatorMatrix &other );
    OperatorMatrix operator*( const OperatorMatrix &other );
    OperatorMatrix elementWiseMul( const OperatorMatrix &other );
    OperatorMatrix concat( const OperatorMatrix &other );
    OperatorMatrix addElement( const std::string &element, int before_or_after = 1 );

    friend std::ostream &operator<<( std::ostream &stream, const OperatorMatrix &mat ) {
        int max_length = 0;
        for ( int i = 0; i < mat.rows(); i++ ) {
            for ( int j = 0; j < mat.cols(); j++ ) {
                max_length = std::max<int>( max_length, mat.at( i, j ).size() );
            }
        }
        max_length += 2;
        for ( int i = 0; i < mat.rows(); i++ ) {
            stream << "[";
            for ( int j = 0; j < mat.cols(); j++ ) {
                stream << fmt::format( "{:<{}}", mat.at( i, j ), max_length );
            }
            stream << "]\n";
        }
        return stream;
    }
};

// Calculates the tensor product of string matrices a and b
// @param &a,&b: Input matrices
// @return Returns a x b where x is the tensor product
static OperatorMatrix tensor( const OperatorMatrix &a, const OperatorMatrix &b ) {
    assert( a.rows() == a.cols() && b.rows() == b.cols() && "Only Square matrices accepted" );
    OperatorMatrix ret = OperatorMatrix( a.rows() * b.rows(), a.cols() * b.cols() );
    for ( int i = 0; i < a.rows(); i++ )
        for ( int j = 0; j < a.cols(); j++ )
            for ( int k = 0; k < b.rows(); k++ )
                for ( int l = 0; l < b.cols(); l++ ) {
                    if ( a.at( i, j ).front() != '0' && b.at( k, l ).front() != '0' )
                        ret.set( i * b.rows() + k, j * b.cols() + l, a( i, j ) + b( k, l ) );
                    else
                        ret.set( i * b.rows() + k, j * b.cols() + l, "0" );
                }
    return ret;
}

static OperatorMatrix tensor( const OperatorMatrix &a, const OperatorMatrix &b, const OperatorMatrix &c ) {
    return tensor( tensor( a, b ), c );
}

// This class is insanely inefficient but thats ok.
class OperatorMatricesText {
   public:
    // Base And Expanded Matrices
    //Operator Matrices
    OperatorMatrix atom_base, photon_base_h, photon_base_v, complete_base;
    OperatorMatrix bare_atom_base, bare_photon_base_h, bare_photon_base_v; // used for tensor product of other bare matrices
    OperatorMatrix photon_create_H, photon_annihilate_H, photon_create_V, photon_annihilate_V;
    OperatorMatrix atom_state_biexciton, atom_state_ground, atom_state_H, atom_state_V, photon_n_H, photon_n_V;
    OperatorMatrix atom_sigmaplus_G_H, atom_sigmaminus_G_H, atom_sigmaplus_H_B, atom_sigmaminus_H_B, atom_sigmaplus_G_V, atom_sigmaminus_G_V, atom_sigmaplus_V_B, atom_sigmaminus_V_B;
    OperatorMatrix atom_inversion_G_H, atom_inversion_G_V, atom_inversion_H_B, atom_inversion_V_B, atom_inversion_G_B;
    OperatorMatrix atom_sigmaminus_G_B;

    // Bare matrices:
    OperatorMatrix bare_photon_create_H, bare_photon_annihilate_H, bare_photon_create_V, bare_photon_annihilate_V;
    OperatorMatrix bare_atom_state_biexciton, bare_atom_state_ground, bare_atom_state_H, bare_atom_state_V, bare_photon_n_H, bare_photon_n_V;
    OperatorMatrix bare_atom_sigmaplus_G_H, bare_atom_sigmaminus_G_H, bare_atom_sigmaplus_H_B, bare_atom_sigmaminus_H_B, bare_atom_sigmaplus_G_V, bare_atom_sigmaminus_G_V, bare_atom_sigmaplus_V_B, bare_atom_sigmaminus_V_B;
    OperatorMatrix bare_atom_inversion_G_H, bare_atom_inversion_G_V, bare_atom_inversion_H_B, bare_atom_inversion_V_B, bare_atom_inversion_G_B;
    OperatorMatrix bare_atom_sigmaminus_G_B;

    // Hamiltons:
    //OperatorMatrix H_0;

    OperatorMatricesText(){};
    void generateOperators( const Parameters &p );
    void generateOperatorsSimple( const Parameters &p );
    void outputOperators( const Parameters &p );
    OperatorMatrix commutator( OperatorMatrix &left, OperatorMatrix &right );
};