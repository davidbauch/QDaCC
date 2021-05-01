#pragma once
// Dependencies
#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "misc/timer.h"
#include "system/parameters.h"

class OperatorMatrices {
   public:
    // Operator Matrix Classes
    Sparse H;
    Sparse H_0;
    Sparse H_I;
    Sparse rho;
    Sparse H_used;
    std::vector<std::string> base;
    //std::vector<double> phononCouplingFactor;
    Dense phononCouplingFactor;

    //Operator Matrices
    Sparse photon_create_H, photon_annihilate_H, photon_create_V, photon_annihilate_V;
    Sparse atom_state_biexciton, atom_state_ground, atom_state_H, atom_state_V, photon_n_H, photon_n_V;
    Sparse atom_sigmaplus_G_H, atom_sigmaminus_G_H, atom_sigmaplus_H_B, atom_sigmaminus_H_B, atom_sigmaplus_G_V, atom_sigmaminus_G_V, atom_sigmaplus_V_B, atom_sigmaminus_V_B;
    Sparse atom_inversion_G_H, atom_inversion_G_V, atom_inversion_H_B, atom_inversion_V_B, atom_inversion_G_B;
    Sparse atom_sigmaminus_G_B;

    // Bare matrices:
    Dense bare_photon_create_H, bare_photon_annihilate_H, bare_photon_create_V, bare_photon_annihilate_V;
    Dense bare_atom_state_biexciton, bare_atom_state_ground, bare_atom_state_H, bare_atom_state_V, bare_photon_n_H, bare_photon_n_V;
    Dense bare_atom_sigmaplus_G_H, bare_atom_sigmaminus_G_H, bare_atom_sigmaplus_H_B, bare_atom_sigmaminus_H_B, bare_atom_sigmaplus_G_V, bare_atom_sigmaminus_G_V, bare_atom_sigmaplus_V_B, bare_atom_sigmaminus_V_B;
    Dense bare_atom_inversion_G_H, bare_atom_inversion_G_V, bare_atom_inversion_H_B, bare_atom_inversion_V_B, bare_atom_inversion_G_B;
    Dense bare_atom_sigmaminus_G_B;

    // Projector Matrices (only few needed)
    Sparse projector_atom_sigmaplus_G_H, projector_atom_sigmaminus_G_H, projector_atom_sigmaplus_H_B, projector_atom_sigmaminus_H_B, projector_atom_sigmaplus_G_V, projector_atom_sigmaminus_G_V, projector_atom_sigmaplus_V_B, projector_atom_sigmaminus_V_B;
    Sparse projector_photon_create_H, projector_photon_annihilate_H, projector_photon_create_V, projector_photon_annihilate_V;

    // Constructor
    OperatorMatrices(){};
    OperatorMatrices( const Parameters &p );

    // Generating operators
    // @param &p: Parameter Class Reference
    // @return Returns True if generating operators was successful
    bool generateOperators( const Parameters &p );

    // Output operators to loglevel 2
    // @param &p: Parameter Class Reference
    void outputOperators( const Parameters &p );

    // ##### Operator functions #####

    // Calculates the tensor product of matrices a and b
    // @param &a,&b: Input matrices
    // @return Returns a x b where x is the tensor product
    template <class M>
    static M tensor( const M &a, const M &b ) {
        assert( a.rows() == a.cols() && b.rows() == b.cols() && "Only Square matrices accepted" );
        M ret = M::Zero( a.cols() * b.cols(), a.rows() * b.rows() );
        for ( int i = 0; i < a.rows(); i++ )
            for ( int j = 0; j < a.cols(); j++ )
                for ( int k = 0; k < b.rows(); k++ )
                    for ( int l = 0; l < b.cols(); l++ ) {
                        ret( i * b.rows() + k, j * b.cols() + l ) = a( i, j ) * b( k, l );
                    }
        return ret;
    }

    // Calculates the tensor product of matrices a, b and c
    // @param &a,&b,&c: Input matrices
    // @return Returns (a x b) x c where x is the tensor product
    template <class M>
    static M tensor( const M &a, const M &b, const M &c ) {
        return tensor<M>( tensor<M>( a, b ), c );
    }

    // Determines the resulting base of the tensor product of two matrices
    // @param &a,&b input string vectors containing the named base
    // @return Returns a string vector containing the combined basis vector names
    static std::vector<std::string> tensor( const std::vector<std::string> &a, const std::vector<std::string> &b ) {
        std::vector<std::string> ret;
        for ( int i = 0; i < (int)a.size(); i++ )
            for ( int k = 0; k < (int)b.size(); k++ ) {
                ret.emplace_back( a.at( i ) + "|" + b.at( k ) );
            }
        return ret;
    }

    // Creates a bosonic creation or annihilation operator matrix
    // @param type: Either OPERATOR_PHOTONIC_CREATE of OPERATOR_PHOTONIC_ANNIHILATE
    // @param maxPhotons: Maximum number of photons; Resulting matrix will have dimension (n+1)x(n+1)
    // @return Returns either creation or annihilation matrix of type M
    template <class M>
    static M create_photonic_operator( int type, int maxPhotons ) {
        M ret = M::Zero( maxPhotons + 1, maxPhotons + 1 );
        for ( int i = 0; i < maxPhotons; i++ ) {
            if ( type == OPERATOR_PHOTONIC_CREATE )
                ret( i + 1, i ) = sqrt( i + 1 );
            if ( type == OPERATOR_PHOTONIC_ANNIHILATE )
                ret( i, i + 1 ) = sqrt( i + 1 );
        }
        return ret;
    }
};