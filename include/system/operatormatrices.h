#pragma once
// Dependencies
#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "misc/timer.h"
#include "system/parameters.h"

class OperatorMatrices {
    // Information Wrapper for single operators:
   public:
    struct matrix_s {
        Dense bra, ket;
        Dense self_hilbert;
        Sparse hilbert;
        Sparse projector;
        int base;
        int direction;
        double energy; // State Energy or Transition Energy
    };

   public:
    // Operator Matrix Classes
    Sparse H;
    Sparse H_0;
    Sparse H_I;
    Sparse rho;
    Sparse H_used;
    std::vector<Sparse> pulse_mat, chirp_mat;
    std::vector<std::string> base;
    std::map<std::string, int> base_index_map; // Maps the index string |a|b|...> onto an integer index
    //std::vector<double> phononCouplingFactor;
    Dense phononCouplingFactor;

    // 3.0 New System Matrices
    std::map<std::string, matrix_s> el_states, ph_states, el_transitions, ph_transitions;
    std::vector<Sparse> lindblad_factors;
    std::vector<Sparse> polaron_factors;
    Dense timetrafo_cachematrix;

    // Constructor
    OperatorMatrices(){};
    OperatorMatrices( Parameters &p );

    // Generating operators
    // @param &p: Parameter Class Reference
    // @return Returns True if generating operators was successful
    bool generateOperators( Parameters &p );

    // Output operators to loglevel 2
    // @param &p: Parameter Class Reference
    void outputOperators( Parameters &p );

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

    template <class M>
    static M tensor( const std::vector<M> &m ) {
        if ( m.size() < 2 )
            return m.front();
        else if ( m.size() == 2 )
            return tensor<M>( m[0], m[1] );
        auto md = m;
        md.pop_back();
        return tensor<M>( tensor<M>( md ), m.back() );
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

    static std::vector<std::string> tensor( const std::vector<std::vector<std::string>> &m ) {
        if ( m.size() < 2 )
            return m.front();
        else if ( m.size() == 2 )
            return tensor( m[0], m[1] );
        auto md = m;
        md.pop_back();
        return tensor( tensor( md ), m.back() );
    }

    // Creates a bosonic creation or annihilation operator matrix
    // @param type: Either OPERATOR_PHOTONIC_CREATE of OPERATOR_PHOTONIC_ANNIHILATE
    // @param maxPhotons: Maximum number of photons; Resulting matrix will have dimension (n+1)x(n+1)
    // @return Returns either creation or annihilation matrix of type M
    template <class M>
    static M create_photonic_operator( const int &type, const int &maxPhotons ) {
        M ret = M::Zero( maxPhotons + 1, maxPhotons + 1 );
        for ( int i = 0; i < maxPhotons; i++ ) {
            if ( type == OPERATOR_PHOTONIC_CREATE )
                ret( i + 1, i ) = sqrt( i + 1 );
            else if ( type == OPERATOR_PHOTONIC_ANNIHILATE )
                ret( i, i + 1 ) = sqrt( i + 1 );
            else if ( type == OPERATOR_PHOTONIC_STATE )
                ret( i + 1, i + 1 ) = i + 1;
        }
        return ret;
    }
};