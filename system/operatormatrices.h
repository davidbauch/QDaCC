#pragma once
// Dependencies
#include "../global.h"
#include "../misc/helperfunctions.h"
#include "../misc/log.h"
#include "../misc/timer.h"

class OperatorMatrices_Parent {
   public:
    // Constructor
    OperatorMatrices_Parent(){};
    OperatorMatrices_Parent( const Parameters &p ){};
    // @overwrite: Generating operators
    virtual bool generateOperators( const Parameters &p ) { return false; };
    // @overwrite: Output operators to loglevel 2:
    virtual void outputOperators( const Parameters &p ){};
    // Operator functions
    template <class M>
    static M tensor( const M &a, const M &b ) {
        assert( a.rows() == a.cols() && b.rows() == b.cols() && "Only Square matrices accepted");
        M ret = M::Zero( a.cols() * b.cols(), a.rows() * b.rows() );
        for ( int i = 0; i < a.rows(); i++ )
            for ( int j = 0; j < a.cols(); j++ )
                for ( int k = 0; k < b.rows(); k++ )
                    for ( int l = 0; l < b.cols(); l++ ) {
                        ret( i * b.rows() + k, j * b.cols() + l ) = a( i, j ) * b( k, l );
                    }
        return ret;
    }
    template <class M>
    static M tensor( const M &a, const M &b, const M &c ) {
        return tensor<M>(tensor<M>(a,b),c);
    }
    static std::vector<std::string> tensor( const std::vector<std::string> &a, const std::vector<std::string> &b ) {
        std::vector<std::string> ret;
        for ( int i = 0; i < (int)a.size(); i++ )
                for ( int k = 0; k < (int)b.size(); k++ ) {
                        ret.emplace_back(a.at( i) + "|" + b.at( k));
                    }
        return ret;
    }
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
    // To be called from child constructor
    void init( const Parameters &p );
    // Basic Operator functions:
};