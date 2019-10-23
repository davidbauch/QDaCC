#pragma once
// Dependencies
#include <stdlib.h>
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
    static M expand_photonic_operator( const M &op, int N ) {
        std::assert( op.cols() == op.rows() && "Rows must be equal to cols" );
        M ret = M::Zero( op.rows() * N, op.rows() * N );
        for ( int i = 0; i < op.cols(); i++ )
            for ( int j = 0; j < op.rows(); j++ )
                for ( int n = 0; n < N; n++ )
                    for ( int m = 0; m < N; m++ )
                        if ( n == m )
                            ret( i * N + n, j * N + m ) = op( i, j );
        return ret;
    }
    template <class M>
    static M expand_atomic_operator( const M &op, int N ) {
        std::assert( op.cols() == op.rows() && "Rows must be equal to cols" );
        M ret = M::Zero( op.rows() * ( N + 1 ), op.rows() * ( N + 1 ) );
        for ( int n = 0; n <= N; n++ )
            for ( int i = 0; i < op.cols(); i++ )
                for ( int j = 0; j < op.rows(); j++ )
                    ret( n * op.cols() + i, n * op.rows() + j ) = op( i, j );
        return ret;
    }
    template <class M>
    static M tensor( const M &a, const M &b ) {
        std::assert( a.rows() == a.cols() && b.rows() == b.cols() );
        M ret = M::Zero( a.cols() * b.cols(), a.rows() * b.rows() );
        for ( int i = 0; i < a.cols(); i++ )
            for ( int j = 0; j < a.rows(); j++ )
                for ( int k = 0; k < b.cols(); k++ )
                    for ( int l = 0; l < b.rows(); l++ ) {
                        ret( i * a.cols() + k, j * a.rows() + l ) = a( i, j ) * b( k, l );
                    }
        return ret;
    }
    template <class M>
    static M create_photonic_operator( const int type, const int maxPhotons ) {
        M ret = M::Zero( maxPhotons + 1, maxPhotons + 1 );
        for ( int i = 0; i < maxPhotons; i++ ) {
            if ( type == OPERATOR_PHOTONIC_CREATE )
                ret( i + 1, i ) = sqrt( i + 1 );
            if ( type == OPERATOR_PHOTONIC_ANNIHILATE )
                ret( i, i + 1 ) = sqrt( i + 1 );
        }
        return ret;
    }
    template <class M>
    static M map_operator( const M &op, const M &map ) {
        std::assert( op.cols() == op.rows() && map.cols() == map.rows() && "Rows must be equal to cols" );
        M ret = M::Zero( op.rows() * ( map.rows() ), op.cols() * ( map.cols() ) );
        for ( int i = 0; i < op.cols(); i++ )
            for ( int j = 0; j < op.rows(); j++ )
                for ( int n = 0; n < map.cols(); n++ )
                    for ( int m = 0; m < map.rows(); m++ )
                        ret( i * map.cols() + n, j * map.rows() + m ) = op( i, j ) * map( n, m );
        return ret;
    }
    // To be called from child constructor
    void init( const Parameters &p );
    // Basic Operator functions:
};