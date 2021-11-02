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
        Dense bra, ket;     // Bra and Ket Operators in self-space
        Dense self_hilbert; // Matrix representation Ket*Bra in self Hilbert space
        Sparse hilbert;     // Matrix representation Ket*Bra in total Hilbert space
        Sparse projector;   // Projector Matrix
        int base;           // Integer index for base
        int direction;      // "Creator" or "Annihilator" style operator
        double energy;      // State Energy or Transition Energy
        matrix_s() : direction( 0 ), base( -1 ), energy( 0 ){};
    };

   public:
    // Operator Matrix Classes
    Sparse H;
    Sparse H_0;
    Sparse H_I;
    Sparse rho;
    Sparse H_used;
    // Basis and maps:
    std::vector<std::string> base;             // Contains all the string bases vectors in total Hilbert space
    std::vector<Dense> base_selfhilbert;       // Contains the individual self-Hilbert bases
    std::vector<Dense> base_hilbert_index;     // Maps the individual self-Hilbert indices onto the corresponding total Hilbert space indices. Used to calculate partial traces. Key is the systems base integer index
    std::map<std::string, int> base_index_map; // Maps the Key index string |a|b|...> onto an integer index
    std::vector<int> phononCouplingIndex;
    std::vector<double> phononCouplingIndexValue;

    // 3.0 New System Matrices
    std::map<std::string, matrix_s> el_states, ph_states, el_transitions, ph_transitions, extra_transitions;
    std::vector<Sparse> lindblad_factors;
    std::vector<Sparse> polaron_factors, polaron_pulse_factors_explicit_time;
    Dense timetrafo_cachematrix;
    std::vector<Sparse> pulse_mat, chirp_mat;
    Dense initial_state_vector_ket;

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