#pragma once
// Dependencies
#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "misc/timer.h"
#include "system/parameters.h"

/**
 * @brief This class stores all the State- and Transition Matrices.
 * It also provides all of the neccessary methods for the creation of these Matrices.
 *
 */
class OperatorMatrices {
   public:
    /**
     * @brief Information Wrapper for single operators:
     *
     */
    struct matrix_s {
        // Bra and Ket Operators in self-space
        Dense bra, ket;
        // Matrix representation Ket*Bra in self Hilbert space
        Dense self_hilbert;
        // Matrix representation Ket*Bra in total Hilbert space
        Sparse hilbert;
        // Projector Matrix. This is a matrix int the total Hilbert space where every nonzero entry in the original hilbert matrix is 1.
        Sparse projector;
        // Integer index for base. Used for easy partial tracing.
        int base;
        // "Creator" or "Annihilator" style operator
        int direction;
        // State Energy or Transition Energy
        double energy;
        // Default Constructor
        matrix_s() : direction( 0 ), base( -1 ), energy( 0 ){};
    };

   public:
    // Operator Matrix Classes
    // Self-Part of the total Hamiltonian
    Sparse H_0;
    // Interaction Part of the total Hamiltonian
    Sparse H_I;
    // Hamilton Operator Used for numerical calculations. This is what the system's dgl_hamilton method will use.
    Sparse H_used;
    // Density Matrix
    Sparse rho;
    // Identity Matrix in the total Hilbert dimension
    Sparse identity;

    // Base Vector mapping the matrix index onto the corresponding state string in total Hilbert space
    std::vector<std::string> base;
    // Contains the individual self-Hilbert bases
    std::vector<Dense> base_selfhilbert;
    // Maps the individual self-Hilbert indices onto the corresponding total Hilbert space indices. Used to calculate partial traces. Key is the systems base integer index
    std::vector<Dense> base_hilbert_index;
    // Maps the Key index string |a|b|...> onto an integer index
    std::map<std::string, int> base_index_map;
    // Maps the total Hilbert space index onto the corresponding phonon group index
    std::vector<int> phonon_hilbert_index_to_group_index;
    // Maps the the phonon group index onto a group of total Hilbert space indices
    std::vector<std::vector<int>> phonon_group_index_to_hilbert_indices;
    // Maps the phonon group index onto their corresponding coupling value
    std::vector<double> phonon_group_index_to_coupling_value;
    // Phonon coupling Matrix for the PME
    Sparse polaron_phonon_coupling_matrix;
    // Cached PME Rates
    std::map<double, Sparse> phi_vector_matrix_cache_u, phi_vector_matrix_cache_g;

    // QDLC 3.0 New System Matrices. Since Version 3.0, the electronic and optical system is not hardcoded into the program anymore. Instead, the state- and transition matrices get generated on startup and stored in these maps.
    std::map<std::string, matrix_s> el_states, ph_states, el_transitions, ph_transitions, extra_transitions;
    // PME precalculated Chi and partial chi / partial t sumands
    std::vector<Sparse> polaron_factors, polaron_pulse_factors_explicit_time;
    // Cache matrix for the timetransformation matrix
    Dense timetrafo_cachematrix;
    // Cache matrices for the pulse and chirp
    std::vector<Sparse> pulse_mat, chirp_mat;
    // The initial State vector
    Dense initial_state_vector_ket;

    // A global output format for Eigen's matrices. Could also be declared to Eigen at compile time.
    Eigen::IOFormat output_format;

    /**
     * @brief Constructor of the OperatorMatrices Class
     *
     */
    OperatorMatrices(){};
    OperatorMatrices( Parameters &p );

    /**
     * @brief Generates the Operators, e.g. creates all the State- and Transition Matrices
     *
     * @return true On Success
     * @return false On Failure
     */
    bool generate_operators( Parameters &p );

    /**
     * @brief Output operators to loglevel 2
     *
     */
    void output_operators( Parameters &p );

    /**
     * @brief Creates a bosonic creation or annihilation operator matrix
     * The Type can be either OPERATOR_PHOTONIC_CREATE or OPERATOR_PHOTONIC_ANNIHILATE
     * The maximum number of photons is given by maxPhotons; The resulting matrix will have dimension (maxPhotons+1)*(maxPhotons+1)
     * @return Either creation or annihilation matrix of type M
     */
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