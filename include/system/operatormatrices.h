#pragma once
// Dependencies
#include "global.h"
#include "misc/helperfunctions.h"
#include "misc/log.h"
#include "misc/timer.h"
#include "system/parameters.h"

namespace QDACC {

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
        MatrixMain bra, ket;
        // Matrix representation Ket*Bra in self Hilbert space
        MatrixMain self_hilbert;
        // Matrix representation Ket*Bra in total Hilbert space
        MatrixMain hilbert;
        // Projector Matrix. This is a matrix int the total Hilbert space where every nonzero entry in the original hilbert matrix is 1.
        MatrixMain projector;
        // Integer index for base. Used for easy partial tracing.
        int base;
        // "Creator" or "Annihilator" style operator
        int direction;
        // State Energy or Transition Energy
        double energy;
        // String identifiers. This is either a state or a transition with a "from" and a "to" state.
        std::string name, name_transposed, from, to;
        // Default Constructor
        matrix_s() : direction( 0 ), base( -1 ), energy( 0 ){};
    };

   public:
    // Operator Matrix Classes
    // Self-Part of the total Hamiltonian
    MatrixMain H_0;
    // Interaction Part of the total Hamiltonian
    MatrixMain H_I;
    // Hamilton Operator Used for numerical calculations. This is what the system's dgl_hamilton method will use.
    MatrixMain H_used;
    // Density Matrix
    MatrixMain rho;
    // Identity Matrix in the total Hilbert dimension
    MatrixMain identity, zero;

    // Base Vector mapping the matrix index onto the corresponding state string in total Hilbert space
    std::vector<std::string> base;
    // Contains the individual self-Hilbert bases
    std::vector<MatrixMain> base_selfhilbert;
    // Maps the individual self-Hilbert indices onto the corresponding total Hilbert space indices. Used to calculate partial traces. Key is the systems base integer index
    std::vector<MatrixMain> base_hilbert_index;
    // Maps the Key index string |a|b|...> onto an integer index
    std::map<std::string, int> base_index_map;
    // Maps the total Hilbert space index onto the corresponding phonon group index
    std::vector<int> phonon_hilbert_index_to_group_index;
    // Maps the the phonon group index onto a group of total Hilbert space indices
    std::vector<std::vector<int>> phonon_group_index_to_hilbert_indices;
    // Maps the phonon group index onto their corresponding coupling value
    std::vector<double> phonon_group_index_to_coupling_value;
    // Phonon coupling Matrix for the PME
    MatrixMain polaron_phonon_coupling_matrix;
    // Cached PME Rates
    std::map<double, MatrixMain> pme_greenfunction_matrix_cache_u, pme_greenfunction_matrix_cache_g;

    // QDACC 3.0 New System Matrices. Since Version 3.0, the electronic and optical system is not hardcoded into the program anymore. Instead, the state- and transition matrices get generated on startup and stored in these maps.
    std::map<std::string, matrix_s> el_states, ph_states, el_transitions, ph_transitions, extra_transitions;
    // PME precalculated Chi and partial chi / partial t sumands
    std::vector<MatrixMain> polaron_factors, polaron_pulse_factors_explicit_time;
    // Cache matrix for the timetransformation matrix
    Dense timetrafo_cachematrix;
    // Cache matrices for the pulse and chirp
    std::vector<MatrixMain> pulse_mat, chirp_mat;
    // The initial State vector
    MatrixMain initial_state_vector_ket;

    // Custom expectation value operators
    std::vector<MatrixMain> numerics_custom_expectation_values_operators;

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
     * @brief Output Matrix using the global output format
     *
     */
    std::string matrixToString(const Dense& mat);
    std::string matrixToString(const Sparse& mat);

    /**
     * @brief Creates a bosonic creation or annihilation operator matrix
     * The Type can be either PhotonicOperator::Create or PhotonicOperator::Annihilate
     * The maximum number of photons is given by maxPhotons; The resulting matrix will have dimension (maxPhotons+1)*(maxPhotons+1)
     * @return Either creation or annihilation matrix of type M
     */
    template <class M>
    static M create_photonic_operator( const QDACC::PhotonicOperator &type, const int &maxPhotons ) {
        M ret( maxPhotons + 1, maxPhotons + 1 );
        ret.setZero();
        for ( int i = 0; i < maxPhotons; i++ ) {
            if ( type == QDACC::PhotonicOperator::Create )
                ret.coeffRef( i + 1, i ) = sqrt( i + 1 );
            else if ( type == QDACC::PhotonicOperator::Annihilate )
                ret.coeffRef( i, i + 1 ) = sqrt( i + 1 );
            else if ( type == QDACC::PhotonicOperator::State )
                ret.coeffRef( i + 1, i + 1 ) = i + 1;
        }
        return ret;
    }
};
} // namespace QDACC