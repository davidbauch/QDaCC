#pragma once

// Includes standard numerical operations such as Kommutator, Lindbladian, ...
class Operatable {
   public:
    // Constructor
    Operatable(){};

    // Functions to determine numerical matrix trace from input types supporting .trace()
    // @param &mat: Input dense matrix
    // @return Matrix trace
    template <typename T, typename M>
    inline T getTrace( const M &mat ) const {
        return mat.trace();
    }
    // Functions to determine numerical matrix trace from sparse input types not supporting .trace()
    // @param &mat: Input sparse matrix
    // @return Matrix trace of type T
    template <typename T>
    inline T getTrace( const Sparse &mat ) const {
        //return getTrace<T>( Dense( mat ) );
        T ret = (T)0.0;
        for ( int k = 0; k < mat.outerSize(); ++k ) {
            for ( Sparse::InnerIterator it( mat, k ); it; ++it ) {
                if ( it.row() == it.col() ) {
                    ret += it.value();
                }
            }
        }
        return ret;
    }

    // Calculates the Kommutator of two matrices of identical type
    // @param &A,&B: Input matrices
    // @return Kommutator of A and B where [A,B] = AB-BA
    template <typename T>
    inline T dgl_kommutator( const T &A, const T &B ) {
        return A * B - B * A;
    }

    // Calculates the Anti-Kommutator of two matrices of identical type
    // @param &A,&B: Input matrices
    // @return Anti-Kommutator of A and B where [A,B]^- = BA-AB
    template <typename T>
    inline T dgl_antikommutator( const T &A, const T &B ) {
        return A * B + B * A;
    }

    // Calculates the Lindbladian of two input matrices where L = 2*A*rho*B - B*A*rho - rho*B*A
    // @param &rho: Input density matrix
    // @param &op: Input matrix O
    // @param &opd: Input matrix O^dagger
    // @return Returns the Lindbladian of rho, O and O^dagger
    template <typename T, typename T2, typename T3>
    inline T dgl_lindblad( const T &rho, const T2 &op, const T3 &opd ) {
        return 2.0 * op * rho * opd - opd * op * rho - rho * opd * op;
    }
};