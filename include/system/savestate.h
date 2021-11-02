#pragma once 

#include "typedef.h"

namespace QDLC {

// Vector for mat/time, tuple
class SaveState {
   public:
    QDLC::Type::Sparse mat;
    double t;
    SaveState( const QDLC::Type::Sparse &mat, const double time ) : mat( mat ), t( time ){};
};
class SaveStateTau {
   public:
    QDLC::Type::Sparse mat1, mat2;
    double t, tau;
    SaveStateTau( const QDLC::Type::Sparse &mat1, const QDLC::Type::Sparse &mat2, const double t, const double tau ) : mat1( mat1 ), mat2( mat2 ), t( t ), tau( tau ){};
    SaveStateTau( const QDLC::Type::Sparse &mat, const double t ) : mat1( mat ), t( t ), tau( 0 ){};
    SaveStateTau(){};
};
class SaveScalar {
   public:
    QDLC::Type::Scalar scalar;
    double t, tau;
    SaveScalar( const QDLC::Type::Scalar &scalar, const double t, const double tau = 0.0 ) : scalar( scalar ), t( t ), tau( tau ){};
};

bool Save_State_sort_t( const SaveStateTau &ss1, const SaveStateTau &ss2 );
bool Save_State_sort_tau( const SaveStateTau &ss1, const SaveStateTau &ss2 );

}