#pragma once

#include "typedef.h"

namespace QDACC {

// Vector for mat/time, tuple
class SaveState {
   public:
    QDACC::Type::MatrixMain mat;
    double t;
    SaveState(){};
    SaveState( const QDACC::Type::MatrixMain &mat, const double time ) : mat( mat ), t( time ){};
};
class SaveStateTau {
   public:
    QDACC::Type::MatrixMain mat1, mat2;
    double t, tau;
    SaveStateTau( const QDACC::Type::MatrixMain &mat1, const QDACC::Type::MatrixMain &mat2, const double t, const double tau = 0.0 ) : mat1( mat1 ), mat2( mat2 ), t( t ), tau( tau ){};
    SaveStateTau( const QDACC::Type::MatrixMain &mat, const double t ) : mat1( mat ), t( t ), tau( 0 ){};
    SaveStateTau(){};
};
class SaveScalar {
   public:
    QDACC::Type::Scalar scalar;
    double t, tau;
    SaveScalar( const QDACC::Type::Scalar &scalar, const double t, const double tau = 0.0 ) : scalar( scalar ), t( t ), tau( tau ){};
};

bool Save_State_sort_t( const SaveStateTau &ss1, const SaveStateTau &ss2 );
bool Save_State_sort_tau( const SaveStateTau &ss1, const SaveStateTau &ss2 );

} // namespace QDACC