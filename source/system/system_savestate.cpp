#include "system/savestate.h"

bool QDACC::Save_State_sort_t( const QDACC::SaveStateTau &ss1, const QDACC::SaveStateTau &ss2 ) {
    return ( ss1.t < ss2.t );
}
bool QDACC::Save_State_sort_tau( const QDACC::SaveStateTau &ss1, const QDACC::SaveStateTau &ss2 ) {
    return ( ss1.tau < ss2.tau );
}