#include "system/savestate.h"

bool QDLC::Save_State_sort_t( const QDLC::SaveStateTau &ss1, const QDLC::SaveStateTau &ss2 ) {
    return ( ss1.t < ss2.t );
}
bool QDLC::Save_State_sort_tau( const QDLC::SaveStateTau &ss1, const QDLC::SaveStateTau &ss2 ) {
    return ( ss1.tau < ss2.tau );
}