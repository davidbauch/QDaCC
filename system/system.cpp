#include "system.h"

void System_Parent::init() {
    terminate_message = global_message_normaltermination;
    // Adjusting inputs:
    Timer &timer_systeminit = createTimer( "System Initialization" );
    logs.level2( "System initialization... " );
    timer_systeminit.start();
    if ( !init_system() ) {
        logs.level2( "System initialization failed! Exitting program...\n" );
        logs.close();
        exit( EXIT_FAILURE );
    }
    logs.level2( "successful. Elapsed time is {}ms\n", timer_systeminit.getWallTime( TIMER_MILLISECONDS ) );
    timer_systeminit.end();
}

bool System_Parent::calculate_spectrum() {
    return parameters.numerics_calculate_spectrum;
}
bool System_Parent::calculate_g2() {
    return parameters.numerics_calculate_g2;
}
bool System_Parent::use_interactionpicture() {
    return parameters.numerics_use_interactionpicture;
}
bool System_Parent::use_rwa() {
    return parameters.numerics_use_rwa;
}
bool System_Parent::output_handlerstrings() {
    return parameters.output_handlerstrings;
}
bool System_Parent::output_operators() {
    return parameters.output_operators;
}
int System_Parent::getSolverRungeKuttaOrder( int dir = DIR_T ) {
    return ( dir == DIR_T ? parameters.numerics_order_t : parameters.numerics_order_tau );
}
int System_Parent::getTimeTransformationAlg() {
    return parameters.numerics_order_timetrafo;
}

// Returns the number of iterations to skip in either t, tau (both the same such that a same sided lattice is given)
int System_Parent::getIterationSkip( ) {
    return parameters.iterations_t_skip;
}