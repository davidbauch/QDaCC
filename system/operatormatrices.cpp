#include "operatormatrices.h"

void OperatorMatrices_Parent::init( const Parameters &p ) {
    Timer &timer_operatormatrices = createTimer( "Operator Matrices", true, false );
    timer_operatormatrices.start();
    logs.level2( "Generating operator matrices... " );
    if ( !generateOperators( p ) ) {
        logs.level2( "Generating operator matrices failed! Exitting program...\n" );
        logs.close();
        exit( EXIT_FAILURE );
    }
    timer_operatormatrices.end();
    logs.level2( "successful. Elapsed time is {}ms\n", timer_operatormatrices.getWallTime( TIMER_MILLISECONDS ) );
}