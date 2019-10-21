#include "parameters.h"

Parameters_Parent::Parameters_Parent() {
    numerics_calculate_spectrum = 0;
    numerics_calculate_g2 = 0;
    numerics_use_interactionpicture = 0;
    numerics_use_rwa = 0;
    numerics_order_timetrafo = TIMETRANSFORMATION_MATRIXEXPONENTIAL;
    numerics_order_t = 4;
    numerics_order_tau = 4;
    output_advanced_log = 0;
    output_handlerstrings = 0;
    output_operators = 0;
    iterations_skips_t = 1;
    iterations_tau_resolution = 1;
    iterations_w_resolution = 1;
    scale_parameters = false;
    scale_value = 1E12;
}

void Parameters_Parent::init( const std::vector<std::string> &arguments ) {
    // Parsing input:
    Timer &timer_parseInput = createTimer( "Parsing parameters" );
    logs.wrapInBar( "Conversion of input variables", LOG_SIZE_FULL, LOG_LEVEL_2, LOG_BAR_0 );
    logs.level2( "\n" );
    logs.level2( "Parsing input variables... " );
    timer_parseInput.start();
    if ( !parseInput( arguments ) ) {
        logs.level2( "Parsing input variables failed! Exitting program...\n" );
        logs.close();
        exit( EXIT_FAILURE );
    }
    timer_parseInput.end();
    logs.level2( "successful. Elapsed time is {}ms\n", timer_parseInput.getWallTime( TIMER_MILLISECONDS ) );

    // Scaling inputs:
    if ( scale_parameters ) {
        logs.level2("Rescaling parameters to {}... ",scale_value);
        scaleInputs( scale_value );
        logs.level2("Done!\n");
    }

    // Adjusting inputs:
    Timer &timer_adjustInput = createTimer( "Adjusting parameters" );
    logs.level2( "Adjusting input variables... " );
    timer_adjustInput.start();
    if ( !adjustInput() ) {
        logs.level2( "Adjusting input variables failed! Exitting program...\n" );
        logs.close();
        exit( EXIT_FAILURE );
    }
    logs.level2( "successful. Elapsed time is {}ms\n", timer_adjustInput.getWallTime( TIMER_MILLISECONDS ) );
    timer_adjustInput.end();
}