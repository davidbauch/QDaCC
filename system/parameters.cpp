#include "parameters.h"

void Parameters_Parent::init(const std::vector<std::string> &arguments) {
    // Parsing input:
    Timer &timer_parseInput = createTimer("Parsing parameters");
    logs.wrapInBar("Conversion of input variables",LOG_SIZE_FULL,LOG_LEVEL_2,LOG_BAR_0); logs.level2("\n");
    logs.level2("Parsing input variables... ");
    timer_parseInput.start();
    if (!parseInput(arguments)) {
        logs.level2("Parsing input variables failed! Exitting program...\n");
        logs.close();
        exit(EXIT_FAILURE);
    }
    timer_parseInput.end();
    logs.level2("successful. Elapsed time is {}ms\n",timer_parseInput.getWallTime(TIMER_MILLISECONDS));

    // Adjusting inputs:
    Timer &timer_adjustInput = createTimer("Adjusting parameters");
    logs.level2("Adjusting input variables... ");
    timer_adjustInput.start();
    if (!adjustInput()) {
        logs.level2("Adjusting input variables failed! Exitting program...\n");
        logs.close();
        exit(EXIT_FAILURE);
    }
    logs.level2("successful. Elapsed time is {}ms\n",timer_adjustInput.getWallTime(TIMER_MILLISECONDS));
    timer_adjustInput.end();
}