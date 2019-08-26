#pragma once
// Dependencies
#include "../misc/helperfunctions.h"
#include "../misc/log.h"
#include "../misc/timer.h"

class Parameters {
    public:
        // Constructor
        Parameters() {};
        Parameters(const std::vector<std::string> &arguments);
        // Log function; Uses log subclass (log.h)
        void log();
        // Help function, output when --help is called
        static void help();
    private:
        // Parsing inputs from string input vector. Uses parseinput functions
        bool parseInput(const std::vector<std::string> &arguments);
        // Adjusting inputs
        bool adjustInput();
};