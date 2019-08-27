#pragma once
// Dependencies
#include <stdlib.h>
#include "../misc/helperfunctions.h"
#include "../misc/log.h"
#include "../misc/timer.h"

class Parameters_Parent {
    // Mandatory Variables
    public:
        std::string subfolder;
        int doSpectrum, doInteractionPicture, doRWA, orderTimeTrafo;
        int orderRKT, orderRKTau;
        int advancedLogging, outputHandlerStrings, outputOperators;
    // Functions
    public:
        // Constructor
        Parameters_Parent() {};
        Parameters_Parent(const std::vector<std::string> &arguments) {};
        // Log function; Uses log subclass (log.h)
        void log();
        // Help function, output when --help is called
        static void help();
    //private:
        // @overwrite: Parsing inputs from string input vector. Uses parseinput functions
        virtual bool parseInput(const std::vector<std::string> &arguments) {return false;};
        // @overwrite: Adjusting inputs
        virtual bool adjustInput() {return false;};
        // To be called from child constructor
        void init(const std::vector<std::string> &arguments);
};