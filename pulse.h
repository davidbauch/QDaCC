#pragma once
#include "global.h"

class Pulse {
    private:
        class Inputs {
            public:    
                std::string pulsetype;
                double pulse_center;
                double pulse_amp;
                double pulse_sigma;
                double pulse_omega;
                Inputs() {};
        };
        // Pulse values
        std::vector<std::complex<double>> pulsearray;
        // Saving the corresponding t doesn't hurt
        std::vector<double> timearray;
        // Helper variables/functions
        Inputs inputs;
        void generate(Parameters &p);
        int size;
        double step;
        std::vector<double> steps;
    public:
        Pulse() {};
        Pulse(Parameters &p);
        void generateFromParam(Parameters &p);
        void fileOutput(std::string filepath, Parameters &p);
        std::complex<double> get(double t) const;
        std::complex<double> get(int i) const;
};