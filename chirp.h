#pragma once
#include "global.h"

class Chirp {
    private:
        class Inputs {
            public:
                std::vector<double> chirp_onofftime;
                std::vector<double> chirp_start;
                std::vector<double> chirp_end;
                std::vector<double> chirp_total;
                Inputs() {}
                Inputs(Parameters &p) {
                    chirp_total.push_back(p.chirp_total);
                    chirp_start.push_back(p.chirp_start);
                    chirp_end.push_back(p.chirp_end);
                    chirp_onofftime.push_back(p.chirp_onofftime);
                }
        };
        // Chirp values
        std::vector<double> chirparray;
        // Saving the corresponding t
        std::vector<double> timearray;
        // Helper variables/functions
        Inputs inputs;
        void generate(Parameters &p);
        int size;
        double step;
        std::vector<double> steps;
    public:
        Chirp() {};
        Chirp(Parameters &p);
        void generateFromInputFile(Parameters &p);
        void generateFromParam(Parameters &p);
        void generateFromSpline(Parameters &p);
        void fileOutput(std::string filepath, Parameters &p);
        double get(double t) const;
        double get(int i) const;
};