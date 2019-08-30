#pragma once
#include "global.h"

class Spectrum {
    private:
        std::vector<std::complex<double>> out;
        MatrixXcd akf_mat;
        std::vector<MatrixXcd> rhos;
        std::vector<double> times;
        int curIt;
    public:
        Spectrum(System &s);
        Spectrum() {};
        void addRho(MatrixXcd &rho, double t);
        void calculateTauDirection(System &s, MatrixXcd &op1, MatrixXcd &op2, ODESolver &solver);
        void calculateSpectrum(System &s);
        int getIterationNumberTau(System &s);
        int getIterationNumberSpectrum(System &s);
        void fileOutput(System &s, std::string filepath);
        int getRhoDim();
        bool queueNow(System &s);
};