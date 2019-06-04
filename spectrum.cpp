#include "spectrum.h"

Spectrum::Spectrum(System &s) {
    int dim = (int)(s.parameters.maxItMainRK/s.parameters.akf_everyXIt) + 10;
    logs.level2("Creating spectrum with dimension: {}... ",dim);
    out.reserve((int)(std::ceil(s.parameters.akf_spec_max_w)+5));
    akf_mat = MatrixXcd::Zero(dim,dim);
    rhos.reserve(dim);
    times.reserve(dim);
    for(int w=0;w<s.parameters.akf_spec_max_w;w++) {
        out.push_back(0);
    }
    logs.level2("done!\n");
}

void Spectrum::addRho(MatrixXcd &rho, double t) {
    rhos.emplace_back(rho);
    times.emplace_back(t);
}

void Spectrum::calculateTauDirection(System &s, MatrixXcd &op_creator, MatrixXcd &op_annihilator, ODESolver &solver){
    Timer &timer = createTimer("RungeKutta-Tau-Loop");
    int totalIterations = getIterationNumberTau(s);
    ProgressBar progressbar = ProgressBar(totalIterations, 60, 0, BAR_HORIZONTAL, true, 0.1, {"|","/","-","\\"}, {".","-","="}, "Done"); 
    timer.start();

    #pragma omp parallel for schedule(dynamic) shared(timer) num_threads(s.parameters.akf_maxThreads) //reduction(+:spectrum.out)
    for (int i = 0; i < getRhoDim(); i++) {
        double t_t = times.at(i);
        //logs.level2("Calculating Tau-contribution at t = {}\n",t_t);
        MatrixXcd rho_tau = s.dgl_calc_rhotau(rhos.at(i),op_annihilator,t_t);
        akf_mat(i,0) = s.dgl_expectationvalue(rho_tau, op_creator, t_t);
        int j = 1; int curIt_tau = 1;
        for (double t_tau=t_t+s.parameters.t_step; t_tau<s.parameters.t_end; t_tau+=s.parameters.t_step) { // t + +s.parameters.t_step
            rho_tau = solver.iterate(rho_tau,s,t_tau, DIR_TAU);
            timer.iterate();
            if (curIt_tau%s.parameters.akf_everyXIt == 0) {
                akf_mat(i,j) = s.dgl_expectationvalue(rho_tau,op_creator,t_tau);
                j++; // equivalent to s.parameters.akf_vecIndex
                curIt_tau = 1;
            } else {
                curIt_tau++;
            }
            outputProgress(s.parameters.outputHandlerStrings,timer,progressbar,totalIterations,"AKF-Tau");
        }
    }
    timer.end();
} 

// Calculates Spectrum from AKF
void Spectrum::calculateSpectrum(System &s){
    Timer &timer = createTimer("Spectrum-Loop");
    int totalIterations = getIterationNumberSpectrum(s);
    ProgressBar progressbar = ProgressBar(totalIterations, 60, 0, BAR_HORIZONTAL, true, 0.1, {"|","/","-","\\"}, {".","-","="}, "Done"); 

    timer.start();
    #pragma omp parallel for schedule(dynamic) shared(timer) num_threads(s.parameters.akf_maxThreads) //reduction(+:spectrum.out)
    for (int spec_w=0; spec_w < s.parameters.akf_spec_max_w; spec_w++){
        std::vector<std::complex<double>> expfunc; 
        expfunc.reserve(getRhoDim());
        for(int spec_tau = 0; spec_tau < (s.parameters.t_end-s.parameters.t_start-0.0)/(s.parameters.t_step*s.parameters.akf_everyXIt); spec_tau++){
            expfunc.emplace_back( std::exp(-1i*s.parameters.akf_spec_possiblew.at(spec_w)*(double)(spec_tau)*s.parameters.t_step*(double)(s.parameters.akf_everyXIt)) );
        }
        for (int i = 0; i < getRhoDim(); i++) {
            double t_t = times.at(i);
            for(int spec_tau = 0; spec_tau < (s.parameters.t_end-s.parameters.t_start-t_t)/(s.parameters.t_step*s.parameters.akf_everyXIt); spec_tau++){
                out.at(spec_w) += expfunc.at(spec_tau)*akf_mat(i,spec_tau);
            }
            outputProgress(s.parameters.outputHandlerStrings,timer,progressbar,totalIterations,"AKF-Spectrum");
        }
        timer.iterate();
    }
    outputProgress(s.parameters.outputHandlerStrings,timer,progressbar,totalIterations,"AKF-Spectrum",PROGRESS_FORCE_OUTPUT);
    timer.end();
}

int Spectrum::getIterationNumberTau(System &s) {
    int num = 0;
    // Tau Direction Iteration steps
    for (int i = 0; i < getRhoDim(); i++) {
        double t_t = times.at(i);
        for (double t_tau=t_t+s.parameters.t_step; t_tau<s.parameters.t_end; t_tau+=s.parameters.t_step) { // t + +s.parameters.t_step
            num++;
        }
    }
    return num;
}
int Spectrum::getIterationNumberSpectrum(System &s) {
    int num = 0;
    // Spectrum steps
    for (int spec_w=0; spec_w < s.parameters.akf_spec_max_w; spec_w++){
        //for (int i = 0; i < getRhoDim(); i++) {
            num++;
        //}
    }
    return num;
}

void Spectrum::fileOutput(System &s, std::string filepath) {
    FILE *spectrumfile = std::fopen(filepath.c_str(), "w");
    if (!spectrumfile) {
        logs.level2("Failed to open outputfile for spectrum!\n");
        return;
    }
    for (int spec_w=0;spec_w<s.parameters.akf_spec_max_w;spec_w++){
        fmt::print(spectrumfile,"{0:.15e}\t{1:.15e}\n",s.parameters.akf_spec_possiblew[spec_w],real(out[spec_w]*s.parameters.t_step*s.parameters.t_step*(double)(s.parameters.akf_everyXIt*s.parameters.akf_everyXIt))); 
    }
    std::fclose(spectrumfile);  
}

int Spectrum::getRhoDim() {
    return rhos.size();
}

/*
timer.start();
#pragma omp parallel for schedule(dynamic) shared(timer) num_threads(s.parameters.akf_maxThreads) //reduction(+:spectrum.out)
for (int spec_w=0; spec_w < s.parameters.akf_spec_max_w; spec_w++){
    for (int i = 0; i < getRhoDim(); i++) {
        double t_t = times.at(i);
        for(int spec_tau = 0; spec_tau < (s.parameters.t_end-s.parameters.t_start-t_t)/(s.parameters.t_step*s.parameters.akf_everyXIt); spec_tau++){
            out.at(spec_w) += std::exp(-1i*s.parameters.akf_spec_possiblew.at(spec_w)*(double)(spec_tau)*s.parameters.t_step*(double)(s.parameters.akf_everyXIt))*akf_mat(i,spec_tau);
        }
        timer.iterate();
        outputProgress(s.parameters.outputHandlerStrings,timer,progressbar,totalIterations,"AKF-Spectrum");
    }
}
*/