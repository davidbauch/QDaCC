#pragma once
#include "pulse.h"

Pulse::Pulse(Parameters &p) {
    int n = (int)( (p.t_end - p.t_start)/p.t_step*6.0 + 5 ); 
    pulsearray.reserve(n);
    timearray.reserve(n);
    step = p.t_step;
    if (p.orderRKT == 5 || p.orderRKTau == 5)
        steps = { 0, 1./5.*p.t_step, 3./10.*p.t_step, 1./2.*p.t_step, 4./5.*p.t_step, 8./9.*p.t_step };
    else 
        steps = {0,0.5*p.t_step};
}

/* Generate array of energy-values corresponding to the Pulse */
void Pulse::generate(Parameters &p) {
    logs.level2("generating type " + inputs.pulsetype + "... ");
    
    if (inputs.pulsetype.compare("gauss_pi") == 0) {
        inputs.pulsetype = "gauss";
    }

    double t;
    for (double t1 = p.t_start; t1 < p.t_end+p.t_step*steps.size(); t1 += step) {
        for (int i = 0; i < steps.size(); i++ ) {
            t = t1 + steps[i];
            std::complex<double> val = 0;   
            if (inputs.pulsetype.compare("cw") == 0 && t>=p.pulse_center)
                val = std::exp(-1i*p.pulse_freq*(t-p.pulse_center));
            else if (inputs.pulsetype.compare("gauss") == 0)
                val = std::exp(-0.5*std::pow((t-p.pulse_center)/p.pulse_sigma,2.) - 1i*p.pulse_freq*(t-p.pulse_center));
            if (std::abs(inputs.pulse_amp*val) > 1E-10)
                pulsearray.push_back(inputs.pulse_amp*val);
            else 
                pulsearray.push_back(0);
            timearray.push_back(t);
        }
    }
 
    pulsearray.shrink_to_fit();
    timearray.shrink_to_fit();
    size = pulsearray.size();
    logs.level2("pulsearray.size() = {}... ", size);
}

/* Create single Pulse from parameters */
void Pulse::generateFromParam(Parameters &p) {
    logs.level2("Creating Pulse from input parameters... ");
    inputs = Inputs();
    inputs.pulse_center = p.pulse_center;
    inputs.pulse_amp    = p.pulse_amp;
    inputs.pulse_omega  = p.pulse_freq;
    inputs.pulse_sigma  = p.pulse_sigma;
    inputs.pulsetype    = p.pulsetype;
    generate(p);
    logs.level2("done!\n");
}

void Pulse::fileOutput(std::string filepath, Parameters &p) {
    FILE *pulsefile = std::fopen(filepath.c_str(),"w");
    if (!pulsefile) {
        logs.level2("Failed to open outputfile for Pulse!\n");
        return;
    }
    for (long unsigned int i = 0; i < timearray.size()-steps.size(); i+=steps.size()) {
    //for (long unsigned int i = 0; i < timearray.size(); i++) {
        std::fprintf(pulsefile,"%.10e\t%.10e\n",timearray.at(i),std::abs(pulsearray.at(i)));
    }
    std::fclose(pulsefile);
}

std::complex<double> Pulse::get(double t) const {
    int i = std::floor(t/step-1)*steps.size();
    while (timearray.at(i) < t) {
        i++;
    }
    if (i < 0 || i >= size) {
        logs.level2("!! Warning: requested pulsevalue at index {} is out of range! pulsearray.size() = {}\n",i,pulsearray.size());
        i = 0;
    }
    //logs.level2("Requested t = {:.10e}, returned t = {:.10e}\n",t,timearray.at(i));
    return pulsearray.at(i);
    //return lerp<std::complex<double>>(pulsearray.at(j),pulsearray.at(j+1),delta);
}

std::complex<double> Pulse::get(int i) const {
    if (i < 0 || i >= size) {
        logs.level2("!! Warning: requested pulsevalue at index {} is out of range!\n",i);
        i = 0;
    }
    return pulsearray.at(i);
}