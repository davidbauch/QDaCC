// MAKRO: LOADSYSTEM(System)
#include "global.h"
#include "system/Exziton2NS/System_Exziton2NS.cpp"
#include "chirp.cpp"
#include "pulse.cpp"
#include "rungekutta.cpp"
#include "spectrum.cpp"


// g++-8 main2.cpp -o ptestCPP_mymatrix3.out -std=c++1y -O3 -DFMT_HEADER_ONLY -fopenmp

/*
    // TODO: // REMOVE: // FIXME: ENTE QUACK QUACK
*/

/* 
    current: test pMac05062019G
 */

// TODO: RWA funktioniert wieder nicht wth
// TODO: 11.06.2019: chirp in timetrafo! dafür beim erstellen des chirps integral ausrechnen (und derivative wenn schon dabei) und ggf mit ausgeben
// dgl_getHamilton soll hamilton abhängig von RWA oder nicht zurückgeben (dann auch kein mist mit der timetrafo)

// last 2 inputs: XY x=loglevem,y=outputhanderlstrings, Z z=workpath
int main(int argc, char* argv[]) {
    logs = Log(std::string(argv[argc-1]) + "logfile.log", std::string(argv[argc-2]).at(0) == '1' );

    // System
    System system = System(argc, argv);
    system.init();
    // Solver
    ODESolver solver = ODESolver(system,false);
    // Spectrum
    Spectrum spectrum;
    if (system.parameters.doSpectrum)
        spectrum = Spectrum(system);
        
    /* Main Loop */
    Timer &rkTimer = createTimer("RungeKutta-Main-Loop");
    ProgressBar progressbar = ProgressBar(system.parameters.maxItTotal, 60, 0, BAR_HORIZONTAL, true, 0.1, {"|","/","-","\\"}, {".","-","="}, "Done"); //{"⬗","⬙","⬖","⬘"}, {"˥","˦","˧","˨","˩","˨","˧","˦"}, {"○","◔","◑","◕","●","◕","◑","◔"}, {"▁","▂","▃","▄","▅","▆","▇","█","▉","▊","▋","▌","▍","▎","▏"}, {"░", "▒", "▓"}, {" ","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"}, {"░", "▒", "▓"},{"▏","▎","▍","▌","▋","▊","▉","█"}
    rkTimer.start();
    
    MatrixXcd rho = system.getRho0();
    if (system.parameters.doSpectrum)
        spectrum.addRho(rho,system.parameters.t_start);
    system.expectationValues(rho,system.parameters.t_start);

    int curIt = 1;
    /* Main Time Loop */
    for (double t_t=system.parameters.t_start+system.parameters.t_step; t_t<system.parameters.t_end; t_t+=system.parameters.t_step) {
        // Runge-Kutta iteration 
        rho = solver.iterate(rho,system,t_t);
        // Save Rho for tau-direction
        if (system.parameters.doSpectrum && curIt%system.parameters.akf_everyXIt == 0) {
            curIt = 1;
            spectrum.addRho(rho,t_t);
        } else {
            curIt++;
        }
        // Expectation Values
        system.expectationValues(rho,t_t);
        // Progress and time output
        rkTimer.iterate();
        outputProgress(system.parameters.outputHandlerStrings,rkTimer,progressbar,system.parameters.maxItTotal,"Normal");
        // Divergent bruder? 
        if (!system.traceValid(rho,t_t)) {
            t_t = system.parameters.t_end+system.parameters.t_step;
        }
    }
    rkTimer.end();

    // Spectrum
    if (system.parameters.doSpectrum) {
        spectrum.calculateTauDirection(system, system.operatorMatrices.photon_create, system.operatorMatrices.photon_annihilate,solver);
        spectrum.calculateSpectrum(system);
        spectrum.fileOutput(system, system.parameters.subfolder + "spectrum.txt");
    }
    
    double finalTime = Timer::summary();
    logs("\nStartcommand: ");
    for (int ii=0;ii<argc;ii++)
        logs("{} ",std::string(argv[ii])); 
    logs("\n\n"+system.terminateMessage+"\n");

    logs.close();
    if (system.parameters.outputHandlerStrings) {
        fmt::print("\n{0} {1:.1f}\n",PREFIX_PERCENT_TIME_FINAL,finalTime);
        fmt::print("{0} Done in {1}\n",PREFIX_SUFFIX, Timer::format(finalTime));
    }
    exit(1);
}