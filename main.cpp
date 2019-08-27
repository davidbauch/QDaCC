// MAKRO: LOADSYSTEM(System)
#include "global.h"
#include "system/Exziton2NS/System_Exziton2NS.cpp"
#include "chirp.cpp"
#include "pulse.cpp"
#include "rungekutta.cpp"
#include "spectrum.cpp"


// g++-8 main2.cpp -o ptestCPP_mymatrix3.out -std=c++1y -O3 -DFMT_HEADER_ONLY -fopenmp
 
//current: test pMac15062019A

//BIGGE TODO:S:
// header files für system
// header files / .cpp files aufteilung vernünftig.
// generell abläufe in konstruktoren auf funktionsaufrufe umändern.


// last 2 inputs: XY x=loglevem,y=outputhanderlstrings, Z z=workpath
int main(int argc, char* argv[]) {

    // Help
    std::vector<std::string> inputs = argv_to_vec(argc,argv);
    if (argc < 3 || vec_find_str("--help",inputs) != -1 || vec_find_str("-help",inputs) != -1 || vec_find_str("-h",inputs) != -1) {
        if (argc < 3)
            fmt::print("Not enough input parameters!\n");
        Parameters::help();
        exit(0);
    }
    logs = Log(std::string(argv[argc-1]) + "logfile.log", vec_find_str("-advLog",inputs) != -1 );

    // System
    System system = System(inputs);
    // Solver
    ODESolver solver = ODESolver(system,false);
    // Spectrum
    Spectrum spectrum;
    if (system.parameters.doSpectrum)
        spectrum = Spectrum(system);
        
    // Main Loop
    Timer &rkTimer = createTimer("RungeKutta-Main-Loop");
    ProgressBar progressbar = ProgressBar(system.parameters.maxItTotal, 60, 0, BAR_VERTICAL, true, 0.1, {" ","▏","▎","▍","▌","▋","▊","▉","█"}); //{"⬗","⬙","⬖","⬘"}, {"˥","˦","˧","˨","˩","˨","˧","˦"}, {"○","◔","◑","◕","●","◕","◑","◔"}, {"▁","▂","▃","▄","▅","▆","▇","█","▉","▊","▋","▌","▍","▎","▏"}, {"░", "▒", "▓"}, {" ","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"}, {"░", "▒", "▓"},{"▏","▎","▍","▌","▋","▊","▉","█"}
    rkTimer.start();
    
    MatrixXcd rho = system.getRho0();
    if (system.parameters.doSpectrum)
        spectrum.addRho(rho,system.parameters.t_start);
    system.expectationValues(rho,system.parameters.t_start);

    int curIt = 1;
    // Main Time Loop
    for (double t_t=system.parameters.t_start+system.parameters.t_step; t_t<system.parameters.t_end; t_t+=system.parameters.t_step) {
        // Runge-Kutta iteration 
        rho = solver.iterate(rho,system,t_t);
        // Save Rho for tau-direction
        // TODO: move to system function system.queueSpectrum() die dann den check übernimmt.
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
        outputProgress(system.parameters.outputHandlerStrings,rkTimer,progressbar,system.parameters.maxItTotal,"T-Direction: ");
        // Divergent bruder? 
        if (!system.traceValid(rho,t_t)) {
            t_t = system.parameters.t_end+system.parameters.t_step;
        }
    }
    rkTimer.end();

    // Spectrum //TODO: in system.calculateSpectrum, damit für 4ns nacher ggf. mehrere spektren berechnet werden können!
    if (system.parameters.doSpectrum) {
        spectrum.calculateTauDirection(system, system.operatorMatrices.photon_create, system.operatorMatrices.photon_annihilate,solver);
        spectrum.calculateSpectrum(system);
        spectrum.fileOutput(system, system.parameters.subfolder + "spectrum.txt");
    }
    system.exit_system();

    double finalTime = Timer::summary();
    logs("\nStartcommand: ");
    for (int ii=0;ii<argc;ii++)
        logs("{} ",std::string(argv[ii])); 
    logs("\n\n"+system.terminate_message+"\n");

    logs.close();
    if (system.parameters.outputHandlerStrings) {
        fmt::print("\n{0} {1:.1f}\n",PREFIX_PERCENT_TIME_FINAL,finalTime);
        fmt::print("{0} Done in {1}\n",PREFIX_SUFFIX, Timer::format(finalTime));
    }
    exit(1);
}