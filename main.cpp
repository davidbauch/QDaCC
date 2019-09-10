#include "global.h"
#include "system/Exziton2NS/System_Exziton2NS.cpp"
#include "chirp.cpp"
#include "pulse.cpp"
#include "solver.cpp"

// g++-8 main2.cpp -o ptestCPP_mymatrix3.out -std=c++1y -O3 -DFMT_HEADER_ONLY -fopenmp

//current: test pMac15062019A

//BIGGE TODO:S:
// header files für system //DONE
// header files / .cpp files aufteilung vernünftig. //DONE
// generell abläufe in konstruktoren auf funktionsaufrufe umändern. more or less //DONE
// BBIGGGG: scaling der parameter!!
// Zugriffe auf parameters auf system.getXXX() abändern! //FIXME
// Abläufe wie t direction, spectrum über systemaufrufe. //DONE

// TODO: ???????????
// alles auf solver bringen: -> spectrum klasse weg! //DONE
// solver.calculate_t_direction() -> vector mit t-direction zeiten und densitymatritzen //DONE
// solver.calculate_tau_direction() //DONE bis auf speichern der matrizen
// solver.calculate_spectrum() -> inputs: t-vector und rho-vector, operatoren b,b+ //DONE
// solver.calculate_g2() -> inputs: t-vector, rho vector, b,b+
// ...

//TODO: output full densitymatrix on demand (-fullDM)
//TODO: alle abgefragten variablen auf funktionsaufrufe system.f() reduzieren
//TODO: array für pulse, wie chirp! (mehrere pulse zulassen, bei chirp sinds mehrere interpolant punkte)
//TODO: inputs für chirp und pulse mit eigener inpuit klasse!
//TODO: spereater skip für tau richtung und dann nochmal spectrum!
//TODO: time evolution des spektrums

// last 2 inputs: XY x=loglevem,y=outputhanderlstrings, Z z=workpath
int main( int argc, char* argv[] ) {
    // Help
    auto inputs = argv_to_vec( argc, argv );
    if ( argc < 3 || vec_find_str( "--help", inputs ) != -1 || vec_find_str( "-help", inputs ) != -1 || vec_find_str( "-h", inputs ) != -1 ) {
        if ( argc < 3 )
            fmt::print( "Not enough input parameters!\n" );
        Parameters::help();
        exit( 0 );
    }
    logs = Log( std::string( argv[argc - 1] ) + "logfile.log", vec_find_str( "-advLog", inputs ) != -1 );

    // System
    System system = System( inputs );
    // Solver
    ODESolver solver = ODESolver( system );
    // Normal Time direction
    solver.calculate_t_direction( system );
    // Spectrum
    if ( system.calculate_spectrum() ) {                                                                                 // obwohl is eig egal. das kann auch hier in der main geändert werden lol. hauptsache rest ist ez to get.
        solver.calculate_g1( system, system.operatorMatrices.photon_create, system.operatorMatrices.photon_annihilate ); //TODO: operatoren aus system.getNextSpectrumOperator() und dann das als vecator
        solver.calculate_spectrum( system );
    }
    if ( system.calculate_g2() ) {

    }
    // Finalizing all calculations
    system.exit_system();

    double finalTime = Timer::summary();
    logs( "\nStartcommand: " );
    for ( int ii = 0; ii < argc; ii++ )
        logs( "{} ", std::string( argv[ii] ) );
    logs( "\n\n" + system.terminate_message + "\n" );

    logs.close();
    if ( system.output_handlerstrings() ) {
        fmt::print( "\n{0} {1:.1f}\n", PREFIX_PERCENT_TIME_FINAL, finalTime );
        fmt::print( "{0} Done in {1}\n", PREFIX_SUFFIX, Timer::format( finalTime ) );
    }
    exit( EXIT_SUCCESS );
}