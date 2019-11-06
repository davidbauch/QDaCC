#define SYSTEM4LS

#include "global.h"
#ifdef SYSTEM2LS
#include "system/Exciton2NS/System_Exciton2NS.cpp"
#endif
#ifdef SYSTEM4LS
#include "system/Exciton4NS/System_Exciton4NS.cpp"
#endif
#include "chirp.cpp"
#include "pulse.cpp"
#include "solver.cpp"

// g++-8 '/Users/davidbauch/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/main.cpp' -o /Users/davidbauch/bin/QDLC-2LS-1.0.out misc/ALGLIB/MAC/*.o -std=c++1y -O3 -DFMT_HEADER_ONLY -fopenmp -Wall
// test
// last input: workpath
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
#ifdef SYSTEM2LS
    if ( system.calculate_spectrum() ) {
        solver.calculate_g1( system, system.operatorMatrices.photon_create, system.operatorMatrices.photon_annihilate );
        solver.calculate_spectrum( system );
    }
    if ( system.calculate_g2() ) {
        solver.calculate_g2_0( system, system.operatorMatrices.photon_create, system.operatorMatrices.photon_annihilate );
    }
#endif
#ifdef SYSTEM4LS
    if ( system.calculate_spectrum() ) {
        if ( system.calculate_spectrum_H() ) {
            solver.calculate_g1( system, system.operatorMatrices.photon_create_H, system.operatorMatrices.photon_annihilate_H );
            solver.calculate_spectrum( system, "spectrum_H.txt" );
        }
        if ( system.calculate_spectrum_V() ) {
            solver.calculate_g1( system, system.operatorMatrices.photon_create_V, system.operatorMatrices.photon_annihilate_V );
            solver.calculate_spectrum( system, "spectrum_V.txt" );
        }
    }
    if ( system.calculate_g2() ) {
        solver.calculate_g2_0( system, system.operatorMatrices.photon_create_H, system.operatorMatrices.photon_annihilate_H, "g2(0)_H.txt" );
        solver.calculate_g2_0( system, system.operatorMatrices.photon_create_V, system.operatorMatrices.photon_annihilate_V, "g2(0)_V.txt" );
    }
#endif

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

/*
2NS:
// Normal Time direction
solver.calculate_t_direction( system );

// Spectrum
if ( system.calculate_spectrum() ) {
    solver.calculate_g1( system, system.operatorMatrices.photon_create, system.operatorMatrices.photon_annihilate );
    solver.calculate_spectrum( system );
}
if ( system.calculate_g2() ) {
    solver.calculate_g2_0( system, system.operatorMatrices.photon_create, system.operatorMatrices.photon_annihilate );
}
*/

/*
4NS
// Normal Time direction
    solver.calculate_t_direction( system );

    // Spectrum
    if ( system.calculate_spectrum() ) {
        solver.calculate_g1( system, system.operatorMatrices.photon_create_H, system.operatorMatrices.photon_annihilate_H );
        solver.calculate_spectrum( system );
    }
    if ( system.calculate_g2() ) {
        solver.calculate_g2_0( system, system.operatorMatrices.photon_create_H, system.operatorMatrices.photon_annihilate_H );
    }
*/