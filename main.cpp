#include "global.h"
#include "system/Exciton4NS/System_Exciton4NS.cpp"
#include "chirp.cpp"
#include "pulse.cpp"
#include "solver.cpp"

// g++-8 '/Users/davidbauch/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/main.cpp' -o /Users/davidbauch/bin/QDLC-2LS-1.0.out misc/ALGLIB/MAC/*.o -std=c++1y -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs -Wall
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

    // Check for Multifile, if true parse all settings
    std::vector<std::vector<std::string>> sets;
    if ( vec_find_str( "--file", inputs ) != -1 ) {
        // Multifile mode: Program will read inputfile from --file and execute all lines that dont start with '#'
        auto params = Parse_Parameters( inputs, {"--file", "--Threads", "-advLog", "-noHandler"}, {1, 1, 1, 1}, "Global Parameters" );
        std::string filename = params.get( 0, "" );
        if ( filename.size() > 1 ) {
            // Catch global parameters
            std::vector<std::string> globalparams;
            for ( unsigned int i = 0; i < inputs.size() - 1; i++ ) {
                if ( inputs.at( i ).compare( "--file" ) != 0 && ( i > 0 && inputs.at( i - 1 ).compare( "--file" ) != 0 ) ) {
                    globalparams.emplace_back( inputs.at( i ) );
                }
            }
            // Parse file and save to sets vector
            std::ifstream file( filename );
            std::string line;
            std::getline( file, line );
            std::string outputname = splitline( line ).at( 1 );
            std::filesystem::create_directories( inputs.back() + outputname );
            std::ofstream fileout( inputs.back() + outputname + "/settings_" + outputname + ".txt", std::ofstream::out );
            std::vector<std::string> set;
            int counter = 0;
            while ( std::getline( file, line ) ) {
                fileout << line << std::endl;
                if ( line.size() > 5 && line.at( 0 ) != '#' && line.at( 0 ) != ' ' && line.at( 0 ) != '\t' ) {
                    set = splitline( line );
                    for ( auto i = set.begin(); i != set.end(); i++ )
                        if ( ( *i ).compare( "#" ) == 0 ) {
                            set = std::vector<std::string>( set.begin(), i );
                            break;
                        }
                }
                if ( set.size() > 0 ) {
                    if ( set.at( 0 ).compare( "#" ) != 0 && vec_find_str( "python3", set ) == -1 ) {
                        if ( set.size() > 1 ) {
                            // Append global parameters. Second calls wont overwrite these, so global parameters ALWAYS overwrite local parameters
                            set.insert( set.begin() + 1, globalparams.begin(), globalparams.end() );
                            // Append final path info
                            set.emplace_back( inputs.back() + outputname + "/" + outputname + "_" + toStr( counter++ ) + "/" );
                            sets.emplace_back( set );
                        }
                    }
                    set.clear();
                }
            }
            fileout.close();
        }
        //for ( auto a : sets ) {
        //    for ( auto b : a )
        //        fmt::print( "{} | ", b );
        //    fmt::print( "\n" );
        //}
        //exit(1);
    } else {
        // Single file mode: Program will only execute passed parameterset
        sets.emplace_back( inputs );
    }

    // Main Program
    for ( auto set : sets ) {
        inputs = set;
        const std::string fp = inputs.back();
        std::filesystem::create_directories( fp );
        // Logfile
        logs = Log( std::string( inputs.back() ) + "logfile.log", vec_find_str( "-advLog", inputs ) != -1 );

        // System
        System system = System( inputs );
        // Solver
        ODESolver solver = ODESolver( system );
        // Normal Time direction
        solver.calculate_t_direction( system );

        // Spectrum
        if ( system.calculate_spectrum() ) {
            if ( system.calculate_spectrum_H() ) {
                solver.calculate_spectrum( system, system.operatorMatrices.photon_create_H, system.operatorMatrices.photon_annihilate_H, "spectrum_H.txt", 1 );
            }
            if ( system.calculate_spectrum_V() ) {
                solver.calculate_spectrum( system, system.operatorMatrices.photon_create_V, system.operatorMatrices.photon_annihilate_V, "spectrum_V.txt", 2 );
            }
        }
        if ( system.calculate_g2() ) {
            solver.calculate_advanced_photon_statistics( system, system.operatorMatrices.photon_create_H, system.operatorMatrices.photon_annihilate_H, system.operatorMatrices.photon_create_V, system.operatorMatrices.photon_annihilate_V, "advanced_photon_statistics.txt" );
        }

        // Finalizing all calculations
        system.exit_system();

        double finalTime = Timer::summary();
        logs( "\nStartcommand: " );
        for ( auto ii : inputs )
            logs( "{} ", ii );
        logs( "\n\n" + system.terminate_message + "\n" );

        logs.close();
        if ( system.output_handlerstrings() ) {
            fmt::print( "\n{0} {1:.1f}\n", PREFIX_PERCENT_TIME_FINAL, finalTime );
            fmt::print( "{0} Done in {1}\n", PREFIX_SUFFIX, Timer::format( finalTime ) );
        }
        Timer::reset();
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