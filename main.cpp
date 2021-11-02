#include "global.h"
#include "chirp.h"
#include "pulse.h"
#include "system/system.h"
#include "solver/solver_ode.h"

int main( int argc, char* argv[] ) {
    // Commandline Arguments:
    QDLC::CommandlineArguments::init( argc, argv );
    // Check for Multifile, if true parse all settings
    std::vector<std::vector<std::string>> sets;
    std::string filename = QDLC::CommandlineArguments::get_parameter( "--file" );
    auto inputs = QDLC::String::argv_to_vec( argc, argv );
    if ( filename.compare( "none" ) != 0 ) {
        // Multifile mode: Program will read inputfile from --file and execute all lines that dont start with '#'
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
        std::string outputname = QDLC::String::splitline( line ).at( 1 );
        std::filesystem::create_directories( inputs.back() + outputname );
        std::ofstream fileout( inputs.back() + outputname + "/settings_" + outputname + ".txt", std::ofstream::out );
        std::vector<std::string> set;
        int counter = 0;
        while ( std::getline( file, line ) ) {
            fileout << line << std::endl;
            if ( line.size() > 5 && line.at( 0 ) != '#' && line.at( 0 ) != ' ' && line.at( 0 ) != '\t' ) {
                set = QDLC::String::splitline( line );
                for ( auto i = set.begin(); i != set.end(); i++ )
                    if ( ( *i ).compare( "#" ) == 0 ) {
                        set = std::vector<std::string>( set.begin(), i );
                        break;
                    }
            }
            if ( set.size() > 0 ) {
                if ( set.at( 0 ).compare( "#" ) != 0 && QDLC::String::vec_find_str( "python3", set ) == -1 ) {
                    if ( set.size() > 1 ) {
                        // Append global parameters. Second calls wont overwrite these, so global parameters ALWAYS overwrite local parameters
                        set.insert( set.begin() + 1, globalparams.begin(), globalparams.end() );
                        // Append final path info
                        set.emplace_back( inputs.back() + outputname + "/" + outputname + "_" + std::to_string( counter++ ) + "/" );
                        sets.emplace_back( set );
                    }
                }
                set.clear();
            }
        }
        fileout.close();
    } else {
        // Single file mode: Program will only execute passed parameterset
        sets.emplace_back( inputs );
    }

    // Main Program
    for ( auto& set : sets ) {
        QDLC::CommandlineArguments::init( set );
        inputs = set;
        const std::string fp = inputs.back();
        std::filesystem::create_directories( fp );
        // Logfile
        int loglevel = ( QDLC::String::vec_find_str( "-advLog", inputs ) != -1 || QDLC::String::vec_find_str( "-L2", inputs ) != -1 ? 2 : ( QDLC::String::vec_find_str( "-L3", inputs ) != -1 ? 3 : 1 ) );
        Log::init( std::string( inputs.back() ) + "logfile.log", loglevel );

        // System
        System system = System( inputs );
        // Solver
        auto solver = QDLC::Numerics::ODESolver( system );
        // Normal Time direction
        solver.calculate_t_direction( system );
        // G1 and G2 statistics
        solver.calculate_advanced_photon_statistics( system );

        //// Spectrum
        //if ( system.parameters.numerics_calculate_spectrum_H ) {
        //    solver.calculate_spectrum( system, system.operatorMatrices.photon_create_H, system.operatorMatrices.photon_annihilate_H, "spectrum_H", 1 );
        //    if ( system.parameters.numerics_output_electronic_emission ) {
        //        solver.calculate_spectrum( system, system.operatorMatrices.atom_sigmaplus_G_H, system.operatorMatrices.atom_sigmaminus_G_H, "electronic_spectrum_H", 0 );
        //    }
        //}
        //if ( system.parameters.numerics_calculate_spectrum_V ) {
        //    solver.calculate_spectrum( system, system.operatorMatrices.photon_create_V, system.operatorMatrices.photon_annihilate_V, "spectrum_V", 2 );
        //    if ( system.parameters.numerics_output_electronic_emission ) {
        //        solver.calculate_spectrum( system, system.operatorMatrices.atom_sigmaplus_G_V, system.operatorMatrices.atom_sigmaminus_G_V, "electronic_spectrum_V", 0 );
        //    }
        //}
        //if ( system.parameters.numerics_calculate_g2 ) {
        //    solver.calculate_advanced_photon_statistics( system, system.operatorMatrices.photon_create_H, system.operatorMatrices.photon_annihilate_H, system.operatorMatrices.photon_create_V, system.operatorMatrices.photon_annihilate_V, "advanced_photon_statistics" );
        //}

        // Finalizing all calculations
        system.exit_system();

        double finalTime = Timers::summary();
        Log::L1( "\nStartcommand: " );
        for ( auto& ii : inputs )
            Log::L1( "{} ", ii );
        Log::L1( "\n\n" + system.terminate_message + "\n" );

        Log::close();
        if ( system.parameters.output_handlerstrings ) {
            fmt::print( "\n{0} {1:.1f}\n", QDLC::Message::Prefix::PERCENT_TIME_FINAL, finalTime );
            fmt::print( "{0} Done in {1}\n", QDLC::Message::Prefix::SUFFIX, Timer::format( finalTime ) );
        }
        Timers::reset();
    }
    exit( EXIT_SUCCESS );
}