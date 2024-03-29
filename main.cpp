#include "global.h"
#include "system/evaluable/chirp.h"
#include "system/evaluable/pulse.h"
#include "system/system.h"
#include "solver/solver_ode.h"

int main( int argc, char* argv[] ) {
    // Commandline Arguments:
    QDACC::CommandlineArguments::init( argc, argv );
    // Check for Multifile, if true parse all settings
    std::vector<std::vector<std::string>> sets;
    std::string filename = QDACC::CommandlineArguments::get_parameter( "--file" );
    auto inputs = QDACC::String::argv_to_vec( argc, argv );
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
        std::string outputname = QDACC::String::splitline( line ).at( 1 );
        std::filesystem::create_directories( inputs.back() + outputname );
        std::ofstream fileout( inputs.back() + outputname + "/settings_" + outputname + ".txt", std::ofstream::out );
        std::vector<std::string> set;
        int counter = 0;
        while ( std::getline( file, line ) ) {
            fileout << line << std::endl;
            if ( line.size() > 5 && line.at( 0 ) != '#' && line.at( 0 ) != ' ' && line.at( 0 ) != '\t' ) {
                set = QDACC::String::splitline( line );
                for ( auto i = set.begin(); i != set.end(); i++ )
                    if ( ( *i ).compare( "#" ) == 0 ) {
                        set = std::vector<std::string>( set.begin(), i );
                        break;
                    }
            }
            if ( set.size() > 0 ) {
                if ( set.at( 0 ).compare( "#" ) != 0 && QDACC::String::vec_find_str( "python3", set ) == -1 ) {
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
        QDACC::CommandlineArguments::init( set );
        inputs = set;
        const std::string fp = inputs.back();
        std::filesystem::create_directories( fp );
        // Logfile
        int loglevel = ( QDACC::String::vec_find_str( "-advLog", inputs ) != -1 || QDACC::String::vec_find_str( "-L2", inputs ) != -1 ? 2 : ( QDACC::String::vec_find_str( "-L3", inputs ) != -1 ? 3 : 1 ) );
        Log::Logger::init( std::string( inputs.back() ) + "logfile.log", loglevel );

        // System
        auto system = QDACC::System( inputs );
        // Solver
        auto solver = QDACC::Numerics::ODESolver( system );

        // TODO:
        // Do Parameter Optimization here and edit system parameters accordingly.

        // Path Integral Visual Output
        if ( system.parameters.output_dict.contains( "path" ) )
            solver.visualize_path( system.operatorMatrices.rho, system );

        // Normal Time direction
        solver.calculate_t_direction( system );

        // G1 and G2 statistics
        solver.calculate_advanced_photon_statistics( system );

        // Output Numerical stuff
        solver.output_numerical_data( system );

        // Finalizing all calculations
        system.exit_system();

        double finalTime = Timers::summary();
        Log::L1( "\nStartcommand: " );
        for ( auto& ii : inputs )
            Log::L1( "{} ", ii );
        Log::L1( "\n" );

        // TODO: remove all handlerstrings, because they are not needed anymore
        if ( system.parameters.output_handlerstrings ) {
            Log::L1( "\n{0} {1:.1f}\n", QDACC::Message::Prefix::PERCENT_TIME_FINAL, finalTime );
            Log::L1( "{0} Done in {1}\n", QDACC::Message::Prefix::SUFFIX, Timer::format( finalTime ) );
        }
        Timers::reset();
        Log::Logger::close();
    }
    exit( EXIT_SUCCESS );
}