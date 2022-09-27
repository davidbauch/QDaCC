#pragma once
#include <omp.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include "misc/ProgressBar.h"
#include "misc/log.h"
#include "typedef.h"

/* Timer Functions */
class Timer {
   private:
    // CPU time
    time_t cpuTimeStarted;
    time_t cpuTimeEnded;
    time_t totalCPUTime;
    // Wall time
    double wallTimeStarted;
    double wallTimeEnded;
    double totalWallTime;
    double totalIterationNum;
    // Runtime
    bool running;
    double outputMod;
    double outputLast;
    std::string name;

   public:
    bool addtoTotalStatistic;
    bool printToSummary;
    Timer();
    Timer( const std::string &_name, bool _addtoTotalStatistic, bool _printToSummary );
    void start();
    void end();
    void iterate( int num = 1 );
    void add( time_t cpu, double wall );
    double getWallTime( double scale = 1.0 );
    double getWallTimeOnce( double scale = 1.0 );
    double getCPUTime( double scale = 1.0 );
    double getCPUTimeOnce( double scale = 1.0 );
    double getAverageIterationTime( double scale = 1.0 );
    double getTotalIterationNumber();
    bool doOutput();
    void setOutputModTime( double s );
    std::string getName();
    static std::string format( double in );
};

class Timers {
   public:
    static constexpr int PROGRESS_FORCE_OUTPUT = 1;
    static constexpr double SECONDS = 1.0;
    static constexpr double MILLISECONDS = 1E3;
    static constexpr double MICROSECONDS = 1E6;
    static constexpr int RUNNING = 2;
    static constexpr int WAITING = 3;

    Timers( Timers & ) = delete;
    static Timers &Get() {
        static Timers instance;
        return instance;
    }
    static Timer &get( const std::string &name ) {
        return Get().Iget( name );
    }
    static Timer &create( const std::string &name, bool addToTotalStatistics = true, bool printToSummary = true ) {
        return Get().Icreate( name, addToTotalStatistics, printToSummary );
    }
    static void outputTimeStrings( Timer &t, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix = "", bool final = false ) {
        return Get().IoutputTimeStrings( t, currentIt, maxItTotal, suffix, final );
    }
    static void outputProgressBar( Timer &t, ProgressBar &p, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix = "", int state = RUNNING ) {
        return Get().IoutputProgressBar( t, p, currentIt, maxItTotal, suffix, state );
    }
    static void outputProgress( Timer &t, ProgressBar &p, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix = "", int state = RUNNING ) {
        return Get().IoutputProgress( t, p, currentIt, maxItTotal, suffix, state );
    }
    static void reset() {
        return Get().Ireset();
    }
    static double summary( bool output = true ) {
        return Get().Isummary( output );
    }
    static bool toggleHandler() {
        return Get().ItoggleHandler();
    }

   private:
    std::vector<Timer> timers;
    bool output_handler = false;
    Timers(){};
    Timer &Iget( const std::string &name ) {
        for ( unsigned int i = 0; i < timers.size(); i++ ) {
            if ( timers.at( i ).getName().compare( name ) == 0 ) {
                return timers.at( i );
            }
        }
    }
    Timer &Icreate( const std::string &name, bool addToTotalStatistic, bool printToSummary ) {
        timers.push_back( { name, addToTotalStatistic, printToSummary } );
        Log::L2( "[Timer] Created timer with name '{}'{}.\n", name, ( addToTotalStatistic ? " which will be added to total statistics" : "" ) );
        return timers.back();
    }
    void IoutputTimeStrings( Timer &t, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix, bool final );
    void IoutputProgressBar( Timer &t, ProgressBar &p, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix, int state );
    void IoutputProgress( Timer &t, ProgressBar &p, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix, int state );
    void Ireset() {
        timers.clear();
    }
    bool ItoggleHandler() {
        output_handler = not output_handler;
        return output_handler;
    }
    double Isummary( bool output );
};