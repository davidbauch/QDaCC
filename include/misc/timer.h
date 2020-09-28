#pragma once
#include "global.h"

#define PROGRESS_FORCE_OUTPUT 1
#define TIMER_SECONDS 1.0
#define TIMER_MILLISECONDS 1E3
#define TIMER_MICROSECONDS 1E6

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
    bool addtoTotalStatistic, printToSummary;

   public:
    Timer();
    Timer( const std::string &_name, bool _addtoTotalStatistic, bool _printToSummary );
    void start();
    void end();
    void iterate( int num = 1 );
    void add( time_t cpu, double wall );
    double getWallTime( int scale = TIMER_SECONDS );
    double getWallTimeOnce( int scale = TIMER_SECONDS );
    double getCPUTime( int scale = TIMER_SECONDS );
    double getCPUTimeOnce( int scale = TIMER_SECONDS );
    double getAverageIterationTime( int scale = TIMER_SECONDS );
    double getTotalIterationNumber();
    bool doOutput();
    void setOutputModTime( double s );
    std::string getName();
    static double summary( bool output = true );
    static std::string format( double in );
    static void reset();
};

extern std::vector<Timer> allTimers;
void outputTimeStrings( Timer &t, const long unsigned int maxItTotal, std::string suffix = "", int final = 0 );
void outputProgressBar( Timer &t, ProgressBar &p, const long unsigned int maxItTotal, std::string suffix = "", int final = 0 );
void outputProgress( int handler, Timer &t, ProgressBar &p, const long unsigned int maxItTotal, std::string suffix = "", int final = 0 );
Timer &createTimer( std::string _name = "Generic timer", bool _addtoTotalStatistic = true, bool _printToSummary = true );