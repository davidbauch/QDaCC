#include "misc/timer.h"

// Constructor
Timer::Timer() : Timer::Timer( "Generic Timer", true, true ){};
Timer::Timer( const std::string &_name, bool _addtoTotalStatistic, bool _printToSummary ) {
    wallTimeStarted = omp_get_wtime();
    wallTimeEnded = wallTimeStarted;
    totalWallTime = 0;
    cpuTimeStarted = clock();
    cpuTimeEnded = cpuTimeStarted;
    totalCPUTime = 0;
    totalIterationNum = 0;
    running = true;
    outputMod = 0.05; // in seconds
    outputLast = 0;
    name = _name;
    addtoTotalStatistic = _addtoTotalStatistic;
    printToSummary = _printToSummary;
}
// Iterate Timer, can be used instead of using start/stop if the range of the timer includes the whole mainloop.
Timer &Timer::iterate( int num ) {
    totalIterationNum += num;
    return *this;
}
// Start Timer
Timer &Timer::start() {
    // Log::L2("Starting timer '{}'... ",name);
    wallTimeStarted = omp_get_wtime();
    cpuTimeStarted = clock();
    running = true;
    return *this;
    // Log::L2("Done!\n");
}
// End Timer
Timer &Timer::end() {
    // Log::L2("Ending timer '{}'... ",name);
    wallTimeEnded = omp_get_wtime();
    cpuTimeEnded = clock();
    running = false;
    totalCPUTime = getCPUTimeOnce();
    totalWallTime = getWallTimeOnce();
    return iterate();
    // Log::L2("Done!\n");
}
// Add time
Timer &Timer::add( time_t cpu, double wall ) {
    totalWallTime += wall;
    totalCPUTime += cpu;
    return iterate();
}
// Elapsed Wall Time between all start,end and between last start/end
double Timer::getWallTime( double scale ) {
    return totalWallTime * scale;
}
// Elapsed Wall Time between start and now
double Timer::getWallTimeOnce( double scale ) {
    if ( running )
        return ( omp_get_wtime() - wallTimeStarted ) * scale;
    return ( wallTimeEnded - wallTimeStarted ) * scale;
}
// Elapsed CPU Time between all start,end and between last start/end
double Timer::getCPUTime( double scale ) {
    return totalCPUTime * scale;
}
double Timer::getCPUTimeOnce( double scale ) {
    if ( running )
        return (double)( clock() - cpuTimeStarted ) / CLOCKS_PER_SEC * scale;
    return (double)( cpuTimeEnded - cpuTimeStarted ) / CLOCKS_PER_SEC * scale;
}
// Average Iteration Number (walltime)
double Timer::getAverageIterationTime( double scale ) {
    if ( running )
        return getWallTimeOnce() / totalIterationNum * scale;
    return totalWallTime / totalIterationNum * scale;
}
// Total Iterations counted
double Timer::getTotalIterationNumber() {
    return totalIterationNum;
}

// Output?
bool Timer::doOutput() {
    if ( getWallTimeOnce() > outputLast + outputMod ) {
        outputLast = getWallTimeOnce();
        return true;
    }
    return false;
}

void Timer::setOutputModTime( double s ) {
    if ( s < 0 ) {
        return;
    }
    outputMod = s;
}

std::string Timer::getName() {
    return name;
}

std::string Timer::format( double in ) {
    int h = std::floor( in / 3600. );
    int m = std::floor( ( in - 3600. * h ) / 60 );
    int s = std::floor( ( in - 3600 * h - 60 * m ) );
    double ms = ( in - std::floor( in ) ) * 1000;
    return fmt::format( "{}{}{}{}{}{}{:.3f}{}", h > 0 ? std::to_string( h ) : "", h > 0 ? "h:" : "", m > 0 ? std::to_string( m ) : "", m > 0 ? "m:" : "", s > 0 ? std::to_string( s ) : "", s > 0 ? "s:" : "", ms, "ms" );
}

// returns total time taken
double Timers::Isummary( bool output ) {
    int len = 0;
    double totalWallTime = 0;
    double totalCPUTime = 0;
    for ( Timer timer : timers ) {
        if ( (int)timer.getName().size() > len )
            len = timer.getName().size();
    }
    if ( output ) {
        Log::Logger::inBar( "Timer" );
        for ( Timer timer : timers ) {
            if ( timer.addtoTotalStatistic ) {
                totalWallTime += timer.getWallTime();
                totalCPUTime += timer.getCPUTime();
            }
            if ( timer.printToSummary ) {
                Log::L1( "{:<{}}: Walltime: {} ", timer.getName(), len, Timer::format( timer.getWallTime() ) );
                if ( timer.getTotalIterationNumber() > 1 )
                    Log::L1( "CPUTime: {} Iterations: {} Average Time per Iteration: {}", ( timer.getCPUTime() != 0 ) ? Timer::format( timer.getCPUTime() ) : "--", timer.getTotalIterationNumber(), Timer::format( timer.getAverageIterationTime() ) );
                Log::L1( "\n" );
            }
        }
        Log::Logger::Bar( 75, Log::Logger::LEVEL_1, Log::Logger::BAR_1 );
        Log::L1( "{:<{}}: Walltime: {} CPUTime: {}\n", "Total", len, Timer::format( totalWallTime ), ( totalCPUTime != 0 ? Timer::format( totalCPUTime ) : "--" ) );
        Log::Logger::Bar();
    }
    return totalWallTime;
}

void Timers::IoutputProgress( Timer &t, ProgressBar &p, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix, int state ) {
    if ( output_handler )
        outputTimeStrings( t, currentIt, maxItTotal, suffix, state );
    else
        outputProgressBar( t, p, currentIt, maxItTotal, suffix, state );
}
void Timers::IoutputProgressBar( Timer &t, ProgressBar &p, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix, int state ) {
    if ( t.doOutput() || state == PROGRESS_FORCE_OUTPUT ) {
        if ( state == WAITING ) {
            p.wait( "", suffix );
        } else if ( state != PROGRESS_FORCE_OUTPUT )
            p.print( currentIt, maxItTotal, fmt::format( "T - {}", Timer::format( ( maxItTotal - currentIt ) * t.getAverageIterationTime() ) ), suffix );
        else {
            p.print( maxItTotal, maxItTotal, fmt::format( "T: {}", Timer::format( t.getWallTime() ) ), suffix );
            fmt::print( "\n" );
        }
    }
}
void Timers::IoutputTimeStrings( Timer &t, const unsigned int currentIt, const unsigned int maxItTotal, const std::string &suffix, bool final ) {
    if ( t.doOutput() || final ) {
        fmt::print( "{0}\t{1:.2f}\n", QDLC::Message::Prefix::PERCENT, ( currentIt / (double)maxItTotal * 100. ) );
        fmt::print( "{0}\t{1:.0f}\n", QDLC::Message::Prefix::PERCENT_TIME, 1.0 * maxItTotal / currentIt * t.getWallTimeOnce() );
        fmt::print( "{0}\t{1}\n", QDLC::Message::Prefix::SUFFIX, suffix );
    }
}