#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <unordered_map>
#include <format>
#include <fstream>
// TODO: fix vector overload
#include <sstream>
#include <algorithm>    

namespace Log {

const static int LEVEL_1 = 1;
const static int LEVEL_2 = 2;
const static int LEVEL_3 = 3;
const static int BAR_SIZE_FULL = 75;
const static int BAR_SIZE_HALF = 40;
const static int BAR_0 = 0;
const static int BAR_1 = 1;
const static int BAR_2 = 2;
const static int BAR_3 = 3;
const static int BAR_4 = 4;
const static int BAR_5 = 5;
const static int BAR_6 = 6;

// Escape Commands
const static std::string RED = "\033[31m";
const static std::string GREEN = "\033[32m";
const static std::string YELLOW = "\033[33m";
const static std::string BLUE = "\033[34m";
const static std::string BOLD = "\033[1m";
const static std::string RESET = "\033[0m";

class Logger {
   public:
    Logger( const Logger & ) = delete;

    static Logger &Get() {
        static Logger instance;
        return instance;
    }
    static void init( const std::string &filepath, int max_log_level ) {
        return Get().Iinit( filepath, max_log_level );
    }
    template <typename... Args>
    static void Warning( const std::string &msg, const Args &...args ) {
        return Get().Ilevel1_log( YELLOW+BOLD+msg+RESET, true, std::make_format_args( args... ) );
    }
    template <typename... Args>
    static void L1( const std::string &msg, const Args &...args ) {
        return Get().Ilevel1_log( msg, true, std::make_format_args( args... ) );
    }
    template <typename... Args>
    static void P1( const std::string &msg, const Args &...args ) {
        return Get().Ilevel1_log( msg, false, std::make_format_args( args... ) );
    }
    template <typename... Args>
    static void L2( const std::string &msg, const Args &...args ) {
        return Get().Ilevel2_log( msg, true, true, std::make_format_args( args... ) );
    }
    template <typename T>
    static void L2( const std::string &msg, const std::vector<T> &vec ) {
        std::stringstream kekw;
        std::ranges::for_each( vec.begin(), vec.end(), [&kekw]( const auto &el ) { kekw << el << ","; } );
        return Get().Ilevel2_log( msg, true, true, std::make_format_args( kekw.str() ) );
    }
    template <typename... Args>
    static void L3( const std::string &msg, const Args &...args ) {
        return Get().Ilevel3_log( msg, true, true, std::make_format_args( args... ) );
    }
    template <typename T>
    static void L3( const std::string &msg, const std::vector<T> &vec ) {
        std::stringstream kekw;
        std::ranges::for_each( vec.begin(), vec.end(), [&kekw]( const auto &el ) { kekw << el << ","; } );
        return Get().Ilevel3_log( msg, true, true, std::make_format_args( kekw.str() ) );
    }

    template <typename... Args>
    static void Error( const std::string &file, const std::string &function, int line, const std::string &msg, const Args &...args ) {
        return Get().Ierror_log( file, function, line, msg, std::make_format_args( args... ) );
    }

    static void Bar( int size = Log::BAR_SIZE_FULL, int level = Log::LEVEL_1, int _bar = Log::BAR_0 ) {
        return Get().Ibar( size, level, _bar );
    }
    static void inBar( const std::string &msg, int size = Log::BAR_SIZE_FULL, int level = Log::LEVEL_1, int _bar = Log::BAR_0 ) {
        return Get().Iinbar( msg, size, level, _bar );
    }
    static void wrapInBar( const std::string &msg, int size = Log::BAR_SIZE_FULL, int level = Log::LEVEL_1, int _barOut = Log::BAR_0, int _barIn = -1 ) {
        return Get().Iwrapinbar( msg, size, level, _barOut, _barIn );
    }
    static void Debug( const std::string &msg = "" ) {
        return Get().Idebug( msg );
    }
    static int max_log_level() {
        return Get().Iget_max_loglevel();
    }
    static void close() {
        return Get().Iclose();
    }
    static void Nothing(){
        // This is a placeholder method for disabled logger levels.
    };

   private:
    int max_loglevel;
    int debug_counter;
    std::unordered_map<std::string, std::string> colormap;
    std::fstream file;
    std::vector<std::string> bars;
    Logger(){};
    std::string repeat( std::string s, int n ) {
        std::string s1 = s;
        for ( int i = 1; i < n; i++ )
            s += s1;
        return s;
    }
    void Iinit( const std::string &filepath, int max_log_level ) {
        file = std::fstream();
        file.open( filepath, std::ios::out );
        if (not file.is_open()){
            std::cout << "Could not open file '" << filepath << "' for logging!\n";
            exit(1);
        }
        //std::setbuf( file, NULL );
        //std::setbuf( stdout, NULL );
        max_loglevel = max_log_level;
        debug_counter = 0;
        bars = { repeat( "=", Log::BAR_SIZE_FULL ),
                 repeat( "-", Log::BAR_SIZE_FULL ),
                 repeat( "o", Log::BAR_SIZE_FULL ),
                 repeat( "#", Log::BAR_SIZE_FULL ),
                 repeat( ":", Log::BAR_SIZE_FULL ),
                 repeat( "+", Log::BAR_SIZE_FULL ),
                 repeat( "*", Log::BAR_SIZE_FULL ) };
        if ( max_loglevel >= 2 ) {
            file << std::format( "Succesfully created logfile '{}'!\n", filepath );
        }
    }
    template <class Args>
    void Ilevel1_log( const std::string &msg, bool to_file, const Args &args ) {
        std::cout << std::vformat( msg, args );
        if ( to_file )
            file << std::vformat( msg, args );
    }
    template <class Args>
    void Ilevel2_log( const std::string &msg, bool use_colormap, bool to_file, const Args &args ) {
        if ( max_loglevel > LEVEL_1 ) {
            if ( use_colormap ) {
                auto p1 = msg.find_first_of( "[" ) + 1;
                auto p2 = msg.find_first_of( "]" );
                if ( p1 == std::string::npos or p2 == std::string::npos )
                    use_colormap = false;
                else {
                    std::string mid = msg.substr( p1, p2 - p1 );
                    if ( not colormap.contains( mid ) )
                        colormap[mid] = "\033[38;5;" + std::to_string( colormap.size() + 2 ) + "m";
                    auto color = colormap[mid];
                    std::cout << std::vformat( "| \033[38;2;100;100;100m" + msg.substr( 0, p1 ) + color + mid + "\033[38;2;100;100;100m" + msg.substr( p2 ) + "\033[0m", args );
                }
            }
            if ( not use_colormap )
                std::cout << std::vformat( "| \033[38;2;100;100;100m" + msg + "\033[0m", args );
            if ( to_file )
                file << std::vformat( " | " + msg, args );
        }
    }
    template <class Args>
    void Ilevel3_log( const std::string &msg, bool use_colormap, bool to_file, const Args &args ) {
        if ( max_loglevel > LEVEL_2 ) {
            if ( use_colormap ) {
                auto p1 = msg.find_first_of( "[" ) + 1;
                auto p2 = msg.find_first_of( "]" );
                if ( p1 == std::string::npos or p2 == std::string::npos )
                    use_colormap = false;
                else {
                    std::string mid = msg.substr( p1, p2 - p1 );
                    if ( not colormap.contains( mid ) )
                        colormap[mid] = "\033[38;5;" + std::to_string( colormap.size() + 2 ) + "m";
                    auto color = colormap[mid];
                    std::cout << std::vformat( " | \033[31m- \033[93m" + msg.substr( 0, p1 ) + color + mid + "\033[38;2;100;100;100m" + msg.substr( p2 ) + "\033[0m", args );
                }
            }
            if ( not use_colormap )
                std::cout << std::vformat( " | \033[31m- \033[93m" + msg + "\033[0m", args );
            if ( to_file )
                file << std::vformat( " | - " + msg, args );
        }
    }
    template <class Args>
    void Ierror_log( const std::string &in_file, const std::string &function, int line, const std::string &msg, const Args &args ) {
        std::string error_message = std::format( "On line {} in function {} in file {} -- ", line, function, in_file );
        std::cout << std::vformat( "\033[31m[ERROR]\033[93m " + error_message + msg + "\033[0m", args );
        file << std::vformat( "[ERROR] " + error_message + msg, args );
    }
    void Ibar( int size, int level, int _bar ) {
        size = std::min( std::max( size, 0 ), Log::BAR_SIZE_FULL );
        std::string ret = bars.at( _bar ).substr( 0, size );
        if ( level == LEVEL_1 ) {
            L1( "{}{}\n", ret, ret );
        } else if ( level == LEVEL_2 ) {
            L2( "{}{}\n", ret, ret );
        } else if ( level == LEVEL_3 ) {
#ifndef LOG_DISABLE_L3
            L3( "{}{}\n", ret, ret );
#endif
        }
    }
    void Iinbar( const std::string &msg, int size, int level, int _bar ) {
        std::string ret = " | " + msg + " | ";
        ret = bars.at( _bar ).substr( 0, std::floor( size - ret.size() / 2.0 ) ) + ret;
        ret += bars.at( _bar ).substr( 0, std::floor( 2 * size - ret.size() ) - 1 );
        if ( (int)ret.size() <= 2 * size )
            ret += bars.at( _bar ).substr( 0, 1 );
        if ( level == LEVEL_1 ) {
            L1( "{}\n", ret );
        } else if ( level == LEVEL_2 ) {
            L2( "{}\n", ret );
        } else if ( level == LEVEL_3 ) {
#ifndef LOG_DISABLE_L3
            L3( "{}\n", ret );
#endif
        }
    }
    void Iwrapinbar( const std::string &msg, int size, int level, int _barOut, int _barIn ) {
        _barIn = ( _barIn == -1 ? _barOut : _barIn );
        Bar( size, level, _barOut );
        inBar( msg, size, level, _barIn );
        Bar( size, level, _barOut );
    }
    void Idebug( const std::string &msg ) {
        L1( "[\033[31mDEBUG\033[0m] {}{}\n", std::format( "({}) ", debug_counter++ ), msg );
    }
    void Iclose() {
        L1( "[END OF LOGFILE]" );
        file.close();
    }
    int Iget_max_loglevel() {
        return max_loglevel;
    }
};

// TODO:
//  #define LOG_ERROR logError(__FILE__, __FUNCTION__, __LINE__)
//  inline void logError(const char* file, const char* func, int line) {cout them}

// Log Macros to disable L2 and especially L3 logging at compile time
//#define LOG_DISABLE_L1
//#define LOG_DISABLE_L2
//#define LOG_DISABLE_L3

#ifdef LOG_DISABLE_L1
#    define L1( fmt, ... ) Logger::Nothing()
#else
#    define L1( fmt, ... ) Logger::L1( fmt, ##__VA_ARGS__ )
#endif

#ifdef LOG_DISABLE_L1
#    define Warning( fmt, ... ) Logger::Nothing()
#else
#    define Warning( fmt, ... ) Logger::Warning( fmt, ##__VA_ARGS__ )
#endif

#ifdef LOG_DISABLE_L2
#    define L2( fmt, ... ) Logger::Nothing()
#else
#    define L2( fmt, ... ) Logger::L2( fmt, ##__VA_ARGS__ )
#endif

#ifdef LOG_DISABLE_L3
#    define L3( fmt, ... ) Logger::Nothing()
#else
#    define L3( fmt, ... ) Logger::L3( fmt, ##__VA_ARGS__ )
#endif

#define E( fmt, ... ) Logger::Error( __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__ )
#define Error( fmt, ... ) Logger::Error( __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__ )

} // namespace Log