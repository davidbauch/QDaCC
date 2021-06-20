#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fmt/core.h> // -DFMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>

class Log {
   public:
    const static int LEVEL_1 = 1;
    const static int LEVEL_2 = 2;
    const static int BAR_SIZE_FULL = 75;
    const static int BAR_SIZE_HALF = 40;
    const static int BAR_0 = 0;
    const static int BAR_1 = 1;
    const static int BAR_2 = 2;
    const static int BAR_3 = 3;
    const static int BAR_4 = 4;
    const static int BAR_5 = 5;
    const static int BAR_6 = 6;

   public:
    Log( const Log & ) = delete;

    static Log &Get() {
        static Log instance;
        return instance;
    }
    static void init( const std::string &filepath, int max_log_level ) {
        return Get().Iinit( filepath, max_log_level );
    }
    template <typename... Args>
    static void L1( const std::string &msg, const Args &...args ) {
        return Get().Ilevel1_log( msg, fmt::make_format_args( args... ) );
    }
    template <typename... Args>
    static void L2( const std::string &msg, const Args &...args ) {
        return Get().Ilevel2_log( msg, fmt::make_format_args( args... ) );
    }
    template <typename... Args>
    static void L3( const std::string &msg, const Args &...args ) {
        return Get().Ilevel3_log( msg, fmt::make_format_args( args... ) );
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

    static void debug( const std::string &msg, int counter = -1 ) {
        return Get().Idebug( msg, counter );
    }
    static void debug( int counter = -1 ) {
        return Get().Idebug( "", counter );
    }
    static void close() {
        return Get().Iclose();
    }

   private:
    int max_loglevel;
    int debug_counter;
    FILE *file;
    std::vector<std::string> bars;
    Log(){};
    std::string repeat( std::string s, int n ) {
        std::string s1 = s;
        for ( int i = 1; i < n; i++ )
            s += s1;
        return s;
    }
    void Iinit( const std::string &filepath, int max_log_level ) {
        file = std::fopen( filepath.c_str(), "w" );
        assert( file );
        std::setbuf( file, NULL );
        std::setbuf( stdout, NULL );
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
            fmt::print( file, "Succesfully created logfile '{}'!\n", filepath );
        }
    }
    template <class Args>
    void Ilevel1_log( const std::string &msg, const Args &args ) {
        fmt::vprint( msg, args );
        fmt::vprint( file, msg, args );
    }
    template <class Args>
    void Ilevel2_log( const std::string &msg, const Args &args ) {
        if ( max_loglevel > LEVEL_1 ) {
            fmt::vprint( " | " + msg, args );
            fmt::vprint( file, " | " + msg, args );
        }
    }
    template <class Args>
    void Ilevel3_log( const std::string &msg, const Args &args ) {
        if ( max_loglevel > LEVEL_2 ) {
            fmt::vprint( " | - " + msg, args );
            fmt::vprint( file, " | - " + msg, args );
        }
    }
    void Ibar( int size, int level, int _bar ) {
        size = std::min( std::max( size, 0 ), Log::BAR_SIZE_FULL );
        std::string ret = bars.at( _bar ).substr( 0, size );
        if ( level == LEVEL_1 ) {
            L1( "{}{}\n", ret, ret );
        } else if ( level == LEVEL_2 ) {
            L2( "{}{}\n", ret, ret );
        }
    }
    void Iinbar( const std::string &msg, int size, int level, int _bar ) {
        std::string ret = " | " + msg + " | ";
        ret = bars.at( _bar ).substr( 0, std::floor( size - ret.size() / 2.0 ) ) + ret + bars.at( _bar ).substr( 0, std::ceil( size - ret.size() / 2.0 ) );
        if ( (int)ret.size() <= 2 * size )
            ret += bars.at( _bar ).substr( 0, 1 );
        if ( level == LEVEL_1 ) {
            L1( "{}\n", ret );
        } else if ( level == LEVEL_2 ) {
            L2( "{}\n", ret );
        }
    }
    void Iwrapinbar( const std::string &msg, int size, int level, int _barOut, int _barIn ) {
        _barIn = ( _barIn == -1 ? _barOut : _barIn );
        Bar( size, level, _barOut );
        inBar( msg, size, level, _barIn );
        Bar( size, level, _barOut );
    }
    void Idebug( const std::string &msg, int counter ) {
        L1( "[DEBUG] {}{}\n", ( counter != -1 ? fmt::format( "({})", debug_counter++ ) : "" ), msg );
    }
    void Iclose() {
        L1( "[END OF LOGFILE]" );
        std::fclose( file );
    }
};