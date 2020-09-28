#pragma once

#include <cmath>
#include <complex>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fmt/core.h> // -DFMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>

#define LOG_LEVEL_2 2
#define LOG_LEVEL_1 0
#define LOG_SIZE_FULL 75
#define LOG_SIZE_HALF 40
#define LOG_BAR_0 0
#define LOG_BAR_1 1
#define LOG_BAR_2 2
#define LOG_BAR_3 3
#define LOG_BAR_4 4
#define LOG_BAR_5 5
#define LOG_BAR_6 6

class Log {
   private:
    FILE *file;
    int deVar = 0;
    bool logLevel1, logLevel2;
    std::vector<std::string> bars;
    std::string repeat( std::string s, int n ) {
        std::string s1 = s;
        for ( int i = 1; i < n; i++ )
            s += s1;
        return s;
    }

   public:
    Log( std::string filepath, bool _level2 = false ) {
        file = std::fopen( filepath.c_str(), "w" );
        assert( file );
        std::setbuf( file, NULL );
        std::setbuf( stdout, NULL );
        logLevel2 = _level2;
        bars = {repeat( "=", LOG_SIZE_FULL ),
                repeat( "-", LOG_SIZE_FULL ),
                repeat( "o", LOG_SIZE_FULL ),
                repeat( "#", LOG_SIZE_FULL ),
                repeat( ":", LOG_SIZE_FULL ),
                repeat( "+", LOG_SIZE_FULL ),
                repeat( "*", LOG_SIZE_FULL )};
        if ( logLevel2 )
            fmt::print( file, "Succesfully created logfile '{}'!\n", filepath );
    }
    Log(){};

    template <typename... Args>
    void operator()( std::string s, const Args &... args ) {
        fmt::vprint( s, fmt::make_format_args( args... ) );
        fmt::vprint( file, s, fmt::make_format_args( args... ) );
    }
    template <typename... Args>
    void level1( std::string s, const Args &... args ) {
        fmt::vprint( s, fmt::make_format_args( args... ) );
        fmt::vprint( file, s, fmt::make_format_args( args... ) );
    }
    template <typename... Args>
    void level2( std::string s, const Args &... args ) const {
        if ( !logLevel2 )
            return;
        fmt::vprint( ": " + s, fmt::make_format_args( args... ) );
        fmt::vprint( file, ": " + s, fmt::make_format_args( args... ) );
    }

    void bar( int size = LOG_SIZE_FULL, int level = LOG_LEVEL_1, int _bar = LOG_BAR_0 ) const {
        size = std::min( std::max( size, 0 ), LOG_SIZE_FULL );
        std::string ret = bars.at( _bar ).substr( 0, size );
        if ( level == 0 ) {
            fmt::print( "{}{}\n", ret, ret );
            fmt::print( file, "{}{}\n", ret, ret );
        } else if ( level == 2 ) {
            level2( "{}{}\n", ret, ret );
        }
    }
    void inBar( std::string s, int size = LOG_SIZE_FULL, int level = LOG_LEVEL_1, int _bar = LOG_BAR_0 ) const {
        std::string ret = " | " + s + " | ";
        ret = bars.at( _bar ).substr( 0, std::floor( size - ret.size() / 2.0 ) ) + ret + bars.at( _bar ).substr( 0, std::ceil( size - ret.size() / 2.0 ) );
        if ( (int)ret.size() <= 2 * size )
            ret += bars.at( _bar ).substr( 0, 1 );
        if ( level == 0 ) {
            fmt::print( "{}\n", ret );
            fmt::print( file, "{}\n", ret );
        } else if ( level == 2 ) {
            level2( "{}\n", ret );
        }
    }

    void wrapInBar( std::string s, int size = LOG_SIZE_FULL, int level = LOG_LEVEL_1, int _barOut = LOG_BAR_0, int _barIn = -1 ) {
        _barIn = ( _barIn == -1 ? _barOut : _barIn );
        bar( size, level, _barOut );
        inBar( s, size, level, _barIn );
        bar( size, level, _barOut );
    }

    void debug( double i = -1, std::string suffix = "" ) {
        level2( "Currently at {} {}\n", ( i == -1 ) ? deVar++ : i, suffix );
    }
    void debug( std::string suffix = "" ) {
        debug( -1, suffix );
    }
    void close() {
        level1( "[END OF LOGFILE]" );
        std::fclose( file );
    }

};