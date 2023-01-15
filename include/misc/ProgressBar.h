#pragma once
#include <cmath>
//#include <complex>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
//#include <stdarg.h>
#include <string>
#include <iostream>
#include <omp.h>      // -fopenmp
#include <fmt/core.h> // -DFMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <vector>
#include <numeric> // accumulate

#define BAR_HORIZONTAL 0
#define BAR_VERTICAL 1

namespace ProgressBarSetting {
const std::vector<std::string> sym1 = { " ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█" };
const std::vector<std::string> sym2 = { " ", ".", ":", "|" };
const std::vector<std::string> sym3 = { " ", "\033[38;2;25;25;25m|\033[0m", "\033[38;2;50;50;50m|\033[0m", "\033[38;2;75;75;75m|\033[0m", "\033[38;2;100;100;100m|\033[0m", "\033[38;2;125;125;125m|\033[0m", "\033[38;2;150;150;150m|\033[0m", "\033[38;2;175;175;175m|\033[0m", "\033[38;2;200;200;200m|\033[0m", "\033[38;2;225;225;225m|\033[0m", "\033[38;2;250;250;250m|\033[0m" };
const std::vector<std::string> sym4 = { " ", "\033[38;2;25;25;25m=\033[0m", "\033[38;2;50;50;50m=\033[0m", "\033[38;2;75;75;75m=\033[0m", "\033[38;2;100;100;100m=\033[0m", "\033[38;2;125;125;125m=\033[0m", "\033[38;2;150;150;150m=\033[0m", "\033[38;2;175;175;175m=\033[0m", "\033[38;2;200;200;200m=\033[0m", "\033[38;2;225;225;225m=\033[0m", "\033[38;2;250;250;250m=\033[0m" };
const std::vector<std::string> sym5 = { " ", "\033[38;2;25;25;25m»\033[0m", "\033[38;2;50;50;50m»\033[0m", "\033[38;2;75;75;75m»\033[0m", "\033[38;2;100;100;100m»\033[0m", "\033[38;2;125;125;125m»\033[0m", "\033[38;2;150;150;150m»\033[0m", "\033[38;2;175;175;175m»\033[0m", "\033[38;2;200;200;200m»\033[0m", "\033[38;2;225;225;225m»\033[0m", "\033[38;2;250;250;250m»\033[0m" };
const std::vector<std::string> sym6 = { " ", "\033[38;2;25;25;25m═\033[0m", "\033[38;2;50;50;50m═\033[0m", "\033[38;2;75;75;75m═\033[0m", "\033[38;2;100;100;100m═\033[0m", "\033[38;2;125;125;125m═\033[0m", "\033[38;2;150;150;150m═\033[0m", "\033[38;2;175;175;175m═\033[0m", "\033[38;2;200;200;200m═\033[0m", "\033[38;2;225;225;225m═\033[0m", "\033[38;2;250;250;250m═\033[0m" };
} // namespace ProgressBarSetting

class ProgressBar {
   private:
    double currentPercent;
    std::string symstr;
    int num0, num1, num2;
    bool half;
    int c = 0;
    double lastSpin = 0;
    std::string addstr( int num, const std::string& str ) {
        std::string ret = "";
        for ( int i = 0; i < num; i++ ) {
            ret += str;
        }
        return ret;
    }

   public:
    int barLength, decimalPoints;
    std::string barEnd;
    std::vector<std::string> spin;
    std::vector<std::string> sym;
    std::string strBarStart = "[";
    std::string strBarEnd = "]";
    bool isSpinning;
    double spinPS;
    int maxSize = 0;
    int type;

    ProgressBar( int _barLength = 60, int _decimalPoints = 0, int _type = BAR_VERTICAL, bool _isSpinning = true, double _spinPS = 0.1, const std::vector<std::string>& _sym = ProgressBarSetting::sym2, const std::vector<std::string>& _spin = { "|", "/", "-", "\\" }, const std::string& _barEnd = "Done" ) : barLength( _barLength ),
                                                                                                                                                                                                                                                                                                                 decimalPoints( _decimalPoints ),
                                                                                                                                                                                                                                                                                                                 barEnd( _barEnd ),
                                                                                                                                                                                                                                                                                                                 spin( _spin ),
                                                                                                                                                                                                                                                                                                                 sym( _sym ),
                                                                                                                                                                                                                                                                                                                 isSpinning( _isSpinning ),
                                                                                                                                                                                                                                                                                                                 spinPS( _spinPS ),
                                                                                                                                                                                                                                                                                                                 type( _type ) {
        symstr = std::accumulate( sym.begin(), sym.end(), std::string( "" ) );
    };
    void calculate( int currentIterations, int maximumIterations ) {
        currentPercent = std::min( ( 1.0 * currentIterations ) / maximumIterations, 0.9999 );
        int size = (int)sym.size() - 1;
        if ( type == 0 ) {
            num0 = std::floor( currentPercent * size ); // min(*,sym.size()-1)
            num2 = std::floor( barLength * size * ( currentPercent - ( 1.0 * num0 ) / size ) );
            num1 = std::max( 0, barLength - num2 );
        } else {
            num0 = std::floor( currentPercent * barLength );
            num1 = std::floor( size * barLength * ( currentPercent - ( 1.0 * num0 ) / barLength ) );
            num1 %= size;
            // num2 = std::max(0,barLength-1-num0); // Redundant
        }
    }
    std::string print( int currentIterations, int maximumIterations, const std::string& barSuffix = "", const std::string& barPrefix = "", bool bold = true ) {
        calculate( currentIterations, maximumIterations );
        std::string ret = "";
        if ( type == 0 )
            ret = addstr( num2, sym.at( ( num0 + 1 ) % sym.size() ) ) + addstr( num1, sym.at( num0 ) );
        else {
            ret = addstr( num0, sym.back() ) + sym.at( num1 ) + addstr( std::max( 0, barLength - num0 - 1 ), sym.front() );
        }
        if ( omp_get_wtime() - lastSpin >= spinPS ) {
            lastSpin = omp_get_wtime();
            c++;
        }
        ret = "\033[2K\033[38;2;255;255;255m" + strBarStart + "\033[0m" + ret + "\033[38;2;255;255;255m" + strBarEnd + ( decimalPoints >= 0 ? fmt::format( " {:.{}f}\%", 1.0 * currentIterations / maximumIterations * 100, decimalPoints ) : "" ) + ( ( isSpinning && currentIterations < maximumIterations ) ? fmt::format( " [{}] ", spin.at( c % spin.size() ) ) : " " ) + ( currentIterations < maximumIterations ? barSuffix : barSuffix + " - " + barEnd ) + " - " + barPrefix + "\033[0m";
        maxSize = ( (int)ret.size() > maxSize ) ? (int)ret.size() : maxSize;
        // fmt::print( "{:<{}}\r", ret, maxSize );
        if ( bold )
            ret = "\033[1m" + ret + "\033[0m";
        fmt::print( "{}\r", ret );
        return ret;
    }
    // Uses the "spin" component to animate a waiting position
    std::string wait( const std::string& barSuffix = "", const std::string& barPrefix = "", bool bold = true ) {
        if ( omp_get_wtime() - lastSpin >= spinPS ) {
            lastSpin = omp_get_wtime();
            c++;
        }
        int position = ( c % ( 2 * sym.size() - 1 ) ) - sym.size() + 1;
        if ( position < 0 )
            position = -position;
        std::string ret = addstr( barLength, sym.at( position ) );
        ret = "\033[2K\033[38;2;255;255;255m" + strBarStart + "\033[0m" + ret + "\033[38;2;255;255;255m" + strBarEnd + " Waiting" + ( isSpinning ? fmt::format( " [{}] ", spin.at( c % spin.size() ) ) : " " ) + barSuffix + " - " + barPrefix + "\033[0m";
        maxSize = ( (int)ret.size() > maxSize ) ? (int)ret.size() : maxSize;
        // fmt::print( "{:<{}}\r", ret, maxSize );
        if ( bold )
            ret = "\033[1m" + ret + "\033[0m";
        fmt::print( "{}\r", ret );
        return ret;
    }
};