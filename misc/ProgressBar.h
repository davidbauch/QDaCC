#pragma once
#include <cmath>
//#include <complex>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
//#include <stdarg.h>
#include <string>
#include <iostream>
#include <omp.h> // -fopenmp
#include <fmt/core.h> // -DFMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>

#define BAR_HORIZONTAL 0
#define BAR_VERTICAL 1

class ProgressBar {
    private:
        double currentPercent;
        int num0,num1,num2;
        bool half;
        int c;
        double lastSpin;
        std::string addstr(int num, std::string str) {
            std::string ret = "";
            for (int i = 0; i < num; i++) {
                ret += str;
            }
            return ret;
        }
    public:
        int maximumIterations;
        int barLength, decimalPoints;
        std::string barEnd;
        std::vector<std::string> spin;
        std::vector<std::string> sym;
        std::string strBarStart;
        std::string strBarEnd;
        bool isSpinning;
        double spinPS;
        int maxSize;
        int type;
        ProgressBar(int _maximumIterations, int _barLength = 20, int _decimalPoints = 0, int _type = BAR_VERTICAL, bool _isSpinning = true, double _spinPS = 0.1, std::vector<std::string> _sym = {".", "-", "="}, std::vector<std::string> _spin = {"|","/","-","\\"}, std::string _barEnd = "Done") {
            maximumIterations   = _maximumIterations;
            barLength           = _barLength;
            decimalPoints       = _decimalPoints;
            isSpinning          = _isSpinning;
            spinPS              = _spinPS;
            lastSpin            = 0;
            sym                 = _sym;
            barEnd              = _barEnd;
            spin                = _spin;
            c                   = 0;
            maxSize             = 0;
            type                = _type;
            strBarStart         = "[";
            strBarEnd           = "]";
        }
        void calculate(int currentIterations) {
            currentPercent  = 1.0*currentIterations/maximumIterations;
            int size        = sym.size() - 1;
            if (type == 0) {
                num0 = std::floor(currentPercent*size); // min(*,sym.size()-1)
                num2 = std::floor( barLength*size*(currentPercent-(1.0*num0)/size) );
                num1 = std::max(0,barLength-num2);
            } else {
                num0 = std::floor(currentPercent*barLength);
                num1 = std::floor( size*barLength*(currentPercent-(1.0*num0)/barLength) ); num1%=size;
                //num2 = std::max(0,barLength-1-num0); // Redundant
            }
        }
        std::string print(int currentIterations, std::string barSuffix = "", std::string barPrefix = "") {
            calculate(currentIterations);
            std::string ret = "";
            if (type == 0)
                ret = addstr(num2, sym.at((num0+1)%sym.size())) + addstr(num1, sym.at(num0));
            else {
                ret = addstr(num0, sym.back()) + sym.at(num1) + addstr(std::max(0,barLength - num0-1), sym.front());
            }
            if (omp_get_wtime() - lastSpin >= spinPS) { 
                lastSpin = omp_get_wtime();
                c++;
            }
            ret = barPrefix + strBarStart + ret + strBarEnd
             + ( decimalPoints >= 0 ? fmt::format(" {:.{}f}\%", 1.0*currentIterations/maximumIterations*100, decimalPoints) : "" )
             + ( (isSpinning && currentIterations < maximumIterations) ? fmt::format(" [{}] ",spin.at(c%spin.size())) : " " )
             + ( currentIterations < maximumIterations ? barSuffix : barEnd );
            maxSize = (ret.size()>maxSize) ? ret.size() : maxSize;
            fmt::print("{:<{}}\r",ret,maxSize);
            return ret;
        }
};