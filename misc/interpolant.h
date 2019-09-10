#pragma once
#include "ALGLIB/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <fmt/core.h> // -DFMT_HEADER_ONLY
#include <fmt/format.h>
#include "ALGLIB/interpolation.h"


class Interpolant {
   private:
    alglib::spline1dinterpolant p;
    std::string X, Y, Z, type;

   public:
    Interpolant() {}
    Interpolant( std::vector<double> &interpolationPointsX, std::vector<double> &interpolationPointsY, std::string t );
    Interpolant( std::vector<double> &interpolationPointsX, std::vector<double> &interpolationPointsY, std::vector<double> &interpolationPointsZ, std::string t );
    void generate( std::string t );
    double evaluate( double x ) const;
    std::vector<double> evaluate( std::vector<double> &xar );
};

// g++ .\interpolant.cpp .\*.o
// set yrange[-1.1:1.1]; plot "out.txt" u 1:2 w l t 'linear',"out.txt" u 1:3 w l t 'cubic',"out.txt" u 1:4 w l t 'hermite',"out.txt" u 1:5 w l t 'dakima',"out.txt" u 1:6 w l t 'catmullrom',"out.txt" u 1:7 w l t 'monotone'