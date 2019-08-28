#pragma once

#include <cmath>
#include <complex>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
//#include <stdarg.h>
#include <string>
#include <fmt/core.h> // -DFMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h> // -fopenmp
#include <iostream>

std::string PREFIX_PERCENT               = "@#PERCENT#@";
std::string PREFIX_PERCENT_TIME          = "@#PERCENTTIME#@";
std::string PREFIX_PERCENT_TIME_FINAL    = "@#PERCENTTIMEFINAL#@";
std::string PREFIX_WARNING               = "@#WARNING#@";
std::string PREFIX_ERROR                 = "@#ERROR#@";
std::string PREFIX_DEBUG                 = "@#DEBUG#@";
std::string PREFIX_OUTPUT                = "@#OUTPUT#@";
std::string PREFIX_SUFFIX                = "@#SUFFIX#@";

#include "misc/log.h"
#include "misc/ProgressBar.h"
#include "misc/helperfunctions.h"
#include "misc/timer.h"

using Eigen::MatrixXcd;
using Eigen::MatrixXd; 

using namespace std::complex_literals;

std::string global_message_normaltermination = "[Programm Terminated normally]";
std::string global_message_error_divergent = "[FATAL ERROR: PROCESS DIVERGED, PROGRAM TERMINATED]";
std::string global_message_error_wrong_number_input = "[FATAL ERROR: WRONG NUMBER OF INPUT PARAMETERS, PROGRAM TERMINATED]";

#define DIR_T 1
#define DIR_TAU 2

#define TIMETRANSFORMATION_PRECALCULATED 0
#define TIMETRANSFORMATION_MATRIXEXPONENTIAL 1
#define TIMETRANSFORMATION_MATRIXEXPONENTIAL_PRECALCULATED 2

#define OPERATOR_PHOTONIC_CREATE 0
#define OPERATOR_PHOTONIC_ANNIHILATE 1