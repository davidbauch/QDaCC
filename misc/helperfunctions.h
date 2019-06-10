#pragma once
#include "../global.h"

#ifndef M_PI
#define M_PI 3.1415926535
#endif

int delta(int i, int j) {
    return (i==j) ? 1 : 0;
}

double Hz_to_eV(double hz){
    return hz*6.582119516885722624E-16;
}

double eV_to_Hz(double eV){
    return eV/6.582119516885722624E-16;
}

double Hz_to_wavelength(double hz){
    return 299792458*2*M_PI/hz*1E9;
}

double rabiFrequency(double deltaE, double g, double n) {
    return std::sqrt( std::pow(deltaE,2) + 4.*std::pow(g,2)*n );
}

// is target in arr of length len?
// returns index where target was found, -1 if target is not in arr
int in(char **arr, int len, char *target) {
  int i;
  for(i = 0; i < len; i++) {
    if(std::strncmp(arr[i], target, strlen(target)) == 0) {
      return i;
    }
  }
  return -1;
}

int instr(const std::string &arr, const std::string tofind, int start = 0) {
    bool found = true;
    for (int i = start; i < arr.size()+1-tofind.size(); i++) {
        found = true;
        for (int j = 0; j < tofind.size(); j++) {
            //fmt::print("comparing {} and {}... ",arr.at(i+j),tofind.at(j));
            if (tofind.at(j) != arr.at(i+j)) {
                found = false;
                j = tofind.size();
            }
        }
        if (found) {
            //fmt::print("found at index {}\n",i);
            return i;
        }
    }
    return -1;
}

// Converts an input string ( like e.g. [1,2,3,4] ) into a vector<string>
std::vector<std::string> str_to_vec(std::string input) {
    std::vector<std::string> ret;
    int s = 1; // Starting index
    int e = 1; // End
    while ((e = instr(input,",",s)) != -1) {
        ret.emplace_back( input.substr(s,e-s) );
        s = e+1;
    }
    ret.emplace_back( input.substr(s,input.size()-s-1) );
    //for (std::string el : ret)
    //    std::cout << el << "\n";
    return ret;
}

#define toStr(x) std::to_string(x)

template <typename T>
T lerp(T a, T b, double c) {
    if (c<0 || c > 1) {
        //logs("Warning: Invalid lerp delta! delta = {:e}\n",c);
        return a;
    }
    return (1.0-c)*a + c*b;
}

bool is_number(const std::string& s)
{
    //for (int i = 0; i < s.size(); i++) {
    //    if (!(std::isdigit(c) || c == 'e' || c == 'E' || c == '-' || c == '+' || c == '.'))
    //}
    return !s.empty() && std::find_if(s.begin(), 
        s.end(), [](char c) { return (!(std::isdigit(c) || c == 'e' || c == 'E' || c == '-' || c == '+' || c == '.')); }) == s.end();
}

// Valid conversions from ns,ps,Hz,eV,meV,mueV to corresponding SI unit (not scaled)
template <typename T>
T convertParam(const std::string input) {
    double value = 0;
    double conversion = 1;
    int index;
    if ( -1 != (index = instr(input,"eV"))) {
        // Found 'eV' as unit (energy), now check for scaling
        if (input.at(index-1) == 'm') {
            // meV
            logs.level2("from meV to Hz...");
            value = eV_to_Hz( std::stod(input.substr(0,index-1)) );
            conversion = 1E-3;
        } else if (input.at(index-2) == 'm' && input.at(index-1) == 'u') {
            // mueV
            logs.level2("from mueV to Hz...");
            value = eV_to_Hz( std::stod(input.substr(0,index-2)) );
            conversion = 1E-6;
        } else if (is_number(input.substr(index-1,1))) {
            // eV
            logs.level2("from eV to Hz...");
            value = eV_to_Hz( std::stod(input.substr(0,index)) );
            conversion = 1.0;
        } else {
            logs.level2("Conversion of input '{}' from eV failed!\n",input);
            return (T)0.0;
        }
    } else if (-1 != (index = instr(input,"s"))) { 
        // Found 's' as unit (time)
        //fmt::print("\n {} {} {} {}\n",index, input.at(index-1)=='n',input.compare(index-1,1,"n") ,input.at(index-1));
        if (input.at(index-1) == 'n') {
            // ns
            logs.level2("from ns to s...");
            value = std::stod(input.substr(0,index-1));
            conversion = 1E-9;
        } else if (input.at(index-1) == 'p') {
            // ps
            logs.level2("from ps to s...");
            value = std::stod(input.substr(0,index-1));
            conversion = 1E-12; //fmt::print("{} {} ... ", value, conversion);
        } else if (is_number(input.substr(index-1,1))) {
            // s
            logs.level2("from s to s...");
            value = std::stod(input.substr(0,index));
            conversion = 1.0;
        } else {
            logs.level2("Conversion from input '{}' from time failed!\n");
            return (T)0.0;
        }
    } else if ( -1 != (index = instr(input,"Hz"))) {
        // Found 'Hz' as unit (Frequency)
            logs.level2("from Hz to Hz...");
        if (is_number(input.substr(index-1,1))) {
            value = std::stod(input.substr(0,index-1));
            conversion = 1.0;
        } else {
            logs.level2("Conversion from input '{}' from frequency failed!\n");
            return (T)0.0;
        }
    } else if (is_number(input)){ 
        // Assuming Frequency input
            logs.level2("no conversion...");
        value = std::stod(input);
    } else {
        // Input type is unknown
        logs.level2("Input Type of input '{}' is unkown!\n",input);
            return (T)0.0;
    }
    logs.level2("done, final value = {}!\n",value*conversion);
    return (T)(value*conversion);
}

template <typename T>
std::vector<T> convertParam(const std::vector<std::string> input) {
    std::vector<T> ret;
    for (std::string in : input) {
        ret.emplace_back( convertParam<T>(in) );
    }
    return ret;
}

template <typename T>
T getNextInput(const std::vector<std::string> &arguments, const std::string name, int &index) {
    logs.level2("Trying to convert input named '{}' at index '{}' to '{}'...",name,index,arguments.at(index));
    return convertParam<T>(arguments.at(index++));
}
std::string getNextInputString(const std::vector<std::string> &arguments, const std::string name, int &index) {
    logs.level2("Trying to convert input named '{}' at index '{}' to '{}'... done\n",name,index,arguments.at(index));
    return arguments.at(index++);
}
template <typename T>
std::vector<T> getNextInputVector(const std::vector<std::string> &arguments, const std::string name, int &index) {
    logs.level2("Trying to convert input named '{}' at index '{}' to '{}'...",name,index,arguments.at(index));
    return convertParam<T>(str_to_vec(arguments.at(index++)));
}
/*template <typename T>
std::enable_if_t<std::is_same<std::string, T>::value, std::string> getNextInput(const std::vector<std::string> &arguments, const std::string name, int &index) {
    logs.level2("Trying to convert input number named '{}' at index '{}' to '{}'... ",name,index,arguments.at(index));
    return arguments.at(index++);
}*/