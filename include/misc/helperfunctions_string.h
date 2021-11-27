#pragma once
#include <string>
#include <sstream>
#include <vector>

namespace QDLC {

namespace String {

// Extracts the value between 2 lit in a string
std::string strip( std::string input, char lit = '\'' );

// Extracts the value between 2 lit in a string in a vector of strings
std::vector<std::string> strip( std::vector<std::string> input, char lit = '\'' );

// Length norm an input string to a totallength. 
std::string trail( std::string input, int totallength, std::string pre = " ", int dir = 0 );
std::string tail( std::string input, int totallength, std::string pre = " " );

// Returns true if input starts with find
bool startswith( std::string input, std::string find );

// Splits std::string at every "lit" and returns std::vector of std::string 
std::vector<std::string> split( std::string input, std::string lit = " " );

// Finds tofind in arr. Returns position if found, -1 if not.
int instr( const std::string &arr, const std::string tofind, int start = 0 );

// Converts an input string ( like e.g. [1,2,3,4] ) into a vector<string>
std::vector<std::string> str_to_vec( std::string input = "[]" );

// Finds a string in vector and returns index of first found element
int vec_find_str( std::string toFind, const std::vector<std::string> &input );

// Split input string at spli
std::vector<std::string> splitline( const std::string &input = "", const char splitter = ' ' );

std::vector<std::string> argv_to_vec( int argc, char **argv );

} // namespace String

} // namespace QDLC