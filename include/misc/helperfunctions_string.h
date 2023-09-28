#pragma once
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

namespace QDACC {

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

// Splits std::string once at "lit" and returns a pair of the split elements
std::pair<std::string, std::string> split_pair( const std::string &input, const std::string &lit = " " );

// Splits std::string at every "lit", reverses the vector and returns a joint string again
std::string split_and_reverse( const std::string &input, const std::string &lit = " " );

// Finds tofind in arr. Returns position if found, -1 if not.
int instr( const std::string &arr, const std::string tofind, int start = 0 );

// Converts an input string ( like e.g. [1,2,3,4] ) into a vector<string>
std::vector<std::string> str_to_vec( std::string input = "[]" );

// Finds a string in vector and returns index of first found element
int vec_find_str( const std::string &toFind, const std::vector<std::string> &input );

// Split input string at spli
std::vector<std::string> splitline( const std::string &input = "", const char splitter = ' ' );

std::vector<std::string> argv_to_vec( int argc, char **argv );

// Replaces part of string. Generates a new string
std::string replace( std::string input, const std::string &toreplace, const std::string &replaceto, bool ignore_capitals = false );

// Converts a String to lower case. Generates a new string
std::string to_lower( std::string str );

// Adds a prefix and a suffix around a word. Generates a new string
std::string add_prefix_and_suffix( std::string input, const std::string &tochange, const std::string &prefix = "", const std::string &suffix = "", bool ignore_capitals = false );

} // namespace String

} // namespace QDACC