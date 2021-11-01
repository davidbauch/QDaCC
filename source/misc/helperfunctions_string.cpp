#include "misc/helperfunctions_string.h"

std::string QDLC::String::strip( std::string input, char lit ){
    int a = input.find_first_of( lit ) + 1;
    int e = input.find_last_of( lit );
    // maybe: return "" if length is 0
    return input.substr( a, e - a );
}

std::vector<std::string> QDLC::String::strip( std::vector<std::string> input, char lit ){
    std::vector<std::string> ret;
    for ( const auto& str : input ) ret.emplace_back( strip( str ) );
    return ret;
}

std::string QDLC::String::trail( std::string input, int totallength, std::string pre, int dir ){
    if ( pre.size() == 0 )
        return input;
    int remaining = ( totallength - input.size() ) / pre.size();
    int trailing = ( totallength - input.size() ) % pre.size();
    for ( int i = 0; i < remaining; i++ )
        if ( dir == 0 )
            input = pre + input;
        else
            input = input + pre;
    if ( dir == 0 )
        input = pre.substr( 0, trailing ) + input;
    else
        input = input + pre.substr( 0, trailing );
    return input;
}

std::string QDLC::String::tail( std::string input, int totallength, std::string pre ){
    return trail( input, totallength, pre, 1 );
}

bool QDLC::String::startswith( std::string input, std::string find ){
    return input.substr( 0, find.size() ).compare( find ) == 0;
}

std::vector<std::string> QDLC::String::split( std::string input, std::string lit ){
    std::vector<std::string> ret;
    long unsigned int i = 0;
    long unsigned int start = 0;
    while ( i < input.size() ) {
        if ( input.at( i ) == lit.at( 0 ) ) {
            bool legit = true;
            for ( long unsigned int j = 1; j < lit.size() && i + j < input.size(); j++ )
                if ( input.at( i + j ) != lit.at( j ) ) legit = false;
            if ( legit ) {
                ret.emplace_back( input.substr( start, i - start ) );
                start = i + lit.size();
            }
        }
        i++;
    }
    ret.emplace_back( input.substr( start ) );
    return ret;
}

int QDLC::String::instr( const std::string &arr, const std::string tofind, int start ) {
    bool found = true;
    for ( int i = start; i < (int)arr.size() + 1 - (int)tofind.size(); i++ ) {
        found = true;
        for ( int j = 0; j < (int)tofind.size(); j++ ) {
            //fmt::print("comparing {} and {}... ",arr.at(i+j),tofind.at(j));
            if ( tofind.at( j ) != arr.at( i + j ) ) {
                found = false;
                j = (int)tofind.size();
            }
        }
        if ( found ) {
            //fmt::print("found at index {}\n",i);
            return i;
        }
    }
    return -1;
}

std::vector<std::string> QDLC::String::str_to_vec( std::string input ) {
    std::vector<std::string> ret;
    if ( (int)input.size() < 3 ) {
        ret = {};
        return ret;
    }
    int s = 1; // Starting index
    int e = 1; // End
    while ( ( e = QDLC::String::instr( input, ",", s ) ) != -1 ) {
        ret.emplace_back( input.substr( s, e - s ) );
        s = e + 1;
    }
    ret.emplace_back( input.substr( s, input.size() - s - 1 ) );
    //for (std::string el : ret)
    //    std::cout << el << "\n";
    return ret;
}

int QDLC::String::vec_find_str( std::string toFind, const std::vector<std::string> &input ) {
    for ( int i = 0; i < (int)input.size(); i++ ) {
        if ( input.at( i ).compare( toFind ) == 0 )
            return i;
    }
    return -1;
}

std::vector<std::string> QDLC::String::splitline( const std::string &input, const char splitter ) {
    std::string token, subtoken;
    std::vector<std::string> set;
    std::istringstream iss( input );
    while ( std::getline( iss, subtoken, '\t' ) ) {
        std::istringstream subiss( subtoken );
        while ( std::getline( subiss, token, splitter ) )
            if ( token.size() > 0 ) {
                set.push_back( token );
            }
    }
    return set;
}

std::vector<std::string> QDLC::String::argv_to_vec( int argc, char **argv ) {
    std::vector<std::string> ret;
    ret.reserve( argc );
    for ( int i = 0; i < argc; i++ )
        ret.push_back( std::string( argv[i] ) );
    return ret;
}