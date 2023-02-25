#include "misc/helperfunctions_string.h"

std::string QDLC::String::strip( std::string input, char lit ) {
    int a = input.find_first_of( lit ) + 1;
    int e = input.find_last_of( lit );
    // maybe: return "" if length is 0
    return input.substr( a, e - a );
}

std::vector<std::string> QDLC::String::strip( std::vector<std::string> input, char lit ) {
    std::vector<std::string> ret;
    for ( const auto &str : input ) ret.emplace_back( strip( str ) );
    return ret;
}

std::string QDLC::String::trail( std::string input, int totallength, std::string pre, int dir ) {
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

std::string QDLC::String::tail( std::string input, int totallength, std::string pre ) {
    return trail( input, totallength, pre, 1 );
}

bool QDLC::String::startswith( std::string input, std::string find ) {
    return input.substr( 0, find.size() ).compare( find ) == 0;
}

std::vector<std::string> QDLC::String::split( std::string input, std::string lit ) {
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

std::string QDLC::String::split_and_reverse( const std::string &input, const std::string &lit ) {
    auto vec = split( input, lit );
    std::ranges::reverse( vec );
    std::string ret = "";
    for ( int i = 0; i < vec.size(); i++ ) {
        ret += ( i == vec.size() - 1 ) ? vec[i] : vec[i] + lit;
    }
    return ret;
}

std::pair<std::string, std::string> QDLC::String::split_pair( const std::string &input, const std::string &lit ) {
    auto pos = input.find( lit );
    return std::make_pair( input.substr( 0, pos ), input.substr( pos + lit.size() ) );
}

int QDLC::String::instr( const std::string &arr, const std::string tofind, int start ) {
    auto pos = arr.find(tofind, start);
    if (pos == std::string::npos)
        return -1;
    return pos;
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
    // for (std::string el : ret)
    //     std::cout << el << "\n";
    return ret;
}

int QDLC::String::vec_find_str( const std::string &toFind, const std::vector<std::string> &input ) {
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
        ret.emplace_back( std::string( argv[i] ) );
    return ret;
}

std::string QDLC::String::replace( std::string input, const std::string &toreplace, const std::string &replaceto, bool ignore_capitals ) {
    auto pos = input.find( toreplace );
    if ( ignore_capitals and pos == std::string::npos )
        pos = QDLC::String::to_lower( input ).find( QDLC::String::to_lower( toreplace ) );
    if ( pos == std::string::npos ) return input;
    return input.substr( 0, pos ) + replaceto + input.substr( pos + toreplace.size(), input.size() );
}

std::string QDLC::String::add_prefix_and_suffix( std::string input, const std::string &tochange, const std::string &prefix, const std::string &suffix, bool ignore_capitals ) {
    auto pos = input.find( tochange );
    if ( ignore_capitals and pos == std::string::npos )
        pos = QDLC::String::to_lower( input ).find( QDLC::String::to_lower( tochange ) );
    if ( pos == std::string::npos ) return input;
    return input.substr( 0, pos ) + prefix + input.substr( pos, tochange.size() ) + suffix + input.substr( pos + tochange.size(), input.size() );
}

std::string QDLC::String::to_lower( std::string str ) {
    std::ranges::transform( str.begin(), str.end(), str.begin(), []( unsigned char c ) { return std::tolower( c ); } );
    return str;
}