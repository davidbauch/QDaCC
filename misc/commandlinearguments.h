#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

// V 1.1 - Subkeys of groups can also be read individually, e.g. --g1 arg arg with argkeys arg1,arg2 result in --arg1 arg being also a valid parameter that can be read. If the parameter does not have a subkey, it cannot be individually set
// V 1.1.1 - Init length formatting now better
// V 1.2 - Better help function, --help arg can now search for arg, no more sorting for now
// V 1.2.1 - Now initialized with global init function. The global object 'cla' is now the standard CommandlineArguments class and can be accessed via the CommandlineArgument namespace.
// V 1.2.2 - Constructor now also accepts argv in vecor format, and commandlinearguments vector get function

namespace CommandlineArguments {

// Extracts the value between 2 lit in a string
std::string strip( std::string input, char lit = '\'' ) {
    int a = input.find_first_of( lit ) + 1;
    int e = input.find_last_of( lit );
    // maybe: return "" if length is 0
    return input.substr( a, e - a );
}
std::vector<std::string> strip( std::vector<std::string> input, char lit = '\'' ) {
    std::vector<std::string> ret;
    for ( const auto& str : input ) ret.emplace_back( strip( str ) );
    return ret;
}
std::string trail( std::string input, int totallength, std::string pre = " ", int dir = 0 ) {
    if (pre.size() == 0)
        return input;
    int remaining = (totallength - input.size())/pre.size();
    int trailing = (totallength - input.size())%pre.size();
    for ( int i = 0; i < remaining; i++ )
        if ( dir == 0 )
            input = pre + input;
        else
            input = input + pre;
    if ( dir == 0 )
            input = pre.substr(0,trailing) + input;
        else
            input = input + pre.substr(0,trailing);
    return input;
}
std::string tail( std::string input, int totallength, std::string pre = " " ) { return trail( input, totallength, pre, 1 ); }

bool startswith( std::string input, std::string find ) { return input.substr( 0, find.size() ).compare( find ) == 0; }

std::vector<std::string> splitString( std::string input, std::string lit = " " ) {
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

class CommandlineArguments {
   private:
    std::vector<std::string> configfile;
    std::vector<std::string> commandlinearguments;
    char cla_splitter = '=';
    int cla_width = 50;
    std::string cla_add_param = "add";
    std::string cla_remove_param = "rem"; //TODO
    std::string cla_edit_param = "mod"; //TODO
    std::string version = "1.2.2";

    class Parameter {
       public:
        std::string subkey;
        std::string value;
        std::string datatype; // muss beim einlesen auf jeden fall gecheckt werden ob eingabe dem richtigen dt entspricht
        std::string description;
        Parameter(){};
        Parameter( std::string _subkey, std::string _value, std::string _datatype, std::string _description ) : subkey( _subkey ), value( _value ), datatype( _datatype ), description( _description ) {}
        Parameter( const Parameter& other ) : subkey( other.subkey ), value( other.value ), datatype( other.datatype ), description( other.description ){};

        // Explicit parameter Conversions:
        int toInt() const {
            if ( isempty() ) std::cerr << "Warning, this parameter is empty!" << std::endl;
            try {
                return (int)std::stoi( value );
            } catch ( const std::exception& ) {
                std::cerr << "Couldn't convert '" << value << "' to int!" << std::endl;
                return 0;
            }
        };
        float toFloat() const {
            if ( isempty() ) std::cerr << "Warning, this parameter is empty!" << std::endl;
            try {
                return (float)std::stof( value );
            } catch ( const std::exception& ) {
                std::cerr << "Couldn't convert '" << value << "' to float!" << std::endl;
                return 0;
            }
        };
        double toDouble() const {
            if ( isempty() ) std::cerr << "Warning, this parameter is empty!" << std::endl;
            try {
                return (double)std::stod( value );
            } catch ( const std::exception& ) {
                std::cerr << "Couldn't convert '" << value << "' to double!" << std::endl;
                return 0;
            }
        };
        std::string toString() const {
            if ( isempty() ) std::cerr << "Warning, this parameter is empty!" << std::endl;
            return value;
        };
        bool toBool() const {
            if ( isempty() ) std::cerr << "Warning, this parameter is empty!" << std::endl;
            return ( value.compare( "true" ) == 0 ? true : false );
        };
        // Implicit paramter Conversions:
        operator int() const { return toInt(); }
        operator float() const { return toFloat(); }
        operator double() const { return toDouble(); }
        operator bool() const { return toBool(); }
        operator std::string() const { return toString(); }

        // Checks if key used comparse to this datastructure
        bool validkey( std::string kk ) const {
            if ( subkey.length() == 0 ) return false;
            if ( subkey.compare( kk ) == 0 ) return true;
            return false;
        }

        // Checks if this is the empty element
        bool isempty() const { return subkey.compare( "empty" ) == 0 && value.compare( "0" ) == 0; }

        // Checks if conversion is possible
        bool canConvert( std::string input ) {
            if ( datatype.compare( "bool" ) == 0 ) return true;
            if ( datatype.compare( "array" ) == 0 && input.front() == '[' && input.back() == ']' ) return true; // TODO: was hiermit
            if ( datatype.compare( "string" ) == 0 ) return true;
            if ( datatype.compare( "int" ) == 0 ) {
                try {
                    std::stoi( input );
                    return true;
                } catch ( const std::exception& ) {
                    return false;
                }
            }
            if ( datatype.compare( "float" ) == 0 ) {
                try {
                    std::stof( input );
                    return true;
                } catch ( const std::exception& ) {
                    return false;
                }
            }
            if ( datatype.compare( "double" ) == 0 ) {
                try {
                    std::stod( input );
                    return true;
                } catch ( const std::exception& ) {
                    return false;
                }
            }
            return false;
        }

        void identify() const { std::cout << "Parameterclass for " << subkey << ", value = " << value << ", datatype = " << datatype << ", description = " << description << std::endl; }
    };

    Parameter param_empty = Parameter( "empty", "0", "0", "0" );
    Parameter param_bool = Parameter( "", "false", "bool", "" );

    class Datastructure {
       public:
        std::vector<std::string> key;
        std::string description;
        std::string type;
        std::string group;
        std::vector<Parameter> parameter;
        Datastructure(){};
        Datastructure( const Datastructure& other ) : key( other.key ), description( other.description ), type( other.type ), parameter( other.parameter ), group(other.group){};
        Datastructure( std::vector<std::string> key, std::string description, std::string type, std::vector<Parameter> parameter, std::string group ) : key(key), description(description), type(type), parameter(parameter), group(group){};

        // Checks if key used comparse to this datastructure. if use_hyphens is True, it is assumed that kk is passed as with either - or -- prefix
        bool validkey( std::string kk, bool use_hyphens = false ) const {
            for ( const auto& k : key ) {
                if ( kk.compare( (use_hyphens ? (type.compare("1hyphen") == 0 ? "-" : "--") : "") + k ) == 0 ) return true;
            }
            return false;
        }

        // Looks for keyword filter
        bool validfilter( std::string filter ) const{
            if (validkey(filter) || key.at(0).find(filter) != std::string::npos) return true;
            if (group.compare(filter) == 0 || group.find(filter) != std::string::npos) return true;
            for (const auto& word : splitString(description)) {
                if (word.compare(filter) == 0 || word.find(filter) != std::string::npos) return true;
            }
            for (const auto& param : parameter) {
                if (param.validkey(filter) || param.subkey.find(filter) != std::string::npos) return true;
                if (param.datatype.compare(filter) == 0 || param.datatype.find(filter) != std::string::npos) return true;
                for (const auto& word : splitString(param.description)) {
                    if (word.compare(filter) == 0 || word.find(filter) != std::string::npos) return true;
                }
            }
            return false;
        }

        void identify( int len1, int len2, std::string div = "|     ", std::ostream& out = std::cout ) const {
            int lines_needed = std::max( (int)( description.size() / len2 ) + 1, (int)key.size() );
            std::vector<std::string> lines;
            std::string cur;
            std::string datatypes = "<";
            for ( const auto& param : parameter ) {
                datatypes += param.datatype + ",";
            }
            datatypes.back() = '>';
            for ( int i = 0; i < lines_needed; i++ ) {
                std::string keyline = tail( ( i < key.size() ? ( ( type.compare( "2hyphen" ) == 0 ? "--" : " -" ) + key.at( i ) ) : "" ) + ( i == 0 ? " " + datatypes : "" ), len1 );
                std::string desline = div + tail( ( i * len2 < description.size() ? ( description.at( i * len2 ) == ' ' ? description.substr( i * len2 + 1, len2 ) : description.substr( i * len2, len2 ) ) : "" ), len2 );
                lines.emplace_back( keyline + desline );
            }
            for ( const auto& param : parameter ) {
                if ( param.datatype.compare( "bool" ) != 0 ) {
                    std::string prefix = " | <" + param.datatype + "> Standard is " + param.value + ". ";
                    int newlen = len2 - prefix.size();
                    bool doneFirstOutput = false;
                    for ( int i = 0; i < (int)( param.description.size() / newlen ) + 1; i++ ) {
                        std::string desline = tail( (!doneFirstOutput ? (param.subkey.length() > 0 ? "  --" : "")+param.subkey : ""), len1 ) + div + ( i == 0 ? prefix : trail( "", prefix.size() ) ) + param.description.substr( i * newlen, newlen );
                        lines.emplace_back( desline );
                        doneFirstOutput = true;
                    }
                }
            }
            for ( const auto& s : lines ) {
                out << s << std::endl;
            }
        }
    };

    std::vector<Datastructure> cla_datastructures;

    bool readConfigFile( std::string filepath ) {
        std::ifstream in( ( filepath + ".settings.cla" ).c_str() );
        if ( !in ) {
            std::cerr << "Cannot open the File : " << filepath << std::endl;
            if ( !generateEmptyConfigfile( filepath ) ) return false;
        }
        std::string str;
        while ( std::getline( in, str ) ) {
            if ( str.size() > 0 && str.at( 0 ) != '#' ) {
                str = str.substr( str.find_first_not_of( "\t " ) );
                configfile.push_back( str );
            }
        }
        in.close();
        return true;
    }

    // Reads all CLA Settings from input file, such as splitter, width, etc
    bool readCLA_SettingsFromConfig() {
        for ( auto line : configfile ) {
            if ( startswith( line, "CLA_SPLIT = " ) ) {
                cla_splitter = *strip( line ).c_str();
            } else if ( startswith( line, "CLA_WIDTH = " ) ) {
                cla_width = std::stoi( strip( line ) );
            } else if ( startswith( line, "CLA_ADD_PARAMETER = " ) ) {
                cla_add_param = strip( line );
            } else if ( startswith( line, "CLA_REMOVE_PARAMETER = " ) ) {
                cla_remove_param = strip( line );
            } else if ( startswith( line, "CLA_EDIT_PARAMETER = " ) ) {
                cla_edit_param = strip( line );
            }
        }
        // std::cout << cla_splitter << ", " << cla_width << ", " << cla_add_param
        // << ", " << cla_remove_param << ", " << cla_edit_param << std::endl;
        return true;
    }

    // Reads one CLA Parameter
    Parameter readCLA_ParameterFromConfig( long unsigned int& i ) {
        Parameter ret = Parameter();
        while ( configfile.at( i ).at( 0 ) != '}' ) {
            if ( startswith( configfile.at( i ), "subkey = " ) )
                ret.subkey = strip( configfile.at( i ) );
            else if ( startswith( configfile.at( i ), "value = " ) )
                ret.value = strip( configfile.at( i ) );
            else if ( startswith( configfile.at( i ), "datatype = " ) )
                ret.datatype = strip( configfile.at( i ) );
            else if ( startswith( configfile.at( i ), "description = " ) )
                ret.description = strip( configfile.at( i ) );
            i++;
        }
        return ret;
    }

    // Reads all the CLA Arguments. If an argument was found, read it and save it
    // to Arguments vector
    bool readCLA_ArgumentsFromConfig() {
        int cla_arguments_found = 0;
        bool currently_in_cla_argument = false;
        Datastructure temp;
        for ( long unsigned int i = 0; i < configfile.size(); i++ ) {
            if ( startswith( configfile.at( i ), "CLA_ARGUMENT = {" ) && !currently_in_cla_argument ) {
                currently_in_cla_argument = true;
                temp = Datastructure();
            }
            if ( currently_in_cla_argument ) {
                if ( startswith( configfile.at( i ), "type = " ) )
                    temp.type = strip( configfile.at( i ) );
                else if ( startswith( configfile.at( i ), "key = " ) )
                    temp.key = strip( splitString( configfile.at( i ), "or" ) );
                else if ( startswith( configfile.at( i ), "description = " ) )
                    temp.description = strip( configfile.at( i ) );
                else if ( startswith( configfile.at( i ), "group = " ) )
                    temp.group = strip( configfile.at( i ) );
                else if ( startswith( configfile.at( i ), "CLA_PARAMETER = {" ) )
                    temp.parameter.emplace_back( readCLA_ParameterFromConfig( i ) );
                else if ( configfile.at( i ).at( 0 ) == '}' ) {
                    cla_datastructures.emplace_back( temp );
                    currently_in_cla_argument = false;
                }
            }
        }
        //map_index();
        //std::sort(cla_datastructures.begin(), cla_datastructures.end(), [](const Datastructure& a, const Datastructure& b)->bool {std::cout << "is " << a.group_index << " greater than " << b.group_index << "?" << std::endl; return a.group_index > b.group_index;});
        return true;
    }

    bool help(std::string filter = "") {
        std::cout << "Help Commandline Arguments v." << version << std::endl;
        int detlen = 50;
        for (const auto& ds : cla_datastructures) {
            if (filter.size() > 0 && !ds.validfilter(filter)) continue;
            int detlen1 = 5 + ds.key.at(0).length();
            for (const auto& p : ds.parameter) {
                detlen1 += 1 + p.datatype.length();
            }
            detlen = std::max(detlen,detlen1);
        }
        int len1 = std::max( (int)( cla_width / 3 ), detlen );
        std::string div = "|     ";
        int len2 = cla_width - len1 - div.size();
        std::cout << tail( "", cla_width, "=" ) << std::endl;
        std::string current_group = "";
        bool first = true;
        if (cla_datastructures.size() > 0)
            current_group = cla_datastructures.at(0).group;
        for ( const auto& datastructure : cla_datastructures ) {
            if (filter.size() > 0 && !datastructure.validfilter(filter)) continue;
            if (datastructure.group.compare(current_group) != 0) {
                std::cout << tail( "", cla_width, "-" ) << std::endl;
                current_group = datastructure.group;
            } else if (!first) {
                std::cout << tail("",len1) + "-" + tail("",len2+div.size()-1) << std::endl;//<< tail( "", cla_width, " . " ) << std::endl;
            }
            datastructure.identify( len1, len2, div );
            first = false;
        }
        return true;
    }

    // Parses the cla arguments from commandline and replaces parameters values.
    bool parseCLA_Arguments( std::vector<std::string>& cla_vector, std::string filepath ) {
        for ( long unsigned int i = 0; i < cla_vector.size(); i++ ) {
            auto str = cla_vector.at(i);
            if ( str.compare( "-h" ) == 0 || str.compare( "--help" ) == 0 ) {
                std::string key = "";
                if (str.compare( "--help" ) == 0 && cla_vector.size() > i+1) {
                    key = cla_vector.at(i+1);
                }
                help(key);
                exit( 1 );
            }
        }
        for ( const auto& str : cla_vector ) {
            if ( str.compare( "--" + cla_add_param ) == 0 || str.compare( "-" + cla_add_param ) == 0 ) {
                cla_wizzard_add(filepath);
                exit( 1 );
            }
        }
        // first, look for indizes where arguments start. it is important, that after a "-" occurs, the next char has to be non-numerical, and as such is not a number but a flag.
        int index = 0;
        for ( const auto& str : cla_vector ) {
            if ( str.at( 0 ) == '-' && !std::isdigit( str.at( 1 ) ) ) {
                bool found = false;
                if ( str.size() > 1 && str.at( 1 ) != '-' ) {
                    for ( auto& datastructure : cla_datastructures ) {
                        // 1 hyphen arguments can only have true/false values
                        if ( datastructure.type.compare( "1hyphen" ) == 0 && datastructure.validkey( str.substr( 1 ) ) ) {
                            for ( auto& param : datastructure.parameter ) {
                                param.value = "true";
                            }
                            found = true;
                        }
                    }
                } else if ( str.size() > 2 && str.at( 1 ) == '-' ) {
                    for ( auto& datastructure : cla_datastructures ) {
                        // Look for whole group
                        if ( datastructure.type.compare( "2hyphen" ) == 0 && datastructure.validkey( str.substr( 2 ) ) ) {
                            int subindex = 1;
                            for ( auto& param : datastructure.parameter ) {
                                if ( subindex > datastructure.parameter.size() ) break;
                                auto nextstr = cla_vector.at( index + subindex );
                                if ( param.canConvert( nextstr ) ) {
                                    param.value = nextstr;
                                    subindex++;
                                }
                            }
                            found = true;
                        }
                        // Look for single subparameter keys (as of version 1.1). Not that subkeys cannot be single hyphen.
                        if ( !found && datastructure.type.compare( "2hyphen" ) == 0 ) {
                            for ( auto& param : datastructure.parameter ) {
                                if ( param.validkey( str.substr( 2 ) ) ) {
                                    auto nextstr = cla_vector.at( index + 1 );
                                    if ( param.canConvert( nextstr ) ) {
                                        param.value = nextstr;
                                        found = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if ( !found ) {
                    std::cout << "Unrecognized commandline option '" << str << "'" << std::endl;
                    help();
                    exit( 0 );
                };
            }
            index++;
        }
        return true;
    }

    // Converts argv into vector:
    bool convert_argv( int argc, char* argv[] ) {
        commandlinearguments.clear();
        for ( int i = 0; i < argc; i++ ) {
            commandlinearguments.emplace_back( argv[i] );
        }
        return true;
    }

    bool writeConfigfile( std::string filepath ) {
        if ( configfile.size() > 0 ) {
            std::ofstream out( ( filepath + ".settings.cla" ).c_str() );
            if ( !out ) {
                std::cerr << "Cannot open the File : " << filepath << std::endl;
                return false;
            }
            // Write Settings to configfile:
            out << "# Sections or 'splitts' the commandline argument from its designated value\nCLA_SPLIT = '" << cla_splitter << "'" << std::endl;
            out << "\n# Width of helper function in characters\nCLA_WIDTH = '" << cla_width << "'" << std::endl;
            out << "\n# Single or dual hyphen reserved commands\nCLA_ADD_PARAMETER = '" << cla_add_param << "'\nCLA_REMOVE_PARAMETER = '" << cla_remove_param << "'\nCLA_EDIT_PARAMETER = '" << cla_edit_param << "'" << std::endl;
            out << "\n# All Configured Parameters:" << std::endl;
            // Write Datastructures
            for ( const auto& datastructure : cla_datastructures ) {
                out << "\nCLA_ARGUMENT = {\n\ttype = '" << datastructure.type << "'\n\tkey = ";
                for ( long unsigned int i = 0; i < datastructure.key.size(); i++ ) {
                    out << "'" << datastructure.key.at( i ) << "'";
                    if ( i < datastructure.key.size() - 1 ) out << " or ";
                }
                out << "\n\tdescription = '" << datastructure.description << "'";
                out << "\n\tgroup = '" << datastructure.group << "'";
                for ( const auto& parameter : datastructure.parameter ) {
                    out << "\n\tCLA_PARAMETER = {\n\t\tsubkey = '" << parameter.subkey << "'\n\t\tvalue = '" << parameter.value << "'\n\t\tdatatype = '" << parameter.datatype << "'\n\t\tdescription = '" << parameter.description << "'\n\t}";
                }
                out << "\n}";
            }
            out.close();
        } else {
            std::cerr << "Configfile size is 0!" << std::endl;
        }
        return true;
    }

    bool generateEmptyConfigfile( std::string filepath ) {
        std::cout << "\nGenerate Empty configfile? (j/n)" << std::endl;
        char in = 'n';
        std::cin >> in;
        if ( in == 'j' ) {
            std::ofstream out( ( filepath + ".settings.cla" ).c_str() );
            if ( !out ) {
                std::cerr << "Cannot open the File : " << filepath << std::endl;
                return false;
            }
            // Write Settings to configfile:
            out << "# Sections or 'splitts' the commandline argument from its designated value\nCLA_SPLIT = '='" << std::endl;
            out << "\n# Width of helper function in characters\nCLA_WIDTH = '50'" << std::endl;
            out << "\n# Single or dual hyphen reserved commands\nCLA_ADD_PARAMETER = 'add'\nCLA_REMOVE_PARAMETER = 'rem'\nCLA_EDIT_PARAMETER = 'mod'" << std::endl;
            out.close();
            return true;
        } else
            return false;
    }

    // Handles Add/Remove/Edit of parameters
    bool cla_wizzard_add( std::string filepath ) {
        std::string description,type,group,jn;
        std::vector<std::string> key;
        std::vector<Parameter> parameter;
        bool done = true;
        std::cout << "Enter type (-/--) ";
        std::string temp;
        std::getline(std::cin,temp);
        if (temp.at(0) == '-' && temp.size() > 1 && temp.at(1) == '-')
            type = "2hyphen";
        else type = "1hyphen";
        do {
            std::string temp;
            std::cout << "Enter a key:";
            std::getline(std::cin,temp);
            //TODO: Check if key already exists
            key.emplace_back(temp);
            std::cout << "Do you want to add another alias? [j/n] "; 
            std::getline(std::cin,jn);
            if (jn.size() > 0 && jn.at(0) == 'j') {done = false;} else done=true;
        } while (!done);
        std::cout << "Enter a groupname" << std::endl;
        std::getline(std::cin,group);
        std::cout << "Enter general description" << std::endl;
        std::getline(std::cin,description);
        if (type.compare("1hyphen") == 0 ) {
            parameter.emplace_back(param_bool);
        } else {
            done = true;
            std::string datatype,key,description,value;
            do {
                std::string temp;
                std::cout << "Adding Parameter " << (parameter.size()+1) << "\nEnter a datatype: ";
                std::getline(std::cin,datatype);
                std::cout << "Enter standard value: ";
                std::getline(std::cin,value);
                std::cout << "Enter identification key: ";
                std::getline(std::cin,key);
                std::cout << "Enter description: ";
                std::getline(std::cin,description);
                parameter.emplace_back(Parameter(key,value,datatype,description));
                std::cout << "Do you want to add another parameter? (" + std::to_string(parameter.size()) + " total now) [j/n]";
                std::getline(std::cin,jn);
                if (jn.size() > 0 && jn.at(0) == 'j') {done = false; } else done=true;
            } while (!done);
        }
        cla_datastructures.emplace_back(Datastructure(key,description,type,parameter,group));
        std::cout << "Save? [j/n] ";
        std::getline(std::cin,jn);
        if (jn.size() > 0 && jn.at(0) == 'j') {writeConfigfile( filepath ); std::cout << "New Structure was saved\n";} else {std::cout << "Structure dismissed.\n";}
        return true;
    }

    // TODO:
    // Modifies CLA Datastructure or Parameter member if subkey is not empty
    bool modify( std::string key, std::string subkey, std::string field, std::string newval ) {
        for ( auto& datastructure : cla_datastructures ) {
            if ( datastructure.validkey( key ) )
                // Modify Datastructure of subkey is empty
                if ( subkey.size() == 0 ) {
                } else {
                    for ( auto& param : datastructure.parameter ) {
                        if ( param.validkey( subkey ) ) {
                            param.value = newval;
                        }
                    }
                }
        }
        return true;
    }

   public:
    CommandlineArguments() {};
    CommandlineArguments( int argc, char* argv[], std::string filepath) {
        convert_argv( argc, argv );
        if ( !readConfigFile( filepath ) ) {
            std::cerr << "No Configfile read/generated, exitting..." << std::endl;
            exit( 0 );
        }
        readCLA_SettingsFromConfig();
        readCLA_ArgumentsFromConfig();
        parseCLA_Arguments( commandlinearguments, filepath );
    }
    CommandlineArguments( const std::vector<std::string>& argv, std::string filepath) {
        commandlinearguments = argv;
        if ( !readConfigFile( filepath ) ) {
            std::cerr << "No Configfile read/generated, exitting..." << std::endl;
            exit( 0 );
        }
        readCLA_SettingsFromConfig();
        readCLA_ArgumentsFromConfig();
        parseCLA_Arguments( commandlinearguments, filepath );
    }

    Parameter get( std::string key, std::string subkey = "", bool use_hyphens = true ) const {
        for ( const auto& datastructure : cla_datastructures ) {
            if ( datastructure.validkey( key, use_hyphens ) ) {
                if ( subkey.size() > 0 ) {
                    for ( const auto& param : datastructure.parameter ) {
                        if ( param.validkey( subkey ) ) {
                            return param;
                        }
                    }
                } else {
                    return datastructure.parameter.front();
                }
            }
        }
        std::cout << "Parameter with key " << key << "/" << subkey << " was not found!" << std::endl;
        return param_empty;
    }
    inline Parameter operator()( std::string key, std::string subkey = "" ) const { return get( key, subkey ); };
    std::vector<std::string> get_vector() const {
        return commandlinearguments;
    }
};

CommandlineArguments cla;

void init( int argc, char* argv[], std::string filepath = "" ) {
    cla = CommandlineArguments(argc, argv, filepath);
}
void init( const std::vector<std::string>& argv, std::string filepath = "" ) {
    cla = CommandlineArguments(argv, filepath);
}

} // namespace CommandlineArguments