#include "global.h"

std::string PREFIX_PERCENT = "@#PERCENT#@";
std::string PREFIX_PERCENT_TIME = "@#PERCENTTIME#@";
std::string PREFIX_PERCENT_TIME_FINAL = "@#PERCENTTIMEFINAL#@";
std::string PREFIX_WARNING = "@#WARNING#@";
std::string PREFIX_ERROR = "@#ERROR#@";
std::string PREFIX_DEBUG = "@#DEBUG#@";
std::string PREFIX_OUTPUT = "@#OUTPUT#@";
std::string PREFIX_SUFFIX = "@#SUFFIX#@";

std::string global_message_normaltermination = "[Programm Terminated normally]";
std::string global_message_error_divergent = "[FATAL ERROR: PROCESS DIVERGED, PROGRAM TERMINATED]";
std::string global_message_error_wrong_number_input = "[FATAL ERROR: WRONG NUMBER OF INPUT PARAMETERS, PROGRAM TERMINATED]";

bool Save_State_sort_t(const SaveStateTau &ss1, const SaveStateTau &ss2) {
    return (ss1.t < ss2.t);
}
bool Save_State_sort_tau(const SaveStateTau &ss1, const SaveStateTau &ss2) {
    return (ss1.tau < ss2.tau);
}

//Test:
Log logs = Log();