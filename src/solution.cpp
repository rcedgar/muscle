#include "myutils.h"
#include <string>  // Make sure to include this

void cmd_solution() {
    // Retrieve the option following -solution
    std::string solution_option = opt(solution);  // If opt() returns std::string, store it correctly

    // Print the retrieved option using ProgressLog
    ProgressLog("Solution option: %s", solution_option.c_str());  // Convert to const char* using .c_str()
}
