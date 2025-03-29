#include <cstdio>
#include <iostream>
#include <string>
#include "utils.hpp"

int execute_command(const std::string& cmd, const std::string& filename) {
    // Redirect stdout and stderr to the file
    std::string full_cmd = cmd + " > " + filename + " 2>&1";

    // Use system() to execute the full shell command (supports pipes)
    int status = system(full_cmd.c_str());

    // Check exit status
    if (WIFEXITED(status)) {
        return WEXITSTATUS(status);  // Extract exit code on Linux
    }
    else {
        return -1;
    }
}