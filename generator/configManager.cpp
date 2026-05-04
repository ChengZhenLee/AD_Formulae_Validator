#include "configManager.h"
#include <fstream>
#include <string>
#include <map>
#include <sstream>


// Retrieve values defined in the config file
void ConfigManager::load(const std::string &filename) {
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[key] = value;
        }
    }
}


// Returns the variable dimension size of tangent seeds
size_t ConfigManager::getTangentDim() {
    return std::stoull(config["V"]);
}


// Returns the variable dimension size of adjoint seeds
size_t ConfigManager::getAdjointDim() {
    return std::stoull(config["U"]);
}