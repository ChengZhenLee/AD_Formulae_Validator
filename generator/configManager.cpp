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


std::string ConfigManager::getType() {
    return config["T"];
}


// Returns the shape of tangent seeds
size_t ConfigManager::getTangentShape() {
    return std::stoull(config["V"]);
}


// Returns the shape of adjoint seeds
size_t ConfigManager::getAdjointShape() {
    return std::stoull(config["U"]);
}


// Returns the shape of inputs x
size_t ConfigManager::getXShape() {
    return std::stoull(config["x"]);
}


// Returns the shape of outputs y
size_t ConfigManager::getYShape() {
    return std::stoull(config["y"]);
}


// Returns the name of the primal function
std::string ConfigManager::getPrimalFunctionName() {
    return config["f"];
}