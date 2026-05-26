#include "configManager.h"
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <iostream>


// Retrieve values defined in the config file
void ConfigManager::load(const std::string &filename) {
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "🛑 ERROR: Could not open the config file at path: \"" << filename << "\"\n";
        std::cerr << "Please check if the file exists and if the relative path is correct.\n";
        return;
    }

    while (std::getline(file, line)) {
        // Remove Windows carriage return '\r' if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            // Trim leading/trailing whitespace from key
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);

            // Trim leading/trailing whitespace from value
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

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


// Returns the AD sequence
std::string ConfigManager::getSequence() {
    return config["sequence"];
}