#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

#include <algorithm>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>


class ConfigManager {
    private:
        std::map<std::string, std::string> config;
        ConfigManager() {}

        std::string trim(const std::string& str) {
            if (str.empty()) return "";
            
            // Find first non-whitespace, non-control character
            auto start = std::find_if(str.begin(), str.end(), [](unsigned char ch) {
                return !std::isspace(ch) && !std::iscntrl(ch);
            });
            
            // Find last non-whitespace, non-control character
            auto end = std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) {
                return !std::isspace(ch) && !std::iscntrl(ch);
            }).base();
            
            return (start < end) ? std::string(start, end) : "";
        }

    public:
        static ConfigManager& getInstance() {
            static ConfigManager instance;
            return instance;
        }

        void load(const std::string &filename) {
            std::ifstream file(filename);
            std::string line;

            if (!file.is_open()) {
                std::cerr << "🛑 ERROR: Could not open the config file at path: \"" << filename << "\"\n";
                std::cerr << "Please check if the file exists and if the relative path is correct.\n";
                return;
            }

            while (std::getline(file, line)) {
                // Skip comments and empty lines early
                if (line.empty() || line[0] == '#') continue;

                std::istringstream is_line(line);
                std::string key, value;
                
                if (std::getline(is_line, key, '=')) {
                    // Check if there's actually a value after the '='
                    if (std::getline(is_line, value)) {
                        
                        // Use the ultra-safe trimmer to strip spaces, \r, \n, and tabs
                        key = trim(key);
                        value = trim(value);

                        if (!key.empty()) {
                            config[key] = value;
                        }
                    }
                }
            }
        }

        std::string getType() {
            return config["T"];
        }

        // Returns the shape of tangent seeds
        size_t getTangentShape() {
            return std::stoull(config["V"]);
        }

        // Returns the shape of adjoint seeds
        size_t getAdjointShape() {
            return std::stoull(config["U"]);
        }

        // Returns the shape of inputs x
        size_t getXShape() {
            return std::stoull(config["x"]);
        }

        // Returns the shape of outputs y
        size_t getYShape() {
            return std::stoull(config["y"]);
        }

        // Returns the name of the primal function
        std::string getPrimalFunctionName() {
            return config["f"];
        }

        // Returns the AD sequence, reversed so that order label 1 corresponds to
        // the LAST character of the configured string, matching "X over Y" AD
        // terminology: e.g. "ta" ("tangent over adjoint") means adjoint is order 1
        // (innermost, first-applied) and tangent is order 2 (outer, wraps it).
        std::string getSequence() {
            std::string sequence = config["sequence"];
            std::reverse(sequence.begin(), sequence.end());
            return sequence;
        }
};

#endif