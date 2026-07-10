#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

// Process-wide singleton holding the key=value pairs from
// generator/configs.txt (x, y, V, U, T, f, sequence), with an optional
// runtime override for `sequence` from the --sequence CLI flag. Every
// generator/validator translation unit reads its configuration through
// this single instance rather than passing config values around.
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
                std::cerr << "ERROR: Could not open the config file at path: \"" << filename << "\"\n";
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
                        
                        // Use the safe trimmer to strip spaces, \r, \n, and tabs
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
        // the FIRST character of the configured string. For example, "ta" means
        // "adjoint-over-tangent" mode. The reversed sequence is used for driver
        // generation, parameter generation, seeding and extracting.
        std::string getSequence() {
            std::string sequence = config["sequence"];
            std::reverse(sequence.begin(), sequence.end());
            return sequence;
        }

        // Returns the sequence exactly as configured (no reversal). Only for
        // display/logging, where showing the user's own typed string matters.
        std::string getRawSequence() {
            return config["sequence"];
        }

        // Overrides a config value at runtime, e.g. a --sequence CLI argument.
        void set(const std::string& key, const std::string& value) {
            config[key] = value;
        }

        // Rewrites the "sequence = ..." line in the config file on disk to
        // `sequence`, so a --sequence CLI override becomes the new default
        // for future runs that don't pass one, instead of reverting to
        // whatever configs.txt said before.
        void persistSequence(const std::string& filename, const std::string& sequence) {
            std::ifstream inFile(filename);
            if (!inFile.is_open()) return;

            std::vector<std::string> lines;
            std::string line;
            bool found = false;

            while (std::getline(inFile, line)) {
                std::istringstream is_line(line);
                std::string key;
                if (!found && std::getline(is_line, key, '=') && trim(key) == "sequence") {
                    lines.push_back("sequence = " + sequence);
                    found = true;
                } else {
                    lines.push_back(line);
                }
            }
            inFile.close();

            if (!found) {
                lines.push_back("sequence = " + sequence);
            }

            std::ofstream outFile(filename);
            if (!outFile.is_open()) return;
            for (const auto& l : lines) {
                outFile << l << "\n";
            }
        }
};

#endif