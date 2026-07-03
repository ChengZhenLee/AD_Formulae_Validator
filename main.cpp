#include <iostream>
#include <random>
#include <variant>
#include <fstream>
#include <filesystem>
#include "generator/generator.h"
#include "generator/configManager.h"
#include "generator/structures.h"
#include "generator/readWrite.h"
#include "generator/formulaDriver.hpp"
#include "generator/validator.h"

namespace fs = std::filesystem;

using ParameterContainer = std::variant<
    std::vector<Param<double>>,
    std::vector<Param<float>>,
    std::vector<Param<int>>
>;

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

template<typename T>
void fillWithRandomValues(X_t<T>& x_input) {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    using DistributionType = std::conditional_t<
        std::is_floating_point_v<T>, 
        std::uniform_real_distribution<T>, 
        std::uniform_int_distribution<T>
    >;
    DistributionType dist(static_cast<T>(0), static_cast<T>(100));
    
    for (auto& element : x_input) {
        element = dist(gen);
    }
}

// Struct to track our build states
struct PipelineCache {
    std::string last_sequence = "";
    size_t last_length = 0;

    void load(const std::string& filepath) {
        std::ifstream f(filepath);
        if (f.is_open()) {
            std::getline(f, last_sequence);
            f >> last_length;
        }
    }

    void save(const std::string& filepath) const {
        std::ofstream f(filepath);
        if (f.is_open()) {
            f << last_sequence << "\n" << last_length;
        }
    }
};

// ============================================================================
// PIPELINE STAGES
// ============================================================================

bool handleADDriverCompilation(const std::string& current_sequence, const PipelineCache& cache) {
    bool needs_rebuild = (current_sequence != cache.last_sequence) || !fs::exists("generator/adDrivers.exe");
    
    if (!needs_rebuild) {
        std::cout << "--> Sequence unchanged. Skipping AD driver generation & compilation.\n";
        return true;
    }

    std::cout << "\n[1/6] Generating & Compiling AD drivers (Sequence Changed)...\n";
    try {
        generateADHeader("generator/adDrivers.h");
        generateADDrivers("generator/adDrivers.cpp");
    } catch (const std::exception &e) {
        std::cerr << "Generation failed: " << e.what() << "\n";
        return false;
    }

    std::string compileCmd = "g++ -std=c++20 -g -O0 generator/adDrivers.cpp "
                             "-I ../generator -I ../include/ad -I ../include -o generator/adDrivers.exe";
    if (std::system(compileCmd.c_str()) != 0) {
        std::cerr << "Compilation of AD driver failed.\n";
        return false;
    }
    return true;
}

bool handleHelperDriverCompilation(size_t current_length, const PipelineCache& cache) {
    bool needs_rebuild = (current_length != cache.last_length) || !fs::exists("generator/adHelper.exe");

    if (!needs_rebuild) {
        std::cout << "--> Sequence length unchanged. Skipping Helper driver generation & compilation.\n";
        return true;
    }

    std::cout << "\n[1.5/6] Generating & Compiling helper drivers (Length Changed)...\n";
    try {
        generateHelperHeader("generator/adHelper.h");
        generateHelperDrivers("generator/adHelper.cpp");
    } catch (const std::exception &e) {
        std::cerr << "Helper driver generation failed: " << e.what() << "\n";
        return false;
    }

    std::string helperCompileCmd = "g++ -std=c++20 -g -O0 generator/adHelper.cpp "
                                   "-I ../generator -I ../include/ad -I ../include -o generator/adHelper.exe";
    if (std::system(helperCompileCmd.c_str()) != 0) {
        std::cerr << "Compilation of helper driver failed.\n";
        return false;
    }
    return true;
}

ParameterContainer generateSeeds(const std::string& target_type, ConfigManager& cm, X_t<double>& x_input) {
    std::cout << "\n[2/6] Generating Input & Derivative Seeds...\n";
    ParameterContainer parameters;

    if (target_type == "double") {
        parameters = generateRandomSeeds<double>(cm.getSequence(), x_input);
        writeParameters<double>(std::get<std::vector<Param<double>>>(parameters), "generator/parameters.bin");
        
        auto derivSeeds = generateDerivativeSeeds<double>(cm.getSequence(), x_input);
        writeParameters<double>(derivSeeds, "generator/derivatives_seeds.bin");
    } 
    else if (target_type == "float") {
        X_t<float> x_input_float(x_input.begin(), x_input.end());
        parameters = generateRandomSeeds<float>(cm.getSequence(), x_input_float);
        writeParameters<float>(std::get<std::vector<Param<float>>>(parameters), "generator/parameters.bin");
        
        auto derivSeeds = generateRandomSeeds<float>(cm.getSequence(), x_input_float);
        writeParameters<float>(derivSeeds, "generator/derivatives_seeds.bin");
    } 
    else if (target_type == "int") {
        X_t<int> x_input_int(x_input.begin(), x_input.end());
        parameters = generateRandomSeeds<int>(cm.getSequence(), x_input_int);
        writeParameters<int>(std::get<std::vector<Param<int>>>(parameters), "generator/parameters.bin");
        
        auto derivSeeds = generateRandomSeeds<int>(cm.getSequence(), x_input_int);
        writeParameters<int>(derivSeeds, "generator/derivatives_seeds.bin");
    }
    return parameters;
}

bool runExecutables() {
    std::cout << "\n[4/6] Running compiled drivers...\n";
    
    // Windows vs Posix compatibility handling for execution paths if required, 
    // keeping your original windows dot-slash format here:
    if (std::system(".\\generator\\adHelper.exe") != 0) {
        std::cerr << "Execution of helper driver failed!\n";
        return false;
    }
    if (std::system(".\\generator\\adDrivers.exe") != 0) {
        std::cerr << "Execution of AD driver failed!\n";
        return false;
    }
    return true;
}

ParameterContainer runFormula(const std::string& target_type, const ParameterContainer& parameters) {
    std::cout << "\n[5/6] Running formula driver...\n";
    if (target_type == "double") {
        return runFormulaDriver<double>(std::get<std::vector<Param<double>>>(parameters));
    } else if (target_type == "float") {
        return runFormulaDriver<float>(std::get<std::vector<Param<float>>>(parameters));
    } else {
        return runFormulaDriver<int>(std::get<std::vector<Param<int>>>(parameters));
    }
}

void validate(const std::string& target_type, const ParameterContainer& formulaResults) {
    std::cout << "\n[6/6] Validating results...\n";
    if (target_type == "double") {
        validateParameters<double>(std::get<std::vector<Param<double>>>(formulaResults), "generator/results.bin", "generator/validation_results.txt");
    } else if (target_type == "float") {
        validateParameters<float>(std::get<std::vector<Param<float>>>(formulaResults), "generator/results.bin", "generator/validation_results.txt");
    } else if (target_type == "int") {
        validateParameters<int>(std::get<std::vector<Param<int>>>(formulaResults), "generator/results.bin", "generator/validation_results.txt");
    }
    std::cout << "Validator executed successfully.\n";
}

// ============================================================================
// MAIN PIPELINE EXECUTION
// ============================================================================

int main(int argc, char** argv) {
    std::cout << "=== AD Validator Pipeline Started ===\n";

    ConfigManager::getInstance().load("generator/configs.txt");
    ConfigManager cm = ConfigManager::getInstance();
    
    std::string current_sequence = cm.getSequence();
    size_t current_length = current_sequence.length();
    std::string target_type = cm.getType();

    // Load caching information
    PipelineCache cache;
    const std::string cache_file = "generator/.pipeline_cache";
    cache.load(cache_file);

    // Phase 1: Conditional Compilation
    if (!handleADDriverCompilation(current_sequence, cache)) return 1;
    if (!handleHelperDriverCompilation(current_length, cache)) return 2;

    // Update and save cache right away if compilations succeeded
    cache.last_sequence = current_sequence;
    cache.last_length = current_length;
    cache.save(cache_file);

    // Phase 2: Input Generation
    X_t<double> x_input(cm.getXShape());
    fillWithRandomValues(x_input);
    
    ParameterContainer parameters;
    try {
        parameters = generateSeeds(target_type, cm, x_input);
    } catch (const std::exception &e) {
        std::cerr << "Seed generation failed: " << e.what() << "\n";
        return 3;
    }

    // Phase 3 & 4: Execute drivers
    if (!runExecutables()) return 4;

    // Phase 5 & 6: Formula & Validation
    try {
        ParameterContainer formulaResults = runFormula(target_type, parameters);
        validate(target_type, formulaResults);
    } catch (const std::exception &e) {
        std::cerr << "Pipeline execution/validation failure: " << e.what() << "\n";
        return 5;
    }

    return 0;
}