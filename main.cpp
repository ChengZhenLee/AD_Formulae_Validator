#include <iostream>
#include <random>
#include <variant>
#include "generator/generator.h"
#include "generator/configManager.h"
#include "generator/structures.h"
#include "generator/utils.h"
#include "generator/readWrite.h"
#include "generator/formulaDriver.hpp"


// Helper function to run a system command and capture its text output
std::string runCommand(const std::string& cmd) {
    std::string result;
    char buffer[128];
    // Use _popen on Windows, popen on Linux/macOS
    #ifdef _WIN32
        FILE* pipe = _popen(cmd.c_str(), "r");
    #else
        FILE* pipe = popen(cmd.c_str(), "r");
    #endif

    if (!pipe) return "ERROR";
    
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        result += buffer;
    }

    #ifdef _WIN32
        _pclose(pipe);
    #else
        pclose(pipe);
    #endif
        return result;
}


// Helper function to fill input with random values
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


int main(int argc, char** argv) {
    std::cout << "=== AD Validator Pipeline Started ===\n";

    // ----------------------------------------------------
    // PHASE 1: LOAD CONFIG & GENERATE CODE
    // ----------------------------------------------------
    std::cout << "\n[1/6] Generating AD drivers...\n";
    ConfigManager::getInstance().load("generator/configs.txt");
    ConfigManager cm = ConfigManager::getInstance();

    try {
        generateADHeader("generator/adDrivers.h");
        generateADDrivers("generator/adDrivers.cpp"); // Generates generator/driver_ad.cpp and generator/driver_formula.cpp
    } catch (const std::exception &e) {
        std::cerr << "Generation failed: " << e.what() << "\n";
        return 1;
    }
    std::cout << "Succesfully generated AD drivers\n";

    std::cout << "\n[1.5/6] Generating and compiling helper drivers for derivatives...\n";
    try {
        generateHelperHeader("generator/adHelper.h");
        generateHelperDrivers("generator/adHelper.cpp");
    } catch (const std::exception &e) {
        std::cerr << "Helper driver generation failed: " << e.what() << "\n";
        return 2;
    }
    std::cout << "Succesfully generated helper drivers\n";


    // ----------------------------------------------------
    // PHASE 2: GENERATE AND WRITE PARAMETERS INTO A FILE
    // ----------------------------------------------------
    std::cout << "\n[2/6] Generating Input Seeds...\n";

    // 1. Determine how x_input is populated (Manual vs. Random)
    using ParameterContainer = std::variant<
        std::vector<Param<double>>,
        std::vector<Param<float>>,
        std::vector<Param<int>>
    >;
    ParameterContainer parameters;
    X_t<double> x_input; 
    size_t input_size = cm.getXShape(); // Define or fetch the required size of your input vector

    // Default or fallback: Initialize a random input
    x_input.resize(input_size);
    fillWithRandomValues(x_input);

    // Determine what data type the input is and generate the seeds randomly
    std::string target_type = cm.getType();
    if (target_type == "double") {
        parameters = generateRandomSeeds<double>(cm.getSequence(), x_input);
        writeParameters<double>(std::get<std::vector<Param<double>>>(parameters), "generator/parameters.bin");
    } 
    else if (target_type == "float") {
        // If x_input was double, convert it to float for this branch
        X_t<float> x_input_float(x_input.begin(), x_input.end());
        
        parameters = generateRandomSeeds<float>(cm.getSequence(), x_input_float);
        writeParameters<float>(std::get<std::vector<Param<float>>>(parameters), "generator/parameters.bin");
    } 
    else if (target_type == "int") {
        // Convert x_input to int for this branch
        X_t<int> x_input_int(x_input.begin(), x_input.end());
        
        parameters = generateRandomSeeds<int>(cm.getSequence(), x_input_int);
        writeParameters<int>(std::get<std::vector<Param<int>>>(parameters), "generator/parameters.bin");
    } 
    else {
        std::cerr << "Error: Unsupported data type requested in config: " << target_type << "\n";
    }


    // Generate derivative seeds: create parameters with X and X_k (with identity seeding)
    std::cout << "\n[2.5/6] Generating derivative seeds...\n";
    
    if (target_type == "double") {
        auto derivSeeds = generateDerivativeSeeds<double>(cm.getSequence(), x_input);
        writeParameters<double>(derivSeeds, "generator/derivatives_seeds.bin");
    }
    else if (target_type == "float") {
        X_t<float> x_input_float(x_input.begin(), x_input.end());
        auto derivSeeds = generateRandomSeeds<float>(cm.getSequence(), x_input_float);
        writeParameters<float>(derivSeeds, "generator/derivatives_seeds.bin");
    }
    else if (target_type == "int") {
        X_t<int> x_input_int(x_input.begin(), x_input.end());
        auto derivSeeds = generateRandomSeeds<int>(cm.getSequence(), x_input_int);
        writeParameters<int>(derivSeeds, "generator/derivatives_seeds.bin");
    }


    // ----------------------------------------------------
    // PHASE 3: COMPILE AD DRIVER ONLY
    // ----------------------------------------------------
    // std::cout << "\n[3/6] Compiling dynamic AD driver...\n";
    
    // // Compile only the freshly generated AD source file
    // std::string compileCmd = "g++ -std=c++20 -g -O0 "
    //                      "generator/adDrivers.cpp "
    //                      "-I ../generator "
    //                      "-I ../include/ad "      // To find ad.h
    //                      "-I ../include "          // To find Eigen and nlohmann
    //                      "-o generator/adDrivers.exe";

    // int compileStatus = std::system(compileCmd.c_str());

    // if (compileStatus != 0) {
    //     std::cerr << "Compilation of generated AD driver failed with return code: " << compileStatus << "\n";
    //     return 2;
    // }
    // std::cout << "      AD drivers compilation successful.\n";

    // // Compile the helper driver
    // std::cout << "\n[3.5/6] Compiling helper driver...\n";
    // std::string helperCompileCmd = "g++ -std=c++20 -g -O0 "
    //                                "generator/adHelper.cpp "
    //                                "-I ../generator "
    //                                "-I ../include/ad "
    //                                "-I ../include "
    //                                "-o generator/adHelper.exe";

    // int helperCompileStatus = std::system(helperCompileCmd.c_str());
    // if (helperCompileStatus != 0) {
    //     std::cerr << "Compilation of helper driver failed with return code: " << helperCompileStatus << "\n";
    //     return 2;
    // }
    // std::cout << "      Helper driver compilation successful.\n";


    // // ----------------------------------------------------
    // // PHASE 4: RUN THE DRIVERS
    // // ----------------------------------------------------

    // // Run the helper driver to compute derivatives
    // std::cout << "\n[4/6] Running helper driver to compute derivatives...\n";
    // int helperRunStatus = std::system(".\\generator\\adHelper.exe");
    // if (helperRunStatus != 0) {
    //     std::cerr << "Execution of helper driver failed with return code: " << helperRunStatus <<"\n";
    //     return 3;
    // }
    // std::cout << "      Helper driver executed successfully. Derivatives computed.\n";


    // //Run the compiled AD driver executable
    // std::cout << "\n[4.5/6] Running compiled AD driver...\n";
    // int runStatus = std::system(".\\generator\\adDrivers.exe");
    // if (runStatus != 0) {
    //     std::cerr << "Execution of compiled AD driver failed!\n";
    //     return 3;
    // }
    // std::cout << "      AD driver executed successfully.\n";


    // Run the formula driver on the same parameters
    std::cout << "[5/6] Running formula driver...\n";
    try {
        if (target_type == "double") {
            runFormulaDriver<double>(std::deque<Param<double>>(
                std::get<std::vector<Param<double>>>(parameters).begin(),
                std::get<std::vector<Param<double>>>(parameters).end()
            ));
        } else if (target_type == "float") {
            runFormulaDriver<float>(std::deque<Param<float>>(
                std::get<std::vector<Param<float>>>(parameters).begin(),
                std::get<std::vector<Param<float>>>(parameters).end()
            ));
        } else if (target_type == "int") {
            runFormulaDriver<int>(std::deque<Param<int>>(
                std::get<std::vector<Param<int>>>(parameters).begin(),
                std::get<std::vector<Param<int>>>(parameters).end()
            ));
        }
        std::cout << "      Formula driver executed successfully.\n";
    } catch (const std::exception &e) {
        std::cerr << "Formula driver failed: " << e.what() << "\n";
        return 4;
    }
    


    // ----------------------------------------------------
    // PHASE 5: VALIDATION
    // ----------------------------------------------------
}
