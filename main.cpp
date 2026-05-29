#include <iostream>
#include <random>
#include <variant>
#include "generator/generator.h"
#include "generator/configManager.h"
#include "generator/structures.h"
#include "generator/utils.h"
#include "generator/readWrite.h"


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


// Helper function to get the user's input
template<typename T>
X_t<T> getManualInput(size_t size) {
    X_t<T> input_data(size);
    std::cout << "Enter " << size << " values for the input vector:\n";
    for (size_t i = 0; i < size; ++i) {
        std::cout << "[" << i << "]: ";
        std::cin >> input_data[i];
    }
    return input_data;
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
    ConfigManager::getInstance().load("generator/configs.txt");
    ConfigManager cm = ConfigManager::getInstance();

    try {
        generateADHeader("generator/adDrivers.h");
        generateADDrivers("generator/adDrivers.cpp"); // Generates generator/driver_ad.cpp and generator/driver_formula.cpp
        std::cout << "[1/4] Generated AD drivers and header successfully.\n";
    } catch (const std::exception &e) {
        std::cerr << "Generation failed: " << e.what() << "\n";
        return 1;
    }


    // ----------------------------------------------------
    // PHASE 2: GENERATE AND WRITE PARAMETERS INTO A FILE
    // ----------------------------------------------------
    // 1. Determine how x_input is populated (Manual vs. Random)
    using ParameterContainer = std::variant<
        std::vector<Param<double>>,
        std::vector<Param<float>>,
        std::vector<Param<int>>
    >;
    ParameterContainer parameters;
    X_t<double> x_input; 
    size_t input_size = cm.getXShape(); // Define or fetch the required size of your input vector

    std::cout << "Choose input method:\n[1] Generate Random Vector\n[2] Manually Input Vector\nSelection: ";
    int choice;
    std::cin >> choice;

    if (choice == 2) {
        x_input = getManualInput<double>(input_size);
    } else {
        // Default or fallback: Initialize a random input
        x_input.resize(input_size);
        fillWithRandomValues(x_input);
    }

    // Determine what data type the input is and generate the seeds randomly
    std::string target_type = cm.getType();
    if (target_type == "double") {
        parameters = generateRandomSeeds<double>(cm.getSequence(), x_input);
        writeParameters<double>(std::get<std::vector<Param<double>>>(parameters), "generator/parameters.txt");
    } 
    else if (target_type == "float") {
        // If x_input was double, convert it to float for this branch
        X_t<float> x_input_float(x_input.begin(), x_input.end());
        
        parameters = generateRandomSeeds<float>(cm.getSequence(), x_input_float);
        writeParameters<float>(std::get<std::vector<Param<float>>>(parameters), "generator/parameters.txt");
    } 
    else if (target_type == "int") {
        // Convert x_input to int for this branch
        X_t<int> x_input_int(x_input.begin(), x_input.end());
        
        parameters = generateRandomSeeds<int>(cm.getSequence(), x_input_int);
        writeParameters<int>(std::get<std::vector<Param<int>>>(parameters), "generator/parameters.txt");
    } 
    else {
        std::cerr << "Error: Unsupported data type requested in config: " << target_type << "\n";
    }


    // ----------------------------------------------------
    // PHASE 3: COMPILE AD DRIVER ONLY
    // ----------------------------------------------------
    std::cout << "[2/4] Compiling dynamic AD driver...\n";
    
    // Compile only the freshly generated AD source file
    std::string compileCmd = "g++ -std=c++20 "
                         "generator/adDrivers.cpp "
                         "-I ../generator "
                         "-I ../include/ad "      // To find ad.h
                         "-I ../include "          // To find Eigen and nlohmann
                         "-o generator/adDrivers.exe";

    int compileStatus = std::system(compileCmd.c_str());

    if (compileStatus != 0) {
        std::cerr << "Compilation of generated AD driver failed!\n";
        return 2;
    }
    std::cout << "      Compilation successful.\n";


    // ----------------------------------------------------
    // PHASE 4: RUN BOTH IN PARALLEL
    // ----------------------------------------------------
    


    // ----------------------------------------------------
    // PHASE 5: VALIDATION
    // ----------------------------------------------------
}
