#include "generator.h"
#include "configManager.h"
#include <string>
#include <format>

// Generate the parameters needed for a given order and sequence of AD
std::vector<Param> generateParameters(int order, std::string sequence) {
    std::vector<Param> result = {};
    char curMode = sequence.front();

    std::vector<Param> recResult = generateParameters(--order, sequence.substr(1));

    if (curMode == 't') {
        // Append the new parameters
    } else if (curMode == 'a') {
        // Append the new parameters
    }

    return result;
}


// Generate the complex nested AD types
std::string generateNestedADType(int order, std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string V = std::to_string(cm.getTangentDim());
    std::string U = std::to_string(cm.getAdjointDim());

    std::string result = "T";

    for (int i = sequence.length() - 1; i >= 0; i--) {
        if (sequence[i] == 't') {
            result = std::format("T_t<{},{}>", result, V);
        } else if (sequence[i] == 'a') {
            result = std::format("A_t<{},{}>", result, U);
        }
    }

    return result;
}


// Generate the x variable for complex nested AD types
std::string generateXNestedAdType(int order, std::string sequence) {
    std::string complexType = generateXNestedAdType(order, sequence);
    return std::format("X_t<{}>", complexType);
}

// Generate the y variable for complex nested AD types
std::string generateYNestedAdType(int order, std::string sequence) {
        std::string complexType = generateXNestedAdType(order, sequence);
    return std::format("Y_t<{}>", complexType);
}