#include "generator.h"
#include "configManager.h"
#include <string>
#include <format>


// Generate the parameters needed for a given order and sequence of AD
std::vector<Param> generateParameters(int order, std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t xShape = cm.getXShape();
    size_t yShape = cm.getYShape();
    size_t V = cm.getTangentShape();
    size_t U = cm.getAdjointShape();

    // Base case
    if (order == 0 || sequence.empty()) {
        std::vector<Param> base;
        // The original input X
        base.emplace_back(std::deque<size_t>{xShape}, "X", ParamRole::Input);
        // The original output Y
        base.emplace_back(std::deque<size_t>{yShape}, "Y", ParamRole::Output);
        return base;
    }

    char curMode = sequence[0];

    // Recursive result
    std::vector<Param> result = generateParameters(order - 1, sequence.substr(1));

    size_t currentSize = result.size();
    for (size_t i = 0; i < currentSize; i++) {
        Param p = result[i];

        std::deque<size_t> newShapes = p.tensor.shape;
        ParamRole newRole = p.role;
        std::string newName = std::format("{}_{}", p.name, std::to_string(order));

        if (curMode == 't') {
            newShapes.push_back(V);
        } else if (curMode == 'a') {
            newShapes.push_front(U);
            newRole = (p.role == ParamRole::Input) ? ParamRole::Output : ParamRole::Input;
        }

        result.emplace_back(newShapes, newName, newRole);
    }

    return result;
}


// Generate the complex nested AD types
std::string generateNestedADType(int order, std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string V = std::to_string(cm.getTangentShape());
    std::string U = std::to_string(cm.getAdjointShape());

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