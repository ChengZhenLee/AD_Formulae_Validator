#include "configManager.h"
#include "structures.h"
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
        base.emplace_back(std::map<size_t, size_t>{}, std::deque<size_t>{xShape}, "X", ParamRole::Input);
        // The original output Y
        base.emplace_back(std::map<size_t, size_t>{}, std::deque<size_t>{yShape}, "Y", ParamRole::Output);
        return base;
    }

    char curMode = sequence[0];

    // Recursive result
    std::vector<Param> result = generateParameters(order - 1, sequence.substr(1));

    size_t currentSize = result.size();
    for (size_t i = 0; i < currentSize; i++) {
        Param p = result[i];

        std::map<size_t, size_t> newOShape = p.orderedShape;
        std::deque<size_t> newShapes = p.tensor.shape;
        ParamRole newRole = p.role;
        std::string newName = std::format("{}_{}", p.name, std::to_string(order));

        if (curMode == 't') {
            newShapes.push_back(V);
            newOShape[order] = V;
        } else if (curMode == 'a') {
            newShapes.push_front(U);
            newOShape[order] = U;

            // In adjoint mode, swap inputs to outputs and vice versa
            newRole = (p.role == ParamRole::Input) ? ParamRole::Output : ParamRole::Input;
        }

        result.emplace_back(newOShape, newShapes, newName, newRole);
    }

    return result;
}


// Generate the complex nested AD types
std::string generateNestedADType(int order, std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string V = std::to_string(cm.getTangentShape());
    std::string U = std::to_string(cm.getAdjointShape());

    std::string result = cm.getType();

    for (int i = 0; i < sequence.length(); i++) {
        if (sequence[i] == 't') {
            result = std::format("T_t<{},{}>", result, V);
        } else if (sequence[i] == 'a') {
            result = std::format("A_t<{},{}>", result, U);
        }
    }

    return result;
}


// Generate the x variable for complex nested AD types
std::string generateXNestedADType(int order, std::string sequence) {
    std::string complexType = generateNestedADType(order, sequence);
    return std::format("X_t<{}>", complexType);
}


// Generate the y variable for complex nested AD types
std::string generateYNestedADType(int order, std::string sequence) {
        std::string complexType = generateNestedADType(order, sequence);
    return std::format("Y_t<{}>", complexType);
}


// Nested AD types must have a .value() member
template<typename T>
concept isADNestedType = requires(T x) {
    { x.value() };
};


// Seed the nested AD type
// Note: the initial call in the driver is done with empty coords
template<typename ADNested>
void seedAD(ADNested& x, Param& p, size_t curOrder, std::string& sequence, std::deque<size_t> coords) {
    // Check if it ADNested is a valid AD nested type at compile time
    if constexpr (!isADNestedType<ADNested>) {
        if (curOrder == 0) {
            x = p.tensor.data[p.tensor.getIndex(coords)];
            return;
        }
    }

    char curType = sequence[sequence.length() - curOrder];

    if (p.isActive(curOrder)) {
        size_t curShape = p.getShapeAt(curOrder);

        if (curType == 't') {
            for (size_t i = 0; i < curShape; i++) {
                coords.push_back(i);
                ADSeed(x.tangent(i), p, curOrder - 1, sequence, coords);
                coords.pop_back();
            }

        } else if (curType == 'a') {
            for (size_t i = 0; i < curShape; i++ ) {
                coords.push_front(i);
                ADSeed(x.adjoint(i), p, curOrder - 1, sequence, coords);
                coords.pop_back();
            }
        }
    } else {
        ADSeed(x.value(), p, curOrder - 1, sequence, coords);
    }
}


template<typename ADNested>
void extractAD(ADNested& y, Param& p, size_t curOrder, std::string& sequence, std::deque<size_t> coords) {
    // Check if it ADNested is a valid AD nested type at compile time
    if constexpr (!isADNestedType<ADNested>) {
        if (curOrder == 0) {
            p.tensor.data[p.tensor.getIndex(coords)] = y;
            return;
        }
    }


    char curType = sequence[sequence.length() - curOrder];

    if (p.isActive(curOrder)) {
        size_t curShape = p.getShapeAt(curOrder);

        if (curType == 't') {
            for (size_t i = 0; i < curShape; i++) {
                coords.push_back(i);
                ADSeed(y.tangent(i), p, curOrder - 1, sequence, coords);
                coords.pop_back();
            }

        } else if (curType == 'a') {
            for (size_t i = 0; i < curShape; i++ ) {
                coords.push_front(i);
                ADSeed(y.adjoint(i), p, curOrder - 1, sequence, coords);
                coords.pop_back();
            }
        }
    } else {
        ADSeed(y.value(), p, curOrder - 1, sequence, coords);
    }
}


// Seeds all parameters for a specific order/layer
template<typename ADNested>
void seedADForOrder(ADNested& x, std::vector<Param>& parameters, size_t curOrder, std::string& sequence) {
    for (Param p : parameters) {
        if (p.highestOrder == curOrder && p.isActive(curOrder)) {
            seedAD(x, p, curOrder, sequence, {});
        }
    }
}


// Extracts all values for a specific order/layer
template<typename ADNested>
void extractADForOrder(ADNested& y, std::vector<Param>& parameters, size_t curOrder, std::string& sequence) {
    for (Param p : parameters) {
        if (p.highestOrder == curOrder && p.isActive(curOrder)) {
            extractAD(y, p, curOrder, sequence, {});
        }
    }
}


template<typename ADNested>
void seedPrimal(ADNested& x, std::vector<Param> parameters, std::deque<size_t> coords) {
    if constexpr (!isADNestedType<ADNested>) {
        for (Param p : parameters) {
            if (p.name.length() == 1 && p.name[0] == 'X') {
                x = p.tensor.data[p.tensor.getIndex(coords)];
            }
        }
    }

    seedPrimal(x.value(), parameters, coords);
}


template<typename ADNested>
void extractPrimal(ADNested& y, std::vector<Param> parameters, std::deque<size_t> coords) {
    if constexpr (!isADNestedType<ADNested>) {
        for (Param p : parameters) {
            if (p.name.length() == 1 && p.name[0] == 'Y') {
                p.tensor.data[p.tensor.getIndex(coords)] = y;
            }
        }
    }
}


// Determine the current layer's AD type
std::string getCurrentLayerADType(std::string ADNested, size_t curOrder, std::string sequence) {
    // adjoint over tangent -> 'at' -> order of adjoint is 2 -> 'a' -> A_t<double, U>
    // adjoint over tangent -> 'at' -> order of tangent is 1 -> 'at' -> T_t<A_t<double, U>, V>
    std::string subsequence = sequence.substr(0, sequence.length() - curOrder + 1);
    return generateNestedADType(curOrder, subsequence);
}


// Determine the function name of the current layer
std::string getCurrentLayerFunctionName(size_t curOrder) {
    ConfigManager cm = ConfigManager::getInstance();
    if (curOrder == 0) return cm.getPrimalFunctionName();
    
    return std::format("AD_F_{}", curOrder);
}