#ifndef UTILS_HPP
#define UTILS_HPP

#include "structures.h"
#include "configManager.h"
#include "generator.h"
#include "readWrite.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <format>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <type_traits>
#include "user_function.h"

// Nested AD types must have a .value() member
template<typename T, typename = void>
struct has_value_member : std::false_type {};

template<typename T>
struct has_value_member<T, std::void_t<decltype(std::declval<T>().value())>> : std::true_type {};

template<typename T>
inline constexpr bool isADNestedType = has_value_member<T>::value;

// Generate the parameters needed for a given order and sequence of AD
template <typename T>
std::vector<Param<T>> generateParameters(std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t xShape = cm.getXShape();
    size_t yShape = cm.getYShape();
    size_t V = cm.getTangentShape();
    size_t U = cm.getAdjointShape();
    size_t order = sequence.size();

    // Base case
    if (order == 0 || sequence.empty()) {
        std::vector<Param<T>> base;
        base.emplace_back(std::map<size_t, size_t>{}, std::deque<size_t>{xShape}, "X", ParamRole::Input);
        base.emplace_back(std::map<size_t, size_t>{}, std::deque<size_t>{yShape}, "Y", ParamRole::Output);
        return base;
    }

    char curMode = sequence[order - 1];
    std::vector<Param<T>> result = generateParameters<T>(sequence.substr(0, order - 1));

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
            newRole = (p.role == ParamRole::Input) ? ParamRole::Output : ParamRole::Input;
        }

        result.emplace_back(newOShape, newShapes, newName, newRole);
    }

    return result;
}

// Initialize the seed of every parameter with random data
template<typename T>
std::vector<Param<T>> generateRandomSeeds(std::string sequence, const X_t<T>& x) {
    std::vector<Param<T>> parameters = generateParameters<T>(sequence);
    std::random_device rd;
    std::mt19937 gen(rd());

    using DistributionType = std::conditional_t<
        std::is_floating_point_v<T>,
        std::uniform_real_distribution<T>,
        std::uniform_int_distribution<T>
    >;

    DistributionType dist(
        static_cast<T>(*std::min_element(x.begin(), x.end())),
        static_cast<T>(*std::max_element(x.begin(), x.end()))
    );

    for (Param<T>& p : parameters) {
        if (p.name == "X") {
            p.tensor.data.resize(x.size());
            for (int i = 0; i < x.size(); i++) {
                p.tensor.data[i] = x[i];
            }
        } else if (p.name.rfind("X", 0) == 0 && p.name.length() > 1) {
            for (int i = 0; i < p.tensor.data.size(); i++) {
                p.tensor.data[i] = dist(gen);
            }
        }
    }

    return parameters;
}


// TODO: seed X_1, X_2... with identity matrices
// Change the name of the appropriate output tensors to F_1_2...
template<typename T>
std::vector<Param<T>> generateDerivativeSeeds(std::string sequence, const X_t<T>& x) {
    size_t order = sequence.length();
    std::string derivativeSequence = "";
    for (size_t i = 0; i < order; i++) {
        derivativeSequence += 't';
    }
    std::vector<Param<T>> parameters = generateParameters<T>(derivativeSequence);

    for (Param<T>& p : parameters) {
        // Seed X_1, X_2, X_3... with identity matrices
        if (p.name.rfind("X", 0) == 0 && p.name.length() == 3) {
            // Seed with identity
            size_t rows = p.tensor.shape[0];
            size_t cols = p.tensor.shape.size() > 1 ? p.tensor.shape[1] : rows;

            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    size_t idx = p.tensor.getIndex({i, j});
                    p.tensor.data[idx] = (i == j) ? static_cast<T>(1) :static_cast<T>(0);
                }
            }
        }

        // It's easier to change all the "Y" to "F"
        // But be careful not to change the original y=f(x)
        if (p.name.rfind("Y", 0) == 0 && p.name.length() > 1) {
            p.name[0] = 'F';
        }
    }

    return parameters;
}


// TODO: get the derivatives from generator/derivatives.txt
// the derivatives names all start with F
template<typename T>
std::vector<Param<T>> getDerivatives(std::string sequence, const std::vector<T>& x0, T h) {
    // Read back the results and return the computed derivative parameter values.
    std::vector<Param<T>> results = readParameters<T>("generator/derivatives.txt");
    std::vector<Param<T>> derivatives;
    for (auto& p : results) {
        if (p.name.rfind("F", 0)) {
            derivatives.push_back(p);
        }
    }

    return derivatives;
}

// Seed the nested AD type
template<typename ADNested, typename T>
void seedAD(ADNested& x, const Param<T>& p, size_t curOrder, const std::string& sequence, std::deque<size_t> coords) {
    if constexpr(!isADNestedType<ADNested>) {
        if (curOrder > sequence.length()) {
            x = p.tensor.data[p.tensor.getIndex(coords)];
            return;
        }
    }

    char curType = sequence[sequence.length() - curOrder];
    if (p.isActive(curOrder)) {
        size_t curShape = p.getShapeAt(curOrder);
        if (curType == 't') {
            // COMPILER GUARD: Only compile this loop if x actually has a .tangent() method
            if constexpr (requires { x.tangent(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    std::deque<size_t> nextCoords = coords;
                    nextCoords.push_back(i);
                    seedAD(x.tangent(i), p, curOrder + 1, sequence, nextCoords);
                }
            }
        } else if (curType == 'a') {
            // COMPILER GUARD: Only compile this loop if x actually has an .adjoint() method
            if constexpr (requires { x.adjoint(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    std::deque<size_t> nextCoords = coords;
                    nextCoords.push_front(i);
                    seedAD(x.adjoint(i), p, curOrder + 1, sequence, nextCoords);
                }
            }
        }
    } else {
        if constexpr(requires { x.value(); })
        seedAD(x.value(), p, curOrder + 1, sequence, coords);
    }
}

// Extract the nested AD type
template<typename ADNested, typename T>
void extractAD(ADNested& y, Param<T>& p, size_t curOrder, const std::string& sequence, std::deque<size_t> coords) {
    if constexpr(!isADNestedType<ADNested>) {
        if (curOrder > sequence.length()) {
            p.tensor.data[p.tensor.getIndex(coords)] = y;
            return;
        }
    }

    char curType = sequence[sequence.length() - curOrder];
    if (p.isActive(curOrder)) {
        size_t curShape = p.getShapeAt(curOrder);
        
        if (curType == 't') {
            // COMPILER GUARD: Only compile this loop if y actually has a .tangent() method
            if constexpr(requires { y.tangent(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    std::deque<size_t> nextCoords = coords;
                    nextCoords.push_back(i);
                    extractAD(y.tangent(i), p, curOrder + 1, sequence, nextCoords);
                }
            }
            
        } else if (curType == 'a') {
            // COMPILER GUARD: Only compile this loop if y actually has a .adjoint() method
            if constexpr(requires { y.adjoint(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    std::deque<size_t> nextCoords = coords;
                    nextCoords.push_front(i);
                    extractAD(y.adjoint(i), p, curOrder + 1, sequence, nextCoords);
                }
            }
        }
    } else {
        if constexpr(requires { y.value(); })
            extractAD(y.value(), p, curOrder + 1, sequence, coords);
    }
}

// Seeds all parameters for a specific order/layer
template<typename ADNested, typename T>
void seedADForOrder(ADNested& x, std::vector<Param<T>>& parameters, size_t targetOrder, const std::string& sequence) {
    ConfigManager cm = ConfigManager::getInstance();

    for (Param<T>& p : parameters) {
        if (p.highestOrder == targetOrder && p.isActive(targetOrder) && p.role == ParamRole::Input) {
            for (size_t i = 0; i < cm.getXShape(); i++) {
                seedAD(x[i], p, 1, sequence, {i});
            }
        }
    }
}

// Extracts all values for a specific order/layer
template<typename ADNested, typename T>
void extractADForOrder(ADNested& y, std::vector<Param<T>>& parameters, size_t targetOrder, const std::string& sequence) {
    ConfigManager cm = ConfigManager::getInstance();

    for (Param<T>& p : parameters) {
        if (p.highestOrder == targetOrder && p.isActive(targetOrder) && p.role == ParamRole::Output) {
            for (size_t j = 0; j < cm.getYShape(); j++) {
                extractAD(y[j], p, 1, sequence, {j});
            }
        }
    }
}

// Seed the primal value for the nested type
template<typename ADNested, typename T>
void seedPrimal(ADNested& x, const std::vector<Param<T>>& parameters) {
    const size_t xSize = x.size();

    for (const auto& p : parameters) {
        if (p.name.size() == 1 && p.name[0] == 'X') {
            for (size_t i = 0; i < xSize; ++i) {
                seedElement(x[i], p, i);
            }
        }
    }
}

// Extract the primal value for the nested type
template<typename ADNested, typename T>
void extractPrimal(ADNested& y, std::vector<Param<T>>& parameters) {
    const size_t ySize = y.size();

    for (Param<T>& p : parameters) {
        if (p.name.length() == 1 && p.name[0] == 'Y') {
            for (size_t j = 0; j < ySize; j++) {
                extractElement(y[j], p, j);
            }
        }
    }
}

// Helper functions for seeding and extracting
template<typename ADNested, typename T>
void seedElement(ADNested& element, const Param<T> &p, size_t index) {
    if constexpr(isADNestedType<ADNested>) {
        seedElement(element.value(), p, index);
    } else {
        element = p.tensor.data[p.tensor.getIndex({index})];
    }
}


template<typename ADNested, typename T>
void extractElement(ADNested& element, Param<T> &p, size_t index) {
    if constexpr(isADNestedType<ADNested>) {
        extractElement(element.value(), p, index);
    } else {
        p.tensor.data[p.tensor.getIndex({index})] = element;
    }
}


// Find a specific parameter by name from a list of parameters
template <typename T>
Param<T>& findParamByName(const std::string& targetName, std::deque<Param<T>>& parameters) {
    for (auto& p : parameters) {
        if (p.name == targetName) {
            return p;
        }
    }
    throw std::runtime_error(std::format("Parameter not found: {}", targetName));
}

template <typename T>
const Param<T>& findParamByName(const std::string& targetName, const std::deque<Param<T>>& parameters) {
    for (const auto& p : parameters) {
        if (p.name == targetName) {
            return p;
        }
    }
    throw std::runtime_error(std::format("Parameter not found: {}", targetName));
}

#endif
