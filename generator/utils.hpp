#ifndef UTILS_HPP
#define UTILS_HPP

#include "structures.h"
#include "configManager.h"
#include <algorithm>
#include <cmath>
#include <format>
#include <random>
#include <stdexcept>
#include <type_traits>

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

    char curMode = sequence[0];
    std::vector<Param<T>> result = generateParameters<T>(sequence.substr(1));

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
        } else if (p.name.rfind("X", 0) == 0) {
            for (int i = 0; i < p.tensor.data.size(); i++) {
                p.tensor.data[i] = dist(gen);
            }
        }
    }

    return parameters;
}

// Calculate derivatives using finite differences
template<typename T>
std::vector<Param<T>> getDerivatives(std::string sequence, const std::vector<T>& x0, T h) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t order = sequence.length();
    std::vector<Param<T>> derivatives;
    size_t xShape = cm.getXShape();
    size_t yShape = cm.getYShape();

    if (x0.size() < xShape) {
        throw std::runtime_error("Primal x vector too small for derivative computation");
    }

    auto eval_f = [&](const std::vector<T>& x)->std::vector<T> {
        std::vector<T> y(yShape);
        f<T>(const_cast<std::vector<T>&>(x), y);
        return y;
    };

    if (order >= 1) {
        for (size_t k = 1; k <= order; ++k) {
            std::deque<size_t> shape;
            shape.push_back(yShape);
            for (size_t t = 0; t < k; ++t) shape.push_back(xShape);

            Param<T> pk(std::map<size_t, size_t>{}, shape, std::format("F_{}", k), ParamRole::Input);
            std::vector<size_t> idx(k, 0);

            auto compute_and_store = [&](const std::vector<size_t>& idx) {
                std::vector<T> acc(yShape, static_cast<T>(0));
                size_t combos = 1u << k;
                for (size_t mask = 0; mask < combos; ++mask) {
                    T prod_sign = static_cast<T>(1);
                    std::vector<T> x_shift = x0;
                    for (size_t t = 0; t < k; ++t) {
                        int s = ((mask >> t) & 1) ? 1 : -1;
                        prod_sign *= static_cast<T>(s);
                        x_shift[idx[t]] += static_cast<T>(s) * h;
                    }
                    std::vector<T> y_shift = eval_f(x_shift);
                    for (size_t j = 0; j < yShape; ++j) acc[j] += prod_sign * y_shift[j];
                }

                T denom = std::pow(static_cast<T>(2), static_cast<T>(k)) * std::pow(h, static_cast<T>(k));
                for (size_t j = 0; j < yShape; ++j) acc[j] /= denom;

                for (size_t j = 0; j < yShape; ++j) {
                    std::deque<size_t> coords;
                    coords.push_back(j);
                    for (size_t t = 0; t < k; ++t) coords.push_back(idx[t]);
                    size_t index = pk.tensor.getIndex(coords);
                    pk.tensor.data[index] = acc[j];
                }
            };

            bool done = false;
            while (!done) {
                compute_and_store(idx);
                for (size_t pos = 0; pos < k; ++pos) {
                    if (++idx[pos] < xShape) break;
                    idx[pos] = 0;
                    if (pos == k - 1) done = true;
                }
                if (k == 0) break;
            }

            derivatives.push_back(std::move(pk));
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
                    coords.push_back(i);
                    extractAD(y.tangent(i), p, curOrder + 1, sequence, coords);
                }
            }
            
        } else if (curType == 'a') {
            // COMPILER GUARD: Only compile this loop if y actually has a .adjoint() method
            if constexpr(requires { y.adjoint(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    coords.push_front(i);
                    extractAD(y.adjoint(i), p, curOrder + 1, sequence, coords);
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
    size_t xShape = cm.getXShape();

    for (Param<T>& p : parameters) {
        if (p.highestOrder == targetOrder && p.isActive(targetOrder)) {
            for (size_t i = 0; i < xShape; i++) {
                seedAD(x[i], p, 1, sequence, {i});
            }
        }
    }
}

// Extracts all values for a specific order/layer
template<typename ADNested, typename T>
void extractADForOrder(ADNested& y, std::vector<Param<T>>& parameters, size_t targetOrder, const std::string& sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t yShape = cm.getXShape();

    for (Param<T>& p : parameters) {
        if (p.highestOrder == targetOrder && p.isActive(targetOrder)) {
            for (size_t j = 0; j < yShape; j++) {
                extractAD(y[j], p, 1, sequence, {j});
            }
        }
    }
}

// Seed the primal value for the nested type
template<typename ADNested, typename T>
void seedPrimal(ADNested& x, const std::vector<Param<T>>& parameters) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t xShape = cm.getXShape();

    for (const Param<T>& p : parameters) {
        if (p.name.length() == 1 && p.name[0] == 'Y') {
            for (size_t i = 0; i < xShape; i++) {
                seedElement(x[i], p, i);
            }
        }
    }
}

// Extract the primal value for the nested type
template<typename ADNested, typename T>
void extractPrimal(ADNested& y, std::vector<Param<T>>& parameters) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t yShape = cm.getYShape();

    for (Param<T>& p : parameters) {
        if (p.name.length() == 1 && p.name[0] == 'Y') {
            for (size_t j = 0; j < yShape; j++) {
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
