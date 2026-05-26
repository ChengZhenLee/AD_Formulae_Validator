#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <deque>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <limits>
#include <cmath>
#include "ad.h"


// AD tangent type
template<typename T, int size>
using T_t=ad::tangent_t<T, size>;

// AD adjoint type
template<typename T, int size>
using A_t=ad::adjoint_t<T, size>;

// input type
template<typename T>
using X_t=std::deque<T>;

// output type
template<typename T>
using Y_t=std::deque<T>;

// Serialized Tensor
template<typename T>
struct Tensor {
    std::vector<T> data;
    std::deque<size_t> shape;
    std::vector<size_t> strides;

    // Initializer
    Tensor(std::deque<size_t> dims) : 
        shape(dims), strides(dims.size()) {
        // Assign data to 0.0
        size_t total_size = 1;
        for (auto d : dims) total_size *= d;
        data.assign(total_size, 0.0);

        // Calculate the strides
        size_t currentStride = 1;
        for (int i = dims.size() - 1; i >= 0; i--) {
            strides[i] = currentStride;
            currentStride *= dims[i];
        }
    }

    size_t getIndex(const std::deque<size_t>& coords) const {
        size_t index = 0;
        for (size_t i = 0; i < coords.size(); i++) {
            index += coords[i] * strides[i];
        }
        return index;
    }

    // Addition
    Tensor<T> operator+(const Tensor& other) const {
        if (shape != other.shape) {
            throw std::runtime_error("Tensor shapes must match for addition.");
        }

        Tensor result(shape);
        for (size_t i = 0; i < data.size(); i++) {
            result.data[i] = this->data[i] + other.data[i];
        }
        return result;
    }

    // Matrix/tensor contraction: contract last dim of `a` with first dim of `b`
    // This convenience overload delegates to the general contraction below.
    static Tensor<T> product(const Tensor& a, const Tensor& b) {
        // default: contract a's last axis with b's first axis
        std::vector<size_t> contractA = { a.shape.size() - 1 };
        std::vector<size_t> contractB = { 0 };
        return product(a, b, contractA, contractB);
    }

    // General contraction: contract axes specified in contractA with corresponding axes in contractB
    static Tensor<T> product(const Tensor& a, const Tensor& b, const std::vector<size_t>& contractA, const std::vector<size_t>& contractB) {
        if (contractA.size() != contractB.size()) {
            throw std::invalid_argument("contractA and contractB must have the same length");
        }

        size_t r = contractA.size();
        // validate axes and matching dims
        for (size_t i = 0; i < r; ++i) {
            if (contractA[i] >= a.shape.size() || contractB[i] >= b.shape.size())
                throw std::out_of_range("contraction axis out of range");
            if (a.shape[contractA[i]] != b.shape[contractB[i]])
                throw std::invalid_argument("contracted dimensions must match");
        }

        // Build lists of non-contracted axes
        std::vector<size_t> nonA, nonB;
        for (size_t i = 0; i < a.shape.size(); ++i) {
            if (std::find(contractA.begin(), contractA.end(), i) == contractA.end()) nonA.push_back(i);
        }
        for (size_t i = 0; i < b.shape.size(); ++i) {
            if (std::find(contractB.begin(), contractB.end(), i) == contractB.end()) nonB.push_back(i);
        }

        // Result shape: dims of nonA (in order) followed by dims of nonB (in order)
        std::deque<size_t> resShape;
        for (size_t ax : nonA) resShape.push_back(a.shape[ax]);
        for (size_t ax : nonB) resShape.push_back(b.shape[ax]);

        Tensor result(resShape);

        // Iterators over result coordinates
        std::vector<size_t> resCoords(result.shape.size(), 0);
        bool doneRes = result.shape.empty();
        if (!doneRes) doneRes = false;

        auto increment = [&](std::vector<size_t>& coords, const std::deque<size_t>& shape)->bool {
            for (size_t i = 0; i < coords.size(); ++i) {
                coords[i]++;
                if (coords[i] < shape[i]) return false;
                coords[i] = 0;
            }
            return true; // wrapped around
        };

        // handle empty result (scalar)
        if (result.shape.empty()) {
            // only contracted sum
            std::vector<size_t> contractIdx(r, 0);
            bool doneC = (r == 0);
            T sum = static_cast<T>(0);
            while (!doneC) {
                // build coords for a and b
                std::deque<size_t> coordsA(a.shape.size());
                for (size_t i = 0; i < nonA.size(); ++i) coordsA[nonA[i]] = 0;
                for (size_t i = 0; i < r; ++i) coordsA[contractA[i]] = contractIdx[i];

                std::deque<size_t> coordsB(b.shape.size());
                for (size_t i = 0; i < r; ++i) coordsB[contractB[i]] = contractIdx[i];
                for (size_t i = 0; i < nonB.size(); ++i) coordsB[nonB[i]] = 0;

                sum += a.data[a.getIndex(coordsA)] * b.data[b.getIndex(coordsB)];

                // increment contractIdx
                for (size_t p = 0; p < r; ++p) {
                    contractIdx[p]++;
                    if (contractIdx[p] < a.shape[contractA[p]]) break;
                    contractIdx[p] = 0;
                    if (p == r - 1) doneC = true;
                }
                if (r == 0) break;
            }
            result.data[0] = sum;
            return result;
        }

        // General non-scalar result
        bool finished = false;
        while (!finished) {
            // Build coordsA and coordsB for current resCoords
            std::deque<size_t> coordsA(a.shape.size());
            std::deque<size_t> coordsB(b.shape.size());

            // fill non-contracted axes from resCoords
            size_t pos = 0;
            for (size_t i = 0; i < nonA.size(); ++i) {
                coordsA[nonA[i]] = resCoords[pos++];
            }
            for (size_t i = 0; i < nonB.size(); ++i) {
                coordsB[nonB[i]] = resCoords[pos++];
            }

            // Sum over contracted indices
            std::vector<size_t> contractIdx(r, 0);
            bool doneC = (r == 0);
            T sum = static_cast<T>(0);
            while (!doneC) {
                for (size_t i = 0; i < r; ++i) {
                    coordsA[contractA[i]] = contractIdx[i];
                    coordsB[contractB[i]] = contractIdx[i];
                }

                sum += a.data[a.getIndex(coordsA)] * b.data[b.getIndex(coordsB)];

                // increment contractIdx
                for (size_t p = 0; p < r; ++p) {
                    contractIdx[p]++;
                    if (contractIdx[p] < a.shape[contractA[p]]) break;
                    contractIdx[p] = 0;
                    if (p == r - 1) doneC = true;
                }
                if (r == 0) break;
            }

            // store sum into result
            std::deque<size_t> coordsR(result.shape.size());
            for (size_t i = 0; i < result.shape.size(); ++i) coordsR[i] = resCoords[i];
            result.data[result.getIndex(coordsR)] = sum;

            // increment resCoords
            finished = increment(resCoords, result.shape);
        }

        return result;
    }

    // Compute the result shape for a tensor chain without materializing the data
    static std::deque<size_t> productShape(const std::deque<size_t>& aShape, const std::deque<size_t>& bShape) {
        if (aShape.empty()) return bShape;
        if (bShape.empty()) return aShape;
        if (aShape.back() != bShape.front()) return {};

        std::deque<size_t> resultShape;
        for (size_t i = 0; i < aShape.size() - 1; ++i) resultShape.push_back(aShape[i]);
        for (size_t i = 1; i < bShape.size(); ++i) resultShape.push_back(bShape[i]);
        return resultShape;
    }

    static std::deque<size_t> productShape(const std::deque<Tensor<T>>& chain) {
        if (chain.empty()) return {};
        std::deque<size_t> resultShape = chain.front().shape;
        for (size_t i = 1; i < chain.size(); ++i) {
            resultShape = productShape(resultShape, chain[i].shape);
            if (resultShape.empty()) return {};
        }
        return resultShape;
    }

    // Multiply a list of tensors
    static Tensor<T> product(std::deque<Tensor<T>> tensors) {
        if (tensors.empty()) {
            throw std::invalid_argument("Cannot compute product of an empty tensor list.");
        }

        // Start with the first tensor in the deque
        Tensor<T> result = tensors.front();
        tensors.pop_front();

        // Iteratively multiply the current result by the next tensor in the list
        while (!tensors.empty()) {
            result = product(result, tensors.front());
            tensors.pop_front();
        }

        return result;
    }


    static bool compareTensors(
        const Tensor<T>& a, const Tensor<T>& b,
        double tolerance = std::numeric_limits<T>::epsilon()
    ) {
        if (a.shape != b.shape) return false;

        size_t length = a.data.size();

        for (size_t i = 0; i < length; i++) {
            if (std::abs(static_cast<double>(a.data[i]) - static_cast<double>(b.data[i])) > tolerance)
                return false;
        }

        return true;
    }

    // Validate that a deque of tensors can be multiplied left-to-right
    static bool validateTensorChain(const std::deque<Tensor<T>>& chain) {
        if (chain.size() < 2) return true;
        for (size_t i = 0; i + 1 < chain.size(); ++i) {
            size_t left = chain[i].shape.empty() ? 0 : chain[i].shape.back();
            size_t right = chain[i+1].shape.empty() ? 0 : chain[i+1].shape.front();
            if (left != right) return false;
        }
        return true;
    }
};

enum class ParamRole { Input, Output };

// Input parameters for AD and Formula functions
template<typename T>
struct Param {
    Tensor<T> tensor;
    std::string name;
    ParamRole role;
    std::vector<int> activeOrders = {};
    std::map<size_t, size_t> orderedShape;
    int highestOrder = 0;

    Param(std::map<size_t, size_t> oShape, std::deque<size_t> shape, std::string inputName, ParamRole inputRole) 
        : orderedShape(oShape), tensor(shape), name(inputName), role(inputRole) {
        
        std::stringstream ss(inputName);
        std::string segment;

        // Separate the string by the delimiter '_'
        while (std::getline(ss, segment, '_')) {
            // If the current string is a number and not 'X' or 'Y'
            if (!segment.empty() && std::isdigit(segment[0])) {
                int order = std::stoi(segment);
                activeOrders.push_back(order);

                // Set the highest order
                if (order > highestOrder) { 
                    highestOrder = order; 
                }
            }
        }
    }

    bool isSeed() { 
        return role == ParamRole::Input; 
    }

    bool isActive(int order) {
        return std::find(activeOrders.begin(), activeOrders.end(), order) != activeOrders.end();
    }

    size_t getShapeAt(int order) {
        return orderedShape[order] ? orderedShape[order] : 0;
    }
};


template<typename T>
struct Monomial {
    std::deque<Param<T>> parameters = {};

    Monomial(std::deque<Param<T>> newParameters) : parameters(newParameters) {};
};


template<typename T>
struct Equation {
    Param<T> leftSide;
    std::deque<Monomial<T>> rightSide = {};

    Equation(Param<T> left) : leftSide(left) {};
};

#endif