#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <deque>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <nlohmann/json.hpp>
#include "ad.h"

// Shared data model for the symbolic ("formula driver") side of the
// validator. A differentiated equation system is represented as:
//   Equation   := leftSide (a Param) = sum of Monomial
//   Monomial   := product of Param, contracted over shared index names
//   Param      := a named tensor with Einstein-summation index labels
// This mirrors how the chain rule is applied on paper: each differentiation
// order adds one more contracted index to every term. See
// generator/symbolic/formulaDriver.hpp for how these are built up order-by-order.

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


// A dense, row-major, arbitrary-rank tensor: flat `data` plus `shape`
// (extent per axis) and `strides` (row-major offset per axis) needed to
// convert multi-index coordinates to a flat offset into `data`.
template<typename T>
struct Tensor {
    std::deque<T> data;
    std::deque<size_t> shape;
    std::vector<size_t> strides;

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Tensor, data, shape, strides);

    Tensor() = default;
    
    Tensor(std::deque<size_t> dims) : shape(dims), strides(dims.size()) {
        size_t total_size = 1;
        for (auto d : dims) total_size *= d;
        data.assign(total_size, static_cast<T>(0));

        size_t currentStride = 1;
        for (int i = static_cast<int>(dims.size()) - 1; i >= 0; i--) {
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

    // Einstein-summation contraction of two tensors along paired axes:
    // contractA[k] of `a` is summed against contractB[k] of `b` for every k.
    // The result's axes are the uncontracted axes of `a` (in order) followed
    // by the uncontracted axes of `b` (in order) — e.g. contracting a
    // rank-2 `a` and rank-1 `b` over one axis each reproduces a matrix-vector
    // product.
    static Tensor<T> productGeneralContraction(const Tensor& a, const Tensor& b, const std::vector<size_t>& contractA, const std::vector<size_t>& contractB) {
        if (contractA.size() != contractB.size()) {
            throw std::invalid_argument("contractA and contractB must have the same length");
        }

        size_t r = contractA.size();
        for (size_t i = 0; i < r; ++i) {
            if (contractA[i] >= a.shape.size() || contractB[i] >= b.shape.size())
                throw std::out_of_range("contraction axis out of range");
            if (a.shape[contractA[i]] != b.shape[contractB[i]])
                throw std::invalid_argument("contracted dimensions must match");
        }

        std::vector<size_t> nonA, nonB;
        for (size_t i = 0; i < a.shape.size(); ++i) {
            if (std::find(contractA.begin(), contractA.end(), i) == contractA.end()) nonA.push_back(i);
        }
        for (size_t i = 0; i < b.shape.size(); ++i) {
            if (std::find(contractB.begin(), contractB.end(), i) == contractB.end()) nonB.push_back(i);
        }

        std::deque<size_t> resShape;
        for (size_t ax : nonA) resShape.push_back(a.shape[ax]);
        for (size_t ax : nonB) resShape.push_back(b.shape[ax]);

        Tensor result(resShape);
        std::vector<size_t> resCoords(result.shape.size(), 0);

        auto increment = [](std::vector<size_t>& coords, const std::deque<size_t>& targetShape) -> bool {
            for (size_t i = 0; i < coords.size(); ++i) {
                coords[i]++;
                if (coords[i] < targetShape[i]) return false;
                coords[i] = 0;
            }
            return true;
        };

        bool finished = false;
        while (!finished) {
            std::deque<size_t> coordsA(a.shape.size(), 0);
            std::deque<size_t> coordsB(b.shape.size(), 0);

            size_t pos = 0;
            for (size_t i = 0; i < nonA.size(); ++i) coordsA[nonA[i]] = resCoords[pos++];
            for (size_t i = 0; i < nonB.size(); ++i) coordsB[nonB[i]] = resCoords[pos++];

            std::vector<size_t> contractIdx(r, 0);
            bool doneC = (r == 0);
            T sum = static_cast<T>(0);

            while (!doneC) {
                for (size_t i = 0; i < r; ++i) {
                    coordsA[contractA[i]] = contractIdx[i];
                    coordsB[contractB[i]] = contractIdx[i];
                }

                sum += a.data[a.getIndex(coordsA)] * b.data[b.getIndex(coordsB)];

                for (size_t p = 0; p < r; ++p) {
                    contractIdx[p]++;
                    if (r > 0 && contractIdx[p] < a.shape[contractA[p]]) break;
                    contractIdx[p] = 0;
                    if (p == r - 1) doneC = true;
                }
                if (r == 0) break;
            }

            if (result.shape.empty()) {
                result.data[0] = sum;
                break;
            } else {
                std::deque<size_t> coordsR(result.shape.size());
                for (size_t i = 0; i < result.shape.size(); ++i) coordsR[i] = resCoords[i];
                result.data[result.getIndex(coordsR)] = sum;
            }

            finished = increment(resCoords, result.shape);
        }

        return result;
    }

    static Tensor<T> permuteAxes(const Tensor<T>& source, const std::vector<size_t>& permutation) {
        // 1. Safety Gate: Verify permutation dimensions match source rank
        if (permutation.size() != source.shape.size()) {
            throw std::invalid_argument("Permutation rank must match tensor rank.");
        }

        // 2. Determine the new shape and calculate the original tensor's strides
        std::deque<size_t> newShape(source.shape.size());
        for (size_t i = 0; i < permutation.size(); ++i) {
            newShape[i] = source.shape[permutation[i]];
        }

        std::vector<size_t> oldStrides = source.strides;
        
        // 3. Allocate the new output tensor with the re-ordered shape boundary
        Tensor<T> result(newShape); // Assumes constructor allocates underlying flat data vector
        
        // 4. Multi-dimensional coordinate tracker
        std::vector<size_t> currentCoords(newShape.size(), 0);
        size_t totalElements = source.data.size(); // Total flat elements (e.g., d0 * d1 * d2)

        for (size_t outFlatIndex = 0; outFlatIndex < totalElements; ++outFlatIndex) {
            // Find corresponding coordinate index mapping in the original tensor space
            size_t inFlatIndex = 0;
            for (size_t i = 0; i < permutation.size(); ++i) {
                // Map the output dimension's current coordinate back to the original stride factor
                inFlatIndex += currentCoords[i] * oldStrides[permutation[i]];
            }

            // Copy data point across the memory boundary
            result.data[outFlatIndex] = source.data[inFlatIndex];

            // Increment multi-dimensional coordinates manually (Row-Major odometer pattern)
            for (int i = static_cast<int>(newShape.size()) - 1; i >= 0; --i) {
                currentCoords[i]++;
                if (currentCoords[i] < newShape[i]) {
                    break; // Coordinate successfully advanced
                }
                currentCoords[i] = 0; // Roll-over axis index to 0, carry over to next higher dimension
            }
        }

        return result;
    }
};

enum class ParamRole { Input, Output };

NLOHMANN_JSON_SERIALIZE_ENUM(ParamRole, {
    {ParamRole::Input, "INPUT"},
    {ParamRole::Output, "OUTPUT"},
});

// A named node in the differentiation tree: the primal X/Y, a derivative
// seed, or a derived quantity, together with the tensor data and the
// Einstein-summation index labels.
//
// Names encode the derivative path they were created by, e.g. 
// "X_1_2" reads as "X differentiated at order 1, then again at order 2".
// `activeOrders` is parsed straight out of that name 
// (every "_<n>" suffix) and `highestOrder` is its max. 
// Symbolic derivative results carry the name prefix "F"
// instead of "Y" (see generateDerivativeSeeds()) purely to avoid colliding
// with the true primal Y = f(X).
//
// `indexNames` are the Einstein-summation labels for each tensor axis, in
// the same order as `tensor.shape`: "i" for the input axis, "j" for the
// output axis, "v_k"/"u_k" for the tangent/adjoint seed axis introduced at
// order k.
template<typename T>
struct Param {
    Tensor<T> tensor;
    std::string name;
    ParamRole role;
    std::vector<int> activeOrders = {};
    std::map<size_t, size_t> orderedShape;
    std::deque<std::string> indexNames = {};
    int highestOrder = 0;

    // Macro automatically generates to_json() and from_json()
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Param, tensor, name, role, activeOrders, orderedShape, highestOrder, indexNames)

    Param() = default;

    Param(std::map<size_t, size_t> oShape, std::deque<size_t> shape, std::string inputName, ParamRole inputRole, std::deque<std::string> indexNames)
        : orderedShape(oShape), tensor(shape), name(inputName), role(inputRole), indexNames(indexNames) {

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

    // Whether this Param's name has a "_<order>" suffix for the given order,
    // i.e. whether it was produced by (or seeds) that differentiation step.
    bool isActive(int order) const {
        return std::find(activeOrders.begin(), activeOrders.end(), order) != activeOrders.end();
    }

    // Extent of the axis introduced at the given differentiation order
    // (the tangent width V or adjoint width U configured for that order).
    size_t getShapeAt(int order) const {
        return orderedShape.at(order);
    }
};


// One chain-rule term: a product of Params, implicitly contracted over
// whatever index names they share (see evaluateMonomialSequence in
// formulaDriver.hpp). E.g. the chain rule term dF/dX * dX/dv for a tangent
// step is represented as Monomial{F_param, X_param}, contracted over the
// input index "i" that both share.
template<typename T>
struct Monomial {
    std::deque<Param<T>> parameters = {};

    Monomial(std::deque<Param<T>> newParameters) : parameters(newParameters) {};
};


// One derived equation: leftSide = sum of rightSide monomials. The full
// symbolic derivation for a sequence is a std::deque<Equation<T>>, one
// entry per Param produced so far, accumulated in formulaDriver().
template<typename T>
struct Equation {
    Param<T> leftSide;
    std::deque<Monomial<T>> rightSide = {};

    Equation(Param<T> left) : leftSide(left) {};
};

#endif