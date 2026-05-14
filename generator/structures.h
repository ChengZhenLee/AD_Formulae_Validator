#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <deque>
#include <string>
#include <sstream>
#include <map>


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

    size_t getIndex(const std::deque<size_t>& coords) {
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

    // Matrix multiplication
    static Tensor<T> product(const Tensor& a, const Tensor& b) {
        if (a.shape.back() != b.shape.front()) {
            throw std::invalid_argument("Inner dimensions must match for Tensor");
        }

        // Determine result shape
        std::deque<size_t> resShape;
        for (size_t i = 0; i < a.shape.size() - 1; i++) resShape.push_back(a.shape[i]);
        for (size_t i = 1; i < b.shape.size(); i++) resShape.push_back(b.shape[i]);

        Tensor result(resShape);

        // Calculate size for the loops
        size_t commonDim = a.shape.back();
        size_t numOuterA = a.data.size() / commonDim; // product of (left) outer shapes
        size_t numOuterB = b.data.size() / commonDim; // product of (right) outer shapes

        for (size_t i = 0; i < numOuterA; ++i) {
            for (size_t j = 0; j < numOuterB; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < commonDim; ++k) {
                    // Map flat indices using the strides
                    // a_index = (outer_a_offset) + k * last_stride_of_a
                    // b_index = k * first_stride_of_b + (outer_b_offset)
                    sum += a.data[i * commonDim + k] * 
                           b.data[k * b.strides[0] + j];
                }
                result.data[i * numOuterB + j] = sum;
            }
        }
        return result;
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