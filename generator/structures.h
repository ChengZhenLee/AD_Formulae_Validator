#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <deque>
#include <string>


// Serialized Tensor
struct Tensor {
    std::vector<double> data;
    std::deque<size_t> shape;
    std::vector<size_t> strides;

    // Initializer
    Tensor(std::deque<size_t> dims) : shape(dims) {
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

    size_t getIndex(const std::vector<size_t>& coords) {
        size_t index = 0;
        for (size_t i = 0; i < coords.size(); i++) {
            index += coords[i] * strides[i];
        }
        return index;
    }
};

enum class ParamRole { Input, Output };

// Input parameters for AD and Formula functions
struct Param {
    Tensor tensor;
    std::string name;
    ParamRole role;

    Param(std::deque<size_t> shape, std::string inputName, ParamRole inputRole) 
        : tensor(shape), name(inputName), role(inputRole) {
    }

    bool isSeed() const { return role == ParamRole::Input; }
};

#endif