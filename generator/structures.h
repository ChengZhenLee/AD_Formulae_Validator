#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <string>


// Serialized Tensor
struct Tensor {
    std::vector<double> data;
    std::vector<size_t> shape;
    std::vector<size_t> strides;

    // Initializer
    Tensor(std::vector<size_t> dims) : shape(dims) {
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

// Input parameters for AD and Formula functions
struct Param {
    Tensor tensorShape;
    std::string name;

    Param(std::vector<size_t> dims, std::string inputName) 
        : tensorShape(dims), name(inputName) {
    }
};

#endif