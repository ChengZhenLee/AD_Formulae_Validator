#ifndef VALIDATOR_H
#define VALIDATOR_H

#include "readWrite.h"
#include "structures.h"
#include <fstream>
#include <iostream>

template<typename T>
void validateParameters(std::vector<Param<T>> formulaResults, 
    std::string inputFilename, std::string outputFilename) {
    std::vector<Param<T>> adResults = readParameters<T>(inputFilename);
    std::ofstream outFile(outputFilename);

    if (!outFile.is_open()) {
        std::cerr << "Unable to open file: " << outputFilename << "\n";
    }
    
    for (auto& fp : formulaResults) {
        for (auto& adp : adResults) {
            if (fp.name == adp.name) {
                T tolerance = std::pow(std::numeric_limits<T>::epsilon(), 0.5);
                std::string reportLine = std::format("{}: ", fp.name);
                reportLine += compareTensors(fp.tensor, adp.tensor, tolerance) + "\n";
                outFile << reportLine;
            }
        }
    }
}


// Returns a report string 
template<typename T>
std::string compareTensors(Tensor<T> a, Tensor<T> b, T tolerance) {
    // Compare shape
    if (a.shape.size() != b.shape.size()) return "Incorrect tensor shape";

    for (int i = 0; i < a.shape.size(); i++) {
        if (a.shape[i] != b.shape[i]) return "Incorrect tensor shape";
    }

    // Compare the actual values
    if (a.data.size() != b.data.size()) return "Incorrect tensor data size";

    for (int i = 0; i < a.data.size(); i++) {
        if (abs(a.data[i] - b.data[i]) > tolerance) return "Incorrect data values";
    }

    return "valid";
}

#endif