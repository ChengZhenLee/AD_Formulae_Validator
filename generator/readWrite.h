#ifndef READWRITE_H
#define READWRITE_H


#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>
#include "structures.h"


// Write parameters as binary into an output file
template<typename T>
void writeParameters(const std::vector<Param<T>> &parameters, const std::string &filename) {

    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open  " << filename << "for writing.\n";
        return;
    }

    // Convert the JSON to binary Message Pack object
    nlohmann::json j = parameters;
    std::vector<uint8_t> binary_data = nlohmann::json::to_msgpack(j);
    outFile.write(reinterpret_cast<const char*>(binary_data.data()), binary_data.size());

    outFile.close();
}


// Read the parameters from a file
template<typename T>
std::vector<Param<T>> readParameters(const std::string &filename) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile.is_open()) {
        std::cerr << "Failed to open  " << filename << "for reading.\n";
        return {};
    }

    // Read raw file bytes
    std::vector<uint8_t> binary_data(
        (std::istreambuf_iterator<char>(inFile)), 
        std::istreambuf_iterator<char>());

    // Unpack into json object
    nlohmann::json j = nlohmann::json::from_msgpack(binary_data);

    inFile.close();

    std::vector<Param<T>> result;
    if (!j.is_array()) return result;

    for (const auto &elem : j) {
        // Extract fields
        std::map<size_t, size_t> orderedShape = elem.at("orderedShape").get<std::map<size_t, size_t>>();
        std::deque<size_t> shape = elem.at("tensor").at("shape").get<std::deque<size_t>>();
        std::string name = elem.at("name").get<std::string>();
        ParamRole role = elem.at("role").get<ParamRole>();
        std::deque<std::string> indexNames = elem.at("indexNames").get<std::deque<std::string>>();

        // Construct Param using available constructor
        Param<T> p(orderedShape, shape, name, role, indexNames);

        // Fill tensor data
        if (elem.at("tensor").contains("data")) {
            p.tensor.data = elem.at("tensor").at("data").get<std::deque<T>>();
        }

        // Overwrite activeOrders and highestOrder if present
        if (elem.contains("activeOrders")) p.activeOrders = elem.at("activeOrders").get<std::vector<int>>();
        if (elem.contains("highestOrder")) p.highestOrder = elem.at("highestOrder").get<int>();

        result.push_back(std::move(p));
    }

    return result;
}


#endif