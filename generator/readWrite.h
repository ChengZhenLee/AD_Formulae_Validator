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
    std::cout << "Successfully saved to " << filename << "\n";
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

    return j.get<std::vector<Param<T>>>();
}