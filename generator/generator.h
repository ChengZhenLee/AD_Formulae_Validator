#ifndef GENERATOR_H
#define GENERATOR_H


#include <string>

void generateADHeader(std::string filename);

void generateADDrivers(std::string filename);

std::string generateInterface(
    std::string sequence, 
    std::string XADNested, 
    std::string YADNested
);

std::string generateTangent(
    size_t currentOrder, 
    std::string sequence,
    std::string XADNested,
    std::string YADNested
);

std::string generateAdjoint(
    size_t currentOrder, 
    std::string sequence,
    std::string XADNested,
    std::string YADNested
);

std::string generateMain(
    std::string sequence
);

#endif