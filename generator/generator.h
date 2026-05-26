#ifndef GENERATOR_H
#define GENERATOR_H


#include <string>

void generateADHeader();

void generateADDrivers();

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

#endif