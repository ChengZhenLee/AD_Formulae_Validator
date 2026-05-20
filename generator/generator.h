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
    int currentOrder, 
    std::string sequence,
    std::string XADNested,
    std::string YADNested
);

std::string generateAdjoint(
    int currentOrder, 
    std::string sequence,
    std::string XADNested,
    std::string YADNested
);

#endif