#ifndef GENERATOR_H
#define GENERATOR_H


#include "structures.h"

std::string generateDrivers(std::string sequence);

std::string generateTangent(
    int currentOrder, 
    std::string sequence,
    std::string ADNested
);

std::string generateAdjoint(
    int currentOrder, 
    std::string sequence,
    std::string ADNested
);

#endif