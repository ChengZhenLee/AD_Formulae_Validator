#ifndef GENERATOR_H
#define GENERATOR_H


#include "structures.h"

std::string generateTangent(
    int currentOrder, int totalOrder, 
    std::vector<Param> parameters
);

std::string generateAdjoint(
    int currentOrder, int totalOrder, 
    std::vector<Param> parameters
);

#endif