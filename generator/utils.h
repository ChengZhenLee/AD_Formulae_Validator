#ifndef UTILS_H
#define UTILS_H

#include "generator.h"

std::vector<Param> generateParameters(int order, std::string sequence);

std::string generateNestedADType(int order, std::string sequence);

std::string generateXNestedAdType(int order, std::string sequence);

std::string generateYNestedAdType(int order, std::string sequence);


#endif
