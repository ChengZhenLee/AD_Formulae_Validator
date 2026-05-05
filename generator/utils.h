#ifndef UTILS_H
#define UTILS_H

#include "generator.h"

std::vector<Param> generateParameters(int order, std::string sequence);

std::string generateNestedADType(int order, std::string sequence);

std::string generateXNestedADType(int order, std::string sequence);

std::string generateYNestedADType(int order, std::string sequence);

template<typename ADNested>
void seedAD(ADNested& x, const Param& p, size_t curOrder, const std::string& sequence, std::vector<size_t> coords);

#endif
