#ifndef UTILS_H
#define UTILS_H

#include "generator.h"

std::vector<Param> generateParameters(int order, std::string sequence);

std::string generateNestedADType(int order, std::string sequence);

std::string generateXNestedADType(int order, std::string sequence);

std::string generateYNestedADType(int order, std::string sequence);

template<typename ADNested>
void seedAD(ADNested& x, const Param& p, size_t curOrder, const std::string& sequence, std::vector<size_t> coords);

template<typename ADNested>
void extractAD(ADNested& y, Param& p, size_t curOrder, std::string& sequence, std::vector<size_t> coords);

template<typename ADNested>
void seedADForOrder(ADNested& x, std::vector<Param>& parameters, size_t curOrder, std::string& sequence);

template<typename ADNested>
void extractADForOrder(ADNested& y, std::vector<Param>& parameters, size_t curOrder, std::string& sequence);

template<typename ADNested>
void seedPrimal(ADNested& x, std::vector<Param> parameters, std::vector<size_t> coords);

template<typename ADNested>
void extractPrimal(ADNested& x, std::vector<Param> parameters, std::vector<size_t> coords);

std::string getCurrentLayerADType(std::string ADNested, size_t curOrder, std::string sequence);

std::string getCurrentLayerFunctionName(size_t curOrder);

#endif
