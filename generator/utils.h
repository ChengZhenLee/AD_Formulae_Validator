#ifndef UTILS_H
#define UTILS_H


#include "structures.h"

template<typename T>
std::vector<Param<T>> generateParameters(int order, std::string sequence);

std::string generateNestedADType(int order, std::string sequence);

std::string generateXNestedADType(int order, std::string sequence);

std::string generateYNestedADType(int order, std::string sequence);

template<typename ADNested, typename T>
void seedAD(ADNested& x, const Param<T>& p, size_t curOrder, const std::string& sequence, std::vector<size_t> coords);

template<typename ADNested, typename T>
void extractAD(ADNested& y, Param<T>& p, size_t curOrder, std::string& sequence, std::vector<size_t> coords);

template<typename ADNested, typename T>
void seedADForOrder(ADNested& x, std::vector<Param<T>>& parameters, size_t curOrder, std::string& sequence);

template<typename ADNested, typename T>
void extractADForOrder(ADNested& y, std::vector<Param<T>>& parameters, size_t curOrder, std::string& sequence);

template<typename ADNested, typename T>
void seedPrimal(ADNested& x, std::vector<Param<T>> parameters, std::vector<size_t> coords);

template<typename ADNested, typename T>
void extractPrimal(ADNested& x, std::vector<Param<T>> parameters, std::vector<size_t> coords);

std::string getCurrentLayerADType(std::string ADNested, size_t curOrder, std::string sequence);

std::string getCurrentLayerFunctionName(size_t curOrder);

std::string generateRegisterInputString(size_t curOrder);

std::string generateResetTapeString(std::string sequence);

template<typename T>
Param<T> findParamByName(std::string targetName, std::deque<Param<T>> parameters);

template<typename T>
std::deque<Param<T>> getDerivatives(size_t order);

#endif
