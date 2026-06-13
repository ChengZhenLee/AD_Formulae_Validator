#ifndef UTILS_H
#define UTILS_H

#include "structures.h"
#include <deque>
#include <string>
#include <format>

template<typename T>
std::vector<Param<T>> generateParameters(std::string sequence);

template<typename T>
std::vector<Param<T>> generateRandomSeeds(std::string sequence, const X_t<T>& x);

template<typename T>
std::vector<Param<T>> generateDerivativeSeeds(std::string sequence, const X_t<T>& x);

template<typename ADNested, typename T>
void seedAD(ADNested& x, const Param<T>& p, size_t curOrder, const std::string& sequence, std::deque<size_t> coords);

template<typename ADNested, typename T>
void extractAD(ADNested& y, Param<T>& p, size_t curOrder, const std::string& sequence, std::deque<size_t> coords);

template<typename ADNested, typename T>
void seedADForOrder(ADNested& x, std::vector<Param<T>>& parameters, size_t curOrder, const std::string& sequence);

template<typename ADNested, typename T>
void extractADForOrder(ADNested& y, std::vector<Param<T>>& parameters, size_t curOrder, const std::string& sequence);

template<typename ADNested, typename T>
void seedPrimal(ADNested& x, const std::vector<Param<T>>& parameters);

template<typename ADNested, typename T>
void extractPrimal(ADNested& x, const std::vector<Param<T>>& parameters);

template <typename T>
Param<T>& findParamByName(const std::string& targetName, std::deque<Param<T>>& parameters);

template <typename T>
Param<T>& findParamByName(const std::string& targetName, const std::deque<Param<T>>& parameters);

template<typename T>
std::vector<Param<T>> getDerivatives(std::string sequence, const std::vector<T>& x0, T h = static_cast<T>(1e-6));

inline std::string makeParamName(const std::string& baseName, size_t order) {
    if (order == 0) {
        return baseName;
    }
    return std::format("{}_{}", baseName, order);
}

#include "utils.hpp"

#endif
